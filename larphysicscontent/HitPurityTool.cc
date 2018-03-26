/**
 *  @file   larphysicscontent/HitPurityTool.cc
 *
 *  @brief  Implementation of the hit purity tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/HitPurityTool.h"

#include "Pandora/AlgorithmHeaders.h"

#include "../root/Common.h" // ATTN temporary
#include "TNtuple.h"        // ATTN temporary
#include "TCanvas.h"        // ATTN temporary

using namespace pandora;

namespace lar_physics_content
{

HitPurityTool::HitPurityTool() :
    m_maxImpurityScore(1.f),
    m_valueAverageSearchRadius(20UL),
    m_nearestNeighbourNumber(5UL),
    m_makePlots(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitPurityTool::Run(const Algorithm *const pAlgorithm, LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector, float &excessCaloValue)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
       
    if (trackHitEnergyVector.size() < 2UL)
        return true;

    std::sort(trackHitEnergyVector.begin(), trackHitEnergyVector.end(),
        [](const LArTrackHitValue &lhs, const LArTrackHitValue &rhs)
        {
            return lhs.Coordinate() < rhs.Coordinate();
        });
        
    const float lowerCoordinateBound = trackHitEnergyVector.front().Coordinate();
    const float upperCoordinateBound = trackHitEnergyVector.back().Coordinate();
    const float range = upperCoordinateBound - lowerCoordinateBound;
    const float minProtectedCoordinate = lowerCoordinateBound + range * 0.05;
    const float maxProtectedCoordinate = upperCoordinateBound - range * 0.05;
    
    LArAnalysisParticleHelper::TrackHitValueVector changedTrackHitEnergies;
    const std::size_t nHits(trackHitEnergyVector.size());
    bool somethingChanged = true;
    
    while (somethingChanged)
    {
        somethingChanged = false;
        const FloatVector valueAverages = this->GetValueAverages(trackHitEnergyVector, nHits);
        
        float scaleFactor(0.f), mean(0.f), sigma(0.f);
        
        if (!this->GetStatistics(trackHitEnergyVector, nHits, scaleFactor, mean, sigma))
            return true;
        
        const FloatVector impurityScores = this->CalculateImpurityScores(trackHitEnergyVector, scaleFactor, nHits, mean, sigma);
        
        for (std::size_t i = 0UL; i < nHits; ++i)
        {
            if ((trackHitEnergyVector[i].Coordinate() < minProtectedCoordinate) || 
                (trackHitEnergyVector[i].Coordinate() > maxProtectedCoordinate))
            {
                continue;
            }
            
            const float score = impurityScores[i];
            
            if (score > m_maxImpurityScore)
            {
                if (valueAverages[i] < trackHitEnergyVector[i].CaloValue())
                {
                    changedTrackHitEnergies.push_back(trackHitEnergyVector[i]);
                    excessCaloValue += trackHitEnergyVector[i].CaloValue() - valueAverages[i];
                    trackHitEnergyVector[i].CaloValue(valueAverages[i]);
                    somethingChanged = true;
                }
            }
        }
    }
    
    // ATTN temporary
    if (m_makePlots)
        this->MakePlots(trackHitEnergyVector, changedTrackHitEnergies);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector HitPurityTool::GetValueAverages(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector,
    const std::size_t nHits) const
{
    FloatVector valueAverages(nHits, 0.f);
    
    for (int i = 0UL; i < nHits; ++i)
    {
        float valueAverage(0.f);
        std::size_t nValues(0UL);
        
        for (int j = std::max(0, i - static_cast<int>(m_valueAverageSearchRadius)); j < std::min(nHits, i + m_valueAverageSearchRadius); ++j)
        {
            if (j == i)
                continue;
                
            valueAverage += trackHitEnergyVector[j].CaloValue();
            ++nValues;
        }
        
        valueAverages[i] = valueAverage / static_cast<float>(nValues);
    }
    
    return valueAverages;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitPurityTool::GetStatistics(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits,
    float &scaleFactor, float &mean, float &sigma) const
{
    if (nHits < 3UL)
        return false;
    
    float coordinateRange(0.f), caloValueRange(0.f);
    this->CalculateRanges(trackHitEnergyVector, nHits, coordinateRange, caloValueRange);
    
    if (caloValueRange < std::numeric_limits<float>::epsilon())
        return false;
    
    scaleFactor = coordinateRange / caloValueRange;
    
    mean  = this->CalculateMean(trackHitEnergyVector, nHits, scaleFactor);
    sigma = this->CalculateStandardDeviation(trackHitEnergyVector, nHits, scaleFactor, mean);
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitPurityTool::CalculateRanges(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits,
    float &coordinateRange, float &caloValueRange) const
{
    float minCoordinate(std::numeric_limits<float>::max()), maxCoordinate(std::numeric_limits<float>::min());
    float minCaloValue(std::numeric_limits<float>::max()), maxCaloValue(std::numeric_limits<float>::min());
    
    for (std::size_t i = 0UL; i < nHits; ++i)
    {
        const float coordinate = trackHitEnergyVector[i].Coordinate();
        const float caloValue  = trackHitEnergyVector[i].CaloValue();
        
        if (coordinate < minCoordinate)
            minCoordinate = coordinate;
            
        if (coordinate > maxCoordinate)
            maxCoordinate = coordinate;
            
        if (caloValue < minCaloValue)
            minCaloValue = caloValue;
            
        if (caloValue > maxCaloValue)
            maxCaloValue = caloValue;
    }
    
    coordinateRange = maxCoordinate - minCoordinate;
    caloValueRange  = maxCaloValue - minCaloValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float HitPurityTool::CalculateMean(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits,
    const float scaleFactor) const
{
    float mean(0.f);
    std::size_t numberOfDistances(0UL);
    
    for (std::size_t i = 0UL; i < nHits; ++i)
    {
        const FloatVector distanceVector = this->GetSortedDistanceVector(trackHitEnergyVector, trackHitEnergyVector[i], scaleFactor, nHits);
        
        for (std::size_t j = 1UL; j < std::min((1UL + m_nearestNeighbourNumber), nHits - 1UL); ++j)
        {
            mean += distanceVector[j];
            ++numberOfDistances;
        }
    }
    
    return mean / static_cast<float>(numberOfDistances);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float HitPurityTool::CalculateStandardDeviation(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits,
    const float scaleFactor, const float mean) const
{
    float squaredSummedDeviation(0.f);
    std::size_t numDistances(0UL);
    
    for (std::size_t i = 0UL; i < nHits; ++i)
    {
        const FloatVector distanceVector = this->GetSortedDistanceVector(trackHitEnergyVector, trackHitEnergyVector[i], scaleFactor, nHits);
        
        for (std::size_t j = 1UL; j < std::min((1UL + m_nearestNeighbourNumber), nHits - 1UL); ++j)
        {
            const float deviation(distanceVector[j] - mean);
            squaredSummedDeviation += deviation * deviation;
            ++numDistances;
        }
    }
    
    return std::sqrt(squaredSummedDeviation / (static_cast<float>(numDistances) - 1.f));
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector HitPurityTool::CalculateImpurityScores(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector,
    const float scaleFactor, const std::size_t nHits, const float mean, const float sigma) const
{
    FloatVector impurityScores(nHits, 0.f);
        
    for (std::size_t i = 0UL; i < nHits; ++i)
    {
        const FloatVector distanceVector = this->GetSortedDistanceVector(trackHitEnergyVector, trackHitEnergyVector[i], scaleFactor, nHits);
        
        float impurityScore(0.f);
        std::size_t numberOfDistances(0UL);
        
        for (std::size_t j = 1UL; j < std::min((1UL + m_nearestNeighbourNumber), nHits - 1UL); ++j)
        {
            impurityScore += distanceVector[j];
            ++numberOfDistances;
        }
        
        impurityScores[i] = ((impurityScore / static_cast<float>(numberOfDistances)) - mean) / sigma;
    }
    
    return impurityScores;
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector HitPurityTool::GetSortedDistanceVector(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector,
    const LArTrackHitValue &currentTrackHitValue, const float scaleFactor, const std::size_t nHits) const
{
    FloatVector distanceVector(nHits, 0.f);
        
    for (std::size_t j = 0UL; j < nHits; ++j)
    {
        const float deltaCoordinate = trackHitEnergyVector[j].Coordinate() - currentTrackHitValue.Coordinate();
        const float scaledDeltaCaloValue = scaleFactor * (trackHitEnergyVector[j].CaloValue() - currentTrackHitValue.CaloValue());
        distanceVector[j] = std::sqrt(deltaCoordinate * deltaCoordinate + scaledDeltaCaloValue * scaledDeltaCaloValue);
    }
    
    std::sort(distanceVector.begin(), distanceVector.end(), std::less<float>());
    
    return distanceVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

//ATTN temporary
void HitPurityTool::MakePlots(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergyVector,
    const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergiesChanged) const
{
    static int uniquePlotIdentifier(0);
    
    TNtuple *const pNtuple        = new TNtuple("TrackHitsUncorrected", "TrackHitsUncorrected", "Coordinate:CaloValue");
    TNtuple *const pNtupleChanged = new TNtuple("TrackHitsCorrected", "TrackHitsCorrected", "Coordinate:CaloValue");

    float minCoordinate = std::numeric_limits<float>::max();
    float maxCoordinate = std::numeric_limits<float>::min();

    float minCaloValue = std::numeric_limits<float>::max();
    float maxCaloValue = std::numeric_limits<float>::min();

    for (const LArTrackHitValue &trackHitEnergy : trackHitEnergyVector)
    {
        if (trackHitEnergy.Coordinate() < minCoordinate)
            minCoordinate = trackHitEnergy.Coordinate();

         if (trackHitEnergy.Coordinate() > maxCoordinate)
            maxCoordinate = trackHitEnergy.Coordinate();

        if (trackHitEnergy.CaloValue() < minCaloValue)
            minCaloValue = trackHitEnergy.CaloValue();

        if (trackHitEnergy.CaloValue() > maxCaloValue)
            maxCaloValue = trackHitEnergy.CaloValue();

        pNtuple->Fill(trackHitEnergy.Coordinate(), trackHitEnergy.CaloValue());
    }

    for (const LArTrackHitValue &trackHitEnergy : trackHitEnergiesChanged)
    {
        if (trackHitEnergy.Coordinate() < minCoordinate)
            minCoordinate = trackHitEnergy.Coordinate();

         if (trackHitEnergy.Coordinate() > maxCoordinate)
            maxCoordinate = trackHitEnergy.Coordinate();

        if (trackHitEnergy.CaloValue() < minCaloValue)
            minCaloValue = trackHitEnergy.CaloValue();

        if (trackHitEnergy.CaloValue() > maxCaloValue)
            maxCaloValue = trackHitEnergy.CaloValue();
            
        pNtupleChanged->Fill(trackHitEnergy.Coordinate(), trackHitEnergy.CaloValue());
    }

    TCanvas *pCanvas = nullptr;

    {
        auto plotSettings        = g_defaultPlotSettings2D;
        plotSettings.plotType   = SCATTER;
        plotSettings.pointColor = kRed;
        plotSettings.xMin = minCoordinate;
        plotSettings.xMax = maxCoordinate;
        plotSettings.yMin = 0.f;
        plotSettings.yMax = maxCaloValue;
        strcpy(plotSettings.title,  "Selection of hits to Birks-correct");
        strcpy(plotSettings.xTitle, "Projected 3D track coordinate (cm)");
        strcpy(plotSettings.yTitle, "#frac{dQ}{dx} (ADC/cm)");

        pCanvas = PlotNtuple2D(pNtuple, "Coordinate", "CaloValue", std::to_string(uniquePlotIdentifier++).c_str(),
                               plotSettings);

        plotSettings.plotType   = SAME;
        plotSettings.newCanvas  = false;
        plotSettings.pointColor = kBlue;

        PlotNtuple2D(pNtupleChanged, "Coordinate", "CaloValue", std::to_string(uniquePlotIdentifier++).c_str(),
                     plotSettings);
    }

    PandoraMonitoringApi::Pause(this->GetPandora());

    if (pCanvas)
        pCanvas->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitPurityTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxImpurityScore", m_maxImpurityScore));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValueAverageSearchRadius", m_valueAverageSearchRadius));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NearestNeighbourNumber", m_nearestNeighbourNumber));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MakePlots", m_makePlots));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
