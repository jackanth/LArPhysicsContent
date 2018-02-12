/**
 *  @file   LArPhysicsContent/src/BirksHitSelectionTool.cxx
 *
 *  @brief  Implementation of the Birks' hit selection tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "BirksHitSelectionTool.h"

#include "../root/Common.h" // ATTN temporary

#include "TNtuple.h"
#include "TCanvas.h"

using namespace pandora;

namespace lar_physics_content
{

BirksHitSelectionTool::BirksHitSelectionTool() :
    m_birksSelectionBinWidth(5.f),
    m_birksSelectionSearchXWidthFraction(0.1f),
    m_birksSelectionSearchYWidth(0.1f),
    m_birksSelectionMinDensityFraction(0.25f),
    m_makePlots(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BirksHitSelectionTool::Run(const Algorithm *const pAlgorithm, LArTrackHitEnergy::Vector &trackHitEnergyVector, int &uniquePlotIdentifier)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    this->DecideWhichTrackHitsToCorrect(trackHitEnergyVector, uniquePlotIdentifier);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BirksHitSelectionTool::DecideWhichTrackHitsToCorrect(LArTrackHitEnergy::Vector &trackHitEnergyVector, int &uniquePlotIdentifier) const
{
    if (trackHitEnergyVector.empty())
    {
        std::cout << "BirksHitSelectionTool: track hit energy vector was empty" << std::endl;
        return;
    }
    
    // Get the extremal coordinates and the matrix of 2D inter-point distances.
    float minCoordinate = 0.f, maxCoordinate = 0.f, minUncorrectedEnergy = 0.f, maxUncorrectedEnergy = 0.f;
    std::tie(minCoordinate, maxCoordinate, minUncorrectedEnergy, maxUncorrectedEnergy) = this->GetExtremalBirksSelectionValues(trackHitEnergyVector);
    
    const PointDistanceMatrix matrixOfPointDistances = this->CreateMatrixOfPointDistances(trackHitEnergyVector);
     
    // Split the plot into a number of bins.
    const float coordinateRange = maxCoordinate - minCoordinate;
    
    if (coordinateRange < std::numeric_limits<float>::epsilon())
    {
        std::cout << "BirksHitSelectionTool: coordinate range was too small" << std::endl;
        return;
    }
    
    const int numBins = static_cast<int>(std::ceil(coordinateRange / this->m_birksSelectionBinWidth));
    
    LArTrackHitEnergy::Vector allTrackHitEnergiesChanged;
    bool changeMade = true;
    const float searchRadiusX = (maxCoordinate - minCoordinate) * this->m_birksSelectionSearchXWidthFraction;
    
    if (this->m_makePlots)
        this->ProduceHitSelectionPlots(trackHitEnergyVector, allTrackHitEnergiesChanged, false, true, true,
                                                    uniquePlotIdentifier);
    
    // Iteratively get rid of outliers to the distribution.
    while (changeMade)
    {
        FloatVector averageBinPointDensity, averageBinEnergy;
        std::tie(averageBinPointDensity, averageBinEnergy) = this->GetBinAverages(trackHitEnergyVector, searchRadiusX,
                                                                                  matrixOfPointDistances, numBins, minCoordinate);
        
        LArTrackHitEnergy::Vector trackHitEnergiesChanged;
        changeMade = false;
        
        this->MakeBirksSelectionChanges(trackHitEnergyVector, searchRadiusX, matrixOfPointDistances, numBins, changeMade,
                                                     minCoordinate, averageBinPointDensity, averageBinEnergy, trackHitEnergiesChanged,
                                                     allTrackHitEnergiesChanged);

        if (this->m_makePlots)
        {
            this->ProduceHitSelectionPlots(trackHitEnergyVector, trackHitEnergiesChanged, false, true, true, 
                                                        uniquePlotIdentifier);
            this->ProduceHitSelectionPlots(trackHitEnergyVector, trackHitEnergiesChanged, true, true, true, 
                                                        uniquePlotIdentifier);
        }
    }
    
    if (this->m_makePlots)
    {
        this->ProduceHitSelectionPlots(trackHitEnergyVector, allTrackHitEnergiesChanged, false, true, true,
                                                    uniquePlotIdentifier);
        this->ProduceHitSelectionPlots(trackHitEnergyVector, allTrackHitEnergiesChanged, true, true, true,
                                                    uniquePlotIdentifier);
        this->ProduceHitSelectionPlots(trackHitEnergyVector, allTrackHitEnergiesChanged, false, true, false,
                                                    uniquePlotIdentifier);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
std::tuple<float, float, float, float> BirksHitSelectionTool::GetExtremalBirksSelectionValues(
                                                                                   const LArTrackHitEnergy::Vector &trackHitEnergyVector) const
{
    // Find the min and max coordinate and uncorrected energy.
    float minCoordinate = std::numeric_limits<float>::max();
    float maxCoordinate = std::numeric_limits<float>::min();
    
    float minUncorrectedEnergy = std::numeric_limits<float>::max();
    float maxUncorrectedEnergy = std::numeric_limits<float>::min();
    
    for (const LArTrackHitEnergy &trackHitEnergy : trackHitEnergyVector)
    {
        if (trackHitEnergy.Coordinate() < minCoordinate)
            minCoordinate = trackHitEnergy.Coordinate();
            
         if (trackHitEnergy.Coordinate() > maxCoordinate)
            maxCoordinate = trackHitEnergy.Coordinate();   
            
        if (trackHitEnergy.UncorrectedEnergy() < minUncorrectedEnergy)
            minUncorrectedEnergy = trackHitEnergy.UncorrectedEnergy();
            
         if (trackHitEnergy.UncorrectedEnergy() > maxUncorrectedEnergy)
            maxUncorrectedEnergy = trackHitEnergy.UncorrectedEnergy();  
    }
    
    return std::make_tuple(minCoordinate, maxCoordinate, minUncorrectedEnergy, maxUncorrectedEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

BirksHitSelectionTool::PointDistanceMatrix BirksHitSelectionTool::CreateMatrixOfPointDistances(
                                                                                   const LArTrackHitEnergy::Vector &trackHitEnergyVector) const
{
    // Create a matrix of inter-point distances.
    PointDistanceMatrix matrixOfPointDistances;
    matrixOfPointDistances.reserve(trackHitEnergyVector.size());
    
    for (const LArTrackHitEnergy &trackHitEnergy_i : trackHitEnergyVector)
    {
        PointDistanceVector distanceVector;
        distanceVector.reserve(trackHitEnergyVector.size());
        
        for (const LArTrackHitEnergy &trackHitEnergy_j : trackHitEnergyVector)
        {
            const float xDistance = std::fabs(trackHitEnergy_j.Coordinate() - trackHitEnergy_i.Coordinate());
            const float yDistance = std::fabs(trackHitEnergy_j.UncorrectedEnergy() - trackHitEnergy_i.UncorrectedEnergy());
            
            distanceVector.emplace_back(xDistance, yDistance);
        }
        
        matrixOfPointDistances.push_back(std::move(distanceVector));
    }
    
    return matrixOfPointDistances;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<FloatVector, FloatVector> BirksHitSelectionTool::GetBinAverages(const LArTrackHitEnergy::Vector &trackHitEnergyVector, 
                                                                       const float searchRadiusX,
                                                                       const PointDistanceMatrix &matrixOfPointDistances,
                                                                       const int numBins, const float minCoordinate) const
{
    FloatVector averageBinPointDensity = FloatVector(numBins, 0.f);
    FloatVector averageBinEnergy       = FloatVector(numBins, 0.f);
    
    for (int iBin = 0; iBin < numBins; ++iBin)
    {
        const float xLower = minCoordinate + this->m_birksSelectionBinWidth * static_cast<float>(iBin);
        const float xUpper = minCoordinate + this->m_birksSelectionBinWidth * static_cast<float>(iBin + 1);
        
        std::tie(averageBinPointDensity[iBin], averageBinEnergy[iBin]) = this->GetBinAverage(trackHitEnergyVector, xLower, xUpper,
                                                                                             searchRadiusX, matrixOfPointDistances);
    }
    
    return std::make_tuple(averageBinPointDensity, averageBinEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float> BirksHitSelectionTool::GetBinAverage(const LArTrackHitEnergy::Vector &trackHitEnergyVector, const float xLower,
                                                          const float xUpper, const float searchRadiusX,
                                                          const PointDistanceMatrix &matrixOfPointDistances) const
{
    int   totalLocalPoints  = 0;
    int   numberOfHitsInBin = 0;
    float totalBinEnergy    = 0.f;
    
    for (int i = 0; i < trackHitEnergyVector.size(); ++i)
    {
        const LArTrackHitEnergy &trackHitEnergy_i = trackHitEnergyVector.at(i);
        
        if (!trackHitEnergy_i.ApplyCorrection())
            continue;
        
        if (trackHitEnergy_i.Coordinate() > xLower && trackHitEnergy_i.Coordinate() < xUpper)
        {
            int numberOfLocalHits = 0;
            
            for (int j = 0; j < trackHitEnergyVector.size(); ++j)
            {
                if (!trackHitEnergyVector.at(j).ApplyCorrection())
                    continue;
                
                const auto &distancePair = matrixOfPointDistances[i][j];

                if (distancePair.first < searchRadiusX && distancePair.second < this->m_birksSelectionSearchYWidth)
                    ++numberOfLocalHits;
            }
            
            totalBinEnergy   += trackHitEnergy_i.UncorrectedEnergy();
            totalLocalPoints += numberOfLocalHits;
            ++numberOfHitsInBin;
        }
    }
    
    return std::make_tuple(static_cast<float>(totalLocalPoints) / static_cast<float>(numberOfHitsInBin),
                           totalBinEnergy / static_cast<float>(numberOfHitsInBin));
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void BirksHitSelectionTool::MakeBirksSelectionChanges(LArTrackHitEnergy::Vector &trackHitEnergyVector, const float searchRadiusX,
                                                  const PointDistanceMatrix &matrixOfPointDistances, const int numBins,
                                                  bool &changeMade, const float minCoordinate, const FloatVector &averageBinPointDensity,
                                                  const FloatVector &averageBinEnergy, LArTrackHitEnergy::Vector &trackHitEnergiesChanged,
                                                  LArTrackHitEnergy::Vector &allTrackHitEnergiesChanged) const
{
    for (int iBin = 0; iBin < numBins; ++iBin)
    {    
        const float xLower = minCoordinate + this->m_birksSelectionBinWidth * static_cast<float>(iBin);
        const float xUpper = minCoordinate + this->m_birksSelectionBinWidth * static_cast<float>(iBin + 1);
        
        for (int i = 0; i < trackHitEnergyVector.size(); ++i)
        {
            LArTrackHitEnergy &trackHitEnergy_i = trackHitEnergyVector.at(i);
            
            if (!trackHitEnergy_i.ApplyCorrection())
                continue;
            
            if (trackHitEnergy_i.Coordinate() > xLower && trackHitEnergy_i.Coordinate() < xUpper)
            {
                int numberOfLocalHits = 0;
                
                for (int j = 0; j < trackHitEnergyVector.size(); ++j)
                {
                    if (!trackHitEnergyVector.at(j).ApplyCorrection())
                        continue;
                    
                    const auto &distancePair = matrixOfPointDistances[i][j];

                    if (distancePair.first < searchRadiusX && distancePair.second < this->m_birksSelectionSearchYWidth)
                        ++numberOfLocalHits;
                }
                
                if (static_cast<float>(numberOfLocalHits) < this->m_birksSelectionMinDensityFraction * averageBinPointDensity[iBin] &&
                    trackHitEnergy_i.UncorrectedEnergy() > averageBinEnergy[iBin])
                {
                    trackHitEnergy_i.ApplyCorrection(false);
                    trackHitEnergiesChanged.push_back(trackHitEnergy_i);
                    allTrackHitEnergiesChanged.push_back(trackHitEnergy_i);
                    changeMade = true;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BirksHitSelectionTool::ProduceHitSelectionPlots(const LArTrackHitEnergy::Vector &trackHitEnergyVector,
                                                 const LArTrackHitEnergy::Vector &trackHitEnergiesChanged, const bool birksPointsOnly,
                                                 const bool pause, const bool plotUncorrectedEnergiesOnly, int &uniquePlotIdentifier) const
{
    TNtuple *const pNtuple        = new TNtuple("TrackHitsUncorrected", "TrackHitsUncorrected", "Coordinate:UncorrectedEnergy");
    TNtuple *const pNtupleChanged = new TNtuple("TrackHitsCorrected", "TrackHitsCorrected", "Coordinate:UncorrectedEnergy");
    
    float minCoordinate = std::numeric_limits<float>::max();
    float maxCoordinate = std::numeric_limits<float>::min();
    
    float minUncorrectedEnergy = std::numeric_limits<float>::max();
    float maxUncorrectedEnergy = std::numeric_limits<float>::min();
    
    for (const LArTrackHitEnergy &trackHitEnergy : trackHitEnergyVector)
    {
        if (trackHitEnergy.Coordinate() < minCoordinate)
            minCoordinate = trackHitEnergy.Coordinate();
            
         if (trackHitEnergy.Coordinate() > maxCoordinate)
            maxCoordinate = trackHitEnergy.Coordinate();   
            
        if (trackHitEnergy.UncorrectedEnergy() < minUncorrectedEnergy)
            minUncorrectedEnergy = trackHitEnergy.UncorrectedEnergy();
            
        if (trackHitEnergy.UncorrectedEnergy() > maxUncorrectedEnergy)
            maxUncorrectedEnergy = trackHitEnergy.UncorrectedEnergy();  
            
        if (!plotUncorrectedEnergiesOnly && trackHitEnergy.ApplyCorrection())
        {
            if (trackHitEnergy.CorrectedEnergy() < minUncorrectedEnergy)
                minUncorrectedEnergy = trackHitEnergy.CorrectedEnergy();
            
            if (trackHitEnergy.CorrectedEnergy() > maxUncorrectedEnergy)
                maxUncorrectedEnergy = trackHitEnergy.CorrectedEnergy();  
        }

        if (trackHitEnergy.ApplyCorrection())
        {
            pNtuple->Fill(trackHitEnergy.Coordinate(), plotUncorrectedEnergiesOnly ? 
                                                       trackHitEnergy.UncorrectedEnergy() : trackHitEnergy.CorrectedEnergy());
        }
    }
    
    for (const LArTrackHitEnergy &trackHitEnergy : trackHitEnergiesChanged)
        pNtupleChanged->Fill(trackHitEnergy.Coordinate(), trackHitEnergy.UncorrectedEnergy());
    
    TCanvas *pCanvas = nullptr;
    
    {
        auto plotSettings        = g_defaultPlotSettings2D;
        plotSettings.plotType   = SCATTER;
        plotSettings.pointColor = kRed;
        plotSettings.xMin = minCoordinate;
        plotSettings.xMax = maxCoordinate;
        plotSettings.yMin = 0.f;
        plotSettings.yMax = maxUncorrectedEnergy;
        strcpy(plotSettings.title,  "Selection of hits to Birks-correct");
        strcpy(plotSettings.xTitle, "Projected 3D track coordinate (cm)");
        strcpy(plotSettings.yTitle, "#frac{dQ}{dx} (ADC/cm)");
        
        pCanvas = PlotNtuple2D(pNtuple, "Coordinate", "UncorrectedEnergy", std::to_string(uniquePlotIdentifier++).c_str(),
                               plotSettings);
        
        if (!birksPointsOnly)
        {
            plotSettings.plotType   = SAME;
            plotSettings.newCanvas  = false;
            plotSettings.pointColor = kBlue;
            
            PlotNtuple2D(pNtupleChanged, "Coordinate", "UncorrectedEnergy", std::to_string(uniquePlotIdentifier++).c_str(),
                         plotSettings);
        }
    }
    
    if (pause)
    {
        PandoraMonitoringApi::Pause(this->GetPandora());
        
        if (pCanvas)
            pCanvas->Close();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BirksHitSelectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "BirksSelectionBinWidth", this->m_birksSelectionBinWidth));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "BirksSelectionSearchXWidthFraction", this->m_birksSelectionSearchXWidthFraction));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "BirksSelectionSearchYWidth", this->m_birksSelectionSearchYWidth));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "BirksSelectionMinDensityFraction", this->m_birksSelectionMinDensityFraction));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MakePlots", this->m_makePlots));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
