/**
 *  @file   larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.cc
 *
 *  @brief  Implementation of the energy estimator ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.h"
#include "larphysicscontent/LArObjects/LArRootRegistry.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmHeaders.h"

#include "TCanvas.h"
#include "TF1.h"
#include "THStack.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
EnergyEstimatorNtupleTool::EnergyEstimatorNtupleTool() :
    NtupleVariableBaseTool(),
    m_writeEnergiesToNtuple(true),
    m_useParticleId(true),
    m_trainingSetMode(false),
    m_makePlots(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyEstimatorNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteEnergiesToNtuple", m_writeEnergiesToNtuple));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseParticleId", m_useParticleId));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingSetMode", m_trainingSetMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MakePlots", m_makePlots));

    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessEvent(const PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessPrimary(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMcParticle, const MCParticleList *const mcParticleList)
{
    if (m_trainingSetMode)
        return this->ProduceTrainingRecords(pPfo, pfoList, pMcParticle, mcParticleList);

    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessCosmicRay(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMcParticle, const MCParticleList *const mcParticleList)
{
    if (m_trainingSetMode)
        return this->ProduceTrainingRecords(pPfo, pfoList, pMcParticle, mcParticleList);

    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyEstimatorNtupleTool::CaloHitToThreeDDistance(
    const float hitWidth, const ThreeDSlidingFitResult &trackFit, const CartesianVector &threeDPosition, float &threeDDistance) const
{
    CartesianVector fitDirection(0.f, 0.f, 0.f);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArAnalysisHelper::GetFittedDirectionAtThreeDPosition(trackFit, threeDPosition, true, fitDirection));
    const float wirePitch = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W);

    try
    {
        threeDDistance = this->CellToThreeDDistance(hitWidth, wirePitch, fitDirection);
    }

    catch(...)
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::CellToThreeDDistance(const float hitWidth, const float wirePitch, const CartesianVector &fitDirection) const
{
    float polarAngle = 0.f, azimuthalAngle = 0.f;
    std::tie(polarAngle, azimuthalAngle) = this->GetPolarAnglesFromDirection(fitDirection);

    if (polarAngle <= std::numeric_limits<float>::epsilon()) // negligible polar angle, so the wire pitch is the hit separation
        return wirePitch;

    const float cosPhi_sinTheta = std::fabs(std::cos(azimuthalAngle) * std::sin(polarAngle));
    const float sinPhi_sinTheta = std::fabs(std::sin(azimuthalAngle) * std::sin(polarAngle));

    float dx_p = std::numeric_limits<float>::max();
    float dx_w = std::numeric_limits<float>::max();

    if (cosPhi_sinTheta > std::numeric_limits<float>::epsilon())
        dx_p = wirePitch / cosPhi_sinTheta;

    if (sinPhi_sinTheta > std::numeric_limits<float>::epsilon())
        dx_w = hitWidth / sinPhi_sinTheta;

    if ((dx_p < std::numeric_limits<float>::max()) || (dx_w < std::numeric_limits<float>::max()))
        return std::min(dx_p, dx_w);

    // For non-negligible polar angle, this shouldn't happen
    std::cerr << "EnergyEstimatorNtupleTool: Failed to calculate hit's 3D distance" << std::endl;
    throw STATUS_CODE_FAILURE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float> EnergyEstimatorNtupleTool::GetPolarAnglesFromDirection(const CartesianVector &direction) const
{
    const float polarAngle    = std::acos(std::fabs(direction.GetY()));
    const float sinPolarAngle = std::sin(polarAngle);

    if (sinPolarAngle <= std::numeric_limits<float>::epsilon())
        return {0.f, 0.f}; // negligible polar angle means azimuthal angle is undefined

    const float azimuthalAngle = std::asin(std::fabs(direction.GetX() / sinPolarAngle));
    return {polarAngle, azimuthalAngle};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProduceTrainingRecords(
    const ParticleFlowObject *const pPfo, const PfoList &, const MCParticle *const pMcParticle, const MCParticleList *const)
{
    std::vector<LArNtupleRecord> records;

    if (pPfo && pMcParticle)
    {
        LArNtupleRecord::RFloatVector dQdXVector;
        LArNtupleRecord::RFloat       showerAdcCount(0.f);

        for (const ParticleFlowObject *const pDownstreamPfo : this->GetAllDownstreamPfos(pPfo))
        {
            CaloHitList collectionPlaneHits;
            LArPfoHelper::GetCaloHits(pDownstreamPfo, TPC_VIEW_W, collectionPlaneHits);

            if (LArPfoHelper::IsShower(pDownstreamPfo) || !this->GetTrackFit(pDownstreamPfo))
            {
                for (const CaloHit *const pCaloHit : collectionPlaneHits)
                    showerAdcCount += pCaloHit->GetInputEnergy();
            }

            else
            {
                try
                {
                    const LArNtupleHelper::TrackFitSharedPtr &spTrackFit = this->GetTrackFit(pDownstreamPfo);
                    auto [pfodQdX, pfoShowerAdc] = this->SelectNoisyHits(*spTrackFit, collectionPlaneHits, pMcParticle);

                    showerAdcCount += pfoShowerAdc;
                    dQdXVector.insert(dQdXVector.end(), pfodQdX.begin(), pfodQdX.end());
                }

                catch (...)
                {
                    continue;
                }
            }
        }

        records.emplace_back("dQdX", dQdXVector);
        records.emplace_back("ShowerAdcCount", showerAdcCount);
        records.emplace_back("NumVectorEntries", static_cast<LArNtupleRecord::RUInt>(dQdXVector.size()));
    }

    else
    {
        records.emplace_back("dQdX", LArNtupleRecord::RFloatVector());
        records.emplace_back("ShowerAdcCount", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("NumVectorEntries", static_cast<LArNtupleRecord::RUInt>(0U));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<FloatVector, float> EnergyEstimatorNtupleTool::SelectNoisyHits(
    const ThreeDSlidingFitResult &trackFit, CaloHitList caloHitList, const MCParticle *const pMCParticle) const
{
    if (caloHitList.empty())
        return {FloatVector{}, 0.f};

    HitCalorimetryInfoMap hitInfoMap = this->CalculateHitCalorimetryInfo(trackFit, caloHitList);

    if (this->GetNumGoodHits(hitInfoMap) < 5UL)
        return {FloatVector{}, this->SumNoisyCharge(hitInfoMap, true)};

    if (m_makePlots)
    {
        this->PlotHitdQ(caloHitList, hitInfoMap);
        this->PlotHitdX(caloHitList, hitInfoMap);
        this->PlotHitdQdX(caloHitList, hitInfoMap);
        this->PlotTrueMatchedHitdQdX(caloHitList, hitInfoMap, pMCParticle);
    }

    // if (TF1 *pHitGenerationFunction = this->FitHitGenerationFunction(caloHitList, hitInfoMap))
    //    this->GenerateNewHits(trackFit, caloHitList, hitInfoMap, pHitGenerationFunction);

    TF1 *pdQdXFunction = this->FitdQdXFunction(caloHitList, hitInfoMap, 0UL);

    bool        changeMade(true);
    std::size_t iteration(1UL);

    while (changeMade)
    {
        pdQdXFunction = this->FitdQdXFunction(caloHitList, hitInfoMap, iteration++);
        changeMade    = this->IdentifyNoisyHits(caloHitList, hitInfoMap, pdQdXFunction);

        if (this->GetNumGoodHits(hitInfoMap) < 5UL)
            return {FloatVector{}, this->SumNoisyCharge(hitInfoMap, true)};
    }

    this->IdentifyNoisyUnprojectedHits(hitInfoMap);

    const auto [liveLength, deadLength, minCoordinate, maxCoordinate] = this->GetGapParameters(trackFit, caloHitList, hitInfoMap);

    if (liveLength <= std::numeric_limits<float>::epsilon())
    {
        std::cerr << "EnergyEstimatorNtupleTool: The live track length was too small" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    float scaledAdjustedNoisyCharge = this->SumNoisyCharge(hitInfoMap, false) * (liveLength + deadLength) / liveLength;

    FloatVector dQdXVector;
    std::size_t numNonNoisyHits(0UL);

    for (const auto &mapPair : hitInfoMap)
    {
        HitCalorimetryInfo &hitInfo = *mapPair.second;

        if (hitInfo.m_isNoise)
            continue;

        ++numNonNoisyHits;

        if (!hitInfo.m_projectionSuccessful)
        {
            scaledAdjustedNoisyCharge += hitInfo.m_dQ;
            continue;
        }

        dQdXVector.push_back(hitInfo.m_dQ / hitInfo.m_dX);
    }

    if (pdQdXFunction && numNonNoisyHits > 0UL)
    {
        const float lengthPerHit = liveLength / static_cast<float>(numNonNoisyHits);

        for (float coord = minCoordinate; coord <= maxCoordinate; coord += lengthPerHit)
        {
            CartesianVector threeDTrackPosition(0.f, 0.f, 0.f);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, trackFit.GetGlobalFitPosition(coord, threeDTrackPosition));
            const CartesianVector twoDTrackPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), threeDTrackPosition, TPC_VIEW_W);

            if (!LArGeometryHelper::IsInGap(this->GetPandora(), twoDTrackPosition, TPC_VIEW_W, lengthPerHit / 2.f))
                continue;

             dQdXVector.push_back(pdQdXFunction->Eval(coord));
        }
    }

    if (m_makePlots)
        this->MakeColourCodedHitPlot("Identifying noisy hits (final iteration)", caloHitList, hitInfoMap, pdQdXFunction);

    return {dQdXVector, scaledAdjustedNoisyCharge};
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::IdentifyNoisyUnprojectedHits(const HitCalorimetryInfoMap &hitInfoMap) const
{
    float averageNoisedQ(0.f), averageNonNoisedQ(0.f);
    std::size_t numNoiseHits(0UL), numNonNoiseHits(0UL);

    for (const auto &mapPair : hitInfoMap)
    {
        HitCalorimetryInfo &hitInfo = *mapPair.second;

        if (!hitInfo.m_projectionSuccessful)
            continue;

        if (hitInfo.m_isNoise)
        {
            averageNoisedQ += hitInfo.m_dQ;
            ++numNoiseHits;
        }

        else
        {
            averageNonNoisedQ += hitInfo.m_dQ;
            ++numNonNoiseHits;
        }
    }

    if (numNoiseHits == 0UL && numNonNoiseHits == 0UL)
    {
        std::cerr << "EnergyEstimatorNtupleTool: No hits were successfully projected" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (numNoiseHits == 0UL) // nothing is noise
    {
        for (const auto &mapPair : hitInfoMap)
        {
            HitCalorimetryInfo &hitInfo = *mapPair.second;

            if (hitInfo.m_projectionSuccessful)
                continue;

            hitInfo.m_isNoise = false;
        }

        return;
    }

    if (numNonNoiseHits == 0UL) // everything is noise
    {
        for (const auto &mapPair : hitInfoMap)
        {
            HitCalorimetryInfo &hitInfo = *mapPair.second;

            if (hitInfo.m_projectionSuccessful)
                continue;

            hitInfo.m_isNoise = true;
        }

        return;
    }

    // numNonNoiseHits > 0 and numNoiseHits > 0
    averageNoisedQ /= static_cast<float>(numNoiseHits);
    averageNonNoisedQ /= static_cast<float>(numNonNoiseHits);

    const float noiseBoundary = (averageNoisedQ + averageNonNoisedQ) / 2.f;
    const bool noiseIsHigher(averageNoisedQ > averageNonNoisedQ);

    for (const auto &mapPair : hitInfoMap)
    {
        HitCalorimetryInfo &hitInfo = *mapPair.second;

        if (hitInfo.m_projectionSuccessful)
            continue;

        if (hitInfo.m_dQ > noiseBoundary)
            hitInfo.m_isNoise = noiseIsHigher;

        else
            hitInfo.m_isNoise = !noiseIsHigher;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::SumNoisyCharge(const HitCalorimetryInfoMap &hitInfoMap, const bool useAllHits) const
{
    float noisyCharge(0.f);

    for (const auto &mapPair : hitInfoMap)
    {
        HitCalorimetryInfo &hitInfo = *mapPair.second;

        if (!useAllHits && !hitInfo.m_isNoise)
            continue;

        noisyCharge += hitInfo.m_dQ;
    }

    return noisyCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

EnergyEstimatorNtupleTool::HitCalorimetryInfoMap EnergyEstimatorNtupleTool::CalculateHitCalorimetryInfo(
    const ThreeDSlidingFitResult &trackFit, const CaloHitList &caloHitList) const
{
    HitCalorimetryInfoMap hitInfoMap;

    for (const CaloHit *const pCaloHit : caloHitList)
        hitInfoMap.emplace(pCaloHit,
            this->CalculateHitCalorimetryInfo(trackFit, pCaloHit->GetPositionVector(), pCaloHit->GetInputEnergy(), pCaloHit->GetCellSize1()));

    return hitInfoMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

EnergyEstimatorNtupleTool::HitCalorimetryInfoPtr EnergyEstimatorNtupleTool::CalculateHitCalorimetryInfo(
    const ThreeDSlidingFitResult &trackFit, const CartesianVector &twoDPositionVector, const float dQ, const float hitWidth) const
{
    const HitCalorimetryInfoPtr spHitInfo = HitCalorimetryInfoPtr(new HitCalorimetryInfo());

    CartesianVector threeDPosition(0.f, 0.f, 0.f);
    float           projectionError(0.f), dX(0.f);

    if (STATUS_CODE_SUCCESS != LArAnalysisHelper::ProjectTwoDPositionOntoTrackFit(
                                   this->GetPandora(), trackFit, twoDPositionVector, TPC_VIEW_W, true, threeDPosition, projectionError))
    {
        return spHitInfo;
    }

    if (projectionError > 5.f)
        return spHitInfo;

    if (STATUS_CODE_SUCCESS != this->CaloHitToThreeDDistance(hitWidth, trackFit, threeDPosition, dX))
        return spHitInfo;

    if (dX <= std::numeric_limits<float>::epsilon())
        return spHitInfo;

    spHitInfo->m_projectionSuccessful = true;
    spHitInfo->m_threeDPosition       = threeDPosition;
    spHitInfo->m_projectionError      = projectionError;
    spHitInfo->m_coordinate           = trackFit.GetLongitudinalDisplacement(threeDPosition);
    spHitInfo->m_dQ                   = dQ;
    spHitInfo->m_dX                   = dX;

    return spHitInfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TF1 *EnergyEstimatorNtupleTool::FitHitGenerationFunction(const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const
{
    float minCoordinate(0.f), maxCoordinate(0.f);

    if (STATUS_CODE_SUCCESS != this->CalculateCoordinateRange(caloHitList, hitInfoMap, minCoordinate, maxCoordinate))
    {
        std::cerr << "EnergyEstimatorNtupleTool: Failed to calculate coordinate range" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    const float minCentralCoord = minCoordinate + (maxCoordinate - minCoordinate) / 10.f;
    const float maxCentralCoord = maxCoordinate - (maxCoordinate - minCoordinate) / 10.f;
    std::size_t numEntries(0UL);

    LArRootHelper::PlotOptions options;
    options.m_numXBins = static_cast<std::size_t>(std::round(static_cast<float>(caloHitList.size()) / 5.f));
    options.m_title    = "dQ/dX distribution for central track hits";
    options.m_xTitle   = "dQ/dX (integrated ADC/cm)";

    TH1F *pHistogram = this->MakeOneDHitHistogram(
        options, caloHitList, hitInfoMap, [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<float> {
            if (!hitInfo.m_projectionSuccessful)
                return {};

            if ((hitInfo.m_coordinate < minCentralCoord) || (hitInfo.m_coordinate > maxCentralCoord))
                return {};

            ++numEntries;
            return hitInfo.m_dQ / hitInfo.m_dX;
        });

    // Fit a pair of Landaus (new Landau on RHS).
    TF1 *pFunctionL = this->GetTmpRegistry()->CreateWithUniqueName<TF1>("fHitGeneration",
        [&](double *x, double *p) { return p[0] * TMath::Landau(x[0], p[1], p[2]) + p[3] * TMath::Landau(x[0], p[4], p[5]); }, -500.f, 500.f, 6);
    pFunctionL->SetParameter(0, static_cast<double>(numEntries) / 10.);
    pFunctionL->SetParLimits(0, std::min(5., static_cast<double>(numEntries) / 20.), 10000.f);

    pFunctionL->SetParameter(1, 250.f);
    pFunctionL->SetParLimits(1, 150.f, 350.f);

    pFunctionL->SetParameter(2, 50.f);
    pFunctionL->SetParLimits(2, 0.01f, 200.f);

    pFunctionL->SetParameter(3, static_cast<double>(numEntries) / 10.);
    pFunctionL->SetParLimits(3, std::min(5., static_cast<double>(numEntries) / 20.), 10000.f);

    pFunctionL->SetParameter(4, 200.f);
    pFunctionL->SetParLimits(4, 10.f, 1000.f);

    pFunctionL->SetParameter(5, 50.f);
    pFunctionL->SetParLimits(5, 0.01f, 400.f);

    pHistogram->Fit(pFunctionL, "N");

    // Fit a pair of Landau (new Landau on RHS).
    TF1 *pFunctionR = this->GetTmpRegistry()->CreateWithUniqueName<TF1>("fHitGeneration",
        [&](double *x, double *p) { return p[0] * TMath::Landau(x[0], p[1], p[2]) + p[3] * TMath::Landau(x[0], p[4], p[5]); }, -500.f, 500.f, 6);
    pFunctionR->SetParameter(0, static_cast<double>(numEntries) / 10.);
    pFunctionR->SetParLimits(0, std::min(5., static_cast<double>(numEntries) / 20.), 10000.f);

    pFunctionR->SetParameter(1, 250.f);
    pFunctionR->SetParLimits(1, 150.f, 350.f);

    pFunctionR->SetParameter(2, 50.f);
    pFunctionR->SetParLimits(2, 0.01f, 200.f);

    pFunctionR->SetParameter(3, static_cast<double>(numEntries) / 10.);
    pFunctionR->SetParLimits(3, std::min(5., static_cast<double>(numEntries) / 20.), 10000.f);

    pFunctionR->SetParameter(4, 300.f);
    pFunctionR->SetParLimits(4, 250.f, 1000.f);

    pFunctionR->SetParameter(5, 50.f);
    pFunctionR->SetParLimits(5, 0.01f, 400.f);

    pHistogram->Fit(pFunctionR, "N");

    // Fit a single Landau function.
    TF1 *pFunctionSingle = this->GetTmpRegistry()->CreateWithUniqueName<TF1>(
        "fHitGeneration", [&](double *x, double *p) { return p[0] * TMath::Landau(x[0], p[1], p[2]); }, -500.f, 500.f, 6);
    pFunctionSingle->SetParameter(0, static_cast<double>(numEntries) / 10.);
    pFunctionSingle->SetParLimits(0, std::min(5., static_cast<double>(numEntries) / 20.), 10000.f);

    pFunctionSingle->SetParameter(1, 250.f);
    pFunctionSingle->SetParLimits(1, 150.f, 350.f);

    pFunctionSingle->SetParameter(2, 50.f);
    pFunctionSingle->SetParLimits(2, 0.01f, 200.f);

    pHistogram->Fit(pFunctionSingle, "N");

    // Find the best function.
    TF1 *pFunction = pFunctionR;

    if (pFunctionL && pFunctionL->GetChisquare() > 0.f)
    {
        if (!pFunction || pFunction->GetChisquare() == 0.f || (pFunctionL->GetChisquare() < pFunction->GetChisquare()))
            pFunction = pFunctionL;
    }

    if (pFunctionSingle && pFunctionSingle->GetChisquare() > 0.f)
    {
        if (!pFunction || pFunction->GetChisquare() == 0.f || (pFunctionSingle->GetChisquare() < pFunction->GetChisquare()))
            pFunction = pFunctionSingle;
    }

    if (m_makePlots && pFunction)
        this->PlotHistogramAndFunction(pHistogram, pFunction);

    return pFunction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnergyEstimatorNtupleTool::IdentifyNoisyHits(const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdQdXFunction) const
{
    const float stdev            = this->GetFitStandardDeviation(caloHitList, hitInfoMap, pdQdXFunction);
    bool        somethingChanged = false;

    for (const auto &mapPair : hitInfoMap)
    {
        HitCalorimetryInfo &hitInfo = *mapPair.second;

        if (!hitInfo.m_projectionSuccessful || hitInfo.m_isNoise)
            continue;

        const float trueEnergy = hitInfo.m_dQ / hitInfo.m_dX;
        const float fitEnergy  = pdQdXFunction->Eval(hitInfo.m_coordinate);

        if (trueEnergy - fitEnergy > 3.f * stdev)
        {
            hitInfo.m_isNoise = true;
            somethingChanged  = true;
        }
    }

    return somethingChanged;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TF1 *EnergyEstimatorNtupleTool::FitdQdXFunction(const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, const std::size_t iteration) const
{
    TH2F *pHistogram = this->MakeTwoDHitHistogram(LArRootHelper::PlotOptions(), caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful || hitInfo.m_isNoise)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dQ / hitInfo.m_dX);
        });

    TF1 *pFunctionForwards = this->GetTmpRegistry()->CreateWithUniqueName<TF1>("fdQdX", "x < [1] ? min([0]/(1.0 - x / [1]) + [2], 1000.0) : 0.0", 0.f, 1000.f);

    pFunctionForwards->SetParameter(0, 12.5f);
    pFunctionForwards->SetParLimits(0, 0.f, 25.f);

    pFunctionForwards->SetParameter(1, 400.f);
    pFunctionForwards->SetParLimits(1, 20.f, 1000.f);

    pFunctionForwards->SetParameter(2, 200.f);
    pFunctionForwards->SetParLimits(2, 0.f, 1000.f);

    TF1 *pFunctionBackwards = this->GetTmpRegistry()->CreateWithUniqueName<TF1>("fdQdX", "x < - [1] ? min([0]/(1.0 + x / [1]) + [2], 1000.0) : 0.0", 0.f, 1000.f);

    pFunctionBackwards->SetParameter(0, 12.5f);
    pFunctionBackwards->SetParLimits(0, 0.f, 25.f);

    pFunctionBackwards->SetParameter(1, 400.f);
    pFunctionBackwards->SetParLimits(1, 20.f, 1000.f);

    pFunctionBackwards->SetParameter(2, 200.f);
    pFunctionBackwards->SetParLimits(2, 0.f, 1000.f);

    TF1 *pFunction(nullptr);

    try
    {
        const int forwardsFitStatus = static_cast<Int_t>(pHistogram->Fit(pFunctionForwards, "V"));
        const int backwardsFitStatus = static_cast<Int_t>(pHistogram->Fit(pFunctionBackwards, "V"));

        const bool forwardsFitSuccess = (forwardsFitStatus == 0 ) && (pFunctionForwards->GetChisquare() > 0.);
        const bool backwardsFitSuccess = (backwardsFitStatus == 0 ) && (pFunctionBackwards->GetChisquare() > 0.);

        if (!forwardsFitSuccess && !backwardsFitSuccess)
        {
            std::cerr << "EnergyEstimatorNtupleTool: Forwards- and backwards-dQdX fits failed" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        if (!forwardsFitSuccess)
        {
            std::cerr << "EnergyEstimatorNtupleTool: Forwards-dQdX fit failed" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        if (!backwardsFitSuccess)
        {
            std::cerr << "EnergyEstimatorNtupleTool: Backwards-dQdX fit failed" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        // Both fits are viable
        pFunction = (pFunctionForwards->GetChisquare() < pFunctionBackwards->GetChisquare()) ? pFunctionForwards : pFunctionBackwards;
    }

    catch (...)
    {
        pFunction = nullptr;
    }

    if (!pFunction)
    {
        std::cerr << "EnergyEstimatorNtupleTool: Noisy hit fit failed, trying exponential fit" << std::endl;

        try
        {
            const std::string fitName("expo");
            const int fitStatus = static_cast<Int_t>(pHistogram->Fit(fitName.c_str(), "V"));
            TF1 *pExpFunction = pHistogram->GetFunction(fitName.c_str());

            const bool fitSuccess = (fitStatus == 0) && (pExpFunction->GetChisquare() > 0.);

            if (!fitSuccess)
            {
                std::cerr << "EnergyEstimatorNtupleTool: Exponential dQdX fit failed" << std::endl;
                throw STATUS_CODE_FAILURE;
            }

            pFunction = pExpFunction;
        }

        catch (...)
        {
            pFunction = nullptr;
        }
    }

    if (!pFunction)
    {
        std::cerr << "EnergyEstimatorNtupleTool: Exponential noisy hit fit failed, trying linear fit" << std::endl;

        try
        {
            const std::string fitName("pol1");
            const int fitStatus = static_cast<Int_t>(pHistogram->Fit(fitName.c_str(), "V"));
            TF1 *pLinearFunction = pHistogram->GetFunction(fitName.c_str());

            const bool fitSuccess = (fitStatus == 0) && (pLinearFunction->GetChisquare() > 0.);

            if (!fitSuccess)
            {
                std::cerr << "EnergyEstimatorNtupleTool: Linear dQdX fit failed" << std::endl;
                throw STATUS_CODE_FAILURE;
            }

            pFunction = pLinearFunction;
        }

        catch (...)
        {
            pFunction = nullptr;
        }
    }

    if (!pFunction)
    {
        std::cerr << "EnergyEstimatorNtupleTool: All noisy hit fits failed" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (m_makePlots)
        this->MakeColourCodedHitPlot("Identifying noisy hits (iteration " + std::to_string(iteration) + ")", caloHitList, hitInfoMap, pFunction);

    return pFunction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float, float, float> EnergyEstimatorNtupleTool::GetGapParameters(
    const ThreeDSlidingFitResult &trackFit, CaloHitList &caloHitList, HitCalorimetryInfoMap &hitInfoMap) const
{
    // Work out the hit density per unit 3D distance
    std::size_t numDataPoints(0UL);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const auto findIter = hitInfoMap.find(pCaloHit);

        if (findIter == hitInfoMap.end())
        {
            std::cerr << "EnergyEstimatorNtupleTool: Failed to find hit in calorimetry hit info map" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        const HitCalorimetryInfo &hitInfo = *findIter->second;

        if (!hitInfo.m_projectionSuccessful)
            continue;

        ++numDataPoints;
    }

    float minCoordinate(0.f), maxCoordinate(0.f);

    if (STATUS_CODE_SUCCESS != this->CalculateCoordinateRange(caloHitList, hitInfoMap, minCoordinate, maxCoordinate))
    {
        std::cerr << "EnergyEstimatorNtupleTool: Failed to generate new hits because could not get range" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    const float coordinateSpan = maxCoordinate - minCoordinate;
    const float dLength        = coordinateSpan / static_cast<float>(numDataPoints);
    float       liveLength(0.f);

    for (float coord = minCoordinate, maxBound = maxCoordinate + dLength / 2.f; coord < maxBound; coord += dLength)
    {
        CartesianVector threeDTrackPosition(0.f, 0.f, 0.f);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, trackFit.GetGlobalFitPosition(coord, threeDTrackPosition));
        const CartesianVector twoDTrackPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), threeDTrackPosition, TPC_VIEW_W);

        if (!LArGeometryHelper::IsInGap(this->GetPandora(), twoDTrackPosition, TPC_VIEW_W, dLength / 4.f))
            liveLength += dLength;
    }

    if (liveLength <= std::numeric_limits<float>::epsilon())
    {
        std::cerr << "EnergyEstimatorNtupleTool: Failed to generate new hits because the live length was too small" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    const float deadLength = coordinateSpan > liveLength ? coordinateSpan - liveLength : 0.f;
    return {liveLength, deadLength, minCoordinate, maxCoordinate};
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::GenerateNewHits(
    const ThreeDSlidingFitResult &trackFit, CaloHitList &caloHitList, HitCalorimetryInfoMap &hitInfoMap, TF1 *pHitGenerationFunction) const
{
    // Work out the hit density per unit 3D distance
    std::size_t numDataPoints(0UL);
    float       averageWidth(0.f);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const auto findIter = hitInfoMap.find(pCaloHit);

        if (findIter == hitInfoMap.end())
        {
            std::cerr << "EnergyEstimatorNtupleTool: Failed to find hit in calorimetry hit info map" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        const HitCalorimetryInfo &hitInfo = *findIter->second;

        if (!hitInfo.m_projectionSuccessful)
            continue;

        averageWidth += pCaloHit->GetCellSize1();
        ++numDataPoints;
    }

    float minCoordinate(0.f), maxCoordinate(0.f);

    if ((numDataPoints == 0UL) || (STATUS_CODE_SUCCESS != this->CalculateCoordinateRange(caloHitList, hitInfoMap, minCoordinate, maxCoordinate)))
    {
        std::cerr << "EnergyEstimatorNtupleTool: Failed to generate new hits because could not get range" << std::endl;
        return;
    }

    averageWidth /= static_cast<float>(numDataPoints);
    const float dLength = (maxCoordinate - minCoordinate) / static_cast<float>(numDataPoints);
    float       liveLength(0.f);

    for (float coord = minCoordinate, maxBound = maxCoordinate + dLength / 2.f; coord < maxBound; coord += dLength)
    {
        CartesianVector threeDTrackPosition(0.f, 0.f, 0.f);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, trackFit.GetGlobalFitPosition(coord, threeDTrackPosition));
        const CartesianVector twoDTrackPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), threeDTrackPosition, TPC_VIEW_W);

        if (!LArGeometryHelper::IsInGap(this->GetPandora(), twoDTrackPosition, TPC_VIEW_W, dLength / 4.f))
            liveLength += dLength;
    }

    if (liveLength <= std::numeric_limits<float>::epsilon())
    {
        std::cerr << "EnergyEstimatorNtupleTool: Failed to generate new hits because the live length was too small" << std::endl;
        return;
    }

    const float           lengthPerHit = liveLength / static_cast<float>(numDataPoints);
    CaloHitList           newHits;
    HitCalorimetryInfoMap newHitInfoMap;

    for (float coord = minCoordinate, maxBound = maxCoordinate + lengthPerHit / 2.f; coord < maxBound; coord += lengthPerHit)
    {
        CartesianVector threeDTrackPosition(0.f, 0.f, 0.f);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, trackFit.GetGlobalFitPosition(coord, threeDTrackPosition));
        const CartesianVector twoDTrackPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), threeDTrackPosition, TPC_VIEW_W);

        if (!LArGeometryHelper::IsInGap(this->GetPandora(), twoDTrackPosition, TPC_VIEW_W, lengthPerHit / 4.f))
            continue;

        // We're in a gap
        if (const HitCalorimetryInfoPtr spHitInfo = this->CalculateHitCalorimetryInfo(trackFit, twoDTrackPosition, 0.f, averageWidth)) // dummy dQ for now
        {
            float dQ = -1.f;

            while (dQ < 0.f)
                dQ = spHitInfo->m_dX * pHitGenerationFunction->GetRandom();

            spHitInfo->m_dQ     = dQ;
            spHitInfo->m_isFake = true;

            const CaloHit *const pModelHit = caloHitList.front();

            PandoraApi::CaloHit::Parameters parameters;
            parameters.m_positionVector = twoDTrackPosition;
            parameters.m_cellSize1      = averageWidth;
            parameters.m_inputEnergy    = spHitInfo->m_dQ;
            parameters.m_hitType        = TPC_VIEW_W;

            parameters.m_cellSize0               = pModelHit->GetCellSize0();
            parameters.m_expectedDirection       = pModelHit->GetExpectedDirection();
            parameters.m_cellNormalVector        = pModelHit->GetCellNormalVector();
            parameters.m_cellGeometry            = pModelHit->GetCellGeometry();
            parameters.m_cellThickness           = pModelHit->GetCellThickness();
            parameters.m_nCellRadiationLengths   = pModelHit->GetNCellRadiationLengths();
            parameters.m_nCellInteractionLengths = pModelHit->GetNCellInteractionLengths();
            parameters.m_time                    = pModelHit->GetTime();
            parameters.m_mipEquivalentEnergy     = pModelHit->GetMipEquivalentEnergy();
            parameters.m_electromagneticEnergy   = pModelHit->GetElectromagneticEnergy();
            parameters.m_hadronicEnergy          = pModelHit->GetHadronicEnergy();
            parameters.m_isDigital               = pModelHit->IsDigital();
            parameters.m_hitRegion               = pModelHit->GetHitRegion();
            parameters.m_layer                   = pModelHit->GetLayer();
            parameters.m_isInOuterSamplingLayer  = pModelHit->IsInOuterSamplingLayer();
            parameters.m_pParentAddress          = nullptr;

            const CaloHit *pNewHit(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*this->GetAlgorithm(), parameters, pNewHit));

            newHits.push_back(pNewHit);
            newHitInfoMap.emplace(pNewHit, spHitInfo);
        }
    }

    caloHitList.insert(caloHitList.end(), std::make_move_iterator(newHits.begin()), std::make_move_iterator(newHits.end()));
    hitInfoMap.insert(std::make_move_iterator(newHitInfoMap.begin()), std::make_move_iterator(newHitInfoMap.end()));

    if (m_makePlots)
        this->MakeColourCodedHitPlot("Generating new hits", caloHitList, hitInfoMap, nullptr);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyEstimatorNtupleTool::CalculateCoordinateRange(
    const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, float &minCoordinate, float &maxCoordinate) const
{
    // Work out the hit density per unit 3D distance
    float coordMin(std::numeric_limits<float>::max());
    float coordMax(std::numeric_limits<float>::min());
    bool  minCoordinateSet(false);
    bool  maxCoordinateSet(false);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const auto findIter = hitInfoMap.find(pCaloHit);

        if (findIter == hitInfoMap.end())
        {
            std::cerr << "EnergyEstimatorNtupleTool: Failed to find hit in calorimetry hit info map" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        const HitCalorimetryInfo &hitInfo = *findIter->second;

        if (!hitInfo.m_projectionSuccessful)
            continue;

        if (hitInfo.m_coordinate > coordMax)
        {
            coordMax         = hitInfo.m_coordinate;
            maxCoordinateSet = true;
        }

        if (hitInfo.m_coordinate < coordMin)
        {
            coordMin         = hitInfo.m_coordinate;
            minCoordinateSet = true;
        }
    }

    if (!minCoordinateSet || !maxCoordinateSet)
        return STATUS_CODE_FAILURE;

    minCoordinate = coordMin;
    maxCoordinate = coordMax;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TH1F *EnergyEstimatorNtupleTool::MakeOneDHitHistogram(const LArRootHelper::PlotOptions &options, const CaloHitList &caloHitList,
    const HitCalorimetryInfoMap &hitInfoMap, const HitDataGetter<float> &hitDataGetter) const
{
    LArRootHelper::FloatVector values;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const auto findIter = hitInfoMap.find(pCaloHit);

        if (findIter == hitInfoMap.end())
        {
            std::cerr << "EnergyEstimatorNtupleTool: Failed to find CaloHit in calorimetry map" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        if (const auto oData = hitDataGetter(pCaloHit, *findIter->second))
            values.push_back(*oData);
    }

    return LArRootHelper::CreateOneDHistogram(*this->GetTmpRegistry(), values, options);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TH2F *EnergyEstimatorNtupleTool::MakeTwoDHitHistogram(const LArRootHelper::PlotOptions &options, const CaloHitList &caloHitList,
    const HitCalorimetryInfoMap &hitInfoMap, const HitDataGetter<std::pair<float, float>> &hitDataGetter) const
{
    LArRootHelper::FloatVector xValues, yValues;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const auto findIter = hitInfoMap.find(pCaloHit);

        if (findIter == hitInfoMap.end())
        {
            std::cerr << "EnergyEstimatorNtupleTool: Failed to find CaloHit in calorimetry map" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        if (const auto oData = hitDataGetter(pCaloHit, *findIter->second))
        {
            xValues.push_back(oData->first);
            yValues.push_back(oData->second);
        }
    }

    if (xValues.size() != yValues.size())
    {
        std::cerr << "EnergyEstimatorNtupleTool: ..." << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (xValues.empty())
        return nullptr;

    return LArRootHelper::CreateTwoDHistogram(*this->GetTmpRegistry(), xValues, yValues, options);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::GetFitStandardDeviation(const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdQdXFunction) const
{
    float stdev(0.f);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const auto findIter = hitInfoMap.find(pCaloHit);

        if (findIter == hitInfoMap.end())
            continue;

        const HitCalorimetryInfo &hitInfo = *findIter->second;

        if (!hitInfo.m_projectionSuccessful || hitInfo.m_isNoise)
            continue;

        const float fitEnergy  = pdQdXFunction->Eval(hitInfo.m_coordinate);
        const float trueEnergy = hitInfo.m_dQ / hitInfo.m_dX;

        stdev += (fitEnergy - trueEnergy) * (fitEnergy - trueEnergy);
    }

    if (caloHitList.size() > 1UL)
        stdev /= static_cast<float>(caloHitList.size());

    return std::sqrt(stdev);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::MakeColourCodedHitPlot(
    const std::string &title, const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdQdXFunction) const
{
    LArRootHelper::PlotOptions options;
    options.m_markerColour = kBlack;

    TH2F *pHits = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful || hitInfo.m_isNoise || hitInfo.m_isFake)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dQ / hitInfo.m_dX);
        });

    options.m_markerColour = kBlue;

    TH2F *pNoisyHits = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful || !hitInfo.m_isNoise || hitInfo.m_isFake)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dQ / hitInfo.m_dX);
        });

    options.m_markerColour = kRed;

    TH2F *pFakeHits = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful || hitInfo.m_isNoise || !hitInfo.m_isFake)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dQ / hitInfo.m_dX);
        });

    options.m_markerColour = kMagenta;

    TH2F *pFakeNoisyHits = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful || !hitInfo.m_isNoise || !hitInfo.m_isFake)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dQ / hitInfo.m_dX);
        });

    const std::string titleString(title + ";3D track coordinate (cm);dQ/dx (integrated ADC/cm)");

    THStack *pStack = this->GetTmpRegistry()->CreateWithUniqueName<THStack>("hStack", titleString.c_str());

    if (pHits && pHits->GetEntries() > 0)
        pStack->Add(pHits);

    if (pNoisyHits && pNoisyHits->GetEntries() > 0)
        pStack->Add(pNoisyHits);

    if (pFakeHits && pFakeHits->GetEntries() > 0)
        pStack->Add(pFakeHits);

    if (pFakeNoisyHits && pFakeNoisyHits->GetEntries() > 0)
        pStack->Add(pFakeNoisyHits);

    TCanvas *pCanvas = this->GetPlotsRegistry()->CreateWithUniqueName<TCanvas>("ColorCodedHitPlot", "ColorCodedHitPlot", 10, 10, 900, 600);

    this->GetPlotsRegistry()->DoAsRegistry([&]() {
        pStack->Draw("nostack");

        if (pdQdXFunction)
            pdQdXFunction->Draw("same");

        pCanvas->Write();
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::size_t EnergyEstimatorNtupleTool::GetNumGoodHits(const HitCalorimetryInfoMap &hitInfoMap) const
{
    std::size_t numGoodHits(0UL);

    for (const auto &mapEntry : hitInfoMap)
    {
        if (mapEntry.second->m_projectionSuccessful && !mapEntry.second->m_isNoise)
            ++numGoodHits;
    }

    return numGoodHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::PlotHitdQ(const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const
{
    LArRootHelper::PlotOptions options;
    options.m_title        = "Hit dQ along track fit";
    options.m_markerColour = kBlack;
    options.m_xTitle       = "Track coordinate (cm)";
    options.m_yTitle       = "dQ (integrated ADC)";

    TH2F *pHistogram = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dQ);
        });

    TCanvas *pCanvas = this->GetPlotsRegistry()->CreateWithUniqueName<TCanvas>("HitdQPlot", "HitdQPlot", 10, 10, 900, 600);

    this->GetPlotsRegistry()->DoAsRegistry([&]() {
        pHistogram->Draw();
        pCanvas->Write();
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::PlotHitdX(const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const
{
    LArRootHelper::PlotOptions options;
    options.m_title        = "Hit 3D dX along track fit";
    options.m_markerColour = kBlack;
    options.m_xTitle       = "Track coordinate (cm)";
    options.m_yTitle       = "dX (cm)";

    TH2F *pHistogram = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dX);
        });

    TCanvas *pCanvas = this->GetPlotsRegistry()->CreateWithUniqueName<TCanvas>("HitdXPlot", "HitdXPlot", 10, 10, 900, 600);

    this->GetPlotsRegistry()->DoAsRegistry([&]() {
        pHistogram->Draw();
        pCanvas->Write();
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::PlotHitdQdX(const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const
{
    LArRootHelper::PlotOptions options;
    options.m_title        = "Hit dQ/dX along track fit";
    options.m_markerColour = kBlack;
    options.m_xTitle       = "Track coordinate (cm)";
    options.m_yTitle       = "dQ/dX (integrated ADC/cm)";

    TH2F *pHistogram = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful)
                return {};

            return std::make_pair(hitInfo.m_coordinate, hitInfo.m_dQ / hitInfo.m_dX);
        });

    TCanvas *pCanvas = this->GetPlotsRegistry()->CreateWithUniqueName<TCanvas>("HitdQdXPlot", "HitdQdXPlot", 10, 10, 900, 600);

    this->GetPlotsRegistry()->DoAsRegistry([&]() {
        pHistogram->Draw();
        pCanvas->Write();
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::PlotTrueMatchedHitdQdX(
    const CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, const MCParticle *const pMCParticle) const
{
    LArRootHelper::PlotOptions options;
    options.m_title        = "True-match-weighted dQ/dX along track fit";
    options.m_markerColour = kBlack;
    options.m_xTitle       = "Track coordinate (cm)";
    options.m_yTitle       = "dQ/dX (integrated ADC/cm)";

    TH2F *pHistogram = this->MakeTwoDHitHistogram(options, caloHitList, hitInfoMap,
        [&](const CaloHit *const pCaloHit, const HitCalorimetryInfo &hitInfo) -> std::optional<std::pair<float, float>> {
            if (!hitInfo.m_projectionSuccessful)
                return {};

            const auto  weightMap = pCaloHit->GetMCParticleWeightMap();
            const auto  findIter  = weightMap.find(pMCParticle);
            const float weight    = (findIter == weightMap.end()) ? 0.f : findIter->second;

            return std::make_pair(hitInfo.m_coordinate, weight * hitInfo.m_dQ / hitInfo.m_dX);
        });

    TCanvas *pCanvas = this->GetPlotsRegistry()->CreateWithUniqueName<TCanvas>("TrueMatchedHitdQdXPlot", "TrueMatchedHitdQdXPlot", 10, 10, 900, 600);

    this->GetPlotsRegistry()->DoAsRegistry([&]() {
        pHistogram->Draw();
        pCanvas->Write();
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyEstimatorNtupleTool::PlotHistogramAndFunction(TH1 *pHistogram, TF1 *pFunction) const
{
    TCanvas *pCanvas = this->GetPlotsRegistry()->CreateWithUniqueName<TCanvas>("HistAndFunc", "HistAndFunc", 10, 10, 900, 600);

    this->GetPlotsRegistry()->DoAsRegistry([&]() {
        pHistogram->Draw();
        pFunction->Draw("same");
        pCanvas->Write();
    });
}

} // namespace lar_physics_content
