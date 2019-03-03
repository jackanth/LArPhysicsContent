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
#include "TLatex.h"
#include "TTreeReader.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

EnergyEstimatorNtupleTool::HitCalorimetryInfo::HitCalorimetryInfo() :
    m_projectionSuccessful(false),
    m_threeDPosition(0.f, 0.f, 0.f),
    m_projectionError(0.f),
    m_coordinate(0.f),
    m_dQ(0.f),
    m_dX(0.f),
    m_isNoise(false),
    m_isFake(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnergyEstimatorNtupleTool::EnergyEstimatorNtupleTool() :
    NtupleVariableBaseTool{},
    m_braggGradientTrainingMode{false},
    m_makePlots{false},
    m_modboxRho{0.f},
    m_modboxA{0.f},
    m_modboxB{0.f},
    m_modboxEpsilon{0.f},
    m_modboxWion{0.f},
    m_modboxC{0.f},
    m_modboxFactor{0.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyEstimatorNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BraggGradientTrainingMode", m_braggGradientTrainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MakePlots", m_makePlots));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModBoxRho", m_modboxRho));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModBoxA", m_modboxA));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModBoxB", m_modboxB));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModBoxEpsilon", m_modboxEpsilon));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModBoxWion", m_modboxWion));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModBoxC", m_modboxC));

    if (m_modboxC < std::numeric_limits<float>::epsilon())
    {
        std::cerr << "EnergyEstimatorNtupleTool: ModBox C parameter was too small" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (m_modboxB < std::numeric_limits<float>::epsilon())
    {
        std::cerr << "EnergyEstimatorNtupleTool: ModBox B parameter was too small" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    m_modboxFactor = m_modboxRho * m_modboxEpsilon / m_modboxB;

    if (m_modboxFactor < std::numeric_limits<float>::epsilon())
    {
        std::cerr << "EnergyEstimatorNtupleTool: ModBox (rho * epsilon / B) was too small" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessEvent(const PfoList &, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const pNeutrinoPfo, const PfoList &, const std::shared_ptr<LArInteractionValidationInfo> &)
{
    if (m_braggGradientTrainingMode)
        return {};

    std::vector<LArNtupleRecord> records;

    if (pNeutrinoPfo)
    {
        LArNtupleRecord::RFloat energyEstimator(0.f);
        LArNtupleRecord::RUInt  numTrackHits(0U), numTrackHitsLost(0U);

        for (const ParticleFlowObject *const pPrimary : pNeutrinoPfo->GetDaughterPfoList())
        {
            energyEstimator += this->GetPrimaryRecord<LArNtupleRecord::RFloat>("RecoKineticEnergy", pPrimary);
            numTrackHits += this->GetPrimaryRecord<LArNtupleRecord::RUInt>("NumTrackHits", pPrimary);
            numTrackHitsLost += this->GetPrimaryRecord<LArNtupleRecord::RUInt>("NumTrackHitsLost", pPrimary);
        }

        records.emplace_back("RecoKineticEnergy", energyEstimator);
        records.emplace_back("NumTrackHits", numTrackHits);
        records.emplace_back("NumTrackHitsLost", numTrackHitsLost);
    }

    else
    {
        records.emplace_back("RecoKineticEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("NumTrackHits", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("NumTrackHitsLost", static_cast<LArNtupleRecord::RUInt>(0U));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    std::vector<LArNtupleRecord> records;
    const MCParticle *const      pMcParticle = spMcTarget ? spMcTarget->GetMCParticle() : nullptr;

    (void)pfoList;
    if (m_braggGradientTrainingMode)
        return records; // return this->ProduceBraggGradientTrainingRecords(pPfo, pfoList, pMcParticle);

    std::vector<LArNtupleRecord> energyEstimatorRecords = this->GetEnergyEstimatorRecords(pPfo, pMcParticle);
    records.insert(records.end(), std::make_move_iterator(energyEstimatorRecords.begin()), std::make_move_iterator(energyEstimatorRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    std::vector<LArNtupleRecord> records;
    const MCParticle *const      pMcParticle = spMcTarget ? spMcTarget->GetMCParticle() : nullptr;

    if (m_braggGradientTrainingMode)
        return this->ProduceBraggGradientTrainingRecords(pPfo, pfoList, pMcParticle);

    std::vector<LArNtupleRecord> energyEstimatorRecords = this->GetEnergyEstimatorRecords(pPfo, pMcParticle);
    records.insert(records.end(), std::make_move_iterator(energyEstimatorRecords.begin()), std::make_move_iterator(energyEstimatorRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::GetEnergyEstimatorRecords(
    const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) const
{
    std::vector<LArNtupleRecord> records;

    if (pPfo)
    {
        const auto [dQdXVector, dXVector, showerCharge, numHitsLostToErrors] = this->GetHitCalorimetryInfo(pPfo, pMCParticle);
        float energyEstimator                                                = this->EstimateShowerEnergy(showerCharge);

        for (std::size_t i = 0UL, numHits = dQdXVector.size(); i < numHits; ++i)
            energyEstimator += this->EstimateTrackHitEnergy(dQdXVector.at(i), dXVector.at(i));

        records.emplace_back("RecoKineticEnergy", static_cast<LArNtupleRecord::RFloat>(energyEstimator));
        records.emplace_back("NumTrackHits", static_cast<LArNtupleRecord::RUInt>(dQdXVector.size()));
        records.emplace_back("NumTrackHitsLost", static_cast<LArNtupleRecord::RUInt>(numHitsLostToErrors));
    }

    else
    {
        records.emplace_back("RecoKineticEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("NumTrackHits", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("NumTrackHitsLost", static_cast<LArNtupleRecord::RUInt>(0U));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<LArNtupleRecord::RFloatVector, LArNtupleRecord::RFloatVector, LArNtupleRecord::RFloat, LArNtupleRecord::RUInt>
EnergyEstimatorNtupleTool::GetHitCalorimetryInfo(const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) const
{
    LArNtupleRecord::RFloatVector dQdXVector, dXVector;
    LArNtupleRecord::RFloat       showerCharge(0.f);
    LArNtupleRecord::RUInt        numHitsLostToErrors(0UL);

    for (const ParticleFlowObject *const pDownstreamPfo : this->GetAllDownstreamPfos(pPfo))
    {
        CaloHitList collectionPlaneHits;
        LArPfoHelper::GetCaloHits(pDownstreamPfo, TPC_VIEW_W, collectionPlaneHits);

        if (LArPfoHelper::IsShower(pDownstreamPfo) || !this->GetTrackFit(pDownstreamPfo))
        {
            for (const CaloHit *const pCaloHit : collectionPlaneHits)
                showerCharge += pCaloHit->GetInputEnergy();

            continue;
        }

        // It's tracklike and we have a good track fit
        CaloHitList threeDCaloHits;
        LArPfoHelper::GetCaloHits(pDownstreamPfo, TPC_3D, threeDCaloHits);
        CaloHitMap caloHitMap;

        for (const CaloHit *const pThreeDHit : threeDCaloHits)
        {
            if (const CaloHit *const pTwoDHit = reinterpret_cast<const CaloHit *const>(pThreeDHit->GetParentAddress()))
            {
                if (pTwoDHit->GetHitType() != TPC_VIEW_W)
                    continue;

                caloHitMap.emplace(pTwoDHit, pThreeDHit);
            }
        }

        (void)pMCParticle;

        // const bool                                isBackwards = pPfo->GetMomentum().GetDotProduct(pMCParticle->GetMomentum()) < 0.f;
        // const LArNtupleHelper::TrackFitSharedPtr &spTrackFit  = this->GetTrackFit(pDownstreamPfo);
        // const auto hitChargeVector = this->GetdEdxDistribution(*spTrackFit, collectionPlaneHits, caloHitMap, isBackwards, pMCParticle);

        // const auto detector                        = bf::DetectorHelper::GetMicroBooNEDetector();
        // const auto quickPidAlgorithm               = bf::QuickPidAlgorithm{detector};
        // const auto [braggGradient, braggIntercept] = quickPidAlgorithm.CalculateBraggGradient(hitChargeVector);

        // std::cout << "Grad = " << braggGradient << ", intercept = " << braggIntercept << std::endl;
    }

    if (dQdXVector.size() != dXVector.size())
    {
        std::cerr << "EnergyEstimatorNtupleTool: The size of the dQ/dx vector (" << dQdXVector.size()
                  << ") did not match the size of the dx vector (" << dXVector.size() << ")" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    return {dQdXVector, dXVector, showerCharge, numHitsLostToErrors};
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

    catch (...)
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
    throw StatusCodeException(STATUS_CODE_FAILURE);
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

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProduceBraggGradientTrainingRecords(
    const ParticleFlowObject *const pPfo, const PfoList &, const MCParticle *const pMcParticle)
{
    std::vector<LArNtupleRecord> records;

    float braggGradient1                 = 0.f;
    float braggIntercept1                = 0.f;
    float braggGradient2                 = 0.f;
    float braggIntercept2                = 0.f;
    float braggAvgDetectorThickness      = 0.f;
    float pida                           = 0.f;
    float medianFilteredEnergyLossRate   = 0.f;
    float medianUnfilteredEnergyLossRate = 0.f;

    std::vector<float> maxResidualRanges = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f, 18.f, 19.f, 20.f, 22.f, 24.f, 26.f, 30.f, 35.f, 40.f, 50.f, 70.f, 100.f};

    for (const float maxResidualRange : maxResidualRanges)
    {
        const std::string suffix = "_max" + std::to_string(maxResidualRange);

        if (this->GetBraggGradientParameters(pPfo, pMcParticle, braggGradient1, braggIntercept1, braggGradient2, braggIntercept2,
                braggAvgDetectorThickness, pida, medianFilteredEnergyLossRate, medianUnfilteredEnergyLossRate, maxResidualRange))
        {
            records.emplace_back("HasBraggParameters" + suffix, static_cast<LArNtupleRecord::RBool>(true));
            records.emplace_back("BraggGradient1" + suffix, static_cast<LArNtupleRecord::RFloat>(braggGradient1));
            records.emplace_back("BraggIntercept1" + suffix, static_cast<LArNtupleRecord::RFloat>(braggIntercept1));
            records.emplace_back("BraggGradient2" + suffix, static_cast<LArNtupleRecord::RFloat>(braggGradient2));
            records.emplace_back("BraggIntercept2" + suffix, static_cast<LArNtupleRecord::RFloat>(braggIntercept2));
            records.emplace_back("BraggAverageDetectorThickness" + suffix, static_cast<LArNtupleRecord::RFloat>(braggAvgDetectorThickness));
            records.emplace_back("Pida" + suffix, static_cast<LArNtupleRecord::RFloat>(pida));
            records.emplace_back("MedianUnfilteredEnergyLossRate" + suffix, static_cast<LArNtupleRecord::RFloat>(medianUnfilteredEnergyLossRate));
            records.emplace_back("MedianFilteredEnergyLossRate" + suffix, static_cast<LArNtupleRecord::RFloat>(medianFilteredEnergyLossRate));
        }

        else
        {
            records.emplace_back("HasBraggParameters" + suffix, static_cast<LArNtupleRecord::RBool>(false));
            records.emplace_back("BraggGradient1" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("BraggIntercept1" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("BraggGradient2" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("BraggIntercept2" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("BraggAverageDetectorThickness" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("Pida" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("MedianUnfilteredEnergyLossRate" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("MedianFilteredEnergyLossRate" + suffix, static_cast<LArNtupleRecord::RFloat>(0.f));
        }
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnergyEstimatorNtupleTool::GetBraggGradientParameters(const ParticleFlowObject *const pPfo, const MCParticle *const pMcParticle,
    float &firstOrderGradient, float &firstOrderIntercept, float &secondOrderGradient, float &secondOrderIntercept,
    float &averageDetectorThickness, float &pida, float &medianFilteredEnergyLossRate, float &medianUnfilteredEnergyLossRate, const float maxResidualRange) const
{
    // Check PFO eligibility.
    if (!pPfo || !pMcParticle)
        return false;

    if (this->GetAllDownstreamPfos(pPfo).size() != 1UL)
        return false;

    if (LArPfoHelper::IsShower(pPfo))
        return false;

    const auto spTrackFit = this->GetTrackFit(pPfo);

    if (!spTrackFit)
        return false;

    // It's tracklike and we have a good track fit.
    CaloHitList collectionPlaneHits;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, collectionPlaneHits);

    CaloHitList threeDCaloHits;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, threeDCaloHits);
    CaloHitMap caloHitMap;

    for (const CaloHit *const pThreeDHit : threeDCaloHits)
    {
        if (const CaloHit *const pTwoDHit = reinterpret_cast<const CaloHit *const>(pThreeDHit->GetParentAddress()))
        {
            if (pTwoDHit->GetHitType() != TPC_VIEW_W)
                continue;

            caloHitMap.emplace(pTwoDHit, pThreeDHit);
        }
    }

    // Filter the hit charges to get Bragg peak.
    const bool isReconstructedBackwards = pPfo->GetMomentum().GetDotProduct(pMcParticle->GetMomentum()) < 0.f;
    const bool isRecoBackwardsGoing     = pPfo->GetMomentum().GetZ() < 0.f;
    const bool isBackwards              = isReconstructedBackwards != isRecoBackwardsGoing;

    const auto hitChargeVector = this->GetdEdxDistribution(*spTrackFit, collectionPlaneHits, caloHitMap, isBackwards, pMcParticle);

    double                     maxCoordinate = 0.;
    std::vector<bf::HitCharge> filteredHitCharges;
    std::vector<bf::HitCharge> unfilteredHitCharges;

    for (const auto &hitCharge : hitChargeVector)
    {
        if (hitCharge.Coordinate() > maxCoordinate)
            maxCoordinate = hitCharge.Coordinate();
    }

    for (const auto &hitCharge : hitChargeVector)
    {
        if (maxCoordinate - hitCharge.Coordinate() < maxResidualRange && maxCoordinate - hitCharge.Coordinate() > std::numeric_limits<float>::epsilon())
            filteredHitCharges.push_back(hitCharge);

        else
            unfilteredHitCharges.push_back(hitCharge);
    }

    // Get the Bragg parameters and optionally plot them.
    double firstOrderGradientDouble  = 0.;
    double firstOrderInterceptDouble = 0.;
    if (!bf::QuickPidAlgorithm::CalculateBraggGradient(filteredHitCharges, firstOrderGradientDouble, firstOrderInterceptDouble))
        return false;

    firstOrderGradient  = static_cast<float>(firstOrderGradientDouble);
    firstOrderIntercept = static_cast<float>(firstOrderInterceptDouble);

    if (m_makePlots)
    {
        const float trueKineticEnergy = this->GetPrimaryRecord<LArNtupleRecord::RFloat>("mc_KineticEnergy", pPfo);
        std::cerr << "Plotted particle has incident energy " << trueKineticEnergy << std::endl;
        bf::PlotHelper::SetGlobalPlotStyle();

        TCanvas *pCanvas = this->GetTmpRegistry()->CreateWithUniqueName<TCanvas>("BraggGradientPlot", "BraggGradientPlot", 10, 10, 900, 600);
        auto     braggGraph = bf::PlotHelper::GetBraggGradientGraph(filteredHitCharges);

        unsigned int colour = 7UL;
        switch (std::abs(pMcParticle->GetParticleId()))
        {
            case 13:
                colour = 0UL;
                break;
            case 211:
                colour = 1UL;
                break;
            case 321:
                colour = 2UL;
                break;
            case 2212:
                colour = 3UL;
                break;
            default:
                break;
        }

        braggGraph.GetYaxis()->SetTitle("\\mathrm{d}E/\\mathrm{d}x \\text{ (MeV/cm)}");
        braggGraph.GetXaxis()->SetTitle("1 / \\sqrt{R - x}  \\text{ (cm}^{-0.5}\\text{)}");
        braggGraph.SetMarkerColor(bf::PlotHelper::GetSchemeColourLight(colour));
        braggGraph.SetMarkerStyle(7UL);
        braggGraph.SetMinimum(0.);
        braggGraph.Draw("AP");

        const double xMax = braggGraph.GetXaxis()->GetXmax();
        braggGraph.GetXaxis()->SetRangeUser(0., xMax);

        TF1 *pFunction = this->GetTmpRegistry()->CreateWithUniqueName<TF1>("BraggLine", "[0] + [1] * x", 0., xMax);
        pFunction->SetLineColor(bf::PlotHelper::GetSchemeColour(colour));
        pFunction->SetParameter(0, firstOrderIntercept);
        pFunction->SetParameter(1, firstOrderGradient);
        pFunction->Draw("same");

        auto particleSymbol = std::string{};

        switch (pMcParticle->GetParticleId())
        {
            case 13:
                particleSymbol = "\\mu";
                break;
            case 211:
                particleSymbol = "\\pi^+";
                break;
            case -211:
                particleSymbol = "\\pi^-";
                break;
            case 321:
                particleSymbol = "K^+\\";
                break;
            case -321:
                particleSymbol = "K^-\\";
                break;
            case 2212:
                particleSymbol = "p\\";
                break;
            default:
                particleSymbol = "\\text{PDG " + std::to_string(pMcParticle->GetParticleId()) + "}";
                break;
        }

        std::stringstream labelStream;
        labelStream << std::fixed << std::setprecision(2);
        labelStream << "#splitline{gradient = " << firstOrderGradient << "}{intercept = " << firstOrderIntercept << "}";
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.DrawLatexNDC(0.7, 0.7, labelStream.str().c_str());

        TLatex latexLabel;
        latexLabel.SetTextSize(0.06);
        latexLabel.DrawLatexNDC(0.77, 0.8, particleSymbol.c_str());

        if (gROOT->IsBatch())
            pCanvas->SaveAs((std::string{pCanvas->GetName()} + ".eps").c_str());

        else
        {
            bf::PlotHelper::Pause();
            pCanvas->Close();
        }
    }

    // Calculate the average detector thickness.
    float detectorThickness = 0.f;

    for (const auto &hitCharge : filteredHitCharges)
        detectorThickness += static_cast<float>(hitCharge.Extent());

    detectorThickness /= static_cast<float>(filteredHitCharges.size());
    averageDetectorThickness = detectorThickness;

    // Calculate PIDA.
    float pidaValue = 0.f;

    for (const auto &hitCharge : filteredHitCharges)
    {
        const float residualRange = static_cast<float>(maxCoordinate - hitCharge.Coordinate());
        pidaValue += static_cast<float>(hitCharge.EnergyLossRate()) * std::pow(residualRange, 0.42f);
    }

    pidaValue /= static_cast<float>(filteredHitCharges.size());
    pida = pidaValue;

    // Calculate hasBraggPeak.
    std::vector<double> filteredAvgEnergyLossRates;
    std::vector<double> unFilteredAvgEnergyLossRates;

    for (const auto &hitCharge : filteredHitCharges)
    {
        filteredAvgEnergyLossRates.push_back(static_cast<float>(hitCharge.EnergyLossRate()));
    }

    for (const auto &hitCharge : unfilteredHitCharges)
    {
        unFilteredAvgEnergyLossRates.push_back(static_cast<float>(hitCharge.EnergyLossRate()));
    }

    medianFilteredEnergyLossRate   = static_cast<float>(bf::QuickPidHelper::CalculateMedian(filteredAvgEnergyLossRates));
    medianUnfilteredEnergyLossRate = static_cast<float>(bf::QuickPidHelper::CalculateMedian(unFilteredAvgEnergyLossRates));

    // Calculate second order approx.
    double secondOrderGradientDouble  = 0.;
    double secondOrderInterceptDouble = 0.;

    const auto detector    = bf::DetectorHelper::GetMicroBooNEDetector();
    const auto quickPidAlg = bf::QuickPidAlgorithm{detector};
    if (!quickPidAlg.CalculateSecondOrderBraggGradient(filteredHitCharges, secondOrderGradientDouble, secondOrderInterceptDouble))
        return false;

    secondOrderGradient  = static_cast<float>(secondOrderGradientDouble);
    secondOrderIntercept = static_cast<float>(secondOrderInterceptDouble);

    if (m_makePlots)
    {
        // const float trueKineticEnergy = this->GetPrimaryRecord<LArNtupleRecord::RFloat>("mc_KineticEnergy", pPfo);
        // std::cerr << "Plotted particle has incident energy " << trueKineticEnergy << std::endl;
        bf::PlotHelper::SetGlobalPlotStyle();

        TCanvas *pCanvas = this->GetTmpRegistry()->CreateWithUniqueName<TCanvas>("BraggGradientPlot", "BraggGradientPlot", 10, 10, 900, 600);
        auto     braggGraph = quickPidAlg.GetSecondOrderBraggGradientGraph(filteredHitCharges);

        unsigned int colour = 7UL;
        switch (std::abs(pMcParticle->GetParticleId()))
        {
            case 13:
                colour = 0UL;
                break;
            case 211:
                colour = 1UL;
                break;
            case 321:
                colour = 2UL;
                break;
            case 2212:
                colour = 3UL;
                break;
            default:
                break;
        }

        braggGraph.GetYaxis()->SetTitle("Q\\");
        braggGraph.GetXaxis()->SetTitle("R - x  \\text{ (cm)}");
        braggGraph.SetMarkerColor(bf::PlotHelper::GetSchemeColourLight(colour));
        braggGraph.SetMarkerStyle(7UL);
        braggGraph.Draw("AP");

        const double xMax = braggGraph.GetXaxis()->GetXmax();
        braggGraph.GetXaxis()->SetRangeUser(0., xMax);

        TF1 *pFunction = this->GetTmpRegistry()->CreateWithUniqueName<TF1>("BraggLine", "[0] + [1] * x", 0., xMax);
        pFunction->SetLineColor(bf::PlotHelper::GetSchemeColour(colour));
        pFunction->SetParameter(0, secondOrderIntercept);
        pFunction->SetParameter(1, secondOrderGradient);
        pFunction->Draw("same");

        auto particleSymbol = std::string{};

        switch (pMcParticle->GetParticleId())
        {
            case 13:
                particleSymbol = "\\mu";
                break;
            case 211:
                particleSymbol = "\\pi^+";
                break;
            case -211:
                particleSymbol = "\\pi^-";
                break;
            case 321:
                particleSymbol = "K^+\\";
                break;
            case -321:
                particleSymbol = "K^-\\";
                break;
            case 2212:
                particleSymbol = "p\\";
                break;
            default:
                particleSymbol = "\\text{PDG " + std::to_string(pMcParticle->GetParticleId()) + "}";
                break;
        }

        std::stringstream labelStream;
        labelStream << std::scientific << std::setprecision(2);
        labelStream << "#splitline{gradient = " << secondOrderGradient << "}{intercept = " << secondOrderIntercept << "}";
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.DrawLatexNDC(0.7, 0.7, labelStream.str().c_str());

        TLatex latexLabel;
        latexLabel.SetTextSize(0.06);
        latexLabel.DrawLatexNDC(0.77, 0.8, particleSymbol.c_str());

        if (gROOT->IsBatch())
            pCanvas->SaveAs((std::string{pCanvas->GetName()} + ".eps").c_str());

        else
        {
            bf::PlotHelper::Pause();
            pCanvas->Close();
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<bf::HitCharge> EnergyEstimatorNtupleTool::GetdEdxDistribution(const ThreeDSlidingFitResult &trackFit, CaloHitList caloHitList,
    const CaloHitMap &caloHitMap, const bool isBackwards, const MCParticle *const pMcParticle) const
{
    if (caloHitList.empty())
        return {};

    std::vector<HitCalorimetryInfoPtr> hitInfoVector = this->CalculateHitCalorimetryInfo(trackFit, caloHitList, caloHitMap, isBackwards);
    std::vector<bf::HitCharge>         hitChargeVector;
    std::vector<double>                coordinateVector, dQdxVector;

    for (const auto spHitInfo : hitInfoVector)
    {
        if (!spHitInfo->m_projectionSuccessful)
            continue;

        if (spHitInfo->m_dX < std::numeric_limits<float>::epsilon())
            continue;

        const float dQdxUncorrected = spHitInfo->m_dQ / spHitInfo->m_dX;
        const float dQdxCorrected   = this->CorrectChargeDeposition(dQdxUncorrected, spHitInfo->m_threeDPosition);

        const float dEdx = this->ApplyModBoxCorrection(dQdxCorrected);
        hitChargeVector.emplace_back(spHitInfo->m_coordinate, dEdx, spHitInfo->m_dX);

        coordinateVector.emplace_back(spHitInfo->m_coordinate);
        dQdxVector.emplace_back(dQdxCorrected);
    }

    if (m_makePlots && pMcParticle)
    {
        bf::PlotHelper::SetGlobalPlotStyle();

        double maxCoordinate = 0.;

        for (const double coordinate : coordinateVector)
        {
            if (coordinate > maxCoordinate)
                maxCoordinate = coordinate;
        }

        for (double &coordinate : coordinateVector)
            coordinate = maxCoordinate - coordinate;

        TGraph   chargeDistributionGraph{static_cast<Int_t>(coordinateVector.size()), coordinateVector.data(), dQdxVector.data()};
        TCanvas *pCanvas =
            this->GetTmpRegistry()->CreateWithUniqueName<TCanvas>("ChargeDistributionPlot", "ChargeDistributionPlot", 10, 10, 900, 600);

        unsigned int colour = 7UL;
        switch (std::abs(pMcParticle->GetParticleId()))
        {
            case 13:
                colour = 0UL;
                break;
            case 211:
                colour = 1UL;
                break;
            case 321:
                colour = 2UL;
                break;
            case 2212:
                colour = 3UL;
                break;
            default:
                break;
        }

        chargeDistributionGraph.GetYaxis()->SetTitle("\\mathrm{d}Q/\\mathrm{d}x \\text{ (ADC/cm)}");
        chargeDistributionGraph.GetXaxis()->SetTitle("\\text{Residual range (cm})");
        chargeDistributionGraph.SetMarkerColor(bf::PlotHelper::GetSchemeColour(colour));
        chargeDistributionGraph.SetMarkerStyle(7UL);
        chargeDistributionGraph.SetMinimum(0.);
        chargeDistributionGraph.Draw("AP");

        const double xMax = chargeDistributionGraph.GetXaxis()->GetXmax();
        chargeDistributionGraph.GetXaxis()->SetRangeUser(0., xMax);

        auto particleSymbol = std::string{};

        switch (pMcParticle->GetParticleId())
        {
            case 13:
                particleSymbol = "\\mu";
                break;
            case 211:
                particleSymbol = "\\pi^+";
                break;
            case -211:
                particleSymbol = "\\pi^-";
                break;
            case 321:
                particleSymbol = "K^+\\";
                break;
            case -321:
                particleSymbol = "K^-\\";
                break;
            case 2212:
                particleSymbol = "p\\";
                break;
            default:
                particleSymbol = "\\text{PDG " + std::to_string(pMcParticle->GetParticleId()) + "}";
                break;
        }

        TLatex latexLabel;
        latexLabel.SetTextSize(0.06);
        latexLabel.DrawLatexNDC(0.77, 0.8, particleSymbol.c_str());

        if (gROOT->IsBatch())
            pCanvas->SaveAs((std::string{pCanvas->GetName()} + ".eps").c_str());

        else
        {
            bf::PlotHelper::Pause();
            pCanvas->Close();
        }
    }

    return hitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<EnergyEstimatorNtupleTool::HitCalorimetryInfoPtr> EnergyEstimatorNtupleTool::CalculateHitCalorimetryInfo(
    const ThreeDSlidingFitResult &trackFit, const CaloHitList &caloHitList, const CaloHitMap &caloHitMap, const bool isBackwards) const
{
    std::vector<HitCalorimetryInfoPtr> hitInfoVector;

    if (caloHitList.empty())
        return hitInfoVector;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        hitInfoVector.push_back(this->CalculateHitCalorimetryInfo(
            trackFit, pCaloHit->GetPositionVector(), pCaloHit->GetInputEnergy(), pCaloHit->GetCellSize1(), caloHitMap, pCaloHit));
    }

    std::sort(hitInfoVector.begin(), hitInfoVector.end(), [&](const auto &spLhs, const auto &spRhs) {
        return isBackwards ? spLhs->m_coordinate > spRhs->m_coordinate : spLhs->m_coordinate < spRhs->m_coordinate;
    });

    float           coord    = isBackwards ? -hitInfoVector.front()->m_coordinate : hitInfoVector.front()->m_coordinate;
    CartesianVector position = hitInfoVector.front()->m_threeDPosition;
    float           range    = 0.f;

    for (const auto &spHitInfo : hitInfoVector)
    {
        if (!spHitInfo->m_projectionSuccessful)
            continue;

        bool        failed   = false;
        const float newCoord = isBackwards ? -spHitInfo->m_coordinate : spHitInfo->m_coordinate;

        while (coord < newCoord)
        {
            auto trackPosition = CartesianVector{0.f, 0.f, 0.f};
            if (trackFit.GetGlobalFitPosition(isBackwards ? -coord : coord, trackPosition) != STATUS_CODE_SUCCESS)
            {
                failed = true;
                break;
            }

            range += (trackPosition - position).GetMagnitude();
            position = trackPosition;
            coord += 0.0001f;
        }

        if (failed)
        {
            spHitInfo->m_projectionSuccessful = false;
            continue;
        }

        spHitInfo->m_coordinate = range;
    }

    return hitInfoVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

EnergyEstimatorNtupleTool::HitCalorimetryInfoPtr EnergyEstimatorNtupleTool::CalculateHitCalorimetryInfo(const ThreeDSlidingFitResult &trackFit,
    const CartesianVector &twoDPositionVector, const float dQ, const float hitWidth, const CaloHitMap &caloHitMap, const CaloHit *const pCaloHit) const
{
    const HitCalorimetryInfoPtr spHitInfo = HitCalorimetryInfoPtr(new HitCalorimetryInfo());

    CartesianVector threeDHitPosition(0.f, 0.f, 0.f);
    float           threeDHitProjectionError(std::numeric_limits<float>::max());

    const auto findIter = caloHitMap.find(pCaloHit);

    if (findIter != caloHitMap.end())
    {
        threeDHitPosition = findIter->second->GetPositionVector();
        threeDHitProjectionError =
            (LArGeometryHelper::ProjectPosition(this->GetPandora(), threeDHitPosition, TPC_VIEW_W) - twoDPositionVector).GetMagnitude();
    }

    CartesianVector inferredThreeDPosition(0.f, 0.f, 0.f);
    float           inferredThreeDProjectionError(std::numeric_limits<float>::max()), dX(0.f);

    if (STATUS_CODE_SUCCESS != LArAnalysisHelper::ProjectTwoDPositionOntoTrackFit(this->GetPandora(), trackFit, twoDPositionVector,
                                   TPC_VIEW_W, true, inferredThreeDPosition, inferredThreeDProjectionError))
    {
        inferredThreeDProjectionError = std::numeric_limits<float>::max();
    }

    if (inferredThreeDProjectionError > 5.f && threeDHitProjectionError > 5.f)
        return spHitInfo;

    const CartesianVector threeDPosition = (threeDHitProjectionError < inferredThreeDProjectionError) ? threeDHitPosition : inferredThreeDPosition;
    const float           projectionError = std::min(threeDHitProjectionError, inferredThreeDProjectionError);

    if (STATUS_CODE_SUCCESS != this->CaloHitToThreeDDistance(hitWidth, trackFit, threeDPosition, dX))
        return spHitInfo;

    if (dX <= std::numeric_limits<float>::epsilon())
        return spHitInfo;

    spHitInfo->m_projectionSuccessful = true;
    spHitInfo->m_threeDPosition       = threeDPosition;
    spHitInfo->m_projectionError      = projectionError;
    spHitInfo->m_coordinate           = trackFit.GetLongitudinalDisplacement(threeDPosition); // temporary coordinate
    spHitInfo->m_dQ                   = dQ;
    spHitInfo->m_dX                   = dX;

    return spHitInfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::CorrectChargeDeposition(const float dQdxUncorrected, const CartesianVector &threeDPosition) const
{
    const float positionX = threeDPosition.GetX();
    const float positionY = threeDPosition.GetY();
    const float positionZ = threeDPosition.GetZ();

    return dQdxUncorrected * this->GetDriftCoordinateCorrection(positionX) * this->GetYZCoordinateCorrection(positionY, positionZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::GetDriftCoordinateCorrection(const float xPosition) const
{
    const float parameter0 = 1.35970e-1f;
    const float parameter1 = 5.18310e-3f;
    const float parameter2 = -3.45421f;
    const float parameter3 = 1.37225f;
    const float parameter4 = -4.23881e-3f;
    const float parameter5 = -6.83842e-1f;

    return (parameter0 + parameter1 * (xPosition - parameter2) + parameter3 * std::sin(parameter4 * xPosition - parameter5));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::GetYZCoordinateCorrection(const float yPosition, const float zPosition) const
{
    const bool zCut = (zPosition > 0.f) && (zPosition < 400.f);
    const bool yzCut1 = yPosition > -120.f + zPosition * 220.f / 400.f;
    const bool yzCut2 = yPosition < zPosition * 120.f / 250.f;

    if (zCut && yzCut1 && yzCut2)
        return 1.3f;

    return 1.f;
}

} // namespace lar_physics_content
