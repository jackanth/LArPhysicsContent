/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisAlgorithm.cxx
 *
 *  @brief  Implementation of the lar analysis algorithm class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/AnalysisAlgorithm.h"

#include "larphysicscontent/LArHelpers/LArRootHelper.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "TFile.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
AnalysisAlgorithm::AnalysisAlgorithm() :
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minCoordinates(0.f, 0.f, 0.f),
    m_maxCoordinates(0.f, 0.f, 0.f),
    m_birksSelectionMaxdEdX(500.f),
    m_mcParticleListName(),
    m_parametersFile(),
    m_birksFitNtupleName("BirksFit"),
    m_protonEnergyFromRangeNtupleName("EnergyFromRangeProtons"),
    m_pionMuonEnergyFromRangeNtupleName("EnergyFromRangePionsMuons"),
    m_caloHitListName(),
    m_tmvaWeights(),
    m_addMcInformation(true),
    m_birksFitAlpha(0.f),
    m_birksFitBeta(0.f),
    m_birksFitPole(0.f),
    m_protonEnergyFromRangeDataVector(),
    m_pionMuonEnergyFromRangeDataVector(),
    m_pTrackHitEnergyTool(nullptr),
    m_pMcInfoTool(nullptr),
    m_pHitPurityTool(nullptr),
    m_pTmvaReader(nullptr),
    m_tmvaTrackLength(0.f),
    m_tmvaAvgEnergyDeposition(0.f),
    m_uniquePlotIdentifier(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AnalysisAlgorithm::~AnalysisAlgorithm()
{
    if (m_pTmvaReader)
        delete m_pTmvaReader;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject*& pOutputPfo) const
{
    // Work out what kind of PFO we're dealing with.
    bool isNeutrino(false), isPrimaryNeutrinoDaughter(false), isCosmicRay(false);

    if (LArAnalysisParticleHelper::IsNeutrino(pInputPfo))
        isNeutrino = true;

    else if (LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(pInputPfo))
        isPrimaryNeutrinoDaughter = true;

    else if (LArAnalysisParticleHelper::IsCosmicRay(pInputPfo))
        isCosmicRay = true;

    else
        return;

    // Check there is one vertex for this primary PFO and get it.
    if (pInputPfo->GetVertexList().size() != 1)
    {
        std::cout << "AnalysisAlgorithm: could not create LArAnalysisParticle as the number of PFO vertices was " << pInputPfo->GetVertexList().size() << std::endl;
        return;
    }

    const Vertex *const pVertex = pInputPfo->GetVertexList().front();
    bool gotMcInformation = false;
    LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo;

    if (m_addMcInformation)
        gotMcInformation = this->GetMcInformation(pInputPfo, pfoMcInfo, isNeutrino);

    // Common parameter calculations.
    LArAnalysisParticleParameters analysisParticleParameters;

    analysisParticleParameters.m_particleId = pInputPfo->GetParticleId();
    analysisParticleParameters.m_charge     = pInputPfo->GetCharge();
    analysisParticleParameters.m_mass       = pInputPfo->GetMass();
    analysisParticleParameters.m_energy     = pInputPfo->GetEnergy();
    analysisParticleParameters.m_momentum   = pInputPfo->GetMomentum();

    analysisParticleParameters.m_isVertexFiducial = LArAnalysisParticleHelper::IsPointFiducial(pVertex->GetPosition(), m_minCoordinates,
        m_maxCoordinates);

    analysisParticleParameters.m_vertexPosition = pVertex->GetPosition();
    analysisParticleParameters.m_isShower       = LArPfoHelper::IsShower(pInputPfo);

    unsigned numberOfDownstreamParticles(0U);
    this->CountNumberOfDownstreamParticles(pInputPfo, numberOfDownstreamParticles);
    analysisParticleParameters.m_numberOfDownstreamParticles = numberOfDownstreamParticles;

    analysisParticleParameters.m_fiducialHitFraction = LArAnalysisParticleHelper::GetFractionOfFiducialHits(pInputPfo, m_minCoordinates,
        m_maxCoordinates);

    analysisParticleParameters.m_hasMcInfo = gotMcInformation;

    if (gotMcInformation)
    {
        analysisParticleParameters.m_mcEnergy         = pfoMcInfo.m_mcEnergy;
        analysisParticleParameters.m_mcKineticEnergy  = pfoMcInfo.m_mcKineticEnergy;
        analysisParticleParameters.m_mcMass           = pfoMcInfo.m_mcMass;
        analysisParticleParameters.m_mcMomentum       = pfoMcInfo.m_mcMomentum;
        analysisParticleParameters.m_mcVertexPosition = pfoMcInfo.m_mcVertexPosition;

        if (pfoMcInfo.m_mcMomentum.GetMagnitude() > std::numeric_limits<float>::epsilon())
            analysisParticleParameters.m_mcDirectionCosines = pfoMcInfo.m_mcMomentum.GetUnitVector();

        else
        {
            std::cout << "AnalysisAlgorithm: could not get direction from MC momentum as it was too small" << std::endl;
            analysisParticleParameters.m_mcDirectionCosines = CartesianVector(0.f, 0.f, 0.f);
        }

        analysisParticleParameters.m_mcIsVertexFiducial = LArAnalysisParticleHelper::IsPointFiducial(pfoMcInfo.m_mcVertexPosition,
            m_minCoordinates, m_maxCoordinates);

        analysisParticleParameters.m_mcContainmentFraction = pfoMcInfo.m_mcContainmentFraction;
        analysisParticleParameters.m_mcIsShower         = (pfoMcInfo.m_mcType == LArAnalysisParticle::TYPE::SHOWER);
        analysisParticleParameters.m_mcPdgCode          = pfoMcInfo.m_mcPdgCode;
        analysisParticleParameters.m_mcType             = pfoMcInfo.m_mcType;
        analysisParticleParameters.m_mcTypeTree         = pfoMcInfo.m_mcTypeTree;
        analysisParticleParameters.m_pMcMainMCParticle  = pfoMcInfo.m_pMCParticle;
    }

    // Neutrino-specific calculations.
    if (isNeutrino)
    {
        analysisParticleParameters.m_type             = LArAnalysisParticle::TYPE::NEUTRINO;

        float analysisEnergy(0.f), energySourcedFromRange(0.f), energySourcedFromCorrectedTrackCharge(0.f),
            energySourcedFromTrackCharge(0.f), energySourcedFromShowerCharge(0.f);

        CartesianVector analysisMomentum(0.f, 0.f, 0.f);

        unsigned numberOf3dHits(0U), numberOfCollectionPlaneHits(0U);
        LArAnalysisParticle::TypeTree::List primaryTypeTrees;

        for (const ParticleFlowObject *const pPrimary : pInputPfo->GetDaughterPfoList())
        {
            if (const LArAnalysisParticle *const pAnalysisPrimary = dynamic_cast<const LArAnalysisParticle *>(pPrimary))
            {
                analysisEnergy              += pAnalysisPrimary->KineticEnergy();
                analysisMomentum              += pAnalysisPrimary->AnalysisMomentum();
                numberOf3dHits              += pAnalysisPrimary->NumberOf3dHits();
                numberOfCollectionPlaneHits += pAnalysisPrimary->NumberOfCollectionPlaneHits();
                energySourcedFromRange                += pAnalysisPrimary->KineticEnergy() * pAnalysisPrimary->KineticEnergyFromRangeFraction();
                energySourcedFromCorrectedTrackCharge += pAnalysisPrimary->KineticEnergy() * pAnalysisPrimary->KineticEnergyFromCorrectedTrackChargeFraction();
                energySourcedFromTrackCharge          += pAnalysisPrimary->KineticEnergy() * pAnalysisPrimary->KineticEnergyFromUncorrectedTrackChargeFraction();
                energySourcedFromShowerCharge         += pAnalysisPrimary->KineticEnergy() * pAnalysisPrimary->KineticEnergyFromShowerChargeFraction();
                primaryTypeTrees.push_back(pAnalysisPrimary->GetTypeTree());
            }

            else
                std::cout << "AnalysisAlgorithm: primary daughter of neutrino could not be cast an analysis particle" << std::endl;
        }

        analysisParticleParameters.m_kineticEnergy               = analysisEnergy;
        analysisParticleParameters.m_analysisMomentum            = analysisMomentum;
        analysisParticleParameters.m_directionCosines            = analysisMomentum.GetUnitVector();
        analysisParticleParameters.m_numberOf3dHits              = numberOf3dHits;
        analysisParticleParameters.m_numberOfCollectionPlaneHits = numberOfCollectionPlaneHits;
        analysisParticleParameters.m_typeTree = LArAnalysisParticle::TypeTree(LArAnalysisParticle::TYPE::NEUTRINO, primaryTypeTrees);

        if (analysisEnergy > std::numeric_limits<float>::epsilon())
        {
            analysisParticleParameters.m_kineticEnergyFromRangeFraction                  = energySourcedFromRange / analysisEnergy;
            analysisParticleParameters.m_kineticEnergyFromCorrectedTrackChargeFraction   = energySourcedFromCorrectedTrackCharge / analysisEnergy;
            analysisParticleParameters.m_kineticEnergyFromUncorrectedTrackChargeFraction = energySourcedFromTrackCharge / analysisEnergy;
            analysisParticleParameters.m_kineticEnergyFromShowerChargeFraction           = energySourcedFromShowerCharge / analysisEnergy;
        }

        else
        {
            analysisParticleParameters.m_kineticEnergyFromRangeFraction                  = 0.f;
            analysisParticleParameters.m_kineticEnergyFromCorrectedTrackChargeFraction   = 0.f;
            analysisParticleParameters.m_kineticEnergyFromUncorrectedTrackChargeFraction = 0.f;
            analysisParticleParameters.m_kineticEnergyFromShowerChargeFraction           = 0.f;
        }
    }

    // Calculations specific cosmic rays and primary daughters.
    else if (isCosmicRay || isPrimaryNeutrinoDaughter)
    {
        LArAnalysisParticleHelper::FittedTrackInfoMap fittedTrackInfoMap;
        float excessCaloValue(0.f);

        m_pTrackHitEnergyTool->Run(this, pInputPfo, fittedTrackInfoMap, excessCaloValue,
            [&](LArFittedTrackInfo::TrackHitValueVector &trackHitValueVector, float &excessCaloValue) -> bool
            {
                return m_pHitPurityTool->Run(this, trackHitValueVector, excessCaloValue);
            });;

        // Derive a type for all PFOs.
        LArAnalysisParticle::PfoTypeMap particleTypeMap;

        if (isCosmicRay)
            particleTypeMap.emplace(pInputPfo, LArAnalysisParticle::TYPE::COSMIC_RAY);

        else if (isPrimaryNeutrinoDaughter)
            particleTypeMap.emplace(pInputPfo, this->EstimateParticleType(pInputPfo, fittedTrackInfoMap));

        this->RecursivelyAppendParticleTypeMap(pInputPfo, particleTypeMap, fittedTrackInfoMap);

        float particleEnergy(0.f), energySourcedFromRange(0.f), energySourcedFromShowerCharge(0.f),
            energySourcedFromTrackCharge(0.f), energySourcedFromCorrectedTrackCharge(0.f);

        this->EstimateParticleEnergy(pInputPfo, particleTypeMap, fittedTrackInfoMap, particleEnergy,
            energySourcedFromRange, energySourcedFromShowerCharge, energySourcedFromTrackCharge, energySourcedFromCorrectedTrackCharge);

        const LArAnalysisParticle::TypeTree typeTree = this->CreateTypeTree(pInputPfo, particleTypeMap);

        const CartesianVector initialDirection = this->GetDirectionAtVertex(pInputPfo, fittedTrackInfoMap, pVertex, isCosmicRay);

        const int num3DHits = LArAnalysisParticleHelper::GetHitsOfType(pInputPfo, TPC_3D, true).size();
        const int numWHits  = LArAnalysisParticleHelper::GetHitsOfType(pInputPfo, TPC_VIEW_W, true).size();

        const LArAnalysisParticle::TYPE particleType = particleTypeMap.at(pInputPfo);

        analysisParticleParameters.m_kineticEnergy               = particleEnergy;
        analysisParticleParameters.m_directionCosines            = initialDirection;
        analysisParticleParameters.m_analysisMomentum            = initialDirection * particleEnergy;
        analysisParticleParameters.m_numberOf3dHits              = num3DHits;
        analysisParticleParameters.m_numberOfCollectionPlaneHits = numWHits;
        analysisParticleParameters.m_type                        = particleType;
        analysisParticleParameters.m_typeTree                    = typeTree;

        if (particleEnergy > std::numeric_limits<float>::epsilon())
        {
            analysisParticleParameters.m_kineticEnergyFromRangeFraction                  = energySourcedFromRange / particleEnergy;
            analysisParticleParameters.m_kineticEnergyFromCorrectedTrackChargeFraction   = energySourcedFromCorrectedTrackCharge / particleEnergy;
            analysisParticleParameters.m_kineticEnergyFromUncorrectedTrackChargeFraction = energySourcedFromTrackCharge / particleEnergy;
            analysisParticleParameters.m_kineticEnergyFromShowerChargeFraction           = energySourcedFromShowerCharge / particleEnergy;
        }

        else
        {
            analysisParticleParameters.m_kineticEnergyFromRangeFraction                  = 0.f;
            analysisParticleParameters.m_kineticEnergyFromCorrectedTrackChargeFraction   = 0.f;
            analysisParticleParameters.m_kineticEnergyFromUncorrectedTrackChargeFraction = 0.f;
            analysisParticleParameters.m_kineticEnergyFromShowerChargeFraction           = 0.f;
        }
    }

    // Build AnalysisParticle.
    LArAnalysisParticleFactory analysisParticleFactory;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, analysisParticleParameters,
        pOutputPfo, analysisParticleFactory));

    const LArAnalysisParticle *const pLArAnalysisParticle = dynamic_cast<const LArAnalysisParticle*>(pOutputPfo);

    if (!pLArAnalysisParticle)
    {
        std::cout << "AnalysisAlgorithm: failed to cast" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    // Now update vertex and direction
    PandoraContentApi::ParticleFlowObject::Metadata pfoData;
    pfoData.m_momentum = pInputPfo->GetMomentum();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pOutputPfo, pfoData));

    const Vertex *pOutputVertex(NULL);

    // Need an input vertex to provide a track propagation direction
    const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

    PandoraContentApi::Vertex::Parameters vtxParameters;
    vtxParameters.m_position    = pInputVertex->GetPosition();
    vtxParameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    vtxParameters.m_vertexType  = pInputVertex->GetVertexType();

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, vtxParameters, pOutputVertex));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pOutputPfo, pOutputVertex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisAlgorithm::RecursivelyAppendParticleTypeMap(const ParticleFlowObject *const pPfo,
                                                            LArAnalysisParticle::PfoTypeMap &pfoTypeMap,
                                                            const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap) const
{
    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        pfoTypeMap.emplace(pDaughterPfo, this->EstimateParticleType(pDaughterPfo, fittedTrackInfoMap));
        RecursivelyAppendParticleTypeMap(pDaughterPfo, pfoTypeMap, fittedTrackInfoMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticle::TYPE AnalysisAlgorithm::EstimateParticleType(const ParticleFlowObject *const pPfo,
                                                                  const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap) const
{
    if (LArPfoHelper::IsShower(pPfo))
        return LArAnalysisParticle::TYPE::SHOWER;

    const auto trackHitFindIter = fittedTrackInfoMap.find(pPfo);

    if (trackHitFindIter == fittedTrackInfoMap.end())
        return LArAnalysisParticle::TYPE::TRACK;

    const float particleEnergy = this->EstimateTrackEnergyFromCharge(trackHitFindIter->second.HitChargeVector());
    const float particleRange = trackHitFindIter->second.Range();

    if (particleEnergy > 0.f && particleRange > 0.f)
    {
        m_tmvaTrackLength         = particleRange;
        m_tmvaAvgEnergyDeposition = particleEnergy / particleRange;

        const float bdtResponse = m_pTmvaReader->EvaluateMVA("BDT");

        if (bdtResponse > 0.14)
            return LArAnalysisParticle::TYPE::PROTON;

        if (bdtResponse < -0.14)
            return LArAnalysisParticle::TYPE::PION_MUON;
    }

    return LArAnalysisParticle::TYPE::TRACK;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisAlgorithm::EstimateParticleEnergy(const ParticleFlowObject *const pPfo, const LArAnalysisParticle::PfoTypeMap &typeMap,
    const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
    float &particleEnergy, float &energySourcedFromRange, float &energySourcedFromShowerCharge,
    float &energySourcedFromTrackCharge, float &energySourcedFromCorrectedTrackCharge) const
{
    const LArAnalysisParticle::TYPE particleType = typeMap.at(pPfo);

    const auto fitFindIter = fittedTrackInfoMap.find(pPfo);
    const auto fitFound   = (fitFindIter != fittedTrackInfoMap.end());

    switch (particleType)
    {
        // For showers, we assume everything downstream of a showerlike PFO is part of the shower.
        case LArAnalysisParticle::TYPE::SHOWER:
        {
            const float showerEnergy   = this->EstimateShowerEnergy(pPfo);
            particleEnergy                += showerEnergy;
            energySourcedFromShowerCharge += showerEnergy;
            return;
        }

        case LArAnalysisParticle::TYPE::COSMIC_RAY:
        case LArAnalysisParticle::TYPE::PION_MUON:
        case LArAnalysisParticle::TYPE::PROTON:
        {
            const float energyFromCharge = this->EstimateEnergyFromCharge(pPfo);

            if (fitFound)
            {
                const float energyFromRange = this->EstimateTrackEnergyFromRange(fitFindIter->second, particleType);
                particleEnergy         += energyFromRange;
                energySourcedFromRange += energyFromRange;
            }

            else
            {
                particleEnergy               += energyFromCharge; // not including recombination correction (need a fit)
                energySourcedFromTrackCharge += energyFromCharge;
            }

            break;
        }

        case LArAnalysisParticle::TYPE::TRACK:
        {
            if (fitFound)
            {
                const float energyFromCharge = this->EstimateTrackEnergyFromCharge(fitFindIter->second.HitChargeVector()); // incl recombination correction
                particleEnergy                        += energyFromCharge;
                energySourcedFromCorrectedTrackCharge += energyFromCharge;
            }

            else
            {
                const float energyFromCharge = this->EstimateEnergyFromCharge(pPfo); // not including recombination correction (need a fit)
                particleEnergy               += energyFromCharge;
                energySourcedFromTrackCharge += energyFromCharge;
            }

            break;
        }

        default: std::cout << "AnalysisAlgorithm: unknown particle type - could not calculate energy" << std::endl; break;
    }

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        this->EstimateParticleEnergy(pDaughterPfo, typeMap, fittedTrackInfoMap, particleEnergy,
            energySourcedFromRange, energySourcedFromShowerCharge, energySourcedFromTrackCharge, energySourcedFromCorrectedTrackCharge);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisAlgorithm::EstimateShowerEnergy(const ParticleFlowObject *const pPfo) const
{
    float showerEnergy = this->EstimateEnergyFromCharge(pPfo);

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        showerEnergy += this->EstimateShowerEnergy(pDaughterPfo);

    return showerEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisAlgorithm::EstimateEnergyFromCharge(const ParticleFlowObject *const pPfo) const
{
    float energy = 0.f;

    for (const CaloHit *const pCaloHit : LArAnalysisParticleHelper::GetHitsOfType(pPfo, TPC_VIEW_W, false))
        energy += this->CaloHitToScaledEnergy(pCaloHit);

    return energy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisAlgorithm::CaloHitToScaledEnergy(const CaloHit *const pCaloHit) const
{
    const float dQ = pCaloHit->GetInputEnergy();
    return dQ / (m_birksFitAlpha * 1000.f); // GeV
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisAlgorithm::RecombinationCorrectEnergy(const float hitCharge, const float threeDDistance) const
{
    const float dEdx_uncorrected = hitCharge / (m_birksFitAlpha * threeDDistance);               // MeV / cm
    const float dEdx_corrected   = dEdx_uncorrected / (1.f - dEdx_uncorrected / m_birksFitBeta); // MeV / cm
    return dEdx_corrected * threeDDistance / 1000.f; // GeV
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisAlgorithm::EstimateTrackEnergyFromCharge(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector) const
{
    float trackEnergy = 0.f;

    for (const LArTrackHitValue &trackHitEnergy : trackHitEnergyVector)
        trackEnergy += this->RecombinationCorrectEnergy(trackHitEnergy.CaloValue(), trackHitEnergy.ThreeDDistance());

    return trackEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisAlgorithm::EstimateTrackEnergyFromRange(const LArFittedTrackInfo &fittedTrackInfo, const LArAnalysisParticle::TYPE particleType) const
{
    const float trackRange = fittedTrackInfo.Range();
    float trackEnergy = 0.f;

    switch (particleType)
    {
        case LArAnalysisParticle::TYPE::COSMIC_RAY:
        case LArAnalysisParticle::TYPE::PION_MUON:
        case LArAnalysisParticle::TYPE::PROTON:
        {
            bool foundEnergy = false;

            for (const EnergyFromRangeData &energyFromRangeData : (particleType == LArAnalysisParticle::TYPE::PROTON) ?
                                                                  m_protonEnergyFromRangeDataVector :
                                                                  m_pionMuonEnergyFromRangeDataVector)
            {
                if ((trackRange >= energyFromRangeData.m_rangeMin) && (trackRange < energyFromRangeData.m_rangeMax))
                {
                    trackEnergy = energyFromRangeData.m_energy;
                    foundEnergy = true;
                    break;
                }
            }

            if (!foundEnergy)
                trackEnergy = this->EstimateTrackEnergyFromCharge(fittedTrackInfo.HitChargeVector()); // including recombination correction

            break;
        }

        default:
        {
            std::cout << "AnalysisAlgorithm: given a non-known-track particle to calculate energy from range" << std::endl;
            break;
        }
    }

    return trackEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticle::TypeTree AnalysisAlgorithm::CreateTypeTree(const ParticleFlowObject *const pPfo,
                                                                const LArAnalysisParticle::PfoTypeMap &typeMap) const
{
    LArAnalysisParticle::TypeTree::List daughterTypeTrees;

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        daughterTypeTrees.push_back(this->CreateTypeTree(pDaughterPfo, typeMap));

    return {typeMap.at(pPfo), daughterTypeTrees};
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector AnalysisAlgorithm::GetDirectionAtVertex(const ParticleFlowObject *const pPfo,
                                                        const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap, const Vertex *const pVertex,
                                                        const bool isCosmicRay) const
{
    const auto findIter = fittedTrackInfoMap.find(pPfo);

    if (findIter != fittedTrackInfoMap.end())
        return LArAnalysisParticleHelper::GetFittedDirectionAtPosition(findIter->second.Fit(), pVertex->GetPosition(), !isCosmicRay);

    const CaloHitList all3DHits = LArAnalysisParticleHelper::GetHitsOfType(pPfo, TPC_3D, true);

    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::EigenValues  eigenValues{0.f, 0.f, 0.f};
    CartesianVector            centroid{0.f, 0.f, 0.f};
    LArPcaHelper::RunPca(all3DHits, centroid, eigenValues, eigenVectors);

    if (eigenVectors.empty())
    {
        std::cout << "AnalysisAlgorithm: PCA eigenvectors were empty" << std::endl;
        return CartesianVector{0.f, 0.f, 0.f};
    }

    CartesianVector fitDirection = eigenVectors.at(0);
    const CartesianVector vertexToCentroid = centroid - pVertex->GetPosition();

    // We want the fit direction to be mostly aligned with the vertex-to-centroid direction.
    if (vertexToCentroid.GetDotProduct(fitDirection) < 0.f)
        fitDirection *= -1.f;

    return fitDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisAlgorithm::GetMcInformation(const ParticleFlowObject *const pPfo, LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo,
    const bool isNeutrino) const
{
    const MCParticle *pMcMainMCParticle(nullptr);

    if (isNeutrino)
    {
        const MCParticleList *pMCParticleList(nullptr);

        if ((PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS) || !pMCParticleList)
        {
            std::cout << "AnalysisAlgorithm: could not get MC information because no valid MC particle list name was provided" << std::endl;
            return false;
        }

        MCParticleVector trueNeutrinos;
        LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

        if (trueNeutrinos.empty())
            return false;

        if (trueNeutrinos.size() > 1)
        {
            std::sort(trueNeutrinos.begin(), trueNeutrinos.end(), LArMCParticleHelper::SortByMomentum);
            std::cout << "AnalysisAlgorithm: more than one MC neutrino; picking one with largest momentum" << std::endl;
        }

        pMcMainMCParticle = trueNeutrinos.front();
    }

    else
        pMcMainMCParticle = LArAnalysisParticleHelper::GetMainMCParticle(pPfo);

    if (!pMcMainMCParticle)
        return false;

    const CaloHitList *pCaloHitList(nullptr);

    if ((PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS) || !pCaloHitList)
    {
        std::cout << "AnalysisAlgorithm: could not get MC information because no valid CaloHit list name was provided" << std::endl;
        return false;
    }

    m_pMcInfoTool->Run(this, pMcMainMCParticle, pfoMcInfo);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisAlgorithm::CountNumberOfDownstreamParticles(const ParticleFlowObject *const pPfo, unsigned &numberOfParticles) const
{
    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        ++numberOfParticles;
        this->CountNumberOfDownstreamParticles(pDaughterPfo, numberOfParticles);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowMargins", m_fiducialCutLowMargins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighMargins", m_fiducialCutHighMargins));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParametersFile", m_parametersFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BirksFitNtupleName", m_birksFitNtupleName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProtonEnergyFromRangeNtupleName", m_protonEnergyFromRangeNtupleName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PionMuonEnergyFromRangeNtupleName", m_pionMuonEnergyFromRangeNtupleName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BirksSelectionMaxdEdX", m_birksSelectionMaxdEdX));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AddMcInformation", m_addMcInformation));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TmvaWeights", m_tmvaWeights));

    TNtuple *const pBirksNtuple = LArRootHelper::LoadNTupleFromFile(m_parametersFile, m_birksFitNtupleName);

    if (!pBirksNtuple)
    {
        std::cout << "AnalysisAlgorithm: failed to load Birks fit ntuple from file" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    if (pBirksNtuple->GetEntries() != 1)
    {
        std::cout << "AnalysisAlgorithm: Birks fit ntuple did not have exactly one entry: " << pBirksNtuple->GetEntries() << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    float alpha = 0.f, beta = 0.f, dQdxPole = 0.f;

    pBirksNtuple->SetBranchAddress("Alpha",    &alpha);
    pBirksNtuple->SetBranchAddress("Beta",     &beta);
    pBirksNtuple->SetBranchAddress("dQdXPole", &dQdxPole);
    pBirksNtuple->GetEntry(0);

    m_birksFitAlpha = alpha;
    m_birksFitBeta  = beta;
    m_birksFitPole  = dQdxPole;
    m_birksFitPole  = dQdxPole;

    //--------------------------------------------------------------------------------------------------------------------------------------

    TNtuple *const pEnergyFromRangeNtupleProton = LArRootHelper::LoadNTupleFromFile(m_parametersFile, m_protonEnergyFromRangeNtupleName);

    if (!pEnergyFromRangeNtupleProton)
    {
        std::cout << "AnalysisAlgorithm: failed to load proton energy from range ntuple from file" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    const int nEfrProtonEntries = pEnergyFromRangeNtupleProton->GetEntries();

    float efrProtonRangeMin = 0.f, efrProtonRangeMax = 0.f, efrProtonEnergy = 0.f;

    pEnergyFromRangeNtupleProton->SetBranchAddress("RangeMin", &efrProtonRangeMin);
    pEnergyFromRangeNtupleProton->SetBranchAddress("RangeMax", &efrProtonRangeMax);
    pEnergyFromRangeNtupleProton->SetBranchAddress("Energy",   &efrProtonEnergy);

    for (int i = 0; i < nEfrProtonEntries; ++i)
    {
        pEnergyFromRangeNtupleProton->GetEntry(i);
        m_protonEnergyFromRangeDataVector.emplace_back(efrProtonRangeMin, efrProtonRangeMax, efrProtonEnergy);
    }

    //--------------------------------------------------------------------------------------------------------------------------------------

    TNtuple *const pEnergyFromRangeNtuplePionMuon = LArRootHelper::LoadNTupleFromFile(m_parametersFile, m_pionMuonEnergyFromRangeNtupleName);

    if (!pEnergyFromRangeNtuplePionMuon)
    {
        std::cout << "AnalysisAlgorithm: failed to load proton energy from range ntuple from file" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    const int nEfrPionMuonEntries = pEnergyFromRangeNtuplePionMuon->GetEntries();

    float efrPionMuonRangeMin = 0.f, efrPionMuonRangeMax = 0.f, efrPionMuonEnergy = 0.f;

    pEnergyFromRangeNtuplePionMuon->SetBranchAddress("RangeMin", &efrPionMuonRangeMin);
    pEnergyFromRangeNtuplePionMuon->SetBranchAddress("RangeMax", &efrPionMuonRangeMax);
    pEnergyFromRangeNtuplePionMuon->SetBranchAddress("Energy",   &efrPionMuonEnergy);

    for (int i = 0; i < nEfrPionMuonEntries; ++i)
    {
        pEnergyFromRangeNtuplePionMuon->GetEntry(i);
        m_pionMuonEnergyFromRangeDataVector.emplace_back(efrPionMuonRangeMin, efrPionMuonRangeMax, efrPionMuonEnergy);
    }

    TMVA::Tools::Instance();
    m_pTmvaReader = new TMVA::Reader("!Color:Silent");
    m_pTmvaReader->AddVariable("TrackLength", &m_tmvaTrackLength);
    m_pTmvaReader->AddVariable("AvgEnergyDeposition := TrueEnergy/TrackLength", &m_tmvaAvgEnergyDeposition);
    m_pTmvaReader->BookMVA("BDT", m_tmvaWeights.c_str());

    // Use the detector geometry and the margins to get the maximum and minimum fiducial volume coordinates.
    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), m_fiducialCutLowMargins, m_fiducialCutHighMargins, m_minCoordinates,
        m_maxCoordinates);

    AlgorithmTool *pHitPurityAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "HitPurity", pHitPurityAlgorithmTool));

    if (!(m_pHitPurityTool = dynamic_cast<HitPurityTool *>(pHitPurityAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    AlgorithmTool *pMcInfoAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "McInfo", pMcInfoAlgorithmTool));

    if (!(m_pMcInfoTool = dynamic_cast<McInfoTool *>(pMcInfoAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    AlgorithmTool *pTrackHitEnergyAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackHitEnergy", pTrackHitEnergyAlgorithmTool));

    if (!(m_pTrackHitEnergyTool = dynamic_cast<TrackHitEnergyTool *>(pTrackHitEnergyAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}
} // namespace lar_physics_content
