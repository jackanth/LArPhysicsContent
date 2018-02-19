/**
 *  @file LArPhysicsContent/src/AnalysisDataAlgorithm.cxx
 *
 *  @brief Implementation of the LEE analysis data algorithm class.
 *
 *  $Log: $
 */

#include "AnalysisDataAlgorithm.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/PdgTable.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "LArAnalysisParticleHelper.h"
#include "TFile.h"

namespace lar_physics_content
{
AnalysisDataAlgorithm::AnalysisDataAlgorithm() :
    m_fiducialCutLowXMargin(10.f),
    m_fiducialCutHighXMargin(10.f),
    m_fiducialCutLowYMargin(20.f),
    m_fiducialCutHighYMargin(20.f),
    m_fiducialCutLowZMargin(10.f),
    m_fiducialCutHighZMargin(10.f),
    m_trackSlidingFitWindow(25U),
    m_produceBirksFitData(true),
    m_produceEnergyFromRangeData(true),
    m_producePidData(true),
    m_pfoListName(),
    m_rootDataFileName(),
    m_caloHitListName(),
    m_pidDataProtonsTTreeName("PidDataProtons"),
    m_pidDataMuonsPionsTTreeName("PidDataMuonsPions"),
    m_pRootDataFile(nullptr),
    m_pBirksFitDataTree(nullptr),
    m_pEnergyFromRangeProtonDataNtuple(nullptr),
    m_pEnergyFromRangePionMuonDataNtuple(nullptr),
    m_pBirksHitSelectionTool(nullptr),
    m_uniquePlotIdentifier(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AnalysisDataAlgorithm::~AnalysisDataAlgorithm()
{
    if (this->m_pBirksFitDataTree)
        this->m_pBirksFitDataTree->Write();

    if (this->m_pEnergyFromRangeProtonDataNtuple)
        this->m_pEnergyFromRangeProtonDataNtuple->Write();

    if (this->m_pEnergyFromRangePionMuonDataNtuple)
        this->m_pEnergyFromRangePionMuonDataNtuple->Write();

    if (this->m_pRootDataFile)
    {
        if (this->m_pRootDataFile->IsOpen())
            this->m_pRootDataFile->Close();

        delete this->m_pRootDataFile;
    }

    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_pidDataProtonsTTreeName.c_str(), m_rootDataFileName.c_str(), "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_pidDataMuonsPionsTTreeName.c_str(), m_rootDataFileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisDataAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);

    if ((PandoraContentApi::GetList(*this, this->m_pfoListName, pPfoList) != STATUS_CODE_SUCCESS) || !pPfoList)
    {
        std::cout << "AnalysisAlgorithm: failed to get the PFO list" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    const CaloHitList *pCaloHitList(nullptr);

    if ((PandoraContentApi::GetList(*this, this->m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS) || !pCaloHitList)
    {
        std::cout << "AnalysisDataAlgorithm: could not get MC information because no valid CaloHit list name was provided" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    PfoList pfoList;

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(pPfo) || LArAnalysisParticleHelper::IsCosmicRay(pPfo))
            pfoList.push_back(pPfo);
    }

    // Move the below to ReadSettings

    // Use the detector geometry and the margins to get the maximum and minimum fiducial volume coordinates.
    CartesianVector minCoordinates(0.f, 0.f, 0.f), maxCoordinates(0.f, 0.f, 0.f);

    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), m_fiducialCutLowXMargin, m_fiducialCutHighXMargin,
        m_fiducialCutLowYMargin, this->m_fiducialCutHighYMargin, m_fiducialCutLowZMargin, m_fiducialCutHighZMargin,
        minCoordinates, maxCoordinates);

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        if (LArAnalysisParticleHelper::IsCosmicRay(pPfo))
            std::cout << "Considering cosmic ray" << std::endl;

        else if (LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(pPfo))
            std::cout << "Considering primary particle" << std::endl;

        // Try to get the corresponding MC particle. Skip this daughter if no MC particle can be found.
        MCParticleMap mcParticleMap;
        const MCParticle *pMCParticle(nullptr);

        try
        {
            pMCParticle = LArAnalysisParticleHelper::GetMainMCParticle(pPfo);
        }

        catch (...)
        {
            std::cout << "AnalysisDataAlgorithm: failed to get main MC particle" << std::endl;
            continue;
        }

        if (!pMCParticle)
        {
            std::cout << "AnalysisDataAlgorithm: failed to get main MC particle" << std::endl;
            continue;
        }

        mcParticleMap.emplace(pPfo, pMCParticle);

        float mcEnergy(0.f), mcKineticEnergy(0.f), mcMass(0.f), mcContainmentFraction(0.f), mcHitPurity(0.f), mcHitCompleteness(0.f),
            mcCollectionPlaneHitPurity(0.f), mcCollectionPlaneHitCompleteness(0.f);

        LArAnalysisParticle::TypeTree mcTypeTree;
        LArAnalysisParticle::TYPE mcType(LArAnalysisParticle::TYPE::UNKNOWN);
        CartesianVector mcVertexPosition(0.f, 0.f, 0.f), mcMomentum(0.f, 0.f, 0.f);
        int mcPdgCode(0);

        LArAnalysisParticleHelper::GetMcInformation(pMCParticle, mcEnergy, mcKineticEnergy, mcMass, mcTypeTree, mcType, mcVertexPosition, mcMomentum,
        mcPdgCode, mcContainmentFraction, minCoordinates, maxCoordinates);

        LArAnalysisParticleHelper::CalculateHitPurityAndCompleteness(pPfo, pMCParticle, pCaloHitList, false, mcHitPurity, mcHitCompleteness,
            mcCollectionPlaneHitPurity, mcCollectionPlaneHitCompleteness);

        // Quality cuts.
        std::cout << "Containment fraction is    " << 100.f * mcContainmentFraction << "%" << std::endl;
        std::cout << "Hit purity is              " << 100.f * mcHitPurity << "%" << std::endl;
        std::cout << "Hit completeness is        " << 100.f * mcHitCompleteness << "%" << std::endl;
        std::cout << "W view hit purity is       " << 100.f * mcCollectionPlaneHitPurity << "%" << std::endl;
        std::cout << "W view hit completeness is " << 100.f * mcCollectionPlaneHitCompleteness << "%" << std::endl;

        if (mcContainmentFraction < 0.95f || mcHitPurity < 0.95f || mcHitCompleteness < 0.95f || mcCollectionPlaneHitPurity < 0.95f ||
            mcCollectionPlaneHitCompleteness < 0.95f)
        {
            std::cout << "Skipping..." << std::endl;
            continue;
        }

        // For each tracklike PFO, try to perform a 3D sliding linear fit and store it in a map.
        LArAnalysisParticleHelper::TrackFitMap trackFitMap;
        LArAnalysisParticleHelper::RecursivelyAppendTrackFitMap(this->GetPandora(), pPfo, trackFitMap, this->m_trackSlidingFitWindow);

        // For each tracklike PFO, decide which hits we want to Birks-correct.
        LArAnalysisParticleHelper::LArTrackHitEnergyMap trackHitEnergyMap;
        this->RecursivelyAppendLArTrackHitEnergyMap(pPfo, trackHitEnergyMap, trackFitMap);

        // For each PFO, get its main MC particle and store it in a map.
        RecursivelyAppendMCParticleMap(pPfo, mcParticleMap);

        if (this->m_produceBirksFitData)
           this->ProduceBirksFitData(pPfo, trackFitMap, mcParticleMap, trackHitEnergyMap);

        if (this->m_produceEnergyFromRangeData)
        {
            this->RecursivelyProduceEnergyFromRangeData(pPfo, trackFitMap, mcParticleMap,
                                                        this->m_pEnergyFromRangeProtonDataNtuple, {PROTON});

            this->RecursivelyProduceEnergyFromRangeData(pPfo, trackFitMap, mcParticleMap,
                                                        this->m_pEnergyFromRangePionMuonDataNtuple, {PI_PLUS, PI_MINUS, MU_MINUS});
        }

        if (this->m_producePidData)
            this->RecursivelyProducePidData(pPfo, trackFitMap, mcParticleMap, LArAnalysisParticleHelper::IsCosmicRay(pPfo));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisDataAlgorithm::RecursivelyAppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo,
                                                              LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap,
                                                              const LArAnalysisParticleHelper::TrackFitMap &trackFitMap) const
{
    const auto findIter = trackFitMap.find(pPfo);

    if (findIter != trackFitMap.end())
        trackHitEnergyMap.emplace(pPfo, this->AppendLArTrackHitEnergyMap(pPfo, findIter->second));

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        this->RecursivelyAppendLArTrackHitEnergyMap(pDaughterPfo, trackHitEnergyMap, trackFitMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArTrackHitEnergy::Vector AnalysisDataAlgorithm::AppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo,
                                                                     const ThreeDSlidingFitResult &trackFit) const
{
    // Get all hits and order them by projection along the track fit.
    const CaloHitList collectionPlaneHits = LArAnalysisParticleHelper::GetHitsOfType(pPfo, TPC_VIEW_W, false);

    const LArAnalysisParticleHelper::HitProjectionVector orderedHitProjections = LArAnalysisParticleHelper::OrderHitsByProjectionOnToTrackFit(
                                                                                                             collectionPlaneHits, trackFit);

    LArTrackHitEnergy::Vector trackHitEnergyVector;

    for (const LArAnalysisParticleHelper::HitProjectionPair &projectionPair : orderedHitProjections)
    {
        const CaloHit *const pCaloHit = projectionPair.first;
        const float coordinate        = projectionPair.second;
        const float threeDDistance    = LArAnalysisParticleHelper::CaloHitToThreeDDistance(this->GetPandora(), pCaloHit, trackFit);

        trackHitEnergyVector.emplace_back(pCaloHit, coordinate, pCaloHit->GetInputEnergy(), threeDDistance);
    }

    this->m_pBirksHitSelectionTool->Run(this, trackHitEnergyVector, this->m_uniquePlotIdentifier);
    return trackHitEnergyVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisDataAlgorithm::RecursivelyAppendMCParticleMap(const ParticleFlowObject *const pPfo, MCParticleMap &mcParticleMap) const
{
    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        // Try to get the corresponding MC particle. Skip this daughter if no MC particle can be found.
        const MCParticle *pMCParticle = nullptr;

        try
        {
            pMCParticle = LArAnalysisParticleHelper::GetMainMCParticle(pPfo);
        }

        catch (...)
        {
            std::cout << "AnalysisDataAlgorithm: failed to find main MC particle" << std::endl;
        }

        if (pMCParticle)
            mcParticleMap.emplace(pPfo, pMCParticle);

        this->RecursivelyAppendMCParticleMap(pDaughterPfo, mcParticleMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisDataAlgorithm::ProduceBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                                   const MCParticleMap &mcParticleMap,
                                                   const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap) const
{
    const auto findIter = mcParticleMap.find(pPfo);

    if (findIter == mcParticleMap.end())
    {
        std::cout << "AnalysisDataAlgorithm: failed to find primary particle in MC particle map" << std::endl;
        return false;
    }

    const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(findIter->second);

    if (!LArMCParticleHelper::IsVisible(pPrimaryMCParticle))
    {
        std::cout << "AnalysisDataAlgorithm: primary MC particle was invisible so ignoring this primary" << std::endl;
        return false;
    }

    float truePrimaryEnergy = 1000.f * (pPrimaryMCParticle->GetEnergy() - PdgTable::GetParticleMass(pPrimaryMCParticle->GetParticleId()));

    float totalNoBirksAdcIntegral = 0.f;
    FloatVector birksAdcIntegrals, threeDDistances;

    if (!this->GetBirksFitData(pPfo, trackFitMap, trackHitEnergyMap, totalNoBirksAdcIntegral, birksAdcIntegrals, threeDDistances))
        return false;

    this->m_pBirksFitDataTree->Branch("TruePrimaryEnergy",       &truePrimaryEnergy);
    this->m_pBirksFitDataTree->Branch("TotalNoBirksAdcIntegral", &totalNoBirksAdcIntegral);
    this->m_pBirksFitDataTree->Branch("BirksAdcIntegrals",       &birksAdcIntegrals);
    this->m_pBirksFitDataTree->Branch("ThreeDDistances",         &threeDDistances);
    this->m_pBirksFitDataTree->Fill();
    this->m_pBirksFitDataTree->ResetBranchAddresses();

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisDataAlgorithm::GetBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                               const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap,
                                               float &totalNoBirksAdcIntegral, FloatVector &birksAdcIntegrals,
                                               FloatVector &threeDDistances) const
{
    const bool isShower = LArPfoHelper::IsShower(pPfo);

    const auto fitFindIter       = trackFitMap.find(pPfo);
    const auto hitEnergyFindIter = trackHitEnergyMap.find(pPfo);

    const bool fitFound       = (fitFindIter != trackFitMap.end());
    const bool hitEnergyFound = (hitEnergyFindIter != trackHitEnergyMap.end());

    if (fitFound != hitEnergyFound)
        std::cout << "AnalysisDataAlgorithm: fits should be found for tracks iff track hit energies are found" << std::endl;

    if (isShower)
    {
        totalNoBirksAdcIntegral += this->AddUpShowerAdcs(pPfo);
        return true;
    }

    if (!fitFound || !hitEnergyFound)
    {
        std::cout << "AnalysisDataAlgorithm: Could not find fit for tracklike particle so ignoring this primary" << std::endl;
        return false;
    }

    this->GetTrackAdcsAndDistances(hitEnergyFindIter->second, totalNoBirksAdcIntegral, birksAdcIntegrals, threeDDistances);

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        if (!this->GetBirksFitData(pDaughterPfo, trackFitMap, trackHitEnergyMap, totalNoBirksAdcIntegral, birksAdcIntegrals,
                                   threeDDistances))
        {
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisDataAlgorithm::AddUpShowerAdcs(const ParticleFlowObject *const pPfo) const
{
    float showerAdcs = 0.f;

    for (const CaloHit *const pCaloHit : LArAnalysisParticleHelper::GetHitsOfType(pPfo, TPC_VIEW_W, true))
        showerAdcs += pCaloHit->GetInputEnergy();

    return showerAdcs;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisDataAlgorithm::GetTrackAdcsAndDistances(const LArTrackHitEnergy::Vector &trackHitEnergies,
                                                        float &totalNoBirksAdcIntegral, FloatVector &birksAdcIntegrals,
                                                        FloatVector &threeDDistances) const
{
    for (const LArTrackHitEnergy &trackHitEnergy : trackHitEnergies)
    {
        if (trackHitEnergy.ApplyCorrection())
        {
            birksAdcIntegrals.push_back(trackHitEnergy.UncorrectedEnergy());
            threeDDistances.push_back(trackHitEnergy.CorrectedEnergy());
        }

        else
            totalNoBirksAdcIntegral += trackHitEnergy.UncorrectedEnergy();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisDataAlgorithm::RecursivelyProduceEnergyFromRangeData(const ParticleFlowObject *const pPfo,
                                                                     const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                                                     const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap,
                                                                     TNtuple *const pNtuple, const PdgCodeSet &pdgCodeSet) const
{
    const auto trackFitFindIter   = trackFitMap.find(pPfo);
    const auto mcParticleFindIter = mcParticleMap.find(pPfo);

    if ((trackFitFindIter != trackFitMap.end()) && (mcParticleFindIter != mcParticleMap.end()))
    {
        const int pdgCode = mcParticleFindIter->second->GetParticleId();

        if (pdgCodeSet.find(pdgCode) != pdgCodeSet.end())
        {
            const float range      = LArAnalysisParticleHelper::GetParticleRange(pPfo, trackFitFindIter->second);
            const float trueEnergy = this->GetTrueEnergy(mcParticleFindIter->second);

            if (range > 0.f && trueEnergy > 0.f)
                pNtuple->Fill(range, trueEnergy);
        }
    }

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        this->RecursivelyProduceEnergyFromRangeData(pDaughterPfo, trackFitMap, mcParticleMap, pNtuple, pdgCodeSet);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisDataAlgorithm::GetTrueEnergy(const MCParticle *pMCParticle) const
{
    return (pMCParticle->GetEnergy() - PdgTable::GetParticleMass(pMCParticle->GetParticleId()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisDataAlgorithm::RecursivelyProducePidData(const ParticleFlowObject *const pPfo,
                                                         const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                                         const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap, const bool isCosmicRay) const
{
    const auto trackFitFindIter   = trackFitMap.find(pPfo);
    const auto mcParticleFindIter = mcParticleMap.find(pPfo);

    if ((trackFitFindIter != trackFitMap.end()) && (mcParticleFindIter != mcParticleMap.end()))
    {
        switch (mcParticleFindIter->second->GetParticleId())
        {
            case PROTON:
            {
                float range      = LArAnalysisParticleHelper::GetParticleRange(pPfo, trackFitFindIter->second);
                float trueEnergy = this->GetTrueEnergy(mcParticleFindIter->second);
                int   isCosmicRayInt = static_cast<int>(isCosmicRay);

                if (range > 0.f && trueEnergy > 0.f)
                {
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pidDataProtonsTTreeName.c_str(), "TrackLength", range));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pidDataProtonsTTreeName.c_str(), "TrueEnergy", trueEnergy));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pidDataProtonsTTreeName.c_str(), "IsCosmicRay", isCosmicRayInt));
                    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_pidDataProtonsTTreeName.c_str()));
                }

                break;
            }

            case PI_MINUS:
            case PI_PLUS:
            case MU_MINUS:
            {
                float range      = LArAnalysisParticleHelper::GetParticleRange(pPfo, trackFitFindIter->second);
                float trueEnergy = this->GetTrueEnergy(mcParticleFindIter->second);
                int   isCosmicRayInt = static_cast<int>(isCosmicRay);

                if (range > 0.f && trueEnergy > 0.f)
                {
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pidDataMuonsPionsTTreeName.c_str(), "TrackLength", range));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pidDataMuonsPionsTTreeName.c_str(), "TrueEnergy", trueEnergy));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pidDataMuonsPionsTTreeName.c_str(), "IsCosmicRay", isCosmicRayInt));
                    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_pidDataMuonsPionsTTreeName.c_str()));
                }

                break;
            }

            default: break;
        }
    }

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        this->RecursivelyProducePidData(pDaughterPfo, trackFitMap, mcParticleMap, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisDataAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceBirksFitData", this->m_produceBirksFitData));
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceEnergyFromRangeData",
                                                                           this->m_produceEnergyFromRangeData));

     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProducePidData", this->m_producePidData));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName",        this->m_pfoListName));

     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowXMargin", this->m_fiducialCutLowXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighXMargin", this->m_fiducialCutHighXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowYMargin", this->m_fiducialCutLowYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighYMargin", this->m_fiducialCutHighYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowZMargin", this->m_fiducialCutLowZMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighZMargin", this->m_fiducialCutHighZMargin));

     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackSlidingFitWindow",
                                                                           this->m_trackSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootDataFileName", this->m_rootDataFileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", this->m_caloHitListName));


    this->m_pRootDataFile = new TFile(this->m_rootDataFileName.c_str(), "UPDATE");

    if (this->m_produceBirksFitData)
        this->m_pBirksFitDataTree = new TTree("BirksFitData", "BirksFitData");

    if (this->m_produceEnergyFromRangeData)
    {
        this->m_pEnergyFromRangeProtonDataNtuple = new TNtuple("EnergyFromRangeDataProtons", "EnergyFromRangeDataProtons",
                                                               "Range:TrueEnergy");
        this->m_pEnergyFromRangePionMuonDataNtuple = new TNtuple("EnergyFromRangeDataPionsMuons", "EnergyFromRangeDataPionsMuons",
                                                                 "Range:TrueEnergy");
    }

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
        "BirksHitSelection", pAlgorithmTool));

    if (!(this->m_pBirksHitSelectionTool = dynamic_cast<BirksHitSelectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}
} // namespace lar_physics_content
