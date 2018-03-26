/**
 *  @file LArPhysicsContent/src/AnalysisDataAlgorithm.cxx
 *
 *  @brief Implementation of the LEE analysis data algorithm class.
 *
 *  $Log: $
 */

#include "larphysicscontent/AnalysisDataAlgorithm.h"
#include "larphysicscontent/LArAnalysisParticleHelper.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/PdgTable.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

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
    m_mcContainmentFractionLowerBound(0.9f),
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
    m_pHitPurityTool(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AnalysisDataAlgorithm::~AnalysisDataAlgorithm()
{
    if (m_pBirksFitDataTree)
        m_pBirksFitDataTree->Write();

    if (m_pEnergyFromRangeProtonDataNtuple)
        m_pEnergyFromRangeProtonDataNtuple->Write();

    if (m_pEnergyFromRangePionMuonDataNtuple)
        m_pEnergyFromRangePionMuonDataNtuple->Write();

    if (m_pRootDataFile)
    {
        if (m_pRootDataFile->IsOpen())
            m_pRootDataFile->Close();

        delete m_pRootDataFile;
    }

    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_pidDataProtonsTTreeName.c_str(), m_rootDataFileName.c_str(), "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_pidDataMuonsPionsTTreeName.c_str(), m_rootDataFileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisDataAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);

    if ((PandoraContentApi::GetList(*this, m_pfoListName, pPfoList) != STATUS_CODE_SUCCESS) || !pPfoList)
    {
        std::cout << "AnalysisAlgorithm: failed to get the PFO list" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    const CaloHitList *pCaloHitList(nullptr);

    if ((PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS) || !pCaloHitList)
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
        m_fiducialCutLowYMargin, m_fiducialCutHighYMargin, m_fiducialCutLowZMargin, m_fiducialCutHighZMargin,
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

        const LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo = LArAnalysisParticleHelper::GetMcInformation(pMCParticle, minCoordinates, maxCoordinates,
            m_mcContainmentFractionLowerBound);

//        if (mcContainmentFraction < 0.95f || mcHitPurity < 0.95f || mcHitCompleteness < 0.95f || mcCollectionPlaneHitPurity < 0.95f ||
//            mcCollectionPlaneHitCompleteness < 0.95f)
//        {
//            std::cout << "Skipping..." << std::endl;
//            continue;
//        }

        // For each tracklike PFO, try to perform a 3D sliding linear fit and store it in a map.
        LArAnalysisParticleHelper::TrackFitMap trackFitMap;
        LArAnalysisParticleHelper::RecursivelyAppendTrackFitMap(this->GetPandora(), pPfo, trackFitMap, m_trackSlidingFitWindow);

        // For each tracklike PFO, decide which hits we want to Birks-correct.
        LArAnalysisParticleHelper::LArTrackHitEnergyMap trackHitEnergyMap;
        this->RecursivelyAppendLArTrackHitEnergyMap(pPfo, trackHitEnergyMap, trackFitMap);

        // For each PFO, get its main MC particle and store it in a map.
        RecursivelyAppendMCParticleMap(pPfo, mcParticleMap);

        if (m_produceBirksFitData)
           this->ProduceBirksFitData(pPfo, trackFitMap, mcParticleMap, trackHitEnergyMap);

        if (m_produceEnergyFromRangeData)
        {
            this->RecursivelyProduceEnergyFromRangeData(pPfo, trackFitMap, mcParticleMap,
                                                        m_pEnergyFromRangeProtonDataNtuple, {PROTON});

            this->RecursivelyProduceEnergyFromRangeData(pPfo, trackFitMap, mcParticleMap, m_pEnergyFromRangePionMuonDataNtuple, {PI_PLUS, PI_MINUS, MU_MINUS});
        }

        if (m_producePidData)
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

LArAnalysisParticleHelper::TrackHitValueVector AnalysisDataAlgorithm::AppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo,
    const ThreeDSlidingFitResult &trackFit) const
{
    // Get all hits and order them by projection along the track fit.
    const CaloHitList collectionPlaneHits = LArAnalysisParticleHelper::GetHitsOfType(pPfo, TPC_VIEW_W, false);

    const LArAnalysisParticleHelper::HitProjectionVector orderedHitProjections = LArAnalysisParticleHelper::OrderHitsByProjectionOnToTrackFit(
                                                                                                             collectionPlaneHits, trackFit);

    LArAnalysisParticleHelper::TrackHitValueVector trackHitEnergyVector;

    for (const LArAnalysisParticleHelper::HitProjectionPair &projectionPair : orderedHitProjections)
    {
        const CaloHit *const pCaloHit = projectionPair.first;
        const float coordinate        = projectionPair.second;
        const float threeDDistance    = LArAnalysisParticleHelper::CaloHitToThreeDDistance(this->GetPandora(), pCaloHit, trackFit);

        trackHitEnergyVector.emplace_back(pCaloHit, coordinate, threeDDistance, pCaloHit->GetInputEnergy());
    }

    float excessCharge(0.f);
    m_pHitPurityTool->Run(this, trackHitEnergyVector, excessCharge);
    // TODO excess charge
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

    m_pBirksFitDataTree->Branch("TruePrimaryEnergy",       &truePrimaryEnergy);
    m_pBirksFitDataTree->Branch("TotalNoBirksAdcIntegral", &totalNoBirksAdcIntegral);
    m_pBirksFitDataTree->Branch("BirksAdcIntegrals",       &birksAdcIntegrals);
    m_pBirksFitDataTree->Branch("ThreeDDistances",         &threeDDistances);
    m_pBirksFitDataTree->Fill();
    m_pBirksFitDataTree->ResetBranchAddresses();

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

void AnalysisDataAlgorithm::GetTrackAdcsAndDistances(const LArAnalysisParticleHelper::TrackHitValueVector &trackHitEnergies,
    float &totalNoBirksAdcIntegral, FloatVector &birksAdcIntegrals, FloatVector &threeDDistances) const
{
    for (const LArTrackHitValue &trackHitEnergy : trackHitEnergies)
    {
            birksAdcIntegrals.push_back(trackHitEnergy.CaloValue());
            threeDDistances.push_back(trackHitEnergy.ThreeDDistance());
          //  totalNoBirksAdcIntegral += trackHitEnergy.UncorrectedEnergy();
          (void) totalNoBirksAdcIntegral;
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceBirksFitData", m_produceBirksFitData));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceEnergyFromRangeData",
                                                                           m_produceEnergyFromRangeData));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProducePidData", m_producePidData));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName",        m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowXMargin", m_fiducialCutLowXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighXMargin", m_fiducialCutHighXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowYMargin", m_fiducialCutLowYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighYMargin", m_fiducialCutHighYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowZMargin", m_fiducialCutLowZMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighZMargin", m_fiducialCutHighZMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McContainmentFractionLowerBound", m_mcContainmentFractionLowerBound));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackSlidingFitWindow", m_trackSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootDataFileName", m_rootDataFileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));


    m_pRootDataFile = new TFile(m_rootDataFileName.c_str(), "UPDATE");

    if (m_produceBirksFitData)
        m_pBirksFitDataTree = new TTree("BirksFitData", "BirksFitData");

    if (m_produceEnergyFromRangeData)
    {
        m_pEnergyFromRangeProtonDataNtuple = new TNtuple("EnergyFromRangeDataProtons", "EnergyFromRangeDataProtons",
                                                               "Range:TrueEnergy");
        m_pEnergyFromRangePionMuonDataNtuple = new TNtuple("EnergyFromRangeDataPionsMuons", "EnergyFromRangeDataPionsMuons",
                                                                 "Range:TrueEnergy");
    }

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
        "HitPurity", pAlgorithmTool));

    if (!(m_pHitPurityTool = dynamic_cast<HitPurityTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}
} // namespace lar_physics_content
