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
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minCoordinates(0.f, 0.f, 0.f),
    m_maxCoordinates(0.f, 0.f, 0.f),
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
    m_pTrackHitEnergyTool(nullptr),
    m_pMcInfoTool(nullptr),
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

        LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo;
        m_pMcInfoTool->Run(this, pMCParticle, pfoMcInfo);

//        if (mcContainmentFraction < 0.95f || mcHitPurity < 0.95f || mcHitCompleteness < 0.95f || mcCollectionPlaneHitPurity < 0.95f ||
//            mcCollectionPlaneHitCompleteness < 0.95f)
//        {
//            std::cout << "Skipping..." << std::endl;
//            continue;
//        }

        // For each tracklike PFO, try to perform a 3D sliding linear fit and store it in a map.
        LArAnalysisParticleHelper::FittedTrackInfoMap fittedTrackInfoMap;
        float excessCaloValue(0.f);

        m_pTrackHitEnergyTool->Run(this, pPfo, fittedTrackInfoMap, excessCaloValue,
            [&](LArFittedTrackInfo::TrackHitValueVector &trackHitValueVector, float &excessCaloValue) -> bool
            {
                return m_pHitPurityTool->Run(this, trackHitValueVector, excessCaloValue);
            });

        // For each PFO, get its main MC particle and store it in a map.
        RecursivelyAppendMCParticleMap(pPfo, mcParticleMap);

        if (m_produceBirksFitData)
           this->ProduceBirksFitData(pPfo, fittedTrackInfoMap, mcParticleMap);

        if (m_produceEnergyFromRangeData)
        {
            this->RecursivelyProduceEnergyFromRangeData(pPfo, fittedTrackInfoMap, mcParticleMap,
                                                        m_pEnergyFromRangeProtonDataNtuple, {PROTON});

            this->RecursivelyProduceEnergyFromRangeData(pPfo, fittedTrackInfoMap, mcParticleMap, m_pEnergyFromRangePionMuonDataNtuple, {PI_PLUS, PI_MINUS, MU_MINUS});
        }

        if (m_producePidData)
            this->RecursivelyProducePidData(pPfo, fittedTrackInfoMap, mcParticleMap, LArAnalysisParticleHelper::IsCosmicRay(pPfo));
    }

    return STATUS_CODE_SUCCESS;
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

bool AnalysisDataAlgorithm::ProduceBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                                                   const MCParticleMap &mcParticleMap) const
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

    if (!this->GetBirksFitData(pPfo, fittedTrackInfoMap, totalNoBirksAdcIntegral, birksAdcIntegrals, threeDDistances))
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

bool AnalysisDataAlgorithm::GetBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                                               float &totalNoBirksAdcIntegral, FloatVector &birksAdcIntegrals,
                                               FloatVector &threeDDistances) const
{
    if (LArPfoHelper::IsShower(pPfo))
    {
        totalNoBirksAdcIntegral += this->AddUpShowerAdcs(pPfo);
        return true;
    }

    const auto findIter = fittedTrackInfoMap.find(pPfo);
    if (findIter == fittedTrackInfoMap.end())
    {
        std::cout << "AnalysisDataAlgorithm: Could not find fit for tracklike particle so ignoring this primary" << std::endl;
        return false;
    }

    this->GetTrackAdcsAndDistances(findIter->second.HitChargeVector(), totalNoBirksAdcIntegral, birksAdcIntegrals, threeDDistances);

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        if (!this->GetBirksFitData(pDaughterPfo, fittedTrackInfoMap, totalNoBirksAdcIntegral, birksAdcIntegrals, threeDDistances))
            return false;
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

void AnalysisDataAlgorithm::GetTrackAdcsAndDistances(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergies,
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
                                                                     const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                                                                     const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap,
                                                                     TNtuple *const pNtuple, const PdgCodeSet &pdgCodeSet) const
{
    const auto trackFitFindIter   = fittedTrackInfoMap.find(pPfo);
    const auto mcParticleFindIter = mcParticleMap.find(pPfo);

    if ((trackFitFindIter != fittedTrackInfoMap.end()) && (mcParticleFindIter != mcParticleMap.end()))
    {
        const int pdgCode = mcParticleFindIter->second->GetParticleId();

        if (pdgCodeSet.find(pdgCode) != pdgCodeSet.end())
        {
            const float range      = trackFitFindIter->second.Range();
            const float trueEnergy = this->GetTrueEnergy(mcParticleFindIter->second);

            if (range > 0.f && trueEnergy > 0.f)
                pNtuple->Fill(range, trueEnergy);
        }
    }

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        this->RecursivelyProduceEnergyFromRangeData(pDaughterPfo, fittedTrackInfoMap, mcParticleMap, pNtuple, pdgCodeSet);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AnalysisDataAlgorithm::GetTrueEnergy(const MCParticle *pMCParticle) const
{
    return (pMCParticle->GetEnergy() - PdgTable::GetParticleMass(pMCParticle->GetParticleId()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisDataAlgorithm::RecursivelyProducePidData(const ParticleFlowObject *const pPfo,
                                                         const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                                                         const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap, const bool isCosmicRay) const
{
    const auto trackFitFindIter   = fittedTrackInfoMap.find(pPfo);
    const auto mcParticleFindIter = mcParticleMap.find(pPfo);

    if ((trackFitFindIter != fittedTrackInfoMap.end()) && (mcParticleFindIter != mcParticleMap.end()))
    {
        switch (mcParticleFindIter->second->GetParticleId())
        {
            case PROTON:
            {
                float range      = trackFitFindIter->second.Range();
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
                float range      = trackFitFindIter->second.Range();
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
        this->RecursivelyProducePidData(pDaughterPfo, fittedTrackInfoMap, mcParticleMap, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisDataAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowMargins", m_fiducialCutLowMargins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighMargins", m_fiducialCutHighMargins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceBirksFitData", m_produceBirksFitData));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceEnergyFromRangeData",
                                                                           m_produceEnergyFromRangeData));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProducePidData", m_producePidData));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName",        m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootDataFileName", m_rootDataFileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    // Use the detector geometry and the margins to get the maximum and minimum fiducial volume coordinates.
    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), m_fiducialCutLowMargins, m_fiducialCutHighMargins, m_minCoordinates,
        m_maxCoordinates);

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

    return STATUS_CODE_SUCCESS;
}
} // namespace lar_physics_content
