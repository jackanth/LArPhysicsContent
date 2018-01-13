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
#include "DebugDefinitions.h"

#include "TFile.h"

namespace lar_physics_content
{
AnalysisDataAlgorithm::AnalysisDataAlgorithm() :
    m_eventNumber{0U},
    m_produceBirksFitData{true},
    m_produceEnergyFromRangeData{true},
    m_producePidData{true},
    m_pfoListName{},         
    m_fiducialCutXMargin{10.f},
    m_fiducialCutYMargin{20.f},
    m_fiducialCutZMargin{10.f},
    m_trackSlidingFitWindow{25U},
    m_rootDataFileName{},
    m_pRootDataFile{nullptr},
    m_pBirksFitDataTree{nullptr},
    m_pEnergyFromRangeProtonDataNtuple{nullptr},
    m_pEnergyFromRangePionMuonDataNtuple{nullptr},
    m_pPidDataProtonsTTreeName{"PidDataProtons"},
    m_pPidDataMuonsPionsTTreeName{"PidDataMuonsPions"},
    m_uniquePlotIdentifier{0},
    m_pBirksHitSelectionTool(nullptr)
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
    
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_pPidDataProtonsTTreeName.c_str(), m_rootDataFileName.c_str(), "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_pPidDataMuonsPionsTTreeName.c_str(), m_rootDataFileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisDataAlgorithm::Run()
{
    COUT("Event number " << ++this->m_eventNumber);
    
    // Use the detector geometry and the margins to get the maximum and minimum fiducial volume coordinates.
    CartesianVector minCoordinates(0.f, 0.f, 0.f), maxCoordinates(0.f, 0.f, 0.f);
    
    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), this->m_fiducialCutXMargin, this->m_fiducialCutYMargin, this->m_fiducialCutZMargin,
        minCoordinates, maxCoordinates);
    
    // Loop over the reconstructed neutrino PFOs and their primary daughters.
    const PfoList recoNeutrinoList = LArAnalysisParticleHelper::GetRecoNeutrinoList(*this, this->m_pfoListName);
                                                                                  
    for (const ParticleFlowObject *const pPfo : recoNeutrinoList)
    {
        for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        {
            // Check whether fiducial cut is satisfied. If not, ignore this primary PFO.
            if (!LArAnalysisParticleHelper::RecursivelyCheckFiducialCut(pDaughterPfo, minCoordinates, maxCoordinates))
            {
                CERR("Fiducial cut was not satisfied so ignoring this primary");
                continue;
            }
            
            // For each tracklike PFO, try to perform a 3D sliding linear fit and store it in a map.
            LArAnalysisParticleHelper::TrackFitMap trackFitMap;
            LArAnalysisParticleHelper::RecursivelyAppendTrackFitMap(this->GetPandora(), pDaughterPfo, trackFitMap, this->m_trackSlidingFitWindow);
            
            // For each tracklike PFO, decide which hits we want to Birks-correct.
            LArAnalysisParticleHelper::LArTrackHitEnergyMap trackHitEnergyMap;
            this->RecursivelyAppendLArTrackHitEnergyMap(pDaughterPfo, trackHitEnergyMap, trackFitMap);
            
            // For each PFO, get its main MC particle and store it in a map.
            MCParticleMap mcParticleMap;
            RecursivelyAppendMCParticleMap(pDaughterPfo, mcParticleMap);
            
            if (this->m_produceBirksFitData)
               this->ProduceBirksFitData(pDaughterPfo, trackFitMap, mcParticleMap, trackHitEnergyMap);
               
            if (this->m_produceEnergyFromRangeData)
            {
                this->RecursivelyProduceEnergyFromRangeData(pDaughterPfo, trackFitMap, mcParticleMap, 
                                                            this->m_pEnergyFromRangeProtonDataNtuple, {PROTON});
                                                 
                this->RecursivelyProduceEnergyFromRangeData(pDaughterPfo, trackFitMap, mcParticleMap,
                                                            this->m_pEnergyFromRangePionMuonDataNtuple, {PI_PLUS, PI_MINUS, MU_MINUS});
            }
   
            if (this->m_producePidData)
                this->RecursivelyProducePidData(pDaughterPfo, trackFitMap, mcParticleMap);
        }
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
    // Try to get the corresponding MC particle. Skip this daughter if no MC particle can be found.
    const MCParticle *pMCParticle = nullptr;

    try
    {
        pMCParticle = LArMCParticleHelper::GetMainMCParticle(pPfo);
    }

    catch (...)
    {
        CERR("Failed to find main MC particle");
    }

    if (pMCParticle)
        mcParticleMap.emplace(pPfo, pMCParticle);
    
    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        this->RecursivelyAppendMCParticleMap(pDaughterPfo, mcParticleMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisDataAlgorithm::ProduceBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                                   const MCParticleMap &mcParticleMap, 
                                                   const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap) const
{
    const auto findIter = mcParticleMap.find(pPfo);
    
    if (findIter == mcParticleMap.end())
    {
        CERR("Failed to find primary particle in MC particle map");
        return false;
    }
    
    const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(findIter->second);
    
    if (!LArMCParticleHelper::IsVisible(pPrimaryMCParticle))
    {
        CERR("Primary MC particle was invisible so ignoring this primary");
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
        CERR("Fits should be found for tracks iff track hit energies are found");
        
    if (isShower)
    {
        totalNoBirksAdcIntegral += this->AddUpShowerAdcs(pPfo);
        return true;
    }
    
    if (!fitFound || !hitEnergyFound)
    {
        CERR("Could not find fit for tracklike particle so ignoring this primary");
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
                                                         const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap) const
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
                
                if (range > 0.f && trueEnergy > 0.f)
                {
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pPidDataProtonsTTreeName.c_str(), "TrackLength", range));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pPidDataProtonsTTreeName.c_str(), "TrueEnergy", trueEnergy));
                    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_pPidDataProtonsTTreeName.c_str()));
                }
                
                break;
            }
            
            case PI_MINUS:
            case PI_PLUS:
            case MU_MINUS:
            {
                float range      = LArAnalysisParticleHelper::GetParticleRange(pPfo, trackFitFindIter->second);
                float trueEnergy = this->GetTrueEnergy(mcParticleFindIter->second);
                
                if (range > 0.f && trueEnergy > 0.f)
                {
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pPidDataMuonsPionsTTreeName.c_str(), "TrackLength", range));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_pPidDataMuonsPionsTTreeName.c_str(), "TrueEnergy", trueEnergy));
                    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_pPidDataMuonsPionsTTreeName.c_str()));
                }
                
                break;
            }
                
            default: break;
        }
    }
    
    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        this->RecursivelyProducePidData(pDaughterPfo, trackFitMap, mcParticleMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisDataAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceBirksFitData", this->m_produceBirksFitData));
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceEnergyFromRangeData",
                                                                           this->m_produceEnergyFromRangeData));
                                                                           
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProducePidData", this->m_producePidData));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName",        this->m_pfoListName));
                                                                           
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutXMargin",  this->m_fiducialCutXMargin));
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutYMargin",  this->m_fiducialCutYMargin));
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutZMargin",  this->m_fiducialCutZMargin));
    
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackSlidingFitWindow",
                                                                           this->m_trackSlidingFitWindow));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootDataFileName", this->m_rootDataFileName));
          
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
