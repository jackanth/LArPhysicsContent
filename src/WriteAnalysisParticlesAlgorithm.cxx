/**
 *  @file LArPhysicsContent/src/WriteAnalysisParticlesAlgorithm.cxx
 *
 *  @brief Implementation of the write AnalysisParticles algorithm class.
 *
 *  $Log: $
 */

#include "WriteAnalysisParticlesAlgorithm.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
    
WriteAnalysisParticlesAlgorithm::WriteAnalysisParticlesAlgorithm() :
    m_neutrinoPfoListName(),
    m_cosmicRayPfoListName(),
    m_outputFile(),
    m_pOutputTFile(nullptr),
    m_pOutputTree(nullptr),
    m_verbose(false),
    m_treeParameters()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

WriteAnalysisParticlesAlgorithm::~WriteAnalysisParticlesAlgorithm()
{
    if (this->m_pOutputTree)
        this->m_pOutputTree->Write();
        
    if (this->m_pOutputTFile)
    {
        if (this->m_pOutputTFile->IsOpen())
            this->m_pOutputTFile->Close();
            
        delete this->m_pOutputTFile;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode WriteAnalysisParticlesAlgorithm::Run()
{
    // Get the input neutrino PFO list.
    const PfoList *pNeutrinoPfoList(nullptr);

    if ((PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pNeutrinoPfoList) != STATUS_CODE_SUCCESS) || !pNeutrinoPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "WriteAnalysisParticlesAlgorithm: cannot find pfo list " << m_neutrinoPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }
    
    // Get the input cosmic ray PFO list.
  /*  const PfoList *pCosmicRayPfoList(nullptr);

    if ((PandoraContentApi::GetList(*this, m_cosmicRayPfoListName, pCosmicRayPfoList) != STATUS_CODE_SUCCESS) || !pCosmicRayPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "WriteAnalysisParticlesAlgorithm: cannot find pfo list " << m_cosmicRayPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }*/
    
    // Get the unique neutrino, if one exists.
    if (pNeutrinoPfoList->size() > 1)
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: more than one reconstructed neutrino, skipping this event" << std::endl;
        return STATUS_CODE_SUCCESS;
    }
    
    const LArAnalysisParticle *pNeutrinoAnalysisParticle(nullptr);
    
    if (!pNeutrinoPfoList->empty())
    {
        if (!(pNeutrinoAnalysisParticle = dynamic_cast<const LArAnalysisParticle *>(*pNeutrinoPfoList->begin())))
        {
            std::cout << "WriteAnalysisParticlesAlgorithm: neutrino was reconstructed but could not be cast to an analysis particle "
                      << "object" << std::endl;
                      
            return STATUS_CODE_SUCCESS;
        }
    }
    
    // Populate the parameters of the tree for the neutrinos, its daughters, and the cosmic rays.
    this->m_treeParameters = TreeParameters();
    
    if (pNeutrinoAnalysisParticle)
    {
        this->PopulateNeutrinoParameters(*pNeutrinoAnalysisParticle); // if there is no neutrino, the default parameter set is fine
        
        AnalysisParticleList analysisParticleList;
        
        for (const ParticleFlowObject *const pPrimaryPfo : pNeutrinoAnalysisParticle->GetDaughterPfoList())
        {
            if (const LArAnalysisParticle *const pAnalysisParticle = dynamic_cast<const LArAnalysisParticle * const>(pPrimaryPfo))
                analysisParticleList.push_back(pAnalysisParticle);
        }
        
        this->PopulatePrimaryDaughterParameters(analysisParticleList);
    }
    
    //this->PopulateCosmicRayParameters(*pCosmicRayPfoList);
    
    if (this->m_verbose)
        this->PrintTree();
    
    this->m_pOutputTree->Fill();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoParameters(const LArAnalysisParticle &neutrinoAnalysisParticle) const
{
    (void) neutrinoAnalysisParticle;
    
    /*
    m_nu_Exists(false),
    m_nu_IsVertexFiducial(false),
    m_nu_AreAllHitsFiducial(false),
    m_nu_HasMcInfo(false),
    m_nu_Energy(0.f),
    m_nu_EnergyFromChargeOnly(0.f),
    m_nu_LongitudinalEnergy(0.f),
    m_nu_TransverseEnergy(0.f),
    m_nu_VertexX(0.f),
    m_nu_VertexY(0.f),
    m_nu_VertexZ(0.f),
    m_nu_DirectionCosineX(0.f),
    m_nu_DirectionCosineY(0.f),
    m_nu_DirectionCosineZ(0.f),
    m_nu_MomentumX(0.f),
    m_nu_MomentumY(0.f),
    m_nu_MomentumZ(0.f),
    m_nu_NumberOf3dHits(0U),
    m_nu_NumberOfCollectionPlaneHits(0U),
    m_nu_mc_Energy(0.f),
    m_nu_mc_LongitudinalEnergy(0.f),
    m_nu_mc_TransverseEnergy(0.f),
    m_nu_mc_VisibleLongitudinalEnergy(0.f),
    m_nu_mc_VisibleTransverseEnergy(0.f),
    m_nu_mc_VertexX(0.f),
    m_nu_mc_VertexY(0.f),
    m_nu_mc_VertexZ(0.f),
    m_nu_mc_DirectionCosineX(0.f),
    m_nu_mc_DirectionCosineY(0.f),
    m_nu_mc_DirectionCosineZ(0.f),
    m_nu_mc_MomentumX(0.f),
    m_nu_mc_MomentumY(0.f),
    m_nu_mc_MomentumZ(0.f),
    m_nu_mc_IsVertexFiducial(false),
    m_nu_mc_IsContained(false),
    m_nu_mc_IsCorrectlyReconstructed(false),
    m_nu_mc_InteractionType(0),
    m_nu_mc_IsChargedCurrent(false),
    m_nu_mc_VisibleEnergyFraction(0.f),
    m_nu_mc_PdgCode(0),
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulatePrimaryDaughterParameters(const AnalysisParticleList &analysisParticleList) const
{
    for (const LArAnalysisParticle *const pAnalysisParticle : analysisParticleList)
    {
        if (pAnalysisParticle)
            this->AddPrimaryDaughterRecord(*pAnalysisParticle);
    }
    
    // TODO
    
    /*
    m_primary_Number(0U)
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddPrimaryDaughterRecord(const LArAnalysisParticle &analysisParticle) const
{
    // TODO
    (void) analysisParticle;
    
    /*
    m_primary_IsVertexFiducial(),
    m_primary_AreAllHitsFiducial(),
    m_primary_HasMcInfo(),
    m_primary_Energy(),
    m_primary_EnergyFromChargeOnly(),
    m_primary_VertexX(),
    m_primary_VertexY(),
    m_primary_VertexZ(),
    m_primary_DirectionCosineX(),
    m_primary_DirectionCosineY(),
    m_primary_DirectionCosineZ(),
    m_primary_MomentumX(),
    m_primary_MomentumY(),
    m_primary_MomentumZ(),
    m_primary_ParticleType(),
    m_primary_NumberOf3dHits(),
    m_primary_NumberOfCollectionPlaneHits(),
    m_primary_IsShower(),
    m_primary_NumberOfDownstreamParticles(),
    m_primary_mc_Energy(),
    m_primary_mc_VertexX(),
    m_primary_mc_VertexY(),
    m_primary_mc_VertexZ(),
    m_primary_mc_DirectionCosineX(),
    m_primary_mc_DirectionCosineY(),
    m_primary_mc_DirectionCosineZ(),
    m_primary_mc_MomentumX(),
    m_primary_mc_MomentumY(),
    m_primary_mc_MomentumZ(),
    m_primary_mc_IsVertexFiducial(),
    m_primary_mc_IsContained(),
    m_primary_mc_IsCorrectlyReconstructed(),
    m_primary_mc_ParticleType(),
    m_primary_mc_IsShower(),
    m_primary_mc_PdgCode(),
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateCosmicRayParameters(const PfoList &cosmicRayPfoList) const
{
    for (const ParticleFlowObject *const pCosmicRayPfo : cosmicRayPfoList)
    {
        if (pCosmicRayPfo)
            this->AddCosmicRayRecord(*pCosmicRayPfo);
    }
    
    // TO DO
    
    /*
    m_cr_Number(0U)
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddCosmicRayRecord(const ParticleFlowObject &cosmicRayPfo) const
{
    // TODO
    (void) cosmicRayPfo;
    
    /*
    m_cr_IsVertexFiducial(),
    m_cr_AreAllHitsFiducial(),
    m_cr_HasMcInfo(),
    m_cr_Energy(),
    m_cr_EnergyFromChargeOnly(),
    m_cr_VertexX(),
    m_cr_VertexY(),
    m_cr_VertexZ(),
    m_cr_DirectionCosineX(),
    m_cr_DirectionCosineY(),
    m_cr_DirectionCosineZ(),
    m_cr_MomentumX(),
    m_cr_MomentumY(),
    m_cr_MomentumZ(),
    m_cr_NumberOf3dHits(),
    m_cr_NumberOfCollectionPlaneHits(),
    m_cr_mc_Energy(),
    m_cr_mc_VertexX(),
    m_cr_mc_VertexY(),
    m_cr_mc_VertexZ(),
    m_cr_mc_DirectionCosineX(),
    m_cr_mc_DirectionCosineY(),
    m_cr_mc_DirectionCosineZ(),
    m_cr_mc_MomentumX(),
    m_cr_mc_MomentumY(),
    m_cr_mc_MomentumZ(),
    m_cr_mc_IsVertexFiducial(),
    m_cr_mc_IsContained(),
    m_cr_mc_IsCorrectlyReconstructed()
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PrintTree() const
{
    // TODO
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool WriteAnalysisParticlesAlgorithm::IsChargedCurrent(const LArMCParticleHelper::InteractionType interactionType) const
{
    switch (interactionType)
    {
        case LArMCParticleHelper::CCQEL_MU:
        case LArMCParticleHelper::CCQEL_MU_P:
        case LArMCParticleHelper::CCQEL_MU_P_P:
        case LArMCParticleHelper::CCQEL_MU_P_P_P:
        case LArMCParticleHelper::CCQEL_MU_P_P_P_P:
        case LArMCParticleHelper::CCQEL_MU_P_P_P_P_P:
        case LArMCParticleHelper::CCQEL_E:
        case LArMCParticleHelper::CCQEL_E_P:
        case LArMCParticleHelper::CCQEL_E_P_P:
        case LArMCParticleHelper::CCQEL_E_P_P_P:
        case LArMCParticleHelper::CCQEL_E_P_P_P_P:
        case LArMCParticleHelper::CCQEL_E_P_P_P_P_P:
        case LArMCParticleHelper::CCRES_MU:
        case LArMCParticleHelper::CCRES_MU_P:
        case LArMCParticleHelper::CCRES_MU_P_P:
        case LArMCParticleHelper::CCRES_MU_P_P_P:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P_P:
        case LArMCParticleHelper::CCRES_MU_PIPLUS:
        case LArMCParticleHelper::CCRES_MU_P_PIPLUS:
        case LArMCParticleHelper::CCRES_MU_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_MU_P_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_MU_PHOTON:
        case LArMCParticleHelper::CCRES_MU_P_PHOTON:
        case LArMCParticleHelper::CCRES_MU_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_MU_P_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_MU_PIZERO:
        case LArMCParticleHelper::CCRES_MU_P_PIZERO:
        case LArMCParticleHelper::CCRES_MU_P_P_PIZERO:
        case LArMCParticleHelper::CCRES_MU_P_P_P_PIZERO:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P_PIZERO:
        case LArMCParticleHelper::CCRES_MU_P_P_P_P_P_PIZERO:
        case LArMCParticleHelper::CCRES_E:
        case LArMCParticleHelper::CCRES_E_P:
        case LArMCParticleHelper::CCRES_E_P_P:
        case LArMCParticleHelper::CCRES_E_P_P_P:
        case LArMCParticleHelper::CCRES_E_P_P_P_P:
        case LArMCParticleHelper::CCRES_E_P_P_P_P_P:
        case LArMCParticleHelper::CCRES_E_PIPLUS:
        case LArMCParticleHelper::CCRES_E_P_PIPLUS:
        case LArMCParticleHelper::CCRES_E_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_E_P_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_E_P_P_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_E_P_P_P_P_P_PIPLUS:
        case LArMCParticleHelper::CCRES_E_PHOTON:
        case LArMCParticleHelper::CCRES_E_P_PHOTON:
        case LArMCParticleHelper::CCRES_E_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_E_P_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_E_P_P_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_E_P_P_P_P_P_PHOTON:
        case LArMCParticleHelper::CCRES_E_PIZERO:
        case LArMCParticleHelper::CCRES_E_P_PIZERO:
        case LArMCParticleHelper::CCRES_E_P_P_PIZERO:
        case LArMCParticleHelper::CCRES_E_P_P_P_PIZERO:
        case LArMCParticleHelper::CCRES_E_P_P_P_P_PIZERO:
        case LArMCParticleHelper::CCRES_E_P_P_P_P_P_PIZERO:
        case LArMCParticleHelper::CCDIS:
        case LArMCParticleHelper::CCCOH:
            return true;
            
        default: break;
    };
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode WriteAnalysisParticlesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", this->m_neutrinoPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CosmicRayPfoListName", 
        this->m_cosmicRayPfoListName));
        
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", this->m_outputFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Verbose",
        this->m_verbose));
        
    this->m_pOutputTFile = new TFile(this->m_outputFile.c_str(), "UPDATE");
    this->m_pOutputTree = new TTree("PandoraTree", "PandoraTree");
    
    // Neutrino parameters.
    this->m_pOutputTree->Branch("nu_Exists",                           &this->m_treeParameters.m_nu_Exists);
    this->m_pOutputTree->Branch("nu_IsVertexFiducial",                 &this->m_treeParameters.m_nu_IsVertexFiducial);
    this->m_pOutputTree->Branch("nu_AreAllHitsFiducial",               &this->m_treeParameters.m_nu_AreAllHitsFiducial);
    this->m_pOutputTree->Branch("nu_HasMcInfo",                        &this->m_treeParameters.m_nu_HasMcInfo);
    this->m_pOutputTree->Branch("nu_Energy",                           &this->m_treeParameters.m_nu_Energy);
    this->m_pOutputTree->Branch("nu_EnergyFromChargeOnly",             &this->m_treeParameters.m_nu_EnergyFromChargeOnly);
    this->m_pOutputTree->Branch("nu_LongitudinalEnergy",               &this->m_treeParameters.m_nu_LongitudinalEnergy);
    this->m_pOutputTree->Branch("nu_TransverseEnergy",                 &this->m_treeParameters.m_nu_TransverseEnergy);
    this->m_pOutputTree->Branch("nu_VertexX",                          &this->m_treeParameters.m_nu_VertexX);
    this->m_pOutputTree->Branch("nu_VertexY",                          &this->m_treeParameters.m_nu_VertexY);
    this->m_pOutputTree->Branch("nu_VertexZ",                          &this->m_treeParameters.m_nu_VertexZ);
    this->m_pOutputTree->Branch("nu_DirectionCosineX",                 &this->m_treeParameters.m_nu_DirectionCosineX);
    this->m_pOutputTree->Branch("nu_DirectionCosineY",                 &this->m_treeParameters.m_nu_DirectionCosineY);
    this->m_pOutputTree->Branch("nu_DirectionCosineZ",                 &this->m_treeParameters.m_nu_DirectionCosineZ);
    this->m_pOutputTree->Branch("nu_MomentumX",                        &this->m_treeParameters.m_nu_MomentumX);
    this->m_pOutputTree->Branch("nu_MomentumY",                        &this->m_treeParameters.m_nu_MomentumY);
    this->m_pOutputTree->Branch("nu_MomentumZ",                        &this->m_treeParameters.m_nu_MomentumZ);
    this->m_pOutputTree->Branch("nu_NumberOf3dHits",                   &this->m_treeParameters.m_nu_NumberOf3dHits);
    this->m_pOutputTree->Branch("nu_NumberOfCollectionPlaneHits",      &this->m_treeParameters.m_nu_NumberOfCollectionPlaneHits);
    this->m_pOutputTree->Branch("nu_mc_Energy",                        &this->m_treeParameters.m_nu_mc_Energy);
    this->m_pOutputTree->Branch("nu_mc_LongitudinalEnergy",            &this->m_treeParameters.m_nu_mc_LongitudinalEnergy);
    this->m_pOutputTree->Branch("nu_mc_TransverseEnergy",              &this->m_treeParameters.m_nu_mc_TransverseEnergy);
    this->m_pOutputTree->Branch("nu_mc_VisibleLongitudinalEnergy",     &this->m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy);
    this->m_pOutputTree->Branch("nu_mc_VisibleTransverseEnergy",       &this->m_treeParameters.m_nu_mc_VisibleTransverseEnergy);
    this->m_pOutputTree->Branch("nu_mc_VertexX",                       &this->m_treeParameters.m_nu_mc_VertexX);
    this->m_pOutputTree->Branch("nu_mc_VertexY",                       &this->m_treeParameters.m_nu_mc_VertexY);
    this->m_pOutputTree->Branch("nu_mc_VertexZ",                       &this->m_treeParameters.m_nu_mc_VertexZ);
    this->m_pOutputTree->Branch("nu_mc_DirectionCosineX",              &this->m_treeParameters.m_nu_mc_DirectionCosineX);
    this->m_pOutputTree->Branch("nu_mc_DirectionCosineY",              &this->m_treeParameters.m_nu_mc_DirectionCosineY);
    this->m_pOutputTree->Branch("nu_mc_DirectionCosineZ",              &this->m_treeParameters.m_nu_mc_DirectionCosineZ);
    this->m_pOutputTree->Branch("nu_mc_MomentumX",                     &this->m_treeParameters.m_nu_mc_MomentumX);
    this->m_pOutputTree->Branch("nu_mc_MomentumY",                     &this->m_treeParameters.m_nu_mc_MomentumY);
    this->m_pOutputTree->Branch("nu_mc_MomentumZ",                     &this->m_treeParameters.m_nu_mc_MomentumZ);
    this->m_pOutputTree->Branch("nu_mc_IsVertexFiducial",              &this->m_treeParameters.m_nu_mc_IsVertexFiducial);
    this->m_pOutputTree->Branch("nu_mc_IsContained",                   &this->m_treeParameters.m_nu_mc_IsContained);
    this->m_pOutputTree->Branch("nu_mc_IsCorrectlyReconstructed",      &this->m_treeParameters.m_nu_mc_IsCorrectlyReconstructed);
    this->m_pOutputTree->Branch("nu_mc_InteractionType",               &this->m_treeParameters.m_nu_mc_InteractionType);
    this->m_pOutputTree->Branch("nu_mc_IsChargedCurrent",              &this->m_treeParameters.m_nu_mc_IsChargedCurrent);
    this->m_pOutputTree->Branch("nu_mc_VisibleEnergyFraction",         &this->m_treeParameters.m_nu_mc_VisibleEnergyFraction);
    this->m_pOutputTree->Branch("nu_mc_PdgCode",                       &this->m_treeParameters.m_nu_mc_PdgCode);
    
    // Primary neutrino daughter parameters.
    this->m_pOutputTree->Branch("primary_Number",                      &this->m_treeParameters.m_primary_Number);
    
    this->m_pOutputTree->Branch("primary_IsVertexFiducial",            &this->m_treeParameters.m_primary_IsVertexFiducial);
    this->m_pOutputTree->Branch("primary_AreAllHitsFiducial",          &this->m_treeParameters.m_primary_AreAllHitsFiducial);
    this->m_pOutputTree->Branch("primary_HasMcInfo",                   &this->m_treeParameters.m_primary_HasMcInfo);
    this->m_pOutputTree->Branch("primary_Energy",                      &this->m_treeParameters.m_primary_Energy);
    this->m_pOutputTree->Branch("primary_EnergyFromChargeOnly",        &this->m_treeParameters.m_primary_EnergyFromChargeOnly);
    this->m_pOutputTree->Branch("primary_VertexX",                     &this->m_treeParameters.m_primary_VertexX);
    this->m_pOutputTree->Branch("primary_VertexY",                     &this->m_treeParameters.m_primary_VertexY);
    this->m_pOutputTree->Branch("primary_VertexZ",                     &this->m_treeParameters.m_primary_VertexZ);
    this->m_pOutputTree->Branch("primary_DirectionCosineX",            &this->m_treeParameters.m_primary_DirectionCosineX);
    this->m_pOutputTree->Branch("primary_DirectionCosineY",            &this->m_treeParameters.m_primary_DirectionCosineY);
    this->m_pOutputTree->Branch("primary_DirectionCosineZ",            &this->m_treeParameters.m_primary_DirectionCosineZ);
    this->m_pOutputTree->Branch("primary_MomentumX",                   &this->m_treeParameters.m_primary_MomentumX);
    this->m_pOutputTree->Branch("primary_MomentumY",                   &this->m_treeParameters.m_primary_MomentumY);
    this->m_pOutputTree->Branch("primary_MomentumZ",                   &this->m_treeParameters.m_primary_MomentumZ);
    this->m_pOutputTree->Branch("primary_ParticleType",                &this->m_treeParameters.m_primary_ParticleType);
    this->m_pOutputTree->Branch("primary_NumberOf3dHits",              &this->m_treeParameters.m_primary_NumberOf3dHits);
    this->m_pOutputTree->Branch("primary_NumberOfCollectionPlaneHits", &this->m_treeParameters.m_primary_NumberOfCollectionPlaneHits);
    this->m_pOutputTree->Branch("primary_IsShower",                    &this->m_treeParameters.m_primary_IsShower);
    this->m_pOutputTree->Branch("primary_NumberOfDownstreamParticles", &this->m_treeParameters.m_primary_NumberOfDownstreamParticles);
    this->m_pOutputTree->Branch("primary_mc_Energy",                   &this->m_treeParameters.m_primary_mc_Energy);
    this->m_pOutputTree->Branch("primary_mc_VertexX",                  &this->m_treeParameters.m_primary_mc_VertexX);
    this->m_pOutputTree->Branch("primary_mc_VertexY",                  &this->m_treeParameters.m_primary_mc_VertexY);
    this->m_pOutputTree->Branch("primary_mc_VertexZ",                  &this->m_treeParameters.m_primary_mc_VertexZ);
    this->m_pOutputTree->Branch("primary_mc_DirectionCosineX",         &this->m_treeParameters.m_primary_mc_DirectionCosineX);
    this->m_pOutputTree->Branch("primary_mc_DirectionCosineY",         &this->m_treeParameters.m_primary_mc_DirectionCosineY);
    this->m_pOutputTree->Branch("primary_mc_DirectionCosineZ",         &this->m_treeParameters.m_primary_mc_DirectionCosineZ);
    this->m_pOutputTree->Branch("primary_mc_MomentumX",                &this->m_treeParameters.m_primary_mc_MomentumX);
    this->m_pOutputTree->Branch("primary_mc_MomentumY",                &this->m_treeParameters.m_primary_mc_MomentumY);
    this->m_pOutputTree->Branch("primary_mc_MomentumZ",                &this->m_treeParameters.m_primary_mc_MomentumZ);
    this->m_pOutputTree->Branch("primary_mc_IsVertexFiducial",         &this->m_treeParameters.m_primary_mc_IsVertexFiducial);
    this->m_pOutputTree->Branch("primary_mc_IsContained",              &this->m_treeParameters.m_primary_mc_IsContained);
    this->m_pOutputTree->Branch("primary_mc_IsCorrectlyReconstructed", &this->m_treeParameters.m_primary_mc_IsCorrectlyReconstructed);
    this->m_pOutputTree->Branch("primary_mc_ParticleType",             &this->m_treeParameters.m_primary_mc_ParticleType);
    this->m_pOutputTree->Branch("primary_mc_IsShower",                 &this->m_treeParameters.m_primary_mc_IsShower);
    this->m_pOutputTree->Branch("primary_mc_PdgCode",                  &this->m_treeParameters.m_primary_mc_PdgCode);
    
    // Cosmic ray parameters.
    this->m_pOutputTree->Branch("cr_Number",                           &this->m_treeParameters.m_cr_Number);
    this->m_pOutputTree->Branch("cr_IsVertexFiducial",                 &this->m_treeParameters.m_cr_IsVertexFiducial);
    this->m_pOutputTree->Branch("cr_AreAllHitsFiducial",               &this->m_treeParameters.m_cr_AreAllHitsFiducial);
    this->m_pOutputTree->Branch("cr_HasMcInfo",                        &this->m_treeParameters.m_cr_HasMcInfo);
    this->m_pOutputTree->Branch("cr_Energy",                           &this->m_treeParameters.m_cr_Energy);
    this->m_pOutputTree->Branch("cr_EnergyFromChargeOnly",             &this->m_treeParameters.m_cr_EnergyFromChargeOnly);
    this->m_pOutputTree->Branch("cr_VertexX",                          &this->m_treeParameters.m_cr_VertexX);
    this->m_pOutputTree->Branch("cr_VertexY",                          &this->m_treeParameters.m_cr_VertexY);
    this->m_pOutputTree->Branch("cr_VertexZ",                          &this->m_treeParameters.m_cr_VertexZ);
    this->m_pOutputTree->Branch("cr_DirectionCosineX",                 &this->m_treeParameters.m_cr_DirectionCosineX);
    this->m_pOutputTree->Branch("cr_DirectionCosineY",                 &this->m_treeParameters.m_cr_DirectionCosineY);
    this->m_pOutputTree->Branch("cr_DirectionCosineZ",                 &this->m_treeParameters.m_cr_DirectionCosineZ);
    this->m_pOutputTree->Branch("cr_MomentumX",                        &this->m_treeParameters.m_cr_MomentumX);
    this->m_pOutputTree->Branch("cr_MomentumY",                        &this->m_treeParameters.m_cr_MomentumY);
    this->m_pOutputTree->Branch("cr_MomentumZ",                        &this->m_treeParameters.m_cr_MomentumZ);
    this->m_pOutputTree->Branch("cr_NumberOf3dHits",                   &this->m_treeParameters.m_cr_NumberOf3dHits);
    this->m_pOutputTree->Branch("cr_NumberOfCollectionPlaneHits",      &this->m_treeParameters.m_cr_NumberOfCollectionPlaneHits);
    this->m_pOutputTree->Branch("cr_mc_Energy",                        &this->m_treeParameters.m_cr_mc_Energy);
    this->m_pOutputTree->Branch("cr_mc_VertexX",                       &this->m_treeParameters.m_cr_mc_VertexX);
    this->m_pOutputTree->Branch("cr_mc_VertexY",                       &this->m_treeParameters.m_cr_mc_VertexY);
    this->m_pOutputTree->Branch("cr_mc_VertexZ",                       &this->m_treeParameters.m_cr_mc_VertexZ);
    this->m_pOutputTree->Branch("cr_mc_DirectionCosineX",              &this->m_treeParameters.m_cr_mc_DirectionCosineX);
    this->m_pOutputTree->Branch("cr_mc_DirectionCosineY",              &this->m_treeParameters.m_cr_mc_DirectionCosineY);
    this->m_pOutputTree->Branch("cr_mc_DirectionCosineZ",              &this->m_treeParameters.m_cr_mc_DirectionCosineZ);
    this->m_pOutputTree->Branch("cr_mc_MomentumX",                     &this->m_treeParameters.m_cr_mc_MomentumX);
    this->m_pOutputTree->Branch("cr_mc_MomentumY",                     &this->m_treeParameters.m_cr_mc_MomentumY);
    this->m_pOutputTree->Branch("cr_mc_MomentumZ",                     &this->m_treeParameters.m_cr_mc_MomentumZ);
    this->m_pOutputTree->Branch("cr_mc_IsVertexFiducial",              &this->m_treeParameters.m_cr_mc_IsVertexFiducial);
    this->m_pOutputTree->Branch("cr_mc_IsContained",                   &this->m_treeParameters.m_cr_mc_IsContained);
    this->m_pOutputTree->Branch("cr_mc_IsCorrectlyReconstructed",      &this->m_treeParameters.m_cr_mc_IsCorrectlyReconstructed);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TreeParameters::TreeParameters() noexcept :
    m_nu_Exists(false),
    m_nu_IsVertexFiducial(false),
    m_nu_AreAllHitsFiducial(false),
    m_nu_HasMcInfo(false),
    m_nu_Energy(0.f),
    m_nu_EnergyFromChargeOnly(0.f),
    m_nu_LongitudinalEnergy(0.f),
    m_nu_TransverseEnergy(0.f),
    m_nu_VertexX(0.f),
    m_nu_VertexY(0.f),
    m_nu_VertexZ(0.f),
    m_nu_DirectionCosineX(0.f),
    m_nu_DirectionCosineY(0.f),
    m_nu_DirectionCosineZ(0.f),
    m_nu_MomentumX(0.f),
    m_nu_MomentumY(0.f),
    m_nu_MomentumZ(0.f),
    m_nu_NumberOf3dHits(0U),
    m_nu_NumberOfCollectionPlaneHits(0U),
    m_nu_mc_Energy(0.f),
    m_nu_mc_LongitudinalEnergy(0.f),
    m_nu_mc_TransverseEnergy(0.f),
    m_nu_mc_VisibleLongitudinalEnergy(0.f),
    m_nu_mc_VisibleTransverseEnergy(0.f),
    m_nu_mc_VertexX(0.f),
    m_nu_mc_VertexY(0.f),
    m_nu_mc_VertexZ(0.f),
    m_nu_mc_DirectionCosineX(0.f),
    m_nu_mc_DirectionCosineY(0.f),
    m_nu_mc_DirectionCosineZ(0.f),
    m_nu_mc_MomentumX(0.f),
    m_nu_mc_MomentumY(0.f),
    m_nu_mc_MomentumZ(0.f),
    m_nu_mc_IsVertexFiducial(false),
    m_nu_mc_IsContained(false),
    m_nu_mc_IsCorrectlyReconstructed(false),
    m_nu_mc_InteractionType(0),
    m_nu_mc_IsChargedCurrent(false),
    m_nu_mc_VisibleEnergyFraction(0.f),
    m_nu_mc_PdgCode(0),
    m_primary_Number(0U),
    m_primary_IsVertexFiducial(),
    m_primary_AreAllHitsFiducial(),
    m_primary_HasMcInfo(),
    m_primary_Energy(),
    m_primary_EnergyFromChargeOnly(),
    m_primary_VertexX(),
    m_primary_VertexY(),
    m_primary_VertexZ(),
    m_primary_DirectionCosineX(),
    m_primary_DirectionCosineY(),
    m_primary_DirectionCosineZ(),
    m_primary_MomentumX(),
    m_primary_MomentumY(),
    m_primary_MomentumZ(),
    m_primary_ParticleType(),
    m_primary_NumberOf3dHits(),
    m_primary_NumberOfCollectionPlaneHits(),
    m_primary_IsShower(),
    m_primary_NumberOfDownstreamParticles(),
    m_primary_mc_Energy(),
    m_primary_mc_VertexX(),
    m_primary_mc_VertexY(),
    m_primary_mc_VertexZ(),
    m_primary_mc_DirectionCosineX(),
    m_primary_mc_DirectionCosineY(),
    m_primary_mc_DirectionCosineZ(),
    m_primary_mc_MomentumX(),
    m_primary_mc_MomentumY(),
    m_primary_mc_MomentumZ(),
    m_primary_mc_IsVertexFiducial(),
    m_primary_mc_IsContained(),
    m_primary_mc_IsCorrectlyReconstructed(),
    m_primary_mc_ParticleType(),
    m_primary_mc_IsShower(),
    m_primary_mc_PdgCode(),
    m_cr_Number(0U),
    m_cr_IsVertexFiducial(),
    m_cr_AreAllHitsFiducial(),
    m_cr_HasMcInfo(),
    m_cr_Energy(),
    m_cr_EnergyFromChargeOnly(),
    m_cr_VertexX(),
    m_cr_VertexY(),
    m_cr_VertexZ(),
    m_cr_DirectionCosineX(),
    m_cr_DirectionCosineY(),
    m_cr_DirectionCosineZ(),
    m_cr_MomentumX(),
    m_cr_MomentumY(),
    m_cr_MomentumZ(),
    m_cr_NumberOf3dHits(),
    m_cr_NumberOfCollectionPlaneHits(),
    m_cr_mc_Energy(),
    m_cr_mc_VertexX(),
    m_cr_mc_VertexY(),
    m_cr_mc_VertexZ(),
    m_cr_mc_DirectionCosineX(),
    m_cr_mc_DirectionCosineY(),
    m_cr_mc_DirectionCosineZ(),
    m_cr_mc_MomentumX(),
    m_cr_mc_MomentumY(),
    m_cr_mc_MomentumZ(),
    m_cr_mc_IsVertexFiducial(),
    m_cr_mc_IsContained(),
    m_cr_mc_IsCorrectlyReconstructed()
{
}

} // namespace lar_physics_content
