/**
 *  @file LArPhysicsContent/src/WriteAnalysisParticlesAlgorithm.cxx
 *
 *  @brief Implementation of the write AnalysisParticles algorithm class.
 *
 *  $Log: $
 */

#include "WriteAnalysisParticlesAlgorithm.h"

#include "LArAnalysisParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmHeaders.h"

/**
 *  @brief  ...
 */
#define PUSH_TREE_RECORD(treeMember, value) this->m_treeParameters.treeMember.push_back(value)

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
    
WriteAnalysisParticlesAlgorithm::WriteAnalysisParticlesAlgorithm() :
    m_pfoListName(),
    m_outputFile(),
    m_pOutputTFile(nullptr),
    m_pOutputTree(nullptr),
    m_verbose(false),
    m_treeParameters(),
    m_mcParticleListName(),
    m_caloHitListName(),
    m_fiducialCutLowXMargin(10.f),
    m_fiducialCutHighXMargin(10.f),
    m_fiducialCutLowYMargin(20.f),
    m_fiducialCutHighYMargin(20.f),
    m_fiducialCutLowZMargin(10.f),
    m_fiducialCutHighZMargin(10.f),
    m_minCoordinates(0.f, 0.f, 0.f),
    m_maxCoordinates(0.f, 0.f, 0.f),
    m_mcContainmentFractionLowerBound(0.9f),
    m_fiducialHitFractionLowerBound(0.9f),
    m_mcOnlyParticleContainmentCut(0.001f),
    m_mcOnlyParticleEnergyCut(0.0001f)
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
    // Get the input PFO list.
    const PfoList *pPfoList(nullptr);

    if ((PandoraContentApi::GetList(*this, m_pfoListName, pPfoList) != STATUS_CODE_SUCCESS) || !pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "WriteAnalysisParticlesAlgorithm: cannot find pfo list " << m_pfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }
    
    // Get the MC particle list if it exists.
    const MCParticleList *pMCParticleList(nullptr);
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);
    
    // Get the CaloHit list if it exists.
    const CaloHitList *pCaloHitList(nullptr);
    PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList);
    
    this->m_treeParameters = TreeParameters();

    MCPrimaryMap coveredMCPrimaries;
    
    for (const ParticleFlowObject * const pPfo : *pPfoList)
    {
        if (const LArAnalysisParticle *const pAnalysisParticle = dynamic_cast<const LArAnalysisParticle *>(pPfo))
        {
            if (pAnalysisParticle->HasMcInfo() && pAnalysisParticle->McMainMCParticle())
                coveredMCPrimaries.emplace(pAnalysisParticle->McMainMCParticle(), pAnalysisParticle);
        }
    }
    
    for (const ParticleFlowObject * const pPfo : *pPfoList)
    {
        if (const LArAnalysisParticle *const pAnalysisParticle = dynamic_cast<const LArAnalysisParticle *>(pPfo))
        {
            if (LArAnalysisParticleHelper::IsNeutrino(pPfo))  
            {
                if (this->m_treeParameters.m_nu_WasReconstructed)
                {
                    std::cout << "WriteAnalysisParticlesAlgorithm: multiple neutrinos found - only recording one" << std::endl;
                    continue;
                }

                this->PopulateNeutrinoParameters(*pAnalysisParticle, pMCParticleList, pCaloHitList);
                this->m_treeParameters.m_nu_WasReconstructed = true;
            }
            
            else if (LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(pPfo))
            {
                this->AddPrimaryDaughterRecord(*pAnalysisParticle, coveredMCPrimaries);
                ++this->m_treeParameters.m_primary_Number;
            }
            
            else if (LArAnalysisParticleHelper::IsCosmicRay(pPfo))
            {
                this->AddCosmicRayRecord(*pAnalysisParticle, coveredMCPrimaries);
                ++this->m_treeParameters.m_cr_Number;
            }
            
            else
            {
                std::cout << "WriteAnalysisParticlesAlgorithm: analysis particle was not a cosmic ray, neutrino or primary daughter" << std::endl;
                continue;
            }
        }
    }
    
    if (pMCParticleList)
    {
        MCParticleSet mcPrimarySet;
            
        for (const MCParticle *const pMCParticle : *pMCParticleList)
        {
            try
            {
                if ((pMCParticle->GetParentList().size() == 0) && LArMCParticleHelper::IsNeutrino(pMCParticle))
                    mcPrimarySet.insert(pMCParticle);
                
                else
                    mcPrimarySet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle));
            }
            
            catch (...)
            {
                continue;
            }
        }
        
        for (const MCParticle *const pMCPrimary : mcPrimarySet)
        {
            if (coveredMCPrimaries.find(pMCPrimary) != coveredMCPrimaries.end())
                continue;
   
            if (!LArMCParticleHelper::IsNeutrino(pMCPrimary) && !LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) &&
                (pMCPrimary->GetParticleId() != MU_MINUS))
            {
                continue;
            }
            
            float mcEnergy = 0.f, mcKineticEnergy = 0.f, mcMass = 0.f, mcContainmentFraction = 0.f;
            LArAnalysisParticle::TypeTree mcTypeTree;
            LArAnalysisParticle::TYPE mcType(LArAnalysisParticle::TYPE::UNKNOWN);
            CartesianVector mcVertexPosition(0.f, 0.f, 0.f), mcMomentum(0.f, 0.f, 0.f);
            int mcPdgCode(0);
            
            if (!LArAnalysisParticleHelper::GetMcInformation(pMCPrimary, mcEnergy, mcKineticEnergy, mcMass, mcTypeTree, mcType, mcVertexPosition, mcMomentum, mcPdgCode, 
                mcContainmentFraction, this->m_minCoordinates, this->m_maxCoordinates))
            {
                std::cout << "WriteAnalysisParticlesAlgorithm: failed to get MC information for non-reconstructed MC particle" << std::endl;
                continue;
            }
            
            if (mcEnergy < m_mcOnlyParticleEnergyCut)
                continue;
                
            if (mcContainmentFraction < m_mcOnlyParticleContainmentCut)
                continue;

            const bool mcIsVertexFiducial = LArAnalysisParticleHelper::IsPointFiducial(mcVertexPosition, this->m_minCoordinates,
                this->m_maxCoordinates);
                
            CartesianVector mcDirectionCosines(0.f, 0.f, 0.f);
            
            if (mcMomentum.GetMagnitude() < std::numeric_limits<float>::epsilon())
                std::cout << "AnalysisAlgorithm: could not get direction from MC momentum as it was too small" << std::endl;
            
            else
                mcDirectionCosines = mcMomentum.GetUnitVector();
                
            const bool mcIsShower = (mcType == LArAnalysisParticle::TYPE::SHOWER);

            if (LArMCParticleHelper::IsNeutrino(pMCPrimary))
            {
                if (this->m_treeParameters.m_nu_HasMcInfo)
                {
                    std::cout << "WriteAnalysisParticlesAlgorithm: found neutrino that was not covered by MC particles but neutrino MC properties already exist" << std::endl;
                    continue;
                }
                
                this->m_treeParameters.m_nu_HasMcInfo = true;
                
                this->PopulateNeutrinoMcParameters(pMCPrimary, mcEnergy, mcVertexPosition, mcDirectionCosines, mcMomentum, 
                    mcIsVertexFiducial, mcContainmentFraction, mcPdgCode, pMCParticleList, pCaloHitList, 0.f, 0.f, 0.f, 0.f, mcTypeTree);
            }
                
            else if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary))
            {
                ++this->m_treeParameters.m_primary_Number;
                this->AddMcOnlyPrimaryDaughterRecord(pMCPrimary, mcEnergy, mcKineticEnergy, mcMass, mcVertexPosition, mcDirectionCosines, mcMomentum, 
                    mcIsVertexFiducial, mcContainmentFraction, mcType, mcIsShower, mcPdgCode, mcTypeTree);
            }
            
            else if (pMCPrimary->GetParticleId() == MU_MINUS)
            {
                ++this->m_treeParameters.m_cr_Number;
                this->AddMcOnlyCosmicRayRecord(pMCPrimary, mcEnergy, mcKineticEnergy, mcMass, mcVertexPosition, mcDirectionCosines, mcMomentum,
                    mcIsVertexFiducial, mcContainmentFraction, mcType, mcIsShower, mcPdgCode, mcTypeTree);
            }
        }
    }
    
    if (this->m_verbose)
        this->DumpTree();
    
    this->m_pOutputTree->Fill();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoParameters(const LArAnalysisParticle &neutrinoAnalysisParticle, 
    const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList) const
{
    // m_nu_WasReconstructed is dealt with by the calling method.
    
    this->m_treeParameters.m_nu_IsVertexFiducial            = neutrinoAnalysisParticle.IsVertexFiducial();
    this->m_treeParameters.m_nu_IsContained                 = (neutrinoAnalysisParticle.FiducialHitFraction() >= this->m_fiducialHitFractionLowerBound);
    this->m_treeParameters.m_nu_FiducialHitFraction         = neutrinoAnalysisParticle.FiducialHitFraction();
    this->m_treeParameters.m_nu_HasMcInfo                   = neutrinoAnalysisParticle.HasMcInfo();
    this->m_treeParameters.m_nu_VisibleEnergy               = neutrinoAnalysisParticle.KineticEnergy();
    this->m_treeParameters.m_nu_VisibleEnergyFracFromRange                  = neutrinoAnalysisParticle.KineticEnergyFromRangeFraction();
    this->m_treeParameters.m_nu_VisibleEnergyFracFromCorrectedTrackCharge   = neutrinoAnalysisParticle.KineticEnergyFromCorrectedTrackChargeFraction();
    this->m_treeParameters.m_nu_VisibleEnergyFracFromUncorrectedTrackCharge = neutrinoAnalysisParticle.KineticEnergyFromUncorrectedTrackChargeFraction();
    this->m_treeParameters.m_nu_VisibleEnergyFracFromShowerCharge           = neutrinoAnalysisParticle.KineticEnergyFromShowerChargeFraction();
    this->m_treeParameters.m_nu_VertexX                     = neutrinoAnalysisParticle.VertexPosition().GetX();
    this->m_treeParameters.m_nu_VertexY                     = neutrinoAnalysisParticle.VertexPosition().GetY();
    this->m_treeParameters.m_nu_VertexZ                     = neutrinoAnalysisParticle.VertexPosition().GetZ();
    this->m_treeParameters.m_nu_DirectionCosineX            = neutrinoAnalysisParticle.DirectionCosines().GetX();
    this->m_treeParameters.m_nu_DirectionCosineY            = neutrinoAnalysisParticle.DirectionCosines().GetY();
    this->m_treeParameters.m_nu_DirectionCosineZ            = neutrinoAnalysisParticle.DirectionCosines().GetZ();
    this->m_treeParameters.m_nu_TypeTree                    = LArAnalysisParticle::TypeTreeAsString(neutrinoAnalysisParticle.GetTypeTree());
    this->m_treeParameters.m_nu_NumberOf3dHits              = neutrinoAnalysisParticle.NumberOf3dHits();
    this->m_treeParameters.m_nu_NumberOfCollectionPlaneHits = neutrinoAnalysisParticle.NumberOfCollectionPlaneHits();
    this->m_treeParameters.m_nu_NumberOfDownstreamParticles = neutrinoAnalysisParticle.NumberOfDownstreamParticles();
    
    CartesianVector visiblePseudoMomentum(0.f, 0.f, 0.f); // the summed reconstructed directional KEs of the primaries
    
    for (const ParticleFlowObject *const pDaughterPfo : neutrinoAnalysisParticle.GetDaughterPfoList())
    {
        if (const LArAnalysisParticle *const pDaughterAnalysisParticle = dynamic_cast<const LArAnalysisParticle *>(pDaughterPfo))
            visiblePseudoMomentum += pDaughterAnalysisParticle->DirectionCosines() * pDaughterAnalysisParticle->KineticEnergy();
    }
    
    const CartesianVector zAxis(0.f, 0.f, 1.f);
    this->m_treeParameters.m_nu_VisibleTransverseEnergy   = visiblePseudoMomentum.GetDotProduct(zAxis);
    this->m_treeParameters.m_nu_VisibleLongitudinalEnergy = visiblePseudoMomentum.GetCrossProduct(zAxis).GetMagnitude();
    
    if (neutrinoAnalysisParticle.HasMcInfo())
    {
        this->PopulateNeutrinoMcParameters(neutrinoAnalysisParticle.McMainMCParticle(), neutrinoAnalysisParticle.McEnergy(),
            neutrinoAnalysisParticle.McVertexPosition(), neutrinoAnalysisParticle.McDirectionCosines(), 
            neutrinoAnalysisParticle.McMomentum(), neutrinoAnalysisParticle.McIsVertexFiducial(), neutrinoAnalysisParticle.McContainmentFraction(), 
            neutrinoAnalysisParticle.McPdgCode(), pMCParticleList, pCaloHitList, neutrinoAnalysisParticle.McHitPurity(), 
            neutrinoAnalysisParticle.McHitCompleteness(), neutrinoAnalysisParticle.McCollectionPlaneHitPurity(), 
            neutrinoAnalysisParticle.McCollectionPlaneHitCompleteness(), neutrinoAnalysisParticle.GetMcTypeTree());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoMcParameters(const MCParticle *const pMainMcParticle, const float mcEnergy,
    const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum,
    const bool mcIsVertexFiducial, const float mcContainmentFraction, const int mcPdgCode, const MCParticleList *const pMCParticleList, 
    const CaloHitList *const pCaloHitList, const float mcHitPurity, const float mcHitCompleteness, const float mcCollectionPlaneHitPurity,
    const float mcCollectionPlaneHitCompleteness, const LArAnalysisParticle::TypeTree mcTypeTree) const
{
    if (pMainMcParticle)
        this->m_treeParameters.m_nu_mc_McParticleUid = reinterpret_cast<std::uint64_t>(pMainMcParticle->GetUid());
        
    else // it can stay at its default value
        std::cout << "WriteAnalysisParticlesAlgorithm: could not populate MC particle UID" << std::endl;
    
    this->m_treeParameters.m_nu_mc_Energy           = mcEnergy;
    this->m_treeParameters.m_nu_mc_VertexX          = mcVertexPosition.GetX();
    this->m_treeParameters.m_nu_mc_VertexY          = mcVertexPosition.GetY();
    this->m_treeParameters.m_nu_mc_VertexZ          = mcVertexPosition.GetZ();
    this->m_treeParameters.m_nu_mc_DirectionCosineX = mcDirectionCosines.GetX();
    this->m_treeParameters.m_nu_mc_DirectionCosineY = mcDirectionCosines.GetY();
    this->m_treeParameters.m_nu_mc_DirectionCosineZ = mcDirectionCosines.GetZ();
    this->m_treeParameters.m_nu_mc_Momentum         = mcMomentum.GetMagnitude();
    this->m_treeParameters.m_nu_mc_MomentumX        = mcMomentum.GetX();
    this->m_treeParameters.m_nu_mc_MomentumY        = mcMomentum.GetY();
    this->m_treeParameters.m_nu_mc_MomentumZ        = mcMomentum.GetZ();
    this->m_treeParameters.m_nu_mc_IsVertexFiducial = mcIsVertexFiducial;
    this->m_treeParameters.m_nu_mc_IsContained      = (mcContainmentFraction >= m_mcContainmentFractionLowerBound);
    this->m_treeParameters.m_nu_mc_ContainmentFraction = mcContainmentFraction;
    this->m_treeParameters.m_nu_mc_TypeTree         = LArAnalysisParticle::TypeTreeAsString(mcTypeTree);
    this->m_treeParameters.m_nu_mc_PdgCode          = mcPdgCode;
    this->m_treeParameters.m_nu_mc_HitPurity        = mcHitPurity;
    this->m_treeParameters.m_nu_mc_HitCompleteness  = mcHitCompleteness;
    this->m_treeParameters.m_nu_mc_CollectionPlaneHitPurity        = mcCollectionPlaneHitPurity;
    this->m_treeParameters.m_nu_mc_CollectionPlaneHitCompleteness  = mcCollectionPlaneHitCompleteness;
    
    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: no valid MC particle list name provided but neutrino analysis particle has MC info" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    
    if (!pCaloHitList || pCaloHitList->empty())
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: no valid CaloHit list name provided but neutrino analysis particle has MC info" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    
    MCParticleSet mcPrimarySet;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        try
        {
            const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
            
            if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pPrimaryMCParticle))
                mcPrimarySet.insert(pPrimaryMCParticle);
        }
        
        catch (...)
        {
            continue;
        }
    }

    CartesianVector visibleMomentum(0.f, 0.f, 0.f);
    float visibleEnergy(0.f);
    
    for (const MCParticle *const pMCPrimary : mcPrimarySet)
    {
        const float primaryVisibleEnergy = pMCPrimary->GetEnergy() - PdgTable::GetParticleMass(pMCPrimary->GetParticleId());
        visibleEnergy += primaryVisibleEnergy;
        
        if (pMCPrimary->GetMomentum().GetMagnitude() > std::numeric_limits<float>::epsilon())
            visibleMomentum += pMCPrimary->GetMomentum().GetUnitVector() * primaryVisibleEnergy;
    }
    
    const LArInteractionTypeHelper::InteractionType interactionType = this->GetInteractionType(pMCParticleList, pCaloHitList);
    
    // To align with the non-MC definition, we define the longitudinal/transverse components with respect to the z-direction.
    const CartesianVector zDirectionVector(0.f, 0.f, 1.f);

    this->m_treeParameters.m_nu_mc_LongitudinalEnergy        = mcMomentum.GetDotProduct(zDirectionVector);
    this->m_treeParameters.m_nu_mc_TransverseEnergy          = mcMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    this->m_treeParameters.m_nu_mc_VisibleEnergy             = visibleEnergy;
    this->m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy = visibleMomentum.GetDotProduct(zDirectionVector);
    this->m_treeParameters.m_nu_mc_VisibleTransverseEnergy   = visibleMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    this->m_treeParameters.m_nu_mc_InteractionType           = LArInteractionTypeHelper::ToString(interactionType);
    this->m_treeParameters.m_nu_mc_IsChargedCurrent          = this->IsChargedCurrent(interactionType);
    
    float visibleEnergyFraction(0.f);
    
    if (mcEnergy > std::numeric_limits<float>::epsilon())
        visibleEnergyFraction = visibleEnergy / mcEnergy;

    this->m_treeParameters.m_nu_mc_VisibleEnergyFraction = visibleEnergyFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionType WriteAnalysisParticlesAlgorithm::GetInteractionType(const MCParticleList *const pMCParticleList,
    const CaloHitList *const pCaloHitList) const
{
    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);
    
    if (mcNeutrinoVector.empty())
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: could not get interaction type because there were no reconstructed neutrinos" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    
    if (mcNeutrinoVector.size() > 1)
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: could not get interaction type because there was more than one reconstructed neutrino" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    
    const MCParticle *const pMCNeutrino = mcNeutrinoVector.front();

    const LArMCParticleHelper::PrimaryParameters parameters;
    
    if (const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle *>(pMCNeutrino))
        return LArInteractionTypeHelper::GetInteractionType(pLArMCNeutrino, pMCParticleList, pCaloHitList, parameters);

    std::cout << "WriteAnalysisParticlesAlgorithm: could not cast neutrino to LArMCParticle type" << std::endl;
    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddPrimaryDaughterRecord(const LArAnalysisParticle &primaryAnalysisParticle,
    const MCPrimaryMap &coveredMCPrimaries) const
{
    // m_primary_Number is dealt with by the calling method.
    
    PUSH_TREE_RECORD(m_primary_WasReconstructed,            true);
    PUSH_TREE_RECORD(m_primary_IsVertexFiducial,            primaryAnalysisParticle.IsVertexFiducial());
    PUSH_TREE_RECORD(m_primary_IsContained,                (primaryAnalysisParticle.FiducialHitFraction() >= this->m_fiducialHitFractionLowerBound));
    PUSH_TREE_RECORD(m_primary_FiducialHitFraction,         primaryAnalysisParticle.FiducialHitFraction());
    PUSH_TREE_RECORD(m_primary_HasMcInfo,                   primaryAnalysisParticle.HasMcInfo());
    PUSH_TREE_RECORD(m_primary_KineticEnergy,               primaryAnalysisParticle.KineticEnergy());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromRange,                  primaryAnalysisParticle.KineticEnergyFromRangeFraction());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromCorrectedTrackCharge,   primaryAnalysisParticle.KineticEnergyFromCorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromUncorrectedTrackCharge, primaryAnalysisParticle.KineticEnergyFromUncorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromShowerCharge,           primaryAnalysisParticle.KineticEnergyFromShowerChargeFraction());
    PUSH_TREE_RECORD(m_primary_VertexX,                     primaryAnalysisParticle.VertexPosition().GetX());
    PUSH_TREE_RECORD(m_primary_VertexY,                     primaryAnalysisParticle.VertexPosition().GetY());
    PUSH_TREE_RECORD(m_primary_VertexZ,                     primaryAnalysisParticle.VertexPosition().GetZ());
    PUSH_TREE_RECORD(m_primary_DirectionCosineX,            primaryAnalysisParticle.DirectionCosines().GetX());
    PUSH_TREE_RECORD(m_primary_DirectionCosineY,            primaryAnalysisParticle.DirectionCosines().GetY());
    PUSH_TREE_RECORD(m_primary_DirectionCosineZ,            primaryAnalysisParticle.DirectionCosines().GetZ());
    PUSH_TREE_RECORD(m_primary_TypeTree,                    LArAnalysisParticle::TypeTreeAsString(primaryAnalysisParticle.GetTypeTree()));
    PUSH_TREE_RECORD(m_primary_IsShower,                    primaryAnalysisParticle.IsShower());
    PUSH_TREE_RECORD(m_primary_IsTrack,                    !primaryAnalysisParticle.IsShower());
    PUSH_TREE_RECORD(m_primary_IsProton,                   (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::PROTON));
    PUSH_TREE_RECORD(m_primary_IsPionOrMuon,               (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::PION_MUON) || (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_primary_NumberOf3dHits,              primaryAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_primary_NumberOfCollectionPlaneHits, primaryAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_primary_NumberOfDownstreamParticles, primaryAnalysisParticle.NumberOfDownstreamParticles());
    
    if (primaryAnalysisParticle.HasMcInfo())
    {
        bool particleSplitByReco(false);
        
        if (primaryAnalysisParticle.McMainMCParticle())
        {
            particleSplitByReco = (coveredMCPrimaries.count(primaryAnalysisParticle.McMainMCParticle()) > 1U);
            PUSH_TREE_RECORD(m_primary_mc_McParticleUid, reinterpret_cast<std::uint64_t>(primaryAnalysisParticle.McMainMCParticle()->GetUid()));
        }
        
        else
        {
            std::cout << "WriteAnalysisParticlesAlgorithm: primary neutrino daughter had MC info but no associated main MC particle" << std::endl;
            PUSH_TREE_RECORD(m_primary_mc_McParticleUid, 0ULL);
        }
        
        PUSH_TREE_RECORD(m_primary_mc_IsParticleSplitByReco,    particleSplitByReco);
        PUSH_TREE_RECORD(m_primary_mc_Energy,                   primaryAnalysisParticle.McEnergy());
        PUSH_TREE_RECORD(m_primary_mc_KineticEnergy,            primaryAnalysisParticle.McKineticEnergy());
        PUSH_TREE_RECORD(m_primary_mc_Mass,                     primaryAnalysisParticle.McMass());
        PUSH_TREE_RECORD(m_primary_mc_VertexX,                  primaryAnalysisParticle.McVertexPosition().GetX());
        PUSH_TREE_RECORD(m_primary_mc_VertexY,                  primaryAnalysisParticle.McVertexPosition().GetY());
        PUSH_TREE_RECORD(m_primary_mc_VertexZ,                  primaryAnalysisParticle.McVertexPosition().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,         primaryAnalysisParticle.McDirectionCosines().GetX());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,         primaryAnalysisParticle.McDirectionCosines().GetY());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,         primaryAnalysisParticle.McDirectionCosines().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_Momentum,                 primaryAnalysisParticle.McMomentum().GetMagnitude());
        PUSH_TREE_RECORD(m_primary_mc_MomentumX,                primaryAnalysisParticle.McMomentum().GetX());
        PUSH_TREE_RECORD(m_primary_mc_MomentumY,                primaryAnalysisParticle.McMomentum().GetY());
        PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                primaryAnalysisParticle.McMomentum().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,         primaryAnalysisParticle.McIsVertexFiducial());
        PUSH_TREE_RECORD(m_primary_mc_IsContained,             (primaryAnalysisParticle.McContainmentFraction() >= this->m_mcContainmentFractionLowerBound));
        PUSH_TREE_RECORD(m_primary_mc_ContainmentFraction,      primaryAnalysisParticle.McContainmentFraction());
        PUSH_TREE_RECORD(m_primary_mc_TypeTree,                 LArAnalysisParticle::TypeTreeAsString(primaryAnalysisParticle.GetMcTypeTree()));
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                 primaryAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_primary_mc_IsTrack,                 !primaryAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_primary_mc_IsProton,                (primaryAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PROTON));
        PUSH_TREE_RECORD(m_primary_mc_IsPionOrMuon,            (primaryAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PION_MUON));
        PUSH_TREE_RECORD(m_primary_mc_IsCosmicRay,             (primaryAnalysisParticle.McType() == LArAnalysisParticle::TYPE::COSMIC_RAY));
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  primaryAnalysisParticle.McPdgCode());
        PUSH_TREE_RECORD(m_primary_mc_HitPurity,                primaryAnalysisParticle.McHitPurity());
        PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,          primaryAnalysisParticle.McHitCompleteness());
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitPurity,       primaryAnalysisParticle.McCollectionPlaneHitPurity());
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitCompleteness, primaryAnalysisParticle.McCollectionPlaneHitCompleteness());
    }
    
    else
    {
        PUSH_TREE_RECORD(m_primary_mc_IsParticleSplitByReco,    false);
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid,            0ULL);
        PUSH_TREE_RECORD(m_primary_mc_Energy,                   0.f);
        PUSH_TREE_RECORD(m_primary_mc_KineticEnergy,            0.f);
        PUSH_TREE_RECORD(m_primary_mc_Mass,                     0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexX,                  0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexY,                  0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexZ,                  0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,         0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,         0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,         0.f);
        PUSH_TREE_RECORD(m_primary_mc_Momentum,                 0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumX,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumY,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,         false);
        PUSH_TREE_RECORD(m_primary_mc_IsContained,              false);
        PUSH_TREE_RECORD(m_primary_mc_ContainmentFraction,      0.f);
        PUSH_TREE_RECORD(m_primary_mc_TypeTree,                 "-");
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                 false);
        PUSH_TREE_RECORD(m_primary_mc_IsTrack,                  false);
        PUSH_TREE_RECORD(m_primary_mc_IsProton,                 false);
        PUSH_TREE_RECORD(m_primary_mc_IsPionOrMuon,             false);
        PUSH_TREE_RECORD(m_primary_mc_IsCosmicRay,              false);
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  0);
        PUSH_TREE_RECORD(m_primary_mc_HitPurity,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,          0.f);
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitPurity,       0.f);
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitCompleteness, 0.f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyPrimaryDaughterRecord(const MCParticle *const pMainMcParticle, const float mcEnergy, const float mcKineticEnergy, const float mcMass,
    const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum,
    const bool mcIsVertexFiducial, const float mcContainmentFraction, const LArAnalysisParticle::TYPE mcType, const bool mcIsShower, 
    const int mcPdgCode, const LArAnalysisParticle::TypeTree mcTypeTree) const
{
    PUSH_TREE_RECORD(m_primary_WasReconstructed,            false);
    PUSH_TREE_RECORD(m_primary_IsVertexFiducial,            false);
    PUSH_TREE_RECORD(m_primary_IsContained,                 false);
    PUSH_TREE_RECORD(m_primary_FiducialHitFraction,         0.f);
    PUSH_TREE_RECORD(m_primary_HasMcInfo,                   true);
    PUSH_TREE_RECORD(m_primary_KineticEnergy,               0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromRange,                  0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromCorrectedTrackCharge,   0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromUncorrectedTrackCharge, 0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromShowerCharge,           0.f);
    PUSH_TREE_RECORD(m_primary_VertexX,                     0.f);
    PUSH_TREE_RECORD(m_primary_VertexY,                     0.f);
    PUSH_TREE_RECORD(m_primary_VertexZ,                     0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineX,            0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineY,            0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineZ,            0.f);
    PUSH_TREE_RECORD(m_primary_TypeTree,                    "-");
    PUSH_TREE_RECORD(m_primary_IsShower,                    false);
    PUSH_TREE_RECORD(m_primary_IsTrack,                     false);
    PUSH_TREE_RECORD(m_primary_IsProton,                    false);
    PUSH_TREE_RECORD(m_primary_IsPionOrMuon,                false);
    PUSH_TREE_RECORD(m_primary_NumberOf3dHits,              0U);
    PUSH_TREE_RECORD(m_primary_NumberOfCollectionPlaneHits, 0U);
    PUSH_TREE_RECORD(m_primary_NumberOfDownstreamParticles, 0U);
    
    if (pMainMcParticle)
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, reinterpret_cast<std::uint64_t>(pMainMcParticle->GetUid()));
    
    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: primary neutrino daughter had MC info but no associated main MC particle" << std::endl;
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, 0ULL);
    }
    
    PUSH_TREE_RECORD(m_primary_mc_IsParticleSplitByReco,    false); // not covered by any reconstructed particle
    PUSH_TREE_RECORD(m_primary_mc_Energy,                   mcEnergy);
    PUSH_TREE_RECORD(m_primary_mc_KineticEnergy,            mcKineticEnergy);
    PUSH_TREE_RECORD(m_primary_mc_Mass,                     mcMass);
    PUSH_TREE_RECORD(m_primary_mc_VertexX,                  mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_primary_mc_VertexY,                  mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_primary_mc_VertexZ,                  mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,         mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,         mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,         mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_Momentum,                 mcMomentum.GetMagnitude());
    PUSH_TREE_RECORD(m_primary_mc_MomentumX,                mcMomentum.GetX());
    PUSH_TREE_RECORD(m_primary_mc_MomentumY,                mcMomentum.GetY());
    PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,         mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_primary_mc_IsContained,             (mcContainmentFraction >= this->m_mcContainmentFractionLowerBound));
    PUSH_TREE_RECORD(m_primary_mc_ContainmentFraction,      mcContainmentFraction);
    PUSH_TREE_RECORD(m_primary_mc_TypeTree,                 LArAnalysisParticle::TypeTreeAsString(mcTypeTree));
    PUSH_TREE_RECORD(m_primary_mc_IsShower,                 mcIsShower);
    PUSH_TREE_RECORD(m_primary_mc_IsTrack,                 !mcIsShower);
    PUSH_TREE_RECORD(m_primary_mc_IsProton,                (mcType == LArAnalysisParticle::TYPE::PROTON));
    PUSH_TREE_RECORD(m_primary_mc_IsPionOrMuon,            (mcType == LArAnalysisParticle::TYPE::PION_MUON));
    PUSH_TREE_RECORD(m_primary_mc_IsCosmicRay,             (mcType == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  mcPdgCode);
    PUSH_TREE_RECORD(m_primary_mc_HitPurity,                0.f); // undefined
    PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,          0.f); // undefined
    PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitPurity,       0.f); // undefined
    PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitCompleteness, 0.f); // undefined
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddCosmicRayRecord(const LArAnalysisParticle &cosmicRayAnalysisParticle,
    const MCPrimaryMap &coveredMCPrimaries) const
{
    // m_cr_Number is dealt with by the calling method.
    
    PUSH_TREE_RECORD(m_cr_WasReconstructed,            true);
    PUSH_TREE_RECORD(m_cr_IsVertexFiducial,            cosmicRayAnalysisParticle.IsVertexFiducial());
    PUSH_TREE_RECORD(m_cr_IsContained,                (cosmicRayAnalysisParticle.FiducialHitFraction() >= this->m_fiducialHitFractionLowerBound));
    PUSH_TREE_RECORD(m_cr_FiducialHitFraction,         cosmicRayAnalysisParticle.FiducialHitFraction());
    PUSH_TREE_RECORD(m_cr_HasMcInfo,                   cosmicRayAnalysisParticle.HasMcInfo());
    PUSH_TREE_RECORD(m_cr_KineticEnergy,               cosmicRayAnalysisParticle.KineticEnergy());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromRange,                  cosmicRayAnalysisParticle.KineticEnergyFromRangeFraction());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromCorrectedTrackCharge,   cosmicRayAnalysisParticle.KineticEnergyFromCorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromUncorrectedTrackCharge, cosmicRayAnalysisParticle.KineticEnergyFromUncorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromShowerCharge,           cosmicRayAnalysisParticle.KineticEnergyFromShowerChargeFraction());
    PUSH_TREE_RECORD(m_cr_VertexX,                     cosmicRayAnalysisParticle.VertexPosition().GetX());
    PUSH_TREE_RECORD(m_cr_VertexY,                     cosmicRayAnalysisParticle.VertexPosition().GetY());
    PUSH_TREE_RECORD(m_cr_VertexZ,                     cosmicRayAnalysisParticle.VertexPosition().GetZ());
    PUSH_TREE_RECORD(m_cr_DirectionCosineX,            cosmicRayAnalysisParticle.DirectionCosines().GetX());
    PUSH_TREE_RECORD(m_cr_DirectionCosineY,            cosmicRayAnalysisParticle.DirectionCosines().GetY());
    PUSH_TREE_RECORD(m_cr_DirectionCosineZ,            cosmicRayAnalysisParticle.DirectionCosines().GetZ());
    PUSH_TREE_RECORD(m_cr_TypeTree,                    LArAnalysisParticle::TypeTreeAsString(cosmicRayAnalysisParticle.GetTypeTree()));
    PUSH_TREE_RECORD(m_cr_NumberOf3dHits,              cosmicRayAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_cr_NumberOfCollectionPlaneHits, cosmicRayAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_cr_NumberOfDownstreamParticles, cosmicRayAnalysisParticle.NumberOfDownstreamParticles());
    
    if (cosmicRayAnalysisParticle.HasMcInfo())
    {
        bool particleSplitByReco(false);
        
        if (cosmicRayAnalysisParticle.McMainMCParticle())
        {
            particleSplitByReco = (coveredMCPrimaries.count(cosmicRayAnalysisParticle.McMainMCParticle()) > 1U);
            PUSH_TREE_RECORD(m_cr_mc_McParticleUid, reinterpret_cast<std::uint64_t>(cosmicRayAnalysisParticle.McMainMCParticle()->GetUid()));
        }
        
        else
        {
            std::cout << "WriteAnalysisParticlesAlgorithm: cosmic ray had MC info but no associated main MC particle" << std::endl;
            PUSH_TREE_RECORD(m_cr_mc_McParticleUid, 0ULL);
        }
        
        PUSH_TREE_RECORD(m_cr_mc_IsParticleSplitByReco,    particleSplitByReco);
        PUSH_TREE_RECORD(m_cr_mc_Energy,                   cosmicRayAnalysisParticle.McEnergy());
        PUSH_TREE_RECORD(m_cr_mc_KineticEnergy,            cosmicRayAnalysisParticle.McKineticEnergy());
        PUSH_TREE_RECORD(m_cr_mc_Mass,                     cosmicRayAnalysisParticle.McMass());
        PUSH_TREE_RECORD(m_cr_mc_VertexX,                  cosmicRayAnalysisParticle.McVertexPosition().GetX());
        PUSH_TREE_RECORD(m_cr_mc_VertexY,                  cosmicRayAnalysisParticle.McVertexPosition().GetY());
        PUSH_TREE_RECORD(m_cr_mc_VertexZ,                  cosmicRayAnalysisParticle.McVertexPosition().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,         cosmicRayAnalysisParticle.McDirectionCosines().GetX());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,         cosmicRayAnalysisParticle.McDirectionCosines().GetY());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,         cosmicRayAnalysisParticle.McDirectionCosines().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_Momentum,                 cosmicRayAnalysisParticle.McMomentum().GetMagnitude());
        PUSH_TREE_RECORD(m_cr_mc_MomentumX,                cosmicRayAnalysisParticle.McMomentum().GetX());
        PUSH_TREE_RECORD(m_cr_mc_MomentumY,                cosmicRayAnalysisParticle.McMomentum().GetY());
        PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                cosmicRayAnalysisParticle.McMomentum().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,         cosmicRayAnalysisParticle.McIsVertexFiducial());
        PUSH_TREE_RECORD(m_cr_mc_IsContained,             (cosmicRayAnalysisParticle.McContainmentFraction() >= this->m_mcContainmentFractionLowerBound));
        PUSH_TREE_RECORD(m_cr_mc_ContainmentFraction,      cosmicRayAnalysisParticle.McContainmentFraction());
        PUSH_TREE_RECORD(m_cr_mc_TypeTree,                  LArAnalysisParticle::TypeTreeAsString(cosmicRayAnalysisParticle.GetMcTypeTree()));
        PUSH_TREE_RECORD(m_cr_mc_IsShower,                 cosmicRayAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_cr_mc_IsTrack,                 !cosmicRayAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_cr_mc_IsProton,                (cosmicRayAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PROTON));
        PUSH_TREE_RECORD(m_cr_mc_IsPionOrMuon,            (cosmicRayAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PION_MUON));
        PUSH_TREE_RECORD(m_cr_mc_IsCosmicRay,             (cosmicRayAnalysisParticle.McType() == LArAnalysisParticle::TYPE::COSMIC_RAY));
        PUSH_TREE_RECORD(m_cr_mc_PdgCode,                  cosmicRayAnalysisParticle.McPdgCode());
        PUSH_TREE_RECORD(m_cr_mc_HitPurity,                cosmicRayAnalysisParticle.McHitPurity());
        PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,          cosmicRayAnalysisParticle.McHitCompleteness());
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitPurity,       cosmicRayAnalysisParticle.McCollectionPlaneHitPurity());
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitCompleteness, cosmicRayAnalysisParticle.McCollectionPlaneHitCompleteness());
    }
    
    else
    {
        PUSH_TREE_RECORD(m_cr_mc_IsParticleSplitByReco,    false);
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid,            0ULL);
        PUSH_TREE_RECORD(m_cr_mc_Energy,                   0.f);
        PUSH_TREE_RECORD(m_cr_mc_KineticEnergy,            0.f);
        PUSH_TREE_RECORD(m_cr_mc_Mass,                     0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexX,                  0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexY,                  0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexZ,                  0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,         0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,         0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,         0.f);
        PUSH_TREE_RECORD(m_cr_mc_Momentum,                 0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumX,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumY,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,         false);
        PUSH_TREE_RECORD(m_cr_mc_IsContained,              false);
        PUSH_TREE_RECORD(m_cr_mc_ContainmentFraction,      0.f);
        PUSH_TREE_RECORD(m_cr_mc_TypeTree,                 "-");
        PUSH_TREE_RECORD(m_cr_mc_IsShower,                 false);
        PUSH_TREE_RECORD(m_cr_mc_IsTrack,                  false);
        PUSH_TREE_RECORD(m_cr_mc_IsProton,                 false);
        PUSH_TREE_RECORD(m_cr_mc_IsPionOrMuon,             false);
        PUSH_TREE_RECORD(m_cr_mc_IsCosmicRay,              false);
        PUSH_TREE_RECORD(m_cr_mc_PdgCode,                  0);
        PUSH_TREE_RECORD(m_cr_mc_HitPurity,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,          0.f);
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitPurity,       0.f);
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitCompleteness, 0.f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyCosmicRayRecord(const MCParticle *const pMainMcParticle, const float mcEnergy, const float mcKineticEnergy, const float mcMass,
    const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum,
    const bool mcIsVertexFiducial, const float mcContainmentFraction, const LArAnalysisParticle::TYPE mcType, const bool mcIsShower,
    const int mcPdgCode, const LArAnalysisParticle::TypeTree mcTypeTree) const
{
    // m_cr_Number is dealt with by the calling method.
    
    PUSH_TREE_RECORD(m_cr_WasReconstructed,            false);
    PUSH_TREE_RECORD(m_cr_IsVertexFiducial,            false);
    PUSH_TREE_RECORD(m_cr_IsContained,                 false);
    PUSH_TREE_RECORD(m_cr_FiducialHitFraction,         0.f);
    PUSH_TREE_RECORD(m_cr_HasMcInfo,                   true);
    PUSH_TREE_RECORD(m_cr_KineticEnergy,               0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromRange,                  0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromCorrectedTrackCharge,   0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromUncorrectedTrackCharge, 0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromShowerCharge,           0.f);
    PUSH_TREE_RECORD(m_cr_VertexX,                     0.f);
    PUSH_TREE_RECORD(m_cr_VertexY,                     0.f);
    PUSH_TREE_RECORD(m_cr_VertexZ,                     0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineX,            0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineY,            0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineZ,            0.f);
    PUSH_TREE_RECORD(m_cr_TypeTree,                    "-");
    PUSH_TREE_RECORD(m_cr_NumberOf3dHits,              0U);
    PUSH_TREE_RECORD(m_cr_NumberOfCollectionPlaneHits, 0U);
    PUSH_TREE_RECORD(m_cr_NumberOfDownstreamParticles, 0U);
    
    if (pMainMcParticle)
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid, reinterpret_cast<std::uint64_t>(pMainMcParticle->GetUid()));
    
    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: cosmic ray had MC info but no associated main MC particle" << std::endl;
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid, 0ULL);
    }
    
    PUSH_TREE_RECORD(m_cr_mc_IsParticleSplitByReco,    false); // not covered by any reconstructed particle
    PUSH_TREE_RECORD(m_cr_mc_Energy,                   mcEnergy);
    PUSH_TREE_RECORD(m_cr_mc_KineticEnergy,            mcKineticEnergy);
    PUSH_TREE_RECORD(m_cr_mc_Mass,                     mcMass);
    PUSH_TREE_RECORD(m_cr_mc_VertexX,                  mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_cr_mc_VertexY,                  mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_cr_mc_VertexZ,                  mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,         mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,         mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,         mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_Momentum,                 mcMomentum.GetMagnitude());
    PUSH_TREE_RECORD(m_cr_mc_MomentumX,                mcMomentum.GetX());
    PUSH_TREE_RECORD(m_cr_mc_MomentumY,                mcMomentum.GetY());
    PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,         mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_cr_mc_IsContained,             (mcContainmentFraction >= m_mcContainmentFractionLowerBound));
    PUSH_TREE_RECORD(m_cr_mc_ContainmentFraction,      mcContainmentFraction);
    PUSH_TREE_RECORD(m_cr_mc_TypeTree,                  LArAnalysisParticle::TypeTreeAsString(mcTypeTree));
    PUSH_TREE_RECORD(m_cr_mc_IsShower,                 mcIsShower);
    PUSH_TREE_RECORD(m_cr_mc_IsTrack,                 !mcIsShower);
    PUSH_TREE_RECORD(m_cr_mc_IsProton,                (mcType == LArAnalysisParticle::TYPE::PROTON));
    PUSH_TREE_RECORD(m_cr_mc_IsPionOrMuon,            (mcType == LArAnalysisParticle::TYPE::PION_MUON) || (mcType == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_cr_mc_IsCosmicRay,             (mcType == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_cr_mc_PdgCode,                  mcPdgCode);
    PUSH_TREE_RECORD(m_cr_mc_HitPurity,                0.f); // undefined
    PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,          0.f); // undefined
    PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitPurity,       0.f); // undefined
    PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitCompleteness, 0.f); // undefined
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::DumpTree() const
{
    const std::string nuLabel = " - [nu]        ";
    
    std::cout << "Pandora Tree dump:\n";
    std::cout << nuLabel << "Was reconstructed:                " << std::boolalpha << this->m_treeParameters.m_nu_WasReconstructed << std::noboolalpha << '\n';
    std::cout << nuLabel << "Is vertex fiducial:               " << std::boolalpha << this->m_treeParameters.m_nu_IsVertexFiducial << std::noboolalpha << '\n';
    std::cout << nuLabel << "Is contained:                     " << std::boolalpha << this->m_treeParameters.m_nu_IsContained << std::noboolalpha << '\n';
    std::cout << nuLabel << "Fiducial hit fraction:            " << 100.f * this->m_treeParameters.m_nu_FiducialHitFraction << "%\n";
    std::cout << nuLabel << "Has MC info:                      " << std::boolalpha << this->m_treeParameters.m_nu_HasMcInfo << std::noboolalpha << '\n';
    std::cout << nuLabel << "Visible energy:                   " << 1000.f * this->m_treeParameters.m_nu_VisibleEnergy << " MeV\n";
    std::cout << nuLabel << "Longitudinal visible energy:      " << 1000.f * this->m_treeParameters.m_nu_VisibleLongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "Transverse visible energy:        " << 1000.f * this->m_treeParameters.m_nu_VisibleTransverseEnergy << " MeV\n";
    std::cout << nuLabel << "E frac from range:                " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromRange << "%\n";
    std::cout << nuLabel << "E frac from cor. track charge:    " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromCorrectedTrackCharge << "%\n";
    std::cout << nuLabel << "E frac from uncor. track charge:  " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromUncorrectedTrackCharge << "%\n";
    std::cout << nuLabel << "E frac from shower charge:        " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromShowerCharge << "%\n";
    std::cout << nuLabel << "Vertex x:                         " << this->m_treeParameters.m_nu_VertexX << " cm\n";
    std::cout << nuLabel << "Vertex y:                         " << this->m_treeParameters.m_nu_VertexY << " cm\n";
    std::cout << nuLabel << "Vertex z:                         " << this->m_treeParameters.m_nu_VertexZ << " cm\n";
    std::cout << nuLabel << "Dir cosine x:                     " << this->m_treeParameters.m_nu_DirectionCosineX << '\n';
    std::cout << nuLabel << "Dir cosine y:                     " << this->m_treeParameters.m_nu_DirectionCosineY << '\n';
    std::cout << nuLabel << "Dir cosine z:                     " << this->m_treeParameters.m_nu_DirectionCosineZ << '\n';
    std::cout << nuLabel << "Type tree:                        " << this->m_treeParameters.m_nu_TypeTree << '\n';
    std::cout << nuLabel << "Number of 3D hits:                " << this->m_treeParameters.m_nu_NumberOf3dHits << '\n';
    std::cout << nuLabel << "Number of coll plane hits:        " << this->m_treeParameters.m_nu_NumberOfCollectionPlaneHits << '\n';
    std::cout << nuLabel << "Number of downstream pfos:        " << this->m_treeParameters.m_nu_NumberOfDownstreamParticles << '\n';
    std::cout << nuLabel << "MC particle UID:                  " << this->m_treeParameters.m_nu_mc_McParticleUid << '\n';
    std::cout << nuLabel << "MC energy:                        " << 1000.f * this->m_treeParameters.m_nu_mc_Energy << " MeV\n";
    std::cout << nuLabel << "MC longitudinal energy:           " << 1000.f * this->m_treeParameters.m_nu_mc_LongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "MC transverse energy:             " << 1000.f * this->m_treeParameters.m_nu_mc_TransverseEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible energy:                " << 1000.f * this->m_treeParameters.m_nu_mc_VisibleEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible long energy:           " << 1000.f * this->m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible trans energy:          " << 1000.f * this->m_treeParameters.m_nu_mc_VisibleTransverseEnergy << " MeV\n";
    std::cout << nuLabel << "MC vertex x:                      " << this->m_treeParameters.m_nu_mc_VertexX << " cm\n";
    std::cout << nuLabel << "MC vertex y:                      " << this->m_treeParameters.m_nu_mc_VertexY << " cm\n";
    std::cout << nuLabel << "MC vertex z:                      " << this->m_treeParameters.m_nu_mc_VertexZ << " cm\n";
    std::cout << nuLabel << "MC dir cosine x:                  " << this->m_treeParameters.m_nu_mc_DirectionCosineX << '\n';
    std::cout << nuLabel << "MC dir cosine y:                  " << this->m_treeParameters.m_nu_mc_DirectionCosineY << '\n';
    std::cout << nuLabel << "MC dir cosine z:                  " << this->m_treeParameters.m_nu_mc_DirectionCosineZ << '\n';
    std::cout << nuLabel << "MC momentum:                      " << 1000.f * this->m_treeParameters.m_nu_mc_Momentum << " MeV/c\n";
    std::cout << nuLabel << "MC momentum x:                    " << 1000.f * this->m_treeParameters.m_nu_mc_MomentumX << " MeV/c\n";
    std::cout << nuLabel << "MC momentum y:                    " << 1000.f * this->m_treeParameters.m_nu_mc_MomentumY << " MeV/c\n";
    std::cout << nuLabel << "MC momentum z:                    " << 1000.f * this->m_treeParameters.m_nu_mc_MomentumZ << " MeV/c\n";
    std::cout << nuLabel << "MC is vertex fiducial:            " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsVertexFiducial << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC is contained:                  " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsContained << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC containment fraction:          " << 100.f * this->m_treeParameters.m_nu_mc_ContainmentFraction << "%\n";
    std::cout << nuLabel << "MC type tree:                     " << this->m_treeParameters.m_nu_mc_TypeTree << '\n';
    std::cout << nuLabel << "MC interaction type:              " << this->m_treeParameters.m_nu_mc_InteractionType << '\n';
    std::cout << nuLabel << "MC is charged-current:            " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsChargedCurrent << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC visible energy frac:           " << 100.f * this->m_treeParameters.m_nu_mc_VisibleEnergyFraction << "%\n";
    std::cout << nuLabel << "MC PDG code:                      " << this->m_treeParameters.m_nu_mc_PdgCode << '\n';
    std::cout << nuLabel << "MC hit purity:                    " << 100.f * this->m_treeParameters.m_nu_mc_HitPurity << "%\n";
    std::cout << nuLabel << "MC hit completeness:              " << 100.f * this->m_treeParameters.m_nu_mc_HitCompleteness << "%\n";
    std::cout << nuLabel << "MC hit purity (coll plane):       " << 100.f * this->m_treeParameters.m_nu_mc_CollectionPlaneHitPurity << "%\n";
    std::cout << nuLabel << "MC hit completeness (coll plane): " << 100.f * this->m_treeParameters.m_nu_mc_CollectionPlaneHitCompleteness << "%\n";
    
    std::cout << " - [primary]   Number of primaries:              " << this->m_treeParameters.m_primary_Number << '\n';
    
    for (int i = 0; i < this->m_treeParameters.m_primary_Number; ++i)
    {
        const std::string label = " - [primary " + std::to_string(i) + "] ";
        
        std::cout << label << "Was reconstructed:                " << std::boolalpha << this->m_treeParameters.m_primary_WasReconstructed.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is vertex fiducial:               " << std::boolalpha << this->m_treeParameters.m_primary_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is contained:                     " << std::boolalpha << this->m_treeParameters.m_primary_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Fiducial hit fraction:            " << 100.f * this->m_treeParameters.m_primary_FiducialHitFraction.at(i) << "%\n";
        std::cout << label << "Has MC info:                      " << std::boolalpha << this->m_treeParameters.m_primary_HasMcInfo.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Kinetic energy:                   " << 1000.f * this->m_treeParameters.m_primary_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "KE frac from range:               " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromRange.at(i) << "%\n";
        std::cout << label << "KE frac from cor. track charge:   " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromCorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from uncor. track charge: " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromUncorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from shower charge:       " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromShowerCharge.at(i) << "%\n";
        std::cout << label << "Vertex x:                         " << this->m_treeParameters.m_primary_VertexX.at(i) << " cm\n";
        std::cout << label << "Vertex y:                         " << this->m_treeParameters.m_primary_VertexY.at(i) << " cm\n";
        std::cout << label << "Vertex z:                         " << this->m_treeParameters.m_primary_VertexZ.at(i) << " cm\n";
        std::cout << label << "Dir cosine x:                     " << this->m_treeParameters.m_primary_DirectionCosineX.at(i) << '\n';
        std::cout << label << "Dir cosine y:                     " << this->m_treeParameters.m_primary_DirectionCosineY.at(i) << '\n';
        std::cout << label << "Dir cosine z:                     " << this->m_treeParameters.m_primary_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "Type tree:                        " << this->m_treeParameters.m_primary_TypeTree.at(i) << '\n';
        std::cout << label << "Is shower:                        " << std::boolalpha << m_treeParameters.m_primary_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is track:                         " << std::boolalpha << m_treeParameters.m_primary_IsTrack.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is proton:                        " << std::boolalpha << m_treeParameters.m_primary_IsProton.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is pion or muon:                  " << std::boolalpha << m_treeParameters.m_primary_IsPionOrMuon.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Number of 3D hits:                " << this->m_treeParameters.m_primary_NumberOf3dHits.at(i) << '\n';
        std::cout << label << "Number of coll plane hits:        " << this->m_treeParameters.m_primary_NumberOfCollectionPlaneHits.at(i) << '\n';
        std::cout << label << "Number of downstream pfos:        " << this->m_treeParameters.m_primary_NumberOfDownstreamParticles.at(i) << '\n';
        std::cout << label << "MC particle UID:                  " << this->m_treeParameters.m_primary_mc_McParticleUid.at(i) << '\n';
        std::cout << label << "MC is particle split by reco:     " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsParticleSplitByReco.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC energy:                        " << 1000.f * this->m_treeParameters.m_primary_mc_Energy.at(i) << " MeV\n";
        std::cout << label << "MC kinetic energy:                " << 1000.f * this->m_treeParameters.m_primary_mc_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "MC mass:                          " << 1000.f * this->m_treeParameters.m_primary_mc_Mass.at(i) << " MeV/c^2\n";
        std::cout << label << "MC vertex x:                      " << this->m_treeParameters.m_primary_mc_VertexX.at(i) << " cm\n";
        std::cout << label << "MC vertex y:                      " << this->m_treeParameters.m_primary_mc_VertexY.at(i) << " cm\n";
        std::cout << label << "MC vertex z:                      " << this->m_treeParameters.m_primary_mc_VertexZ.at(i) << " cm\n";
        std::cout << label << "MC dir cosine x:                  " << this->m_treeParameters.m_primary_mc_DirectionCosineX.at(i) << '\n';
        std::cout << label << "MC dir cosine y:                  " << this->m_treeParameters.m_primary_mc_DirectionCosineY.at(i) << '\n';
        std::cout << label << "MC dir cosine z:                  " << this->m_treeParameters.m_primary_mc_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "MC momentum:                      " << 1000.f * this->m_treeParameters.m_primary_mc_Momentum.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum x:                    " << 1000.f * this->m_treeParameters.m_primary_mc_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum y:                    " << 1000.f * this->m_treeParameters.m_primary_mc_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum z:                    " << 1000.f * this->m_treeParameters.m_primary_mc_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "MC is vertex fiducial:            " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is contained:                  " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC containment fraction:          " << 100.f * this->m_treeParameters.m_primary_mc_ContainmentFraction.at(i) << "%\n";
        std::cout << label << "MC type tree:                     " << this->m_treeParameters.m_primary_mc_TypeTree.at(i) << '\n';
        std::cout << label << "MC is shower:                     " << std::boolalpha << m_treeParameters.m_primary_mc_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is track:                      " << std::boolalpha << m_treeParameters.m_primary_mc_IsTrack.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is proton:                     " << std::boolalpha << m_treeParameters.m_primary_mc_IsProton.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is pion or muon:               " << std::boolalpha << m_treeParameters.m_primary_mc_IsPionOrMuon.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is cosmic ray:                 " << std::boolalpha << m_treeParameters.m_primary_mc_IsCosmicRay.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC PDG code:                      " << this->m_treeParameters.m_primary_mc_PdgCode.at(i) << '\n';
        std::cout << label << "MC hit purity:                    " << 100.f * this->m_treeParameters.m_primary_mc_HitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness:              " << 100.f * this->m_treeParameters.m_primary_mc_HitCompleteness.at(i) << "%\n";
        std::cout << label << "MC hit purity (coll plane):       " << 100.f * this->m_treeParameters.m_primary_mc_CollectionPlaneHitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness (coll plane): " << 100.f * this->m_treeParameters.m_primary_mc_CollectionPlaneHitCompleteness.at(i) << "%\n";
    }
    
    std::cout << " - [cr]        Number of cosmic rays:            " << this->m_treeParameters.m_cr_Number << '\n';
    
    for (int i = 0; i < this->m_treeParameters.m_cr_Number; ++i)
    {
        const std::string label = " - [cr " + std::to_string(i) + "]      ";
        
        std::cout << label << "Was reconstructed:                " << std::boolalpha << this->m_treeParameters.m_cr_WasReconstructed.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is vertex fiducial:               " << std::boolalpha << this->m_treeParameters.m_cr_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is contained:                     " << std::boolalpha << this->m_treeParameters.m_cr_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Fiducial hit fraction:            " << 100.f * this->m_treeParameters.m_cr_FiducialHitFraction.at(i) << "%\n";
        std::cout << label << "Has MC info:                      " << std::boolalpha << this->m_treeParameters.m_cr_HasMcInfo.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Kinetic energy:                   " << 1000.f * this->m_treeParameters.m_cr_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "KE frac from range:               " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromRange.at(i) << "%\n";
        std::cout << label << "KE frac from cor. track charge:   " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromCorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from uncor. track charge: " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromUncorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from shower charge:       " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromShowerCharge.at(i) << "%\n";
        std::cout << label << "Vertex x:                         " << this->m_treeParameters.m_cr_VertexX.at(i) << " cm\n";
        std::cout << label << "Vertex y:                         " << this->m_treeParameters.m_cr_VertexY.at(i) << " cm\n";
        std::cout << label << "Vertex z:                         " << this->m_treeParameters.m_cr_VertexZ.at(i) << " cm\n";
        std::cout << label << "Dir cosine x:                     " << this->m_treeParameters.m_cr_DirectionCosineX.at(i) << '\n';
        std::cout << label << "Dir cosine y:                     " << this->m_treeParameters.m_cr_DirectionCosineY.at(i) << '\n';
        std::cout << label << "Dir cosine z:                     " << this->m_treeParameters.m_cr_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "Type tree:                        " << this->m_treeParameters.m_cr_TypeTree.at(i) << '\n';
        std::cout << label << "Number of 3D hits:                " << this->m_treeParameters.m_cr_NumberOf3dHits.at(i) << '\n';
        std::cout << label << "Number of coll plane hits:        " << this->m_treeParameters.m_cr_NumberOfCollectionPlaneHits.at(i) << '\n';
        std::cout << label << "Number of downstream pfos:        " << this->m_treeParameters.m_cr_NumberOfDownstreamParticles.at(i) << '\n';
        std::cout << label << "MC particle UID:                  " << this->m_treeParameters.m_cr_mc_McParticleUid.at(i) << '\n';
        std::cout << label << "MC is particle split by reco:     " << std::boolalpha << this->m_treeParameters.m_cr_mc_IsParticleSplitByReco.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC energy:                        " << 1000.f * this->m_treeParameters.m_cr_mc_Energy.at(i) << " MeV\n";
        std::cout << label << "MC kinetic energy:                " << 1000.f * this->m_treeParameters.m_cr_mc_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "MC mass:                          " << 1000.f * this->m_treeParameters.m_cr_mc_Mass.at(i) << " MeV/c^2\n";
        std::cout << label << "MC vertex x:                      " << this->m_treeParameters.m_cr_mc_VertexX.at(i) << " cm\n";
        std::cout << label << "MC vertex y:                      " << this->m_treeParameters.m_cr_mc_VertexY.at(i) << " cm\n";
        std::cout << label << "MC vertex z:                      " << this->m_treeParameters.m_cr_mc_VertexZ.at(i) << " cm\n";
        std::cout << label << "MC dir cosine x:                  " << this->m_treeParameters.m_cr_mc_DirectionCosineX.at(i) << '\n';
        std::cout << label << "MC dir cosine y:                  " << this->m_treeParameters.m_cr_mc_DirectionCosineY.at(i) << '\n';
        std::cout << label << "MC dir cosine z:                  " << this->m_treeParameters.m_cr_mc_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "MC momentum:                      " << 1000.f * this->m_treeParameters.m_cr_mc_Momentum.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum x:                    " << 1000.f * this->m_treeParameters.m_cr_mc_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum y:                    " << 1000.f * this->m_treeParameters.m_cr_mc_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum z:                    " << 1000.f * this->m_treeParameters.m_cr_mc_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "MC is vertex fiducial:            " << std::boolalpha << this->m_treeParameters.m_cr_mc_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is contained:                  " << std::boolalpha << this->m_treeParameters.m_cr_mc_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC containment fraction:          " << 100.f * this->m_treeParameters.m_cr_mc_ContainmentFraction.at(i) << "%\n";
        std::cout << label << "MC type tree:                     " << this->m_treeParameters.m_cr_mc_TypeTree.at(i) << '\n';
        std::cout << label << "MC is shower:                     " << std::boolalpha << m_treeParameters.m_cr_mc_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is track:                      " << std::boolalpha << m_treeParameters.m_cr_mc_IsTrack.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is proton:                     " << std::boolalpha << m_treeParameters.m_cr_mc_IsProton.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is pion or muon:               " << std::boolalpha << m_treeParameters.m_cr_mc_IsPionOrMuon.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is cosmic ray:                 " << std::boolalpha << m_treeParameters.m_cr_mc_IsCosmicRay.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC PDG code:                      " << this->m_treeParameters.m_cr_mc_PdgCode.at(i) << '\n';
        std::cout << label << "MC hit purity:                    " << 100.f * this->m_treeParameters.m_cr_mc_HitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness:              " << 100.f * this->m_treeParameters.m_cr_mc_HitCompleteness.at(i) << "%\n";
        std::cout << label << "MC hit purity (coll plane):       " << 100.f * this->m_treeParameters.m_cr_mc_CollectionPlaneHitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness (coll plane): " << 100.f * this->m_treeParameters.m_cr_mc_CollectionPlaneHitCompleteness.at(i) << "%\n";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool WriteAnalysisParticlesAlgorithm::IsChargedCurrent(const LArInteractionTypeHelper::InteractionType interactionType) const
{
    switch (interactionType)
    {
        case LArInteractionTypeHelper::CCQEL_MU:
        case LArInteractionTypeHelper::CCQEL_MU_P:
        case LArInteractionTypeHelper::CCQEL_MU_P_P:
        case LArInteractionTypeHelper::CCQEL_MU_P_P_P:
        case LArInteractionTypeHelper::CCQEL_MU_P_P_P_P:
        case LArInteractionTypeHelper::CCQEL_MU_P_P_P_P_P:
        case LArInteractionTypeHelper::CCQEL_E:
        case LArInteractionTypeHelper::CCQEL_E_P:
        case LArInteractionTypeHelper::CCQEL_E_P_P:
        case LArInteractionTypeHelper::CCQEL_E_P_P_P:
        case LArInteractionTypeHelper::CCQEL_E_P_P_P_P:
        case LArInteractionTypeHelper::CCQEL_E_P_P_P_P_P:
        case LArInteractionTypeHelper::CCRES_MU:
        case LArInteractionTypeHelper::CCRES_MU_P:
        case LArInteractionTypeHelper::CCRES_MU_P_P:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P_P:
        case LArInteractionTypeHelper::CCRES_MU_PIPLUS:
        case LArInteractionTypeHelper::CCRES_MU_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_MU_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_MU_PHOTON:
        case LArInteractionTypeHelper::CCRES_MU_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_MU_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_MU_PIZERO:
        case LArInteractionTypeHelper::CCRES_MU_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_MU_P_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_MU_P_P_P_P_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_E:
        case LArInteractionTypeHelper::CCRES_E_P:
        case LArInteractionTypeHelper::CCRES_E_P_P:
        case LArInteractionTypeHelper::CCRES_E_P_P_P:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P_P:
        case LArInteractionTypeHelper::CCRES_E_PIPLUS:
        case LArInteractionTypeHelper::CCRES_E_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_E_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P_P_PIPLUS:
        case LArInteractionTypeHelper::CCRES_E_PHOTON:
        case LArInteractionTypeHelper::CCRES_E_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_E_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P_P_PHOTON:
        case LArInteractionTypeHelper::CCRES_E_PIZERO:
        case LArInteractionTypeHelper::CCRES_E_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_E_P_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P_PIZERO:
        case LArInteractionTypeHelper::CCRES_E_P_P_P_P_P_PIZERO:
        case LArInteractionTypeHelper::CCDIS:
        case LArInteractionTypeHelper::CCCOH:
            return true;
            
        default: break;
    };
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode WriteAnalysisParticlesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", this->m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", 
        this->m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", 
        this->m_caloHitListName));    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", this->m_outputFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Verbose",
        this->m_verbose));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowXMargin", this->m_fiducialCutLowXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighXMargin", this->m_fiducialCutHighXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowYMargin", this->m_fiducialCutLowYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighYMargin", this->m_fiducialCutHighYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowZMargin", this->m_fiducialCutLowZMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighZMargin", this->m_fiducialCutHighZMargin));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McContainmentFractionLowerBound", this->m_mcContainmentFractionLowerBound));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialHitFractionLowerBound", this->m_fiducialHitFractionLowerBound));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McOnlyParticleContainmentCut", this->m_mcOnlyParticleContainmentCut));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McOnlyParticleEnergyCut", this->m_mcOnlyParticleEnergyCut));
    
    // Use the detector geometry and the margins to get the maximum and minimum fiducial volume coordinates.
    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), m_fiducialCutLowXMargin, m_fiducialCutHighXMargin,
        m_fiducialCutLowYMargin, this->m_fiducialCutHighYMargin, m_fiducialCutLowZMargin, m_fiducialCutHighZMargin, this->m_minCoordinates, this->m_maxCoordinates);
        
    this->m_pOutputTFile = new TFile(this->m_outputFile.c_str(), "UPDATE");
    this->m_pOutputTree = new TTree("PandoraTree", "PandoraTree");
    
    // Neutrino parameters.
    this->m_pOutputTree->Branch("nu_WasReconstructed",                 &this->m_treeParameters.m_nu_WasReconstructed);
    this->m_pOutputTree->Branch("nu_IsVertexFiducial",                 &this->m_treeParameters.m_nu_IsVertexFiducial);
    this->m_pOutputTree->Branch("nu_IsContained",                      &this->m_treeParameters.m_nu_IsContained);
    this->m_pOutputTree->Branch("nu_FiducialHitFraction",              &this->m_treeParameters.m_nu_FiducialHitFraction);
    this->m_pOutputTree->Branch("nu_HasMcInfo",                        &this->m_treeParameters.m_nu_HasMcInfo);
    this->m_pOutputTree->Branch("nu_VisibleEnergy",                    &this->m_treeParameters.m_nu_VisibleEnergy);
    this->m_pOutputTree->Branch("nu_VisibleLongitudinalEnergy",                   &this->m_treeParameters.m_nu_VisibleLongitudinalEnergy);
    this->m_pOutputTree->Branch("nu_VisibleTransverseEnergy",                     &this->m_treeParameters.m_nu_VisibleTransverseEnergy);
    this->m_pOutputTree->Branch("nu_VisibleEnergyFracFromRange",                  &this->m_treeParameters.m_nu_VisibleEnergyFracFromRange);
    this->m_pOutputTree->Branch("nu_VisibleEnergyFracFromCorrectedTrackCharge",   &this->m_treeParameters.m_nu_VisibleEnergyFracFromCorrectedTrackCharge);
    this->m_pOutputTree->Branch("nu_VisibleEnergyFracFromUncorrectedTrackCharge", &this->m_treeParameters.m_nu_VisibleEnergyFracFromUncorrectedTrackCharge);
    this->m_pOutputTree->Branch("nu_VisibleEnergyFracFromShowerCharge",           &this->m_treeParameters.m_nu_VisibleEnergyFracFromShowerCharge);
    this->m_pOutputTree->Branch("nu_VertexX",                          &this->m_treeParameters.m_nu_VertexX);
    this->m_pOutputTree->Branch("nu_VertexY",                          &this->m_treeParameters.m_nu_VertexY);
    this->m_pOutputTree->Branch("nu_VertexZ",                          &this->m_treeParameters.m_nu_VertexZ);
    this->m_pOutputTree->Branch("nu_DirectionCosineX",                 &this->m_treeParameters.m_nu_DirectionCosineX);
    this->m_pOutputTree->Branch("nu_DirectionCosineY",                 &this->m_treeParameters.m_nu_DirectionCosineY);
    this->m_pOutputTree->Branch("nu_DirectionCosineZ",                 &this->m_treeParameters.m_nu_DirectionCosineZ);
    this->m_pOutputTree->Branch("nu_TypeTree",                         &this->m_treeParameters.m_nu_TypeTree);
    this->m_pOutputTree->Branch("nu_NumberOf3dHits",                   &this->m_treeParameters.m_nu_NumberOf3dHits);
    this->m_pOutputTree->Branch("nu_NumberOfCollectionPlaneHits",      &this->m_treeParameters.m_nu_NumberOfCollectionPlaneHits);
    this->m_pOutputTree->Branch("nu_NumberOfDownstreamParticles",      &this->m_treeParameters.m_nu_NumberOfDownstreamParticles);
    this->m_pOutputTree->Branch("nu_mc_McParticleUid",                 &this->m_treeParameters.m_nu_mc_McParticleUid);
    this->m_pOutputTree->Branch("nu_mc_Energy",                        &this->m_treeParameters.m_nu_mc_Energy);
    this->m_pOutputTree->Branch("nu_mc_LongitudinalEnergy",            &this->m_treeParameters.m_nu_mc_LongitudinalEnergy);
    this->m_pOutputTree->Branch("nu_mc_TransverseEnergy",              &this->m_treeParameters.m_nu_mc_TransverseEnergy);
    this->m_pOutputTree->Branch("nu_mc_VisibleEnergy",                 &this->m_treeParameters.m_nu_mc_VisibleEnergy);
    this->m_pOutputTree->Branch("nu_mc_VisibleLongitudinalEnergy",     &this->m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy);
    this->m_pOutputTree->Branch("nu_mc_VisibleTransverseEnergy",       &this->m_treeParameters.m_nu_mc_VisibleTransverseEnergy);
    this->m_pOutputTree->Branch("nu_mc_VertexX",                       &this->m_treeParameters.m_nu_mc_VertexX);
    this->m_pOutputTree->Branch("nu_mc_VertexY",                       &this->m_treeParameters.m_nu_mc_VertexY);
    this->m_pOutputTree->Branch("nu_mc_VertexZ",                       &this->m_treeParameters.m_nu_mc_VertexZ);
    this->m_pOutputTree->Branch("nu_mc_DirectionCosineX",              &this->m_treeParameters.m_nu_mc_DirectionCosineX);
    this->m_pOutputTree->Branch("nu_mc_DirectionCosineY",              &this->m_treeParameters.m_nu_mc_DirectionCosineY);
    this->m_pOutputTree->Branch("nu_mc_DirectionCosineZ",              &this->m_treeParameters.m_nu_mc_DirectionCosineZ);
    this->m_pOutputTree->Branch("nu_mc_Momentum",                      &this->m_treeParameters.m_nu_mc_Momentum);
    this->m_pOutputTree->Branch("nu_mc_MomentumX",                     &this->m_treeParameters.m_nu_mc_MomentumX);
    this->m_pOutputTree->Branch("nu_mc_MomentumY",                     &this->m_treeParameters.m_nu_mc_MomentumY);
    this->m_pOutputTree->Branch("nu_mc_MomentumZ",                     &this->m_treeParameters.m_nu_mc_MomentumZ);
    this->m_pOutputTree->Branch("nu_mc_IsVertexFiducial",              &this->m_treeParameters.m_nu_mc_IsVertexFiducial);
    this->m_pOutputTree->Branch("nu_mc_IsContained",                   &this->m_treeParameters.m_nu_mc_IsContained);
    this->m_pOutputTree->Branch("nu_mc_ContainmentFraction",           &this->m_treeParameters.m_nu_mc_ContainmentFraction);
    this->m_pOutputTree->Branch("nu_mc_TypeTree",                      &this->m_treeParameters.m_nu_mc_TypeTree);
    this->m_pOutputTree->Branch("nu_mc_InteractionType",               &this->m_treeParameters.m_nu_mc_InteractionType);
    this->m_pOutputTree->Branch("nu_mc_IsChargedCurrent",              &this->m_treeParameters.m_nu_mc_IsChargedCurrent);
    this->m_pOutputTree->Branch("nu_mc_VisibleEnergyFraction",         &this->m_treeParameters.m_nu_mc_VisibleEnergyFraction);
    this->m_pOutputTree->Branch("nu_mc_PdgCode",                       &this->m_treeParameters.m_nu_mc_PdgCode);
    this->m_pOutputTree->Branch("nu_mc_HitPurity",                     &this->m_treeParameters.m_nu_mc_HitPurity);
    this->m_pOutputTree->Branch("nu_mc_HitCompleteness",               &this->m_treeParameters.m_nu_mc_HitCompleteness);
    this->m_pOutputTree->Branch("nu_mc_CollectionPlaneHitPurity",       &this->m_treeParameters.m_nu_mc_CollectionPlaneHitPurity);
    this->m_pOutputTree->Branch("nu_mc_CollectionPlaneHitCompleteness", &this->m_treeParameters.m_nu_mc_CollectionPlaneHitCompleteness);
    
    // Primary neutrino daughter parameters.
    this->m_pOutputTree->Branch("primary_Number",                      &this->m_treeParameters.m_primary_Number);
    
    this->m_pOutputTree->Branch("primary_WasReconstructed",            &this->m_treeParameters.m_primary_WasReconstructed);
    this->m_pOutputTree->Branch("primary_IsVertexFiducial",            &this->m_treeParameters.m_primary_IsVertexFiducial);
    this->m_pOutputTree->Branch("primary_IsContained",                 &this->m_treeParameters.m_primary_IsContained);
    this->m_pOutputTree->Branch("primary_FiducialHitFraction",         &this->m_treeParameters.m_primary_FiducialHitFraction);
    this->m_pOutputTree->Branch("primary_HasMcInfo",                   &this->m_treeParameters.m_primary_HasMcInfo);
    this->m_pOutputTree->Branch("primary_KineticEnergy",               &this->m_treeParameters.m_primary_KineticEnergy);
    this->m_pOutputTree->Branch("primary_KineticEnergyFracFromRange",                  &this->m_treeParameters.m_primary_KineticEnergyFracFromRange);
    this->m_pOutputTree->Branch("primary_KineticEnergyFracFromCorrectedTrackCharge",   &this->m_treeParameters.m_primary_KineticEnergyFracFromCorrectedTrackCharge);
    this->m_pOutputTree->Branch("primary_KineticEnergyFracFromUncorrectedTrackCharge", &this->m_treeParameters.m_primary_KineticEnergyFracFromUncorrectedTrackCharge);
    this->m_pOutputTree->Branch("primary_KineticEnergyFracFromShowerCharge",           &this->m_treeParameters.m_primary_KineticEnergyFracFromShowerCharge);
    this->m_pOutputTree->Branch("primary_VertexX",                     &this->m_treeParameters.m_primary_VertexX);
    this->m_pOutputTree->Branch("primary_VertexY",                     &this->m_treeParameters.m_primary_VertexY);
    this->m_pOutputTree->Branch("primary_VertexZ",                     &this->m_treeParameters.m_primary_VertexZ);
    this->m_pOutputTree->Branch("primary_DirectionCosineX",            &this->m_treeParameters.m_primary_DirectionCosineX);
    this->m_pOutputTree->Branch("primary_DirectionCosineY",            &this->m_treeParameters.m_primary_DirectionCosineY);
    this->m_pOutputTree->Branch("primary_DirectionCosineZ",            &this->m_treeParameters.m_primary_DirectionCosineZ);
    this->m_pOutputTree->Branch("primary_TypeTree",                   &this->m_treeParameters.m_primary_TypeTree);
    this->m_pOutputTree->Branch("primary_IsShower",                    &this->m_treeParameters.m_primary_IsShower);
    this->m_pOutputTree->Branch("primary_IsTrack",                     &this->m_treeParameters.m_primary_IsTrack);
    this->m_pOutputTree->Branch("primary_IsProton",                    &this->m_treeParameters.m_primary_IsProton);
    this->m_pOutputTree->Branch("primary_IsPionOrMuon",                &this->m_treeParameters.m_primary_IsPionOrMuon);
    this->m_pOutputTree->Branch("primary_NumberOf3dHits",              &this->m_treeParameters.m_primary_NumberOf3dHits);
    this->m_pOutputTree->Branch("primary_NumberOfCollectionPlaneHits", &this->m_treeParameters.m_primary_NumberOfCollectionPlaneHits);
    this->m_pOutputTree->Branch("primary_NumberOfDownstreamParticles", &this->m_treeParameters.m_primary_NumberOfDownstreamParticles);
    this->m_pOutputTree->Branch("primary_mc_McParticleUid",            &this->m_treeParameters.m_primary_mc_McParticleUid);
    this->m_pOutputTree->Branch("primary_mc_IsParticleSplitByReco",    &this->m_treeParameters.m_primary_mc_IsParticleSplitByReco);
    this->m_pOutputTree->Branch("primary_mc_Energy",                   &this->m_treeParameters.m_primary_mc_Energy);
    this->m_pOutputTree->Branch("primary_mc_KineticEnergy",            &this->m_treeParameters.m_primary_mc_KineticEnergy);
    this->m_pOutputTree->Branch("primary_mc_Mass",                     &this->m_treeParameters.m_primary_mc_Mass);
    this->m_pOutputTree->Branch("primary_mc_VertexX",                  &this->m_treeParameters.m_primary_mc_VertexX);
    this->m_pOutputTree->Branch("primary_mc_VertexY",                  &this->m_treeParameters.m_primary_mc_VertexY);
    this->m_pOutputTree->Branch("primary_mc_VertexZ",                  &this->m_treeParameters.m_primary_mc_VertexZ);
    this->m_pOutputTree->Branch("primary_mc_DirectionCosineX",         &this->m_treeParameters.m_primary_mc_DirectionCosineX);
    this->m_pOutputTree->Branch("primary_mc_DirectionCosineY",         &this->m_treeParameters.m_primary_mc_DirectionCosineY);
    this->m_pOutputTree->Branch("primary_mc_DirectionCosineZ",         &this->m_treeParameters.m_primary_mc_DirectionCosineZ);
    this->m_pOutputTree->Branch("primary_mc_Momentum",                 &this->m_treeParameters.m_primary_mc_Momentum);
    this->m_pOutputTree->Branch("primary_mc_MomentumX",                &this->m_treeParameters.m_primary_mc_MomentumX);
    this->m_pOutputTree->Branch("primary_mc_MomentumY",                &this->m_treeParameters.m_primary_mc_MomentumY);
    this->m_pOutputTree->Branch("primary_mc_MomentumZ",                &this->m_treeParameters.m_primary_mc_MomentumZ);
    this->m_pOutputTree->Branch("primary_mc_IsVertexFiducial",         &this->m_treeParameters.m_primary_mc_IsVertexFiducial);
    this->m_pOutputTree->Branch("primary_mc_IsContained",              &this->m_treeParameters.m_primary_mc_IsContained);
    this->m_pOutputTree->Branch("primary_mc_ContainmentFraction",      &this->m_treeParameters.m_primary_mc_ContainmentFraction);
    this->m_pOutputTree->Branch("primary_mc_TypeTree",                 &this->m_treeParameters.m_primary_mc_TypeTree);
    this->m_pOutputTree->Branch("primary_mc_IsShower",                 &this->m_treeParameters.m_primary_mc_IsShower);
    this->m_pOutputTree->Branch("primary_mc_IsTrack",                  &this->m_treeParameters.m_primary_mc_IsTrack);
    this->m_pOutputTree->Branch("primary_mc_IsProton",                 &this->m_treeParameters.m_primary_mc_IsProton);
    this->m_pOutputTree->Branch("primary_mc_IsPionOrMuon",             &this->m_treeParameters.m_primary_mc_IsPionOrMuon);
    this->m_pOutputTree->Branch("primary_mc_IsCosmicRay",              &this->m_treeParameters.m_primary_mc_IsCosmicRay);
    this->m_pOutputTree->Branch("primary_mc_PdgCode",                  &this->m_treeParameters.m_primary_mc_PdgCode);
    this->m_pOutputTree->Branch("primary_mc_HitPurity",                &this->m_treeParameters.m_primary_mc_HitPurity);
    this->m_pOutputTree->Branch("primary_mc_HitCompleteness",          &this->m_treeParameters.m_primary_mc_HitCompleteness);
    this->m_pOutputTree->Branch("primary_mc_CollectionPlaneHitPurity",       &this->m_treeParameters.m_primary_mc_CollectionPlaneHitPurity);
    this->m_pOutputTree->Branch("primary_mc_CollectionPlaneHitCompleteness", &this->m_treeParameters.m_primary_mc_CollectionPlaneHitCompleteness);
    
    // Cosmic ray parameters.
    this->m_pOutputTree->Branch("cr_Number",                           &this->m_treeParameters.m_cr_Number);
    this->m_pOutputTree->Branch("cr_WasReconstructed",                 &this->m_treeParameters.m_cr_WasReconstructed);
    this->m_pOutputTree->Branch("cr_IsVertexFiducial",                 &this->m_treeParameters.m_cr_IsVertexFiducial);
    this->m_pOutputTree->Branch("cr_IsContained",                      &this->m_treeParameters.m_cr_IsContained);
    this->m_pOutputTree->Branch("cr_FiducialHitFraction",              &this->m_treeParameters.m_cr_FiducialHitFraction);
    this->m_pOutputTree->Branch("cr_HasMcInfo",                        &this->m_treeParameters.m_cr_HasMcInfo);
    this->m_pOutputTree->Branch("cr_KineticEnergy",                    &this->m_treeParameters.m_cr_KineticEnergy);
    this->m_pOutputTree->Branch("cr_KineticEnergyFracFromRange",                  &this->m_treeParameters.m_cr_KineticEnergyFracFromRange);
    this->m_pOutputTree->Branch("cr_KineticEnergyFracFromCorrectedTrackCharge",   &this->m_treeParameters.m_cr_KineticEnergyFracFromCorrectedTrackCharge);
    this->m_pOutputTree->Branch("cr_KineticEnergyFracFromUncorrectedTrackCharge", &this->m_treeParameters.m_cr_KineticEnergyFracFromUncorrectedTrackCharge);
    this->m_pOutputTree->Branch("cr_KineticEnergyFracFromShowerCharge",           &this->m_treeParameters.m_cr_KineticEnergyFracFromShowerCharge);
    this->m_pOutputTree->Branch("cr_VertexX",                          &this->m_treeParameters.m_cr_VertexX);
    this->m_pOutputTree->Branch("cr_VertexY",                          &this->m_treeParameters.m_cr_VertexY);
    this->m_pOutputTree->Branch("cr_VertexZ",                          &this->m_treeParameters.m_cr_VertexZ);
    this->m_pOutputTree->Branch("cr_DirectionCosineX",                 &this->m_treeParameters.m_cr_DirectionCosineX);
    this->m_pOutputTree->Branch("cr_DirectionCosineY",                 &this->m_treeParameters.m_cr_DirectionCosineY);
    this->m_pOutputTree->Branch("cr_DirectionCosineZ",                 &this->m_treeParameters.m_cr_DirectionCosineZ);
    this->m_pOutputTree->Branch("cr_TypeTree",                         &this->m_treeParameters.m_cr_TypeTree);
    this->m_pOutputTree->Branch("cr_NumberOf3dHits",                   &this->m_treeParameters.m_cr_NumberOf3dHits);
    this->m_pOutputTree->Branch("cr_NumberOfCollectionPlaneHits",      &this->m_treeParameters.m_cr_NumberOfCollectionPlaneHits);
    this->m_pOutputTree->Branch("cr_NumberOfDownstreamParticles",      &this->m_treeParameters.m_cr_NumberOfDownstreamParticles);
    this->m_pOutputTree->Branch("cr_mc_McParticleUid",                 &this->m_treeParameters.m_cr_mc_McParticleUid);
    this->m_pOutputTree->Branch("cr_mc_IsParticleSplitByReco",         &this->m_treeParameters.m_cr_mc_IsParticleSplitByReco);
    this->m_pOutputTree->Branch("cr_mc_Energy",                        &this->m_treeParameters.m_cr_mc_Energy);
    this->m_pOutputTree->Branch("cr_mc_KineticEnergy",                 &this->m_treeParameters.m_cr_mc_KineticEnergy);
    this->m_pOutputTree->Branch("cr_mc_Mass",                          &this->m_treeParameters.m_cr_mc_Mass);
    this->m_pOutputTree->Branch("cr_mc_VertexX",                       &this->m_treeParameters.m_cr_mc_VertexX);
    this->m_pOutputTree->Branch("cr_mc_VertexY",                       &this->m_treeParameters.m_cr_mc_VertexY);
    this->m_pOutputTree->Branch("cr_mc_VertexZ",                       &this->m_treeParameters.m_cr_mc_VertexZ);
    this->m_pOutputTree->Branch("cr_mc_DirectionCosineX",              &this->m_treeParameters.m_cr_mc_DirectionCosineX);
    this->m_pOutputTree->Branch("cr_mc_DirectionCosineY",              &this->m_treeParameters.m_cr_mc_DirectionCosineY);
    this->m_pOutputTree->Branch("cr_mc_DirectionCosineZ",              &this->m_treeParameters.m_cr_mc_DirectionCosineZ);
    this->m_pOutputTree->Branch("cr_mc_Momentum",                      &this->m_treeParameters.m_cr_mc_Momentum);
    this->m_pOutputTree->Branch("cr_mc_MomentumX",                     &this->m_treeParameters.m_cr_mc_MomentumX);
    this->m_pOutputTree->Branch("cr_mc_MomentumY",                     &this->m_treeParameters.m_cr_mc_MomentumY);
    this->m_pOutputTree->Branch("cr_mc_MomentumZ",                     &this->m_treeParameters.m_cr_mc_MomentumZ);
    this->m_pOutputTree->Branch("cr_mc_IsVertexFiducial",              &this->m_treeParameters.m_cr_mc_IsVertexFiducial);
    this->m_pOutputTree->Branch("cr_mc_IsContained",                   &this->m_treeParameters.m_cr_mc_IsContained);
    this->m_pOutputTree->Branch("cr_mc_ContainmentFraction",           &this->m_treeParameters.m_cr_mc_ContainmentFraction);
    this->m_pOutputTree->Branch("cr_mc_TypeTree",                      &this->m_treeParameters.m_cr_mc_TypeTree);
    this->m_pOutputTree->Branch("cr_mc_IsShower",                      &this->m_treeParameters.m_cr_mc_IsShower);
    this->m_pOutputTree->Branch("cr_mc_IsTrack",                       &this->m_treeParameters.m_cr_mc_IsTrack);
    this->m_pOutputTree->Branch("cr_mc_IsProton",                      &this->m_treeParameters.m_cr_mc_IsProton);
    this->m_pOutputTree->Branch("cr_mc_IsPionOrMuon",                  &this->m_treeParameters.m_cr_mc_IsPionOrMuon);
    this->m_pOutputTree->Branch("cr_mc_IsCosmicRay",                   &this->m_treeParameters.m_cr_mc_IsCosmicRay);
    this->m_pOutputTree->Branch("cr_mc_PdgCode",                       &this->m_treeParameters.m_cr_mc_PdgCode);
    this->m_pOutputTree->Branch("cr_mc_HitPurity",                     &this->m_treeParameters.m_cr_mc_HitPurity);
    this->m_pOutputTree->Branch("cr_mc_HitCompleteness",               &this->m_treeParameters.m_cr_mc_HitCompleteness);
    this->m_pOutputTree->Branch("cr_mc_CollectionPlaneHitPurity",       &this->m_treeParameters.m_cr_mc_CollectionPlaneHitPurity);
    this->m_pOutputTree->Branch("cr_mc_CollectionPlaneHitCompleteness", &this->m_treeParameters.m_cr_mc_CollectionPlaneHitCompleteness);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TreeParameters::TreeParameters() noexcept :
    m_nu_WasReconstructed(false),
    m_nu_IsVertexFiducial(false),
    m_nu_IsContained(false),
    m_nu_FiducialHitFraction(0.f),
    m_nu_HasMcInfo(false),
    m_nu_VisibleEnergy(0.f),
    m_nu_VisibleLongitudinalEnergy(0.f),
    m_nu_VisibleTransverseEnergy(0.f),
    m_nu_VisibleEnergyFracFromRange(0.f),
    m_nu_VisibleEnergyFracFromCorrectedTrackCharge(0.f),
    m_nu_VisibleEnergyFracFromUncorrectedTrackCharge(0.f),
    m_nu_VisibleEnergyFracFromShowerCharge(0.f),
    m_nu_VertexX(0.f),
    m_nu_VertexY(0.f),
    m_nu_VertexZ(0.f),
    m_nu_DirectionCosineX(0.f),
    m_nu_DirectionCosineY(0.f),
    m_nu_DirectionCosineZ(0.f),
    m_nu_TypeTree("-"),
    m_nu_NumberOf3dHits(0U),
    m_nu_NumberOfCollectionPlaneHits(0U),
    m_nu_NumberOfDownstreamParticles(0U),
    m_nu_mc_McParticleUid(0ULL),
    m_nu_mc_Energy(0.f),
    m_nu_mc_LongitudinalEnergy(0.f),
    m_nu_mc_TransverseEnergy(0.f),
    m_nu_mc_VisibleEnergy(0.f),
    m_nu_mc_VisibleLongitudinalEnergy(0.f),
    m_nu_mc_VisibleTransverseEnergy(0.f),
    m_nu_mc_VertexX(0.f),
    m_nu_mc_VertexY(0.f),
    m_nu_mc_VertexZ(0.f),
    m_nu_mc_DirectionCosineX(0.f),
    m_nu_mc_DirectionCosineY(0.f),
    m_nu_mc_DirectionCosineZ(0.f),
    m_nu_mc_Momentum(0.f),
    m_nu_mc_MomentumX(0.f),
    m_nu_mc_MomentumY(0.f),
    m_nu_mc_MomentumZ(0.f),
    m_nu_mc_IsVertexFiducial(false),
    m_nu_mc_IsContained(false),
    m_nu_mc_ContainmentFraction(0.f),
    m_nu_mc_TypeTree("-"),
    m_nu_mc_InteractionType(),
    m_nu_mc_IsChargedCurrent(false),
    m_nu_mc_VisibleEnergyFraction(0.f),
    m_nu_mc_PdgCode(0),
    m_nu_mc_HitPurity(0.f),
    m_nu_mc_HitCompleteness(0.f),
    m_nu_mc_CollectionPlaneHitPurity(0.f),
    m_nu_mc_CollectionPlaneHitCompleteness(0.f),
    m_primary_Number(0U),
    m_primary_WasReconstructed(),
    m_primary_IsVertexFiducial(),
    m_primary_IsContained(),
    m_primary_FiducialHitFraction(),
    m_primary_HasMcInfo(),
    m_primary_KineticEnergy(),
    m_primary_KineticEnergyFracFromRange(),
    m_primary_KineticEnergyFracFromCorrectedTrackCharge(),
    m_primary_KineticEnergyFracFromUncorrectedTrackCharge(),
    m_primary_KineticEnergyFracFromShowerCharge(),
    m_primary_VertexX(),
    m_primary_VertexY(),
    m_primary_VertexZ(),
    m_primary_DirectionCosineX(),
    m_primary_DirectionCosineY(),
    m_primary_DirectionCosineZ(),
    m_primary_TypeTree(),
    m_primary_IsShower(),
    m_primary_IsTrack(),
    m_primary_IsProton(),
    m_primary_IsPionOrMuon(),
    m_primary_NumberOf3dHits(),
    m_primary_NumberOfCollectionPlaneHits(),
    m_primary_NumberOfDownstreamParticles(),
    m_primary_mc_McParticleUid(),
    m_primary_mc_IsParticleSplitByReco(),
    m_primary_mc_Energy(),
    m_primary_mc_KineticEnergy(),
    m_primary_mc_Mass(),
    m_primary_mc_VertexX(),
    m_primary_mc_VertexY(),
    m_primary_mc_VertexZ(),
    m_primary_mc_DirectionCosineX(),
    m_primary_mc_DirectionCosineY(),
    m_primary_mc_DirectionCosineZ(),
    m_primary_mc_Momentum(),
    m_primary_mc_MomentumX(),
    m_primary_mc_MomentumY(),
    m_primary_mc_MomentumZ(),
    m_primary_mc_IsVertexFiducial(),
    m_primary_mc_IsContained(),
    m_primary_mc_ContainmentFraction(),
    m_primary_mc_TypeTree(),
    m_primary_mc_IsShower(),
    m_primary_mc_IsTrack(),
    m_primary_mc_IsProton(),
    m_primary_mc_IsPionOrMuon(),
    m_primary_mc_IsCosmicRay(),
    m_primary_mc_PdgCode(),
    m_primary_mc_HitPurity(),
    m_primary_mc_HitCompleteness(),
    m_primary_mc_CollectionPlaneHitPurity(),
    m_primary_mc_CollectionPlaneHitCompleteness(),
    m_cr_Number(0U),
    m_cr_WasReconstructed(),
    m_cr_IsVertexFiducial(),
    m_cr_IsContained(),
    m_cr_FiducialHitFraction(),
    m_cr_HasMcInfo(),
    m_cr_KineticEnergy(),
    m_cr_KineticEnergyFracFromRange(),
    m_cr_KineticEnergyFracFromCorrectedTrackCharge(),
    m_cr_KineticEnergyFracFromUncorrectedTrackCharge(),
    m_cr_KineticEnergyFracFromShowerCharge(),
    m_cr_VertexX(),
    m_cr_VertexY(),
    m_cr_VertexZ(),
    m_cr_DirectionCosineX(),
    m_cr_DirectionCosineY(),
    m_cr_DirectionCosineZ(),
    m_cr_TypeTree(),
    m_cr_NumberOf3dHits(),
    m_cr_NumberOfCollectionPlaneHits(),
    m_cr_NumberOfDownstreamParticles(),
    m_cr_mc_McParticleUid(),
    m_cr_mc_IsParticleSplitByReco(),
    m_cr_mc_Energy(),
    m_cr_mc_KineticEnergy(),
    m_cr_mc_Mass(),
    m_cr_mc_VertexX(),
    m_cr_mc_VertexY(),
    m_cr_mc_VertexZ(),
    m_cr_mc_DirectionCosineX(),
    m_cr_mc_DirectionCosineY(),
    m_cr_mc_DirectionCosineZ(),
    m_cr_mc_Momentum(),
    m_cr_mc_MomentumX(),
    m_cr_mc_MomentumY(),
    m_cr_mc_MomentumZ(),
    m_cr_mc_IsVertexFiducial(),
    m_cr_mc_IsContained(),
    m_cr_mc_ContainmentFraction(),
    m_cr_mc_TypeTree(),
    m_cr_mc_IsShower(),
    m_cr_mc_IsTrack(),
    m_cr_mc_IsProton(),
    m_cr_mc_IsPionOrMuon(),
    m_cr_mc_IsCosmicRay(),
    m_cr_mc_PdgCode(),
    m_cr_mc_HitPurity(),
    m_cr_mc_HitCompleteness(),
    m_cr_mc_CollectionPlaneHitPurity(),
    m_cr_mc_CollectionPlaneHitCompleteness()
{
}

} // namespace lar_physics_content
