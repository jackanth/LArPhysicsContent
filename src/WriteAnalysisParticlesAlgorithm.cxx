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
    m_fiducialCutXMargin(10.f),
    m_fiducialCutYMargin(20.f),
    m_fiducialCutZMargin(10.f),
    m_minCoordinates(0.f, 0.f, 0.f),
    m_maxCoordinates(0.f, 0.f, 0.f),
    m_mcContainmentFractionLowerBound(0.9f)
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

    MCParticleSet coveredMcPrimaries;
    
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
                
                if (pAnalysisParticle->McMainMCParticle())
                    coveredMcPrimaries.insert(pAnalysisParticle->McMainMCParticle());
            }
            
            else if (LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(pPfo))
            {
                this->AddPrimaryDaughterRecord(*pAnalysisParticle);
                ++this->m_treeParameters.m_primary_Number;
                
                if (pAnalysisParticle->McMainMCParticle())
                    coveredMcPrimaries.insert(pAnalysisParticle->McMainMCParticle());
            }
            
            else if (LArAnalysisParticleHelper::IsCosmicRay(pPfo))
            {
                this->AddCosmicRayRecord(*pAnalysisParticle);
                ++this->m_treeParameters.m_cr_Number;
                
                if (pAnalysisParticle->McMainMCParticle())
                    coveredMcPrimaries.insert(pAnalysisParticle->McMainMCParticle());
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
            if (coveredMcPrimaries.find(pMCPrimary) != coveredMcPrimaries.end())
                continue;
   
            if (!LArMCParticleHelper::IsNeutrino(pMCPrimary) && !LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) &&
                (pMCPrimary->GetParticleId() != MU_MINUS))
            {
                continue;
            }
                
            float mcEnergy = 0.f, mcContainmentFraction = 0.f;
            LArAnalysisParticle::TypeTree mcTypeTree;
            LArAnalysisParticle::TYPE mcType(LArAnalysisParticle::TYPE::UNKNOWN);
            CartesianVector mcVertexPosition(0.f, 0.f, 0.f), mcMomentum(0.f, 0.f, 0.f);
            int mcPdgCode(0);
            
            if (!LArAnalysisParticleHelper::GetMcInformation(pMCPrimary, mcEnergy, mcTypeTree, mcType, mcVertexPosition, mcMomentum, mcPdgCode, 
                mcContainmentFraction, this->m_minCoordinates, this->m_maxCoordinates))
            {
                std::cout << "WriteAnalysisParticlesAlgorithm: failed to get MC information for non-reconstructed MC particle" << std::endl;
                continue;
            }

            const bool mcIsVertexFiducial = LArAnalysisParticleHelper::IsPointFiducial(mcVertexPosition, this->m_minCoordinates,
                this->m_maxCoordinates);
                
            const bool mcIsContained = (mcContainmentFraction > this->m_mcContainmentFractionLowerBound);
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
                    mcIsVertexFiducial, mcIsContained, mcPdgCode, pMCParticleList, pCaloHitList, 0.f, 0.f);
            }
                
            else if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary))
            {
                ++this->m_treeParameters.m_primary_Number;
                this->AddMcOnlyPrimaryDaughterRecord(pMCPrimary, mcEnergy, mcVertexPosition, mcDirectionCosines, mcMomentum, 
                    mcIsVertexFiducial, mcIsContained, mcType, mcIsShower, mcPdgCode);
            }
            
            else if (pMCPrimary->GetParticleId() == MU_MINUS)
            {
                ++this->m_treeParameters.m_cr_Number;
                this->AddMcOnlyCosmicRayRecord(pMCPrimary, mcEnergy, mcVertexPosition, mcDirectionCosines, mcMomentum,
                    mcIsVertexFiducial, mcIsContained);
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
    this->m_treeParameters.m_nu_AreAllHitsFiducial          = neutrinoAnalysisParticle.AreAllHitsFiducial();
    this->m_treeParameters.m_nu_HasMcInfo                   = neutrinoAnalysisParticle.HasMcInfo();
    this->m_treeParameters.m_nu_Energy                      = neutrinoAnalysisParticle.AnalysisEnergy();
    this->m_treeParameters.m_nu_EnergyFromChargeOnly        = neutrinoAnalysisParticle.EnergyFromCharge();
    this->m_treeParameters.m_nu_VertexX                     = neutrinoAnalysisParticle.VertexPosition().GetX();
    this->m_treeParameters.m_nu_VertexY                     = neutrinoAnalysisParticle.VertexPosition().GetY();
    this->m_treeParameters.m_nu_VertexZ                     = neutrinoAnalysisParticle.VertexPosition().GetZ();
    this->m_treeParameters.m_nu_DirectionCosineX            = neutrinoAnalysisParticle.DirectionCosines().GetX();
    this->m_treeParameters.m_nu_DirectionCosineY            = neutrinoAnalysisParticle.DirectionCosines().GetY();
    this->m_treeParameters.m_nu_DirectionCosineZ            = neutrinoAnalysisParticle.DirectionCosines().GetZ();
    this->m_treeParameters.m_nu_MomentumX                   = neutrinoAnalysisParticle.AnalysisMomentum().GetX();
    this->m_treeParameters.m_nu_MomentumY                   = neutrinoAnalysisParticle.AnalysisMomentum().GetY();
    this->m_treeParameters.m_nu_MomentumZ                   = neutrinoAnalysisParticle.AnalysisMomentum().GetZ();
    this->m_treeParameters.m_nu_NumberOf3dHits              = neutrinoAnalysisParticle.NumberOf3dHits();
    this->m_treeParameters.m_nu_NumberOfCollectionPlaneHits = neutrinoAnalysisParticle.NumberOfCollectionPlaneHits();
    this->m_treeParameters.m_nu_NumberOfDownstreamParticles = neutrinoAnalysisParticle.NumberOfDownstreamParticles();
    
    // Derived parameters.
    
    // For other particles, we define the direction cosines by normalizing the momentum vector. In this case, since the neutrino isn't visible,
    // we define the neutrino momentum to be the sum of the momenta of its reconstructed daughters, and its direction to be along z.
    const CartesianVector zDirectionVector = CartesianVector{0.f, 0.f, 1.f};
    
    const float longitudinalEnergy = neutrinoAnalysisParticle.AnalysisMomentum().GetDotProduct(zDirectionVector);
    this->m_treeParameters.m_nu_LongitudinalEnergy = longitudinalEnergy;
    
    const CartesianVector transverseMomentum = neutrinoAnalysisParticle.AnalysisMomentum().GetCrossProduct(zDirectionVector);
    this->m_treeParameters.m_nu_TransverseEnergy = transverseMomentum.GetMagnitude();
    
    if (neutrinoAnalysisParticle.HasMcInfo())
    {
        this->PopulateNeutrinoMcParameters(neutrinoAnalysisParticle.McMainMCParticle(), neutrinoAnalysisParticle.McEnergy(),
            neutrinoAnalysisParticle.McVertexPosition(), neutrinoAnalysisParticle.McDirectionCosines(), 
            neutrinoAnalysisParticle.McMomentum(), neutrinoAnalysisParticle.McIsVertexFiducial(), neutrinoAnalysisParticle.McIsContained(), 
            neutrinoAnalysisParticle.McPdgCode(), pMCParticleList, pCaloHitList, neutrinoAnalysisParticle.McHitPurity(), 
            neutrinoAnalysisParticle.McHitCompleteness());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoMcParameters(const MCParticle *const pMainMcParticle, const float mcEnergy,
    const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum,
    const bool mcIsVertexFiducial, const bool mcIsContained, const int mcPdgCode, const MCParticleList *const pMCParticleList, 
    const CaloHitList *const pCaloHitList, const float mcHitPurity, const float mcHitCompleteness) const
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
    this->m_treeParameters.m_nu_mc_MomentumX        = mcMomentum.GetX();
    this->m_treeParameters.m_nu_mc_MomentumY        = mcMomentum.GetY();
    this->m_treeParameters.m_nu_mc_MomentumZ        = mcMomentum.GetZ();
    this->m_treeParameters.m_nu_mc_IsVertexFiducial = mcIsVertexFiducial;
    this->m_treeParameters.m_nu_mc_IsContained      = mcIsContained;
    this->m_treeParameters.m_nu_mc_PdgCode          = mcPdgCode;
    this->m_treeParameters.m_nu_mc_HitPurity        = mcHitPurity;
    this->m_treeParameters.m_nu_mc_HitCompleteness  = mcHitCompleteness;
    
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

        const CartesianVector &mcPrimaryMomentum(pMCPrimary->GetMomentum());
        
        if (mcPrimaryMomentum.GetMagnitude() > std::numeric_limits<float>::epsilon())
            visibleMomentum += mcPrimaryMomentum.GetUnitVector() * primaryVisibleEnergy;
    }
    
    const LArInteractionTypeHelper::InteractionType interactionType = this->GetInteractionType(pMCParticleList, pCaloHitList);
    
    // To align with the non-MC definition, we define the longitudinal/transverse components with respect to the z-direction.
    const CartesianVector zDirectionVector(0.f, 0.f, 1.f);

    this->m_treeParameters.m_nu_mc_LongitudinalEnergy        = mcMomentum.GetDotProduct(zDirectionVector);
    this->m_treeParameters.m_nu_mc_TransverseEnergy          = mcMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    this->m_treeParameters.m_nu_mc_VisibleEnergy             = visibleMomentum.GetMagnitude();
    this->m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy = visibleMomentum.GetZ();
    this->m_treeParameters.m_nu_mc_VisibleTransverseEnergy   = visibleMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    this->m_treeParameters.m_nu_mc_InteractionType           = static_cast<int>(interactionType);
    this->m_treeParameters.m_nu_mc_IsChargedCurrent          = this->IsChargedCurrent(interactionType);
    this->m_treeParameters.m_nu_mc_VisibleEnergyFraction     = visibleEnergy / mcEnergy;
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

void WriteAnalysisParticlesAlgorithm::AddPrimaryDaughterRecord(const LArAnalysisParticle &primaryAnalysisParticle) const
{
    // m_primary_Number is dealt with by the calling method.
    
    PUSH_TREE_RECORD(m_primary_WasReconstructed,            true);
    PUSH_TREE_RECORD(m_primary_IsVertexFiducial,            primaryAnalysisParticle.IsVertexFiducial());
    PUSH_TREE_RECORD(m_primary_AreAllHitsFiducial,          primaryAnalysisParticle.AreAllHitsFiducial());
    PUSH_TREE_RECORD(m_primary_HasMcInfo,                   primaryAnalysisParticle.HasMcInfo());
    PUSH_TREE_RECORD(m_primary_Energy,                      primaryAnalysisParticle.AnalysisEnergy());
    PUSH_TREE_RECORD(m_primary_EnergyFromChargeOnly,        primaryAnalysisParticle.EnergyFromCharge());
    PUSH_TREE_RECORD(m_primary_VertexX,                     primaryAnalysisParticle.VertexPosition().GetX());
    PUSH_TREE_RECORD(m_primary_VertexY,                     primaryAnalysisParticle.VertexPosition().GetY());
    PUSH_TREE_RECORD(m_primary_VertexZ,                     primaryAnalysisParticle.VertexPosition().GetZ());
    PUSH_TREE_RECORD(m_primary_DirectionCosineX,            primaryAnalysisParticle.DirectionCosines().GetX());
    PUSH_TREE_RECORD(m_primary_DirectionCosineY,            primaryAnalysisParticle.DirectionCosines().GetY());
    PUSH_TREE_RECORD(m_primary_DirectionCosineZ,            primaryAnalysisParticle.DirectionCosines().GetZ());
    PUSH_TREE_RECORD(m_primary_MomentumX,                   primaryAnalysisParticle.AnalysisMomentum().GetX());
    PUSH_TREE_RECORD(m_primary_MomentumY,                   primaryAnalysisParticle.AnalysisMomentum().GetY());
    PUSH_TREE_RECORD(m_primary_MomentumZ,                   primaryAnalysisParticle.AnalysisMomentum().GetZ());
    PUSH_TREE_RECORD(m_primary_ParticleType,                static_cast<int>(primaryAnalysisParticle.Type()));
    PUSH_TREE_RECORD(m_primary_NumberOf3dHits,              primaryAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_primary_NumberOfCollectionPlaneHits, primaryAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_primary_IsShower,                    primaryAnalysisParticle.IsShower());
    PUSH_TREE_RECORD(m_primary_NumberOfDownstreamParticles, primaryAnalysisParticle.NumberOfDownstreamParticles());
    
    if (primaryAnalysisParticle.HasMcInfo())
    {
        if (primaryAnalysisParticle.McMainMCParticle())
            PUSH_TREE_RECORD(m_primary_mc_McParticleUid, reinterpret_cast<std::uint64_t>(primaryAnalysisParticle.McMainMCParticle()->GetUid()));
        
        else
        {
            std::cout << "WriteAnalysisParticlesAlgorithm: primary neutrino daughter had MC info but no associated main MC particle" << std::endl;
            PUSH_TREE_RECORD(m_primary_mc_McParticleUid, 0ULL);
        }
        
        PUSH_TREE_RECORD(m_primary_mc_Energy,                   primaryAnalysisParticle.McEnergy());
        PUSH_TREE_RECORD(m_primary_mc_VertexX,                  primaryAnalysisParticle.McVertexPosition().GetX());
        PUSH_TREE_RECORD(m_primary_mc_VertexY,                  primaryAnalysisParticle.McVertexPosition().GetY());
        PUSH_TREE_RECORD(m_primary_mc_VertexZ,                  primaryAnalysisParticle.McVertexPosition().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,         primaryAnalysisParticle.McDirectionCosines().GetX());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,         primaryAnalysisParticle.McDirectionCosines().GetY());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,         primaryAnalysisParticle.McDirectionCosines().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_MomentumX,                primaryAnalysisParticle.McMomentum().GetX());
        PUSH_TREE_RECORD(m_primary_mc_MomentumY,                primaryAnalysisParticle.McMomentum().GetY());
        PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                primaryAnalysisParticle.McMomentum().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,         primaryAnalysisParticle.McIsVertexFiducial());
        PUSH_TREE_RECORD(m_primary_mc_IsContained,              primaryAnalysisParticle.McIsContained());
        PUSH_TREE_RECORD(m_primary_mc_ParticleType,             static_cast<int>(primaryAnalysisParticle.McType()));
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                 primaryAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  primaryAnalysisParticle.McPdgCode());
        PUSH_TREE_RECORD(m_primary_mc_HitPurity,                primaryAnalysisParticle.McHitPurity());
        PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,          primaryAnalysisParticle.McHitCompleteness());
    }
    
    else
    {
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid,            0ULL);
        PUSH_TREE_RECORD(m_primary_mc_Energy,                   0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexX,                  0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexY,                  0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexZ,                  0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,         0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,         0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,         0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumX,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumY,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,         false);
        PUSH_TREE_RECORD(m_primary_mc_IsContained,              false);
        PUSH_TREE_RECORD(m_primary_mc_ParticleType,             static_cast<int>(LArAnalysisParticle::TYPE::UNKNOWN));
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                 false);
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  0);
        PUSH_TREE_RECORD(m_primary_mc_HitPurity,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,          0.f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyPrimaryDaughterRecord(const MCParticle *const pMainMcParticle, const float mcEnergy,
    const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum,
    const bool mcIsVertexFiducial, const bool mcIsContained, const LArAnalysisParticle::TYPE mcType, const bool mcIsShower, 
    const int mcPdgCode) const
{
    PUSH_TREE_RECORD(m_primary_WasReconstructed,            false);
    PUSH_TREE_RECORD(m_primary_IsVertexFiducial,            false);
    PUSH_TREE_RECORD(m_primary_AreAllHitsFiducial,          false);
    PUSH_TREE_RECORD(m_primary_HasMcInfo,                   true);
    PUSH_TREE_RECORD(m_primary_Energy,                      0.f);
    PUSH_TREE_RECORD(m_primary_EnergyFromChargeOnly,        0.f);
    PUSH_TREE_RECORD(m_primary_VertexX,                     0.f);
    PUSH_TREE_RECORD(m_primary_VertexY,                     0.f);
    PUSH_TREE_RECORD(m_primary_VertexZ,                     0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineX,            0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineY,            0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineZ,            0.f);
    PUSH_TREE_RECORD(m_primary_MomentumX,                   0.f);
    PUSH_TREE_RECORD(m_primary_MomentumY,                   0.f);
    PUSH_TREE_RECORD(m_primary_MomentumZ,                   0.f);
    PUSH_TREE_RECORD(m_primary_ParticleType,                static_cast<int>(LArAnalysisParticle::TYPE::UNKNOWN));
    PUSH_TREE_RECORD(m_primary_NumberOf3dHits,              0U);
    PUSH_TREE_RECORD(m_primary_NumberOfCollectionPlaneHits, 0U);
    PUSH_TREE_RECORD(m_primary_IsShower,                    false);
    PUSH_TREE_RECORD(m_primary_NumberOfDownstreamParticles, 0U);
    
    if (pMainMcParticle)
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, reinterpret_cast<std::uint64_t>(pMainMcParticle->GetUid()));
    
    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: primary neutrino daughter had MC info but no associated main MC particle" << std::endl;
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, 0ULL);
    }
    
    PUSH_TREE_RECORD(m_primary_mc_Energy,                   mcEnergy);
    PUSH_TREE_RECORD(m_primary_mc_VertexX,                  mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_primary_mc_VertexY,                  mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_primary_mc_VertexZ,                  mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,         mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,         mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,         mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_MomentumX,                mcMomentum.GetX());
    PUSH_TREE_RECORD(m_primary_mc_MomentumY,                mcMomentum.GetY());
    PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,         mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_primary_mc_IsContained,              mcIsContained);
    PUSH_TREE_RECORD(m_primary_mc_ParticleType,             static_cast<int>(mcType));
    PUSH_TREE_RECORD(m_primary_mc_IsShower,                 mcIsShower);
    PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  mcPdgCode);
    PUSH_TREE_RECORD(m_primary_mc_HitPurity,                0.f);
    PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,          0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddCosmicRayRecord(const LArAnalysisParticle &cosmicRayAnalysisParticle) const
{
    // m_cr_Number is dealt with by the calling method.
    
    PUSH_TREE_RECORD(m_cr_WasReconstructed,            true);
    PUSH_TREE_RECORD(m_cr_IsVertexFiducial,            cosmicRayAnalysisParticle.IsVertexFiducial());
    PUSH_TREE_RECORD(m_cr_AreAllHitsFiducial,          cosmicRayAnalysisParticle.AreAllHitsFiducial());
    PUSH_TREE_RECORD(m_cr_HasMcInfo,                   cosmicRayAnalysisParticle.HasMcInfo());
    PUSH_TREE_RECORD(m_cr_Energy,                      cosmicRayAnalysisParticle.AnalysisEnergy());
    PUSH_TREE_RECORD(m_cr_EnergyFromChargeOnly,        cosmicRayAnalysisParticle.EnergyFromCharge());
    PUSH_TREE_RECORD(m_cr_VertexX,                     cosmicRayAnalysisParticle.VertexPosition().GetX());
    PUSH_TREE_RECORD(m_cr_VertexY,                     cosmicRayAnalysisParticle.VertexPosition().GetY());
    PUSH_TREE_RECORD(m_cr_VertexZ,                     cosmicRayAnalysisParticle.VertexPosition().GetZ());
    PUSH_TREE_RECORD(m_cr_DirectionCosineX,            cosmicRayAnalysisParticle.DirectionCosines().GetX());
    PUSH_TREE_RECORD(m_cr_DirectionCosineY,            cosmicRayAnalysisParticle.DirectionCosines().GetY());
    PUSH_TREE_RECORD(m_cr_DirectionCosineZ,            cosmicRayAnalysisParticle.DirectionCosines().GetZ());
    PUSH_TREE_RECORD(m_cr_MomentumX,                   cosmicRayAnalysisParticle.AnalysisMomentum().GetX());
    PUSH_TREE_RECORD(m_cr_MomentumY,                   cosmicRayAnalysisParticle.AnalysisMomentum().GetY());
    PUSH_TREE_RECORD(m_cr_MomentumZ,                   cosmicRayAnalysisParticle.AnalysisMomentum().GetZ());
    PUSH_TREE_RECORD(m_cr_NumberOf3dHits,              cosmicRayAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_cr_NumberOfCollectionPlaneHits, cosmicRayAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_cr_NumberOfDownstreamParticles, cosmicRayAnalysisParticle.NumberOfDownstreamParticles());
    
    if (cosmicRayAnalysisParticle.HasMcInfo())
    {
        if (cosmicRayAnalysisParticle.McMainMCParticle())
            PUSH_TREE_RECORD(m_cr_mc_McParticleUid, reinterpret_cast<std::uint64_t>(cosmicRayAnalysisParticle.McMainMCParticle()->GetUid()));
        
        else
        {
            std::cout << "WriteAnalysisParticlesAlgorithm: cosmic ray had MC info but no associated main MC particle" << std::endl;
            PUSH_TREE_RECORD(m_cr_mc_McParticleUid, 0ULL);
        }
        
        PUSH_TREE_RECORD(m_cr_mc_Energy,                   cosmicRayAnalysisParticle.McEnergy());
        PUSH_TREE_RECORD(m_cr_mc_VertexX,                  cosmicRayAnalysisParticle.McVertexPosition().GetX());
        PUSH_TREE_RECORD(m_cr_mc_VertexY,                  cosmicRayAnalysisParticle.McVertexPosition().GetY());
        PUSH_TREE_RECORD(m_cr_mc_VertexZ,                  cosmicRayAnalysisParticle.McVertexPosition().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,         cosmicRayAnalysisParticle.McDirectionCosines().GetX());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,         cosmicRayAnalysisParticle.McDirectionCosines().GetY());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,         cosmicRayAnalysisParticle.McDirectionCosines().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_MomentumX,                cosmicRayAnalysisParticle.McMomentum().GetX());
        PUSH_TREE_RECORD(m_cr_mc_MomentumY,                cosmicRayAnalysisParticle.McMomentum().GetY());
        PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                cosmicRayAnalysisParticle.McMomentum().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,         cosmicRayAnalysisParticle.McIsVertexFiducial());
        PUSH_TREE_RECORD(m_cr_mc_IsContained,              cosmicRayAnalysisParticle.McIsContained());
        PUSH_TREE_RECORD(m_cr_mc_HitPurity,                cosmicRayAnalysisParticle.McHitPurity());
        PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,          cosmicRayAnalysisParticle.McHitCompleteness());
    }
    
    else
    {
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid,            0ULL);
        PUSH_TREE_RECORD(m_cr_mc_Energy,                   0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexX,                  0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexY,                  0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexZ,                  0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,         0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,         0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,         0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumX,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumY,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,         false);
        PUSH_TREE_RECORD(m_cr_mc_IsContained,              false);
        PUSH_TREE_RECORD(m_cr_mc_HitPurity,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,          0.f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyCosmicRayRecord(const MCParticle *const pMainMcParticle, const float mcEnergy,
    const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum,
    const bool mcIsVertexFiducial, const bool mcIsContained) const
{
    // m_cr_Number is dealt with by the calling method.
    
    PUSH_TREE_RECORD(m_cr_WasReconstructed,            false);
    PUSH_TREE_RECORD(m_cr_IsVertexFiducial,            false);
    PUSH_TREE_RECORD(m_cr_AreAllHitsFiducial,          false);
    PUSH_TREE_RECORD(m_cr_HasMcInfo,                   true);
    PUSH_TREE_RECORD(m_cr_Energy,                      0.f);
    PUSH_TREE_RECORD(m_cr_EnergyFromChargeOnly,        0.f);
    PUSH_TREE_RECORD(m_cr_VertexX,                     0.f);
    PUSH_TREE_RECORD(m_cr_VertexY,                     0.f);
    PUSH_TREE_RECORD(m_cr_VertexZ,                     0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineX,            0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineY,            0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineZ,            0.f);
    PUSH_TREE_RECORD(m_cr_MomentumX,                   0.f);
    PUSH_TREE_RECORD(m_cr_MomentumY,                   0.f);
    PUSH_TREE_RECORD(m_cr_MomentumZ,                   0.f);
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
    
    PUSH_TREE_RECORD(m_cr_mc_Energy,                   mcEnergy);
    PUSH_TREE_RECORD(m_cr_mc_VertexX,                  mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_cr_mc_VertexY,                  mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_cr_mc_VertexZ,                  mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,         mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,         mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,         mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_MomentumX,                mcMomentum.GetX());
    PUSH_TREE_RECORD(m_cr_mc_MomentumY,                mcMomentum.GetY());
    PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,         mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_cr_mc_IsContained,              mcIsContained);
    PUSH_TREE_RECORD(m_cr_mc_HitPurity,                0.f);
    PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,          0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::DumpTree() const
{
    const std::string nuLabel = " - [nu]        ";
    
    std::cout << "Pandora Tree dump:\n";
    std::cout << nuLabel << "Was reconstructed:         " << std::boolalpha << this->m_treeParameters.m_nu_WasReconstructed << std::noboolalpha << '\n';
    std::cout << nuLabel << "Is vertex fiducial:        " << std::boolalpha << this->m_treeParameters.m_nu_IsVertexFiducial << std::noboolalpha << '\n';
    std::cout << nuLabel << "Are all hits fiducial:     " << std::boolalpha << this->m_treeParameters.m_nu_AreAllHitsFiducial << std::noboolalpha << '\n';
    std::cout << nuLabel << "Has MC info:               " << std::boolalpha << this->m_treeParameters.m_nu_HasMcInfo << std::noboolalpha << '\n';
    std::cout << nuLabel << "Energy:                    " << 1000.f * this->m_treeParameters.m_nu_Energy << " MeV\n";
    std::cout << nuLabel << "Energy from charge only:   " << 1000.f * this->m_treeParameters.m_nu_EnergyFromChargeOnly << " MeV\n";
    std::cout << nuLabel << "Longitudinal energy:       " << 1000.f * this->m_treeParameters.m_nu_LongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "Transverse energy:         " << 1000.f * this->m_treeParameters.m_nu_TransverseEnergy << " MeV\n";
    std::cout << nuLabel << "Vertex x:                  " << this->m_treeParameters.m_nu_VertexX << " cm\n";
    std::cout << nuLabel << "Vertex y:                  " << this->m_treeParameters.m_nu_VertexY << " cm\n";
    std::cout << nuLabel << "Vertex z:                  " << this->m_treeParameters.m_nu_VertexZ << " cm\n";
    std::cout << nuLabel << "Dir cosine x:              " << this->m_treeParameters.m_nu_DirectionCosineX << '\n';
    std::cout << nuLabel << "Dir cosine y:              " << this->m_treeParameters.m_nu_DirectionCosineY << '\n';
    std::cout << nuLabel << "Dir cosine z:              " << this->m_treeParameters.m_nu_DirectionCosineZ << '\n';
    std::cout << nuLabel << "Momentum x:                " << 1000.f * this->m_treeParameters.m_nu_MomentumX << " MeV/c\n";
    std::cout << nuLabel << "Momentum y:                " << 1000.f * this->m_treeParameters.m_nu_MomentumY << " MeV/c\n";
    std::cout << nuLabel << "Momentum z:                " << 1000.f * this->m_treeParameters.m_nu_MomentumZ << " MeV/c\n";
    std::cout << nuLabel << "Number of 3D hits:         " << this->m_treeParameters.m_nu_NumberOf3dHits << '\n';
    std::cout << nuLabel << "Number of coll plane hits: " << this->m_treeParameters.m_nu_NumberOfCollectionPlaneHits << '\n';
    std::cout << nuLabel << "Number of downstream pfos: " << this->m_treeParameters.m_nu_NumberOfDownstreamParticles << '\n';
    std::cout << nuLabel << "MC particle UID:           " << this->m_treeParameters.m_nu_mc_McParticleUid << '\n';
    std::cout << nuLabel << "MC energy:                 " << 1000.f * this->m_treeParameters.m_nu_mc_Energy << " MeV\n";
    std::cout << nuLabel << "MC longitudinal energy:    " << 1000.f * this->m_treeParameters.m_nu_mc_LongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "MC transverse energy:      " << 1000.f * this->m_treeParameters.m_nu_mc_TransverseEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible energy:         " << 1000.f * this->m_treeParameters.m_nu_mc_VisibleEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible long energy:    " << 1000.f * this->m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible trans energy:   " << 1000.f * this->m_treeParameters.m_nu_mc_VisibleTransverseEnergy << " MeV\n";
    std::cout << nuLabel << "MC vertex x:               " << this->m_treeParameters.m_nu_mc_VertexX << " cm\n";
    std::cout << nuLabel << "MC vertex y:               " << this->m_treeParameters.m_nu_mc_VertexY << " cm\n";
    std::cout << nuLabel << "MC vertex z:               " << this->m_treeParameters.m_nu_mc_VertexZ << " cm\n";
    std::cout << nuLabel << "MC dir cosine x:           " << this->m_treeParameters.m_nu_mc_DirectionCosineX << '\n';
    std::cout << nuLabel << "MC dir cosine y:           " << this->m_treeParameters.m_nu_mc_DirectionCosineY << '\n';
    std::cout << nuLabel << "MC dir cosine z:           " << this->m_treeParameters.m_nu_mc_DirectionCosineZ << '\n';
    std::cout << nuLabel << "MC momentum x:             " << 1000.f * this->m_treeParameters.m_nu_mc_MomentumX << " MeV/c\n";
    std::cout << nuLabel << "MC momentum y:             " << 1000.f * this->m_treeParameters.m_nu_mc_MomentumY << " MeV/c\n";
    std::cout << nuLabel << "MC momentum z:             " << 1000.f * this->m_treeParameters.m_nu_mc_MomentumZ << " MeV/c\n";
    std::cout << nuLabel << "MC is vertex fiducial:     " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsVertexFiducial << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC is contained:           " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsContained << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC interaction type:       " << this->m_treeParameters.m_nu_mc_InteractionType << '\n';
    std::cout << nuLabel << "MC is charged-current:     " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsChargedCurrent << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC visible energy frac:    " << 100.f * this->m_treeParameters.m_nu_mc_VisibleEnergyFraction << "%\n";
    std::cout << nuLabel << "MC PDG code:               " << this->m_treeParameters.m_nu_mc_PdgCode << '\n';
    std::cout << nuLabel << "MC hit purity:             " << 100.f * this->m_treeParameters.m_nu_mc_HitPurity << "%\n";
    std::cout << nuLabel << "MC hit completeness:       " << 100.f * this->m_treeParameters.m_nu_mc_HitCompleteness << "%\n";
    
    std::cout << " - [primary]   Number of primaries:       " << this->m_treeParameters.m_primary_Number << '\n';
    
    for (int i = 0; i < this->m_treeParameters.m_primary_Number; ++i)
    {
        const std::string label = " - [primary " + std::to_string(i) + "] ";
        
        std::cout << label << "Was reconstructed:         " << std::boolalpha << this->m_treeParameters.m_primary_WasReconstructed.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is vertex fiducial:        " << std::boolalpha << this->m_treeParameters.m_primary_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Are all hits fiducial:     " << std::boolalpha << this->m_treeParameters.m_primary_AreAllHitsFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Has MC info:               " << std::boolalpha << this->m_treeParameters.m_primary_HasMcInfo.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Energy:                    " << 1000.f * this->m_treeParameters.m_primary_Energy.at(i) << " MeV\n";
        std::cout << label << "Energy from charge only:   " << 1000.f * this->m_treeParameters.m_primary_EnergyFromChargeOnly.at(i) << " MeV\n";
        std::cout << label << "Vertex x:                  " << this->m_treeParameters.m_primary_VertexX.at(i) << " cm\n";
        std::cout << label << "Vertex y:                  " << this->m_treeParameters.m_primary_VertexY.at(i) << " cm\n";
        std::cout << label << "Vertex z:                  " << this->m_treeParameters.m_primary_VertexZ.at(i) << " cm\n";
        std::cout << label << "Dir cosine x:              " << this->m_treeParameters.m_primary_DirectionCosineX.at(i) << '\n';
        std::cout << label << "Dir cosine y:              " << this->m_treeParameters.m_primary_DirectionCosineY.at(i) << '\n';
        std::cout << label << "Dir cosine z:              " << this->m_treeParameters.m_primary_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "Momentum x:                " << 1000.f * this->m_treeParameters.m_primary_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "Momentum y:                " << 1000.f * this->m_treeParameters.m_primary_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "Momentum z:                " << 1000.f * this->m_treeParameters.m_primary_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "Particle type:             " << this->m_treeParameters.m_primary_ParticleType.at(i) << " (" << 
            LArAnalysisParticle::TypeAsString(static_cast<LArAnalysisParticle::TYPE>(this->m_treeParameters.m_primary_ParticleType.at(i))) << ")\n";
        std::cout << label << "Number of 3D hits:         " << this->m_treeParameters.m_primary_NumberOf3dHits.at(i) << '\n';
        std::cout << label << "Number of coll plane hits: " << this->m_treeParameters.m_primary_NumberOfCollectionPlaneHits.at(i) << '\n';
        std::cout << label << "Is shower:                 " << std::boolalpha << this->m_treeParameters.m_primary_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Number of downstream pfos: " << this->m_treeParameters.m_primary_NumberOfDownstreamParticles.at(i) << '\n';
        std::cout << label << "MC particle UID:           " << this->m_treeParameters.m_primary_mc_McParticleUid.at(i) << '\n';
        std::cout << label << "MC energy:                 " << 1000.f * this->m_treeParameters.m_primary_mc_Energy.at(i) << " MeV\n";
        std::cout << label << "MC vertex x:               " << this->m_treeParameters.m_primary_mc_VertexX.at(i) << " cm\n";
        std::cout << label << "MC vertex y:               " << this->m_treeParameters.m_primary_mc_VertexY.at(i) << " cm\n";
        std::cout << label << "MC vertex z:               " << this->m_treeParameters.m_primary_mc_VertexZ.at(i) << " cm\n";
        std::cout << label << "MC dir cosine x:           " << this->m_treeParameters.m_primary_mc_DirectionCosineX.at(i) << '\n';
        std::cout << label << "MC dir cosine y:           " << this->m_treeParameters.m_primary_mc_DirectionCosineY.at(i) << '\n';
        std::cout << label << "MC dir cosine z:           " << this->m_treeParameters.m_primary_mc_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "MC momentum x:             " << 1000.f * this->m_treeParameters.m_primary_mc_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum y:             " << 1000.f * this->m_treeParameters.m_primary_mc_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum z:             " << 1000.f * this->m_treeParameters.m_primary_mc_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "MC is vertex fiducial:     " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is contained:           " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC particle type:          " << this->m_treeParameters.m_primary_mc_ParticleType.at(i) << " (" << 
            LArAnalysisParticle::TypeAsString(static_cast<LArAnalysisParticle::TYPE>(this->m_treeParameters.m_primary_mc_ParticleType.at(i))) << ")\n";
        std::cout << label << "MC is shower:              " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC PDG code:               " << this->m_treeParameters.m_primary_mc_PdgCode.at(i) << '\n';
        std::cout << label << "MC hit purity:             " << 100.f * this->m_treeParameters.m_primary_mc_HitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness:       " << 100.f * this->m_treeParameters.m_primary_mc_HitCompleteness.at(i) << "%\n";
    }
    
    std::cout << " - [cr]        Number of cosmic rays:     " << this->m_treeParameters.m_cr_Number << '\n';
    
    for (int i = 0; i < this->m_treeParameters.m_cr_Number; ++i)
    {
        const std::string label = " - [cr " + std::to_string(i) + "]      ";
        
        std::cout << label << "Was reconstructed:         " << std::boolalpha << this->m_treeParameters.m_cr_WasReconstructed.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is vertex fiducial:        " << std::boolalpha << this->m_treeParameters.m_cr_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Are all hits fiducial:     " << std::boolalpha << this->m_treeParameters.m_cr_AreAllHitsFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Has MC info:               " << std::boolalpha << this->m_treeParameters.m_cr_HasMcInfo.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Energy:                    " << 1000.f * this->m_treeParameters.m_cr_Energy.at(i) << " MeV\n";
        std::cout << label << "Energy from charge only:   " << 1000.f * this->m_treeParameters.m_cr_EnergyFromChargeOnly.at(i) << " MeV\n";
        std::cout << label << "Vertex x:                  " << this->m_treeParameters.m_cr_VertexX.at(i) << " cm\n";
        std::cout << label << "Vertex y:                  " << this->m_treeParameters.m_cr_VertexY.at(i) << " cm\n";
        std::cout << label << "Vertex z:                  " << this->m_treeParameters.m_cr_VertexZ.at(i) << " cm\n";
        std::cout << label << "Dir cosine x:              " << this->m_treeParameters.m_cr_DirectionCosineX.at(i) << '\n';
        std::cout << label << "Dir cosine y:              " << this->m_treeParameters.m_cr_DirectionCosineY.at(i) << '\n';
        std::cout << label << "Dir cosine z:              " << this->m_treeParameters.m_cr_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "Momentum x:                " << 1000.f * this->m_treeParameters.m_cr_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "Momentum y:                " << 1000.f * this->m_treeParameters.m_cr_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "Momentum z:                " << 1000.f * this->m_treeParameters.m_cr_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "Number of 3D hits:         " << this->m_treeParameters.m_cr_NumberOf3dHits.at(i) << '\n';
        std::cout << label << "Number of coll plane hits: " << this->m_treeParameters.m_cr_NumberOfCollectionPlaneHits.at(i) << '\n';
        std::cout << label << "Number of downstream pfos: " << this->m_treeParameters.m_cr_NumberOfDownstreamParticles.at(i) << '\n';
        std::cout << label << "MC particle UID:           " << this->m_treeParameters.m_cr_mc_McParticleUid.at(i) << '\n';
        std::cout << label << "MC energy:                 " << 1000.f * this->m_treeParameters.m_cr_mc_Energy.at(i) << " MeV\n";
        std::cout << label << "MC vertex x:               " << this->m_treeParameters.m_cr_mc_VertexX.at(i) << " cm\n";
        std::cout << label << "MC vertex y:               " << this->m_treeParameters.m_cr_mc_VertexY.at(i) << " cm\n";
        std::cout << label << "MC vertex z:               " << this->m_treeParameters.m_cr_mc_VertexZ.at(i) << " cm\n";
        std::cout << label << "MC dir cosine x:           " << this->m_treeParameters.m_cr_mc_DirectionCosineX.at(i) << '\n';
        std::cout << label << "MC dir cosine y:           " << this->m_treeParameters.m_cr_mc_DirectionCosineY.at(i) << '\n';
        std::cout << label << "MC dir cosine z:           " << this->m_treeParameters.m_cr_mc_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "MC momentum x:             " << 1000.f * this->m_treeParameters.m_cr_mc_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum y:             " << 1000.f * this->m_treeParameters.m_cr_mc_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum z:             " << 1000.f * this->m_treeParameters.m_cr_mc_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "MC is vertex fiducial:     " << std::boolalpha << this->m_treeParameters.m_cr_mc_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is contained:           " << std::boolalpha << this->m_treeParameters.m_cr_mc_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC hit purity:             " << 100.f * this->m_treeParameters.m_cr_mc_HitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness:       " << 100.f * this->m_treeParameters.m_cr_mc_HitCompleteness.at(i) << "%\n";
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
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutXMargin", this->m_fiducialCutXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutYMargin", this->m_fiducialCutYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutZMargin", this->m_fiducialCutZMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McContainmentFractionLowerBound", this->m_mcContainmentFractionLowerBound));
    
    
    // Use the detector geometry and the margins to get the maximum and minimum fiducial volume coordinates.
    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), this->m_fiducialCutXMargin, this->m_fiducialCutYMargin,
        this->m_fiducialCutZMargin, this->m_minCoordinates, this->m_maxCoordinates);
        
    this->m_pOutputTFile = new TFile(this->m_outputFile.c_str(), "UPDATE");
    this->m_pOutputTree = new TTree("PandoraTree", "PandoraTree");
    
    // Neutrino parameters.
    this->m_pOutputTree->Branch("nu_WasReconstructed",                 &this->m_treeParameters.m_nu_WasReconstructed);
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
    this->m_pOutputTree->Branch("nu_mc_MomentumX",                     &this->m_treeParameters.m_nu_mc_MomentumX);
    this->m_pOutputTree->Branch("nu_mc_MomentumY",                     &this->m_treeParameters.m_nu_mc_MomentumY);
    this->m_pOutputTree->Branch("nu_mc_MomentumZ",                     &this->m_treeParameters.m_nu_mc_MomentumZ);
    this->m_pOutputTree->Branch("nu_mc_IsVertexFiducial",              &this->m_treeParameters.m_nu_mc_IsVertexFiducial);
    this->m_pOutputTree->Branch("nu_mc_IsContained",                   &this->m_treeParameters.m_nu_mc_IsContained);
    this->m_pOutputTree->Branch("nu_mc_InteractionType",               &this->m_treeParameters.m_nu_mc_InteractionType);
    this->m_pOutputTree->Branch("nu_mc_IsChargedCurrent",              &this->m_treeParameters.m_nu_mc_IsChargedCurrent);
    this->m_pOutputTree->Branch("nu_mc_VisibleEnergyFraction",         &this->m_treeParameters.m_nu_mc_VisibleEnergyFraction);
    this->m_pOutputTree->Branch("nu_mc_PdgCode",                       &this->m_treeParameters.m_nu_mc_PdgCode);
    this->m_pOutputTree->Branch("nu_mc_HitPurity",                     &this->m_treeParameters.m_nu_mc_HitPurity);
    this->m_pOutputTree->Branch("nu_mc_HitCompleteness",               &this->m_treeParameters.m_nu_mc_HitCompleteness);
    
    // Primary neutrino daughter parameters.
    this->m_pOutputTree->Branch("primary_Number",                      &this->m_treeParameters.m_primary_Number);
    
    this->m_pOutputTree->Branch("primary_WasReconstructed",            &this->m_treeParameters.m_primary_WasReconstructed);
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
    this->m_pOutputTree->Branch("primary_mc_McParticleUid",            &this->m_treeParameters.m_primary_mc_McParticleUid);
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
    this->m_pOutputTree->Branch("primary_mc_ParticleType",             &this->m_treeParameters.m_primary_mc_ParticleType);
    this->m_pOutputTree->Branch("primary_mc_IsShower",                 &this->m_treeParameters.m_primary_mc_IsShower);
    this->m_pOutputTree->Branch("primary_mc_PdgCode",                  &this->m_treeParameters.m_primary_mc_PdgCode);
    this->m_pOutputTree->Branch("primary_mc_HitPurity",                &this->m_treeParameters.m_primary_mc_HitPurity);
    this->m_pOutputTree->Branch("primary_mc_HitCompleteness",          &this->m_treeParameters.m_primary_mc_HitCompleteness);
    
    // Cosmic ray parameters.
    this->m_pOutputTree->Branch("cr_Number",                           &this->m_treeParameters.m_cr_Number);
    this->m_pOutputTree->Branch("cr_WasReconstructed",                 &this->m_treeParameters.m_cr_WasReconstructed);
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
    this->m_pOutputTree->Branch("cr_NumberOfDownstreamParticles",      &this->m_treeParameters.m_cr_NumberOfDownstreamParticles);
    this->m_pOutputTree->Branch("cr_mc_McParticleUid",                 &this->m_treeParameters.m_cr_mc_McParticleUid);
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
    this->m_pOutputTree->Branch("cr_mc_HitPurity",                     &this->m_treeParameters.m_cr_mc_HitPurity);
    this->m_pOutputTree->Branch("cr_mc_HitCompleteness",               &this->m_treeParameters.m_cr_mc_HitCompleteness);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TreeParameters::TreeParameters() noexcept :
    m_nu_WasReconstructed(false),
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
    m_nu_mc_MomentumX(0.f),
    m_nu_mc_MomentumY(0.f),
    m_nu_mc_MomentumZ(0.f),
    m_nu_mc_IsVertexFiducial(false),
    m_nu_mc_IsContained(false),
    m_nu_mc_InteractionType(0),
    m_nu_mc_IsChargedCurrent(false),
    m_nu_mc_VisibleEnergyFraction(0.f),
    m_nu_mc_PdgCode(0),
    m_nu_mc_HitPurity(0.f),
    m_nu_mc_HitCompleteness(0.f),
    m_primary_Number(0U),
    m_primary_WasReconstructed(),
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
    m_primary_mc_McParticleUid(),
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
    m_primary_mc_ParticleType(),
    m_primary_mc_IsShower(),
    m_primary_mc_PdgCode(),
    m_primary_mc_HitPurity(),
    m_primary_mc_HitCompleteness(),
    m_cr_Number(0U),
    m_cr_WasReconstructed(),
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
    m_cr_NumberOfDownstreamParticles(),
    m_cr_mc_McParticleUid(),
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
    m_cr_mc_HitPurity(),
    m_cr_mc_HitCompleteness()
{
}

} // namespace lar_physics_content
