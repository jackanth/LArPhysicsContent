/**
 *  @file LArPhysicsContent/src/WriteAnalysisParticlesAlgorithm.cxx
 *
 *  @brief Implementation of the write AnalysisParticles algorithm class.
 *
 *  $Log: $
 */

#include "WriteAnalysisParticlesAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmHeaders.h"

/**
 *  @brief  ...
 */
#define PUSH_TREE_RECORD(treeMember, value) this->m_treeParameters.treeMember.push_back(value);

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
    // Get the input PFO list.
    const PfoList *pPfoList(nullptr);

    if ((PandoraContentApi::GetList(*this, m_pfoListName, pPfoList) != STATUS_CODE_SUCCESS) || !pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "WriteAnalysisParticlesAlgorithm: cannot find pfo list " << m_pfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }
    
    this->m_treeParameters = TreeParameters();
    
    for (const ParticleFlowObject * const pPfo : *pPfoList)
    {
        if (const LArAnalysisParticle *const pAnalysisParticle = dynamic_cast<const LArAnalysisParticle *>(pPfo))
        {
            if (LArPfoHelper::IsNeutrino(pPfo))  
            {
                if (this->m_treeParameters.m_nu_Exists)
                {
                    std::cout << "WriteAnalysisParticlesAlgorithm: multiple neutrinos found - only recording one" << std::endl;
                    continue;
                }
                    
                this->PopulateNeutrinoParameters(*pAnalysisParticle);
                this->m_treeParameters.m_nu_Exists = true;
            }
                
            else if (!LArPfoHelper::GetParentPfo(pPfo)) // it's a cosmic ray
            {
                this->AddCosmicRayRecord(*pAnalysisParticle);
                ++this->m_treeParameters.m_cr_Number;
            }
            
            else if (LArPfoHelper::GetParentPfo(pPfo) == LArPfoHelper::GetParentNeutrino(pPfo)) // it's a primary daughter
            {
                std::cout << "Found primary" << std::endl;
                this->AddPrimaryDaughterRecord(*pAnalysisParticle);
                ++this->m_treeParameters.m_primary_Number;
            }
        }
    }
    
    if (this->m_verbose)
        this->PrintTree();
    
    this->m_pOutputTree->Fill();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoParameters(const LArAnalysisParticle &neutrinoAnalysisParticle) const
{
    // m_nu_Exists is dealt with by the calling method.
    
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
    
    if (neutrinoAnalysisParticle.HasMcInfo())
    {
        this->m_treeParameters.m_nu_mc_Energy                   = neutrinoAnalysisParticle.McEnergy();
        this->m_treeParameters.m_nu_mc_VertexX                  = neutrinoAnalysisParticle.McVertexPosition().GetX();
        this->m_treeParameters.m_nu_mc_VertexY                  = neutrinoAnalysisParticle.McVertexPosition().GetY();
        this->m_treeParameters.m_nu_mc_VertexZ                  = neutrinoAnalysisParticle.McVertexPosition().GetZ();
        this->m_treeParameters.m_nu_mc_DirectionCosineX         = neutrinoAnalysisParticle.McDirectionCosines().GetX();
        this->m_treeParameters.m_nu_mc_DirectionCosineY         = neutrinoAnalysisParticle.McDirectionCosines().GetY();
        this->m_treeParameters.m_nu_mc_DirectionCosineZ         = neutrinoAnalysisParticle.McDirectionCosines().GetZ();
        this->m_treeParameters.m_nu_mc_MomentumX                = neutrinoAnalysisParticle.McMomentum().GetX();
        this->m_treeParameters.m_nu_mc_MomentumY                = neutrinoAnalysisParticle.McMomentum().GetY();
        this->m_treeParameters.m_nu_mc_MomentumZ                = neutrinoAnalysisParticle.McMomentum().GetZ();
        this->m_treeParameters.m_nu_mc_IsVertexFiducial         = neutrinoAnalysisParticle.McIsVertexFiducial();
        this->m_treeParameters.m_nu_mc_IsContained              = neutrinoAnalysisParticle.McIsContained();
        this->m_treeParameters.m_nu_mc_IsCorrectlyReconstructed = neutrinoAnalysisParticle.McIsCorrectlyReconstructed();
        this->m_treeParameters.m_nu_mc_PdgCode                  = neutrinoAnalysisParticle.McPdgCode();
    }
    
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
        // To align with the non-MC definition, we define the longitudinal/transverse components with respect to the z-direction.
        const float mcLongitudinalEnergy = neutrinoAnalysisParticle.McMomentum().GetDotProduct(zDirectionVector);
        this->m_treeParameters.m_nu_mc_LongitudinalEnergy = mcLongitudinalEnergy;
        
        const CartesianVector mcTransverseMomentum = neutrinoAnalysisParticle.McMomentum().GetCrossProduct(zDirectionVector);
        this->m_treeParameters.m_nu_mc_TransverseEnergy = mcTransverseMomentum.GetMagnitude();
        
        // TODO the remaining parameters need stuff from the event validation algorithm, which is currently being refactored.
        this->m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy = 0.f;
        this->m_treeParameters.m_nu_mc_VisibleTransverseEnergy   = 0.f;
        this->m_treeParameters.m_nu_mc_InteractionType           = 0;
        this->m_treeParameters.m_nu_mc_IsChargedCurrent          = false;
        this->m_treeParameters.m_nu_mc_VisibleEnergyFraction     = 0.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddPrimaryDaughterRecord(const LArAnalysisParticle &primaryAnalysisParticle) const
{
    // m_primary_Number is dealt with by the calling method.
    
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
        PUSH_TREE_RECORD(m_primary_mc_IsCorrectlyReconstructed, primaryAnalysisParticle.McIsCorrectlyReconstructed());
        PUSH_TREE_RECORD(m_primary_mc_ParticleType,             static_cast<int>(primaryAnalysisParticle.McType()));
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                 primaryAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  primaryAnalysisParticle.McPdgCode());
    }
    
    else
    {
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
        PUSH_TREE_RECORD(m_primary_mc_IsCorrectlyReconstructed, false);
        PUSH_TREE_RECORD(m_primary_mc_ParticleType,             static_cast<int>(LArAnalysisParticle::TYPE::UNKNOWN));
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                 false);
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                  0);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddCosmicRayRecord(const LArAnalysisParticle &cosmicRayAnalysisParticle) const
{
    // m_cr_Number is dealt with by the calling method.
    
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
    
    if (cosmicRayAnalysisParticle.HasMcInfo())
    {
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
        PUSH_TREE_RECORD(m_cr_mc_IsCorrectlyReconstructed, cosmicRayAnalysisParticle.McIsCorrectlyReconstructed());
    }
    
    else
    {
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
        PUSH_TREE_RECORD(m_cr_mc_IsCorrectlyReconstructed, false);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PrintTree() const
{
    // TODO
    const std::string nuLabel = " - [nu]        ";
    
    std::cout << "Pandora Tree dump:\n";
    std::cout << nuLabel << "Exists:                    " << std::boolalpha << this->m_treeParameters.m_nu_Exists << std::noboolalpha << '\n';
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
    std::cout << nuLabel << "MC energy:                 " << 1000.f * this->m_treeParameters.m_nu_mc_Energy << " MeV\n";
    std::cout << nuLabel << "MC longitudinal energy:    " << 1000.f * this->m_treeParameters.m_nu_mc_LongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "MC transverse energy:      " << 1000.f * this->m_treeParameters.m_nu_mc_TransverseEnergy << " MeV\n";
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
    std::cout << nuLabel << "MC is reco correct:        " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsCorrectlyReconstructed << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC interaction type:       " << this->m_treeParameters.m_nu_mc_InteractionType << '\n';
    std::cout << nuLabel << "MC is charged-current:     " << std::boolalpha << this->m_treeParameters.m_nu_mc_IsChargedCurrent << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC visible energy frac:    " << 100.f * this->m_treeParameters.m_nu_mc_VisibleEnergyFraction << '\n';
    std::cout << nuLabel << "MC PDG code:               " << this->m_treeParameters.m_nu_mc_PdgCode << '\n';
    
    std::cout << " - [primary]   Number of primaries:       " << this->m_treeParameters.m_primary_Number << '\n';
    
    for (int i = 0; i < this->m_treeParameters.m_primary_Number; ++i)
    {
        const std::string label = " - [primary " + std::to_string(i) + "] ";
        
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
        std::cout << label << "Particle type:             " << this->m_treeParameters.m_primary_ParticleType.at(i) << '\n';
        std::cout << label << "Number of 3D hits:         " << this->m_treeParameters.m_primary_NumberOf3dHits.at(i) << '\n';
        std::cout << label << "Number of coll plane hits: " << this->m_treeParameters.m_primary_NumberOfCollectionPlaneHits.at(i) << '\n';
        std::cout << label << "Is shower:                 " << std::boolalpha << this->m_treeParameters.m_primary_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Number of downstream pfos: " << this->m_treeParameters.m_primary_NumberOfDownstreamParticles.at(i) << '\n';
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
        std::cout << label << "MC is reco correct:        " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsCorrectlyReconstructed.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC particle type:          " << this->m_treeParameters.m_primary_mc_ParticleType.at(i) << '\n';
        std::cout << label << "MC is shower:              " << std::boolalpha << this->m_treeParameters.m_primary_mc_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC PDG code:               " << this->m_treeParameters.m_primary_mc_PdgCode.at(i) << '\n';
    }
    
    std::cout << " - [cr]        Number of cosmic rays:     " << this->m_treeParameters.m_cr_Number << '\n';
    
    for (int i = 0; i < this->m_treeParameters.m_cr_Number; ++i)
    {
        const std::string label = " - [cr " + std::to_string(i) + "]      ";
        
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
        std::cout << label << "MC is reco correct:        " << std::boolalpha << this->m_treeParameters.m_cr_mc_IsCorrectlyReconstructed.at(i) << std::noboolalpha << '\n';
    }
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", this->m_pfoListName));
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
