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
 *  @brief  Macro for pushing a new value to a tree member vector
 */
#define PUSH_TREE_RECORD(treeMember, value) m_treeParameters.treeMember.push_back(value)

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
    if (m_pOutputTree)
        m_pOutputTree->Write();

    if (m_pOutputTFile)
    {
        if (m_pOutputTFile->IsOpen())
            m_pOutputTFile->Close();

        delete m_pOutputTFile;
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

    m_treeParameters = TreeParameters();

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
                if (m_treeParameters.m_nu_WasReconstructed)
                {
                    std::cout << "WriteAnalysisParticlesAlgorithm: multiple neutrinos found - only recording one" << std::endl;
                    continue;
                }

                this->PopulateNeutrinoParameters(*pAnalysisParticle, pMCParticleList);
                m_treeParameters.m_nu_WasReconstructed = true;
            }

            else if (LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(pPfo))
            {
                this->AddPrimaryDaughterRecord(*pAnalysisParticle, coveredMCPrimaries);
                ++m_treeParameters.m_primary_Number;
            }

            else if (LArAnalysisParticleHelper::IsCosmicRay(pPfo))
            {
                this->AddCosmicRayRecord(*pAnalysisParticle, coveredMCPrimaries);
                ++m_treeParameters.m_cr_Number;
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
                mcContainmentFraction, m_minCoordinates, m_maxCoordinates))
            {
                std::cout << "WriteAnalysisParticlesAlgorithm: failed to get MC information for non-reconstructed MC particle" << std::endl;
                continue;
            }

            if (mcEnergy < m_mcOnlyParticleEnergyCut)
                continue;

            if (mcContainmentFraction < m_mcOnlyParticleContainmentCut)
                continue;

            const bool mcIsVertexFiducial = LArAnalysisParticleHelper::IsPointFiducial(mcVertexPosition, m_minCoordinates,
                m_maxCoordinates);

            CartesianVector mcDirectionCosines(0.f, 0.f, 0.f);

            if (mcMomentum.GetMagnitude() < std::numeric_limits<float>::epsilon())
                std::cout << "AnalysisAlgorithm: could not get direction from MC momentum as it was too small" << std::endl;

            else
                mcDirectionCosines = mcMomentum.GetUnitVector();

            const bool mcIsShower = (mcType == LArAnalysisParticle::TYPE::SHOWER);

            if (LArMCParticleHelper::IsNeutrino(pMCPrimary))
            {
                if (m_treeParameters.m_nu_HasMcInfo)
                {
                    std::cout << "WriteAnalysisParticlesAlgorithm: found neutrino that was not covered by MC particles but neutrino MC properties already exist" << std::endl;
                    continue;
                }

                m_treeParameters.m_nu_HasMcInfo = true;

                this->PopulateNeutrinoMcParameters(pMCPrimary, mcEnergy, mcVertexPosition, mcDirectionCosines, mcMomentum,
                    mcIsVertexFiducial, mcContainmentFraction, mcPdgCode, pMCParticleList, 0.f, 0.f, 0.f, 0.f, mcTypeTree);
            }

            else if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary))
            {
                ++m_treeParameters.m_primary_Number;
                this->AddMcOnlyPrimaryDaughterRecord(pMCPrimary, mcEnergy, mcKineticEnergy, mcMass, mcVertexPosition, mcDirectionCosines, mcMomentum,
                    mcIsVertexFiducial, mcContainmentFraction, mcType, mcIsShower, mcPdgCode, mcTypeTree);
            }

            else if (pMCPrimary->GetParticleId() == MU_MINUS)
            {
                ++m_treeParameters.m_cr_Number;
                this->AddMcOnlyCosmicRayRecord(pMCPrimary, mcEnergy, mcKineticEnergy, mcMass, mcVertexPosition, mcDirectionCosines, mcMomentum,
                    mcIsVertexFiducial, mcContainmentFraction, mcType, mcIsShower, mcPdgCode, mcTypeTree);
            }
        }
    }

    if (m_verbose)
        this->DumpTree();

    m_pOutputTree->Fill();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoParameters(const LArAnalysisParticle &neutrinoAnalysisParticle,
    const MCParticleList *const pMCParticleList) const
{
    // m_nu_WasReconstructed is dealt with by the calling method.
    m_treeParameters.m_nu_IsVertexFiducial                            = neutrinoAnalysisParticle.IsVertexFiducial();
    m_treeParameters.m_nu_IsContained                                 = (neutrinoAnalysisParticle.FiducialHitFraction() >= m_fiducialHitFractionLowerBound);
    m_treeParameters.m_nu_FiducialHitFraction                         = neutrinoAnalysisParticle.FiducialHitFraction();
    m_treeParameters.m_nu_HasMcInfo                                   = neutrinoAnalysisParticle.HasMcInfo();
    m_treeParameters.m_nu_VisibleEnergy                               = neutrinoAnalysisParticle.KineticEnergy();
    m_treeParameters.m_nu_VisibleEnergyFracFromRange                  = neutrinoAnalysisParticle.KineticEnergyFromRangeFraction();
    m_treeParameters.m_nu_VisibleEnergyFracFromCorrectedTrackCharge   = neutrinoAnalysisParticle.KineticEnergyFromCorrectedTrackChargeFraction();
    m_treeParameters.m_nu_VisibleEnergyFracFromUncorrectedTrackCharge = neutrinoAnalysisParticle.KineticEnergyFromUncorrectedTrackChargeFraction();
    m_treeParameters.m_nu_VisibleEnergyFracFromShowerCharge           = neutrinoAnalysisParticle.KineticEnergyFromShowerChargeFraction();
    m_treeParameters.m_nu_VertexX                                     = neutrinoAnalysisParticle.VertexPosition().GetX();
    m_treeParameters.m_nu_VertexY                                     = neutrinoAnalysisParticle.VertexPosition().GetY();
    m_treeParameters.m_nu_VertexZ                                     = neutrinoAnalysisParticle.VertexPosition().GetZ();
    m_treeParameters.m_nu_DirectionCosineX                            = neutrinoAnalysisParticle.DirectionCosines().GetX();
    m_treeParameters.m_nu_DirectionCosineY                            = neutrinoAnalysisParticle.DirectionCosines().GetY();
    m_treeParameters.m_nu_DirectionCosineZ                            = neutrinoAnalysisParticle.DirectionCosines().GetZ();
    m_treeParameters.m_nu_TypeTree                                    = LArAnalysisParticleHelper::TypeTreeAsString(neutrinoAnalysisParticle.GetTypeTree());
    m_treeParameters.m_nu_NumberOf3dHits                              = neutrinoAnalysisParticle.NumberOf3dHits();
    m_treeParameters.m_nu_NumberOfCollectionPlaneHits                 = neutrinoAnalysisParticle.NumberOfCollectionPlaneHits();
    m_treeParameters.m_nu_NumberOfDownstreamParticles                 = neutrinoAnalysisParticle.NumberOfDownstreamParticles();

    CartesianVector visiblePseudoMomentum(0.f, 0.f, 0.f); // the summed reconstructed directional KEs of the primaries

    for (const ParticleFlowObject *const pDaughterPfo : neutrinoAnalysisParticle.GetDaughterPfoList())
    {
        if (const LArAnalysisParticle *const pDaughterAnalysisParticle = dynamic_cast<const LArAnalysisParticle *>(pDaughterPfo))
            visiblePseudoMomentum += pDaughterAnalysisParticle->DirectionCosines() * pDaughterAnalysisParticle->KineticEnergy();
    }

    const CartesianVector zAxis(0.f, 0.f, 1.f);
    m_treeParameters.m_nu_VisibleLongitudinalEnergy = visiblePseudoMomentum.GetDotProduct(zAxis);
    m_treeParameters.m_nu_VisibleTransverseEnergy   = visiblePseudoMomentum.GetCrossProduct(zAxis).GetMagnitude();

    if (neutrinoAnalysisParticle.HasMcInfo())
    {
        this->PopulateNeutrinoMcParameters(neutrinoAnalysisParticle.McMainMCParticle(), neutrinoAnalysisParticle.McEnergy(),
            neutrinoAnalysisParticle.McVertexPosition(), neutrinoAnalysisParticle.McDirectionCosines(),
            neutrinoAnalysisParticle.McMomentum(), neutrinoAnalysisParticle.McIsVertexFiducial(), neutrinoAnalysisParticle.McContainmentFraction(),
            neutrinoAnalysisParticle.McPdgCode(), pMCParticleList, neutrinoAnalysisParticle.McHitPurity(),
            neutrinoAnalysisParticle.McHitCompleteness(), neutrinoAnalysisParticle.McCollectionPlaneHitPurity(),
            neutrinoAnalysisParticle.McCollectionPlaneHitCompleteness(), neutrinoAnalysisParticle.GetMcTypeTree());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoMcParameters(const MCParticle *const pMainMcParticle, const float mcEnergy,
    const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines, const CartesianVector &mcMomentum,
    const bool mcIsVertexFiducial, const float mcContainmentFraction, const int mcPdgCode, const MCParticleList *const pMCParticleList,
    const float mcHitPurity, const float mcHitCompleteness, const float mcCollectionPlaneHitPurity,
    const float mcCollectionPlaneHitCompleteness, const LArAnalysisParticle::TypeTree mcTypeTree) const
{
    if (pMainMcParticle)
        m_treeParameters.m_nu_mc_McParticleUid = reinterpret_cast<std::uint64_t>(pMainMcParticle->GetUid());

    else // it can stay at its default value
        std::cout << "WriteAnalysisParticlesAlgorithm: could not populate MC particle UID" << std::endl;

    m_treeParameters.m_nu_mc_Energy                         = mcEnergy;
    m_treeParameters.m_nu_mc_VertexX                        = mcVertexPosition.GetX();
    m_treeParameters.m_nu_mc_VertexY                        = mcVertexPosition.GetY();
    m_treeParameters.m_nu_mc_VertexZ                        = mcVertexPosition.GetZ();
    m_treeParameters.m_nu_mc_DirectionCosineX               = mcDirectionCosines.GetX();
    m_treeParameters.m_nu_mc_DirectionCosineY               = mcDirectionCosines.GetY();
    m_treeParameters.m_nu_mc_DirectionCosineZ               = mcDirectionCosines.GetZ();
    m_treeParameters.m_nu_mc_Momentum                       = mcMomentum.GetMagnitude();
    m_treeParameters.m_nu_mc_MomentumX                      = mcMomentum.GetX();
    m_treeParameters.m_nu_mc_MomentumY                      = mcMomentum.GetY();
    m_treeParameters.m_nu_mc_MomentumZ                      = mcMomentum.GetZ();
    m_treeParameters.m_nu_mc_IsVertexFiducial               = mcIsVertexFiducial;
    m_treeParameters.m_nu_mc_IsContained                    = (mcContainmentFraction >= m_mcContainmentFractionLowerBound);
    m_treeParameters.m_nu_mc_ContainmentFraction            = mcContainmentFraction;
    m_treeParameters.m_nu_mc_TypeTree                       = LArAnalysisParticleHelper::TypeTreeAsString(mcTypeTree);
    m_treeParameters.m_nu_mc_PdgCode                        = mcPdgCode;
    m_treeParameters.m_nu_mc_HitPurity                      = mcHitPurity;
    m_treeParameters.m_nu_mc_HitCompleteness                = mcHitCompleteness;
    m_treeParameters.m_nu_mc_CollectionPlaneHitPurity       = mcCollectionPlaneHitPurity;
    m_treeParameters.m_nu_mc_CollectionPlaneHitCompleteness = mcCollectionPlaneHitCompleteness;

    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: no valid MC particle list name provided but neutrino analysis particle has MC info" << std::endl;
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

    const LArInteractionTypeHelper::InteractionType interactionType = this->GetInteractionType(pMCParticleList);

    // To align with the non-MC definition, we define the longitudinal/transverse components with respect to the z-direction.
    const CartesianVector zDirectionVector(0.f, 0.f, 1.f);

    m_treeParameters.m_nu_mc_LongitudinalEnergy        = mcMomentum.GetDotProduct(zDirectionVector);
    m_treeParameters.m_nu_mc_TransverseEnergy          = mcMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    m_treeParameters.m_nu_mc_VisibleEnergy             = visibleEnergy;
    m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy = visibleMomentum.GetDotProduct(zDirectionVector);
    m_treeParameters.m_nu_mc_VisibleTransverseEnergy   = visibleMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    m_treeParameters.m_nu_mc_InteractionType           = LArInteractionTypeHelper::ToString(interactionType);
    m_treeParameters.m_nu_mc_IsChargedCurrent          = this->IsChargedCurrent(interactionType);

    float visibleEnergyFraction(0.f);

    if (mcEnergy > std::numeric_limits<float>::epsilon())
        visibleEnergyFraction = visibleEnergy / mcEnergy;

    m_treeParameters.m_nu_mc_VisibleEnergyFraction = visibleEnergyFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionType WriteAnalysisParticlesAlgorithm::GetInteractionType(const MCParticleList *const pMCParticleList) const
{
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryVector);
    
    MCParticleList mcPrimaryList(mcPrimaryVector.begin(), mcPrimaryVector.end());
    return LArInteractionTypeHelper::GetInteractionType(mcPrimaryList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddPrimaryDaughterRecord(const LArAnalysisParticle &primaryAnalysisParticle,
    const MCPrimaryMap &coveredMCPrimaries) const
{
    // m_primary_Number is dealt with by the calling method.

    PUSH_TREE_RECORD(m_primary_WasReconstructed,                            true);
    PUSH_TREE_RECORD(m_primary_IsVertexFiducial,                            primaryAnalysisParticle.IsVertexFiducial());
    PUSH_TREE_RECORD(m_primary_IsContained,                                (primaryAnalysisParticle.FiducialHitFraction() >= m_fiducialHitFractionLowerBound));
    PUSH_TREE_RECORD(m_primary_FiducialHitFraction,                         primaryAnalysisParticle.FiducialHitFraction());
    PUSH_TREE_RECORD(m_primary_HasMcInfo,                                   primaryAnalysisParticle.HasMcInfo());
    PUSH_TREE_RECORD(m_primary_KineticEnergy,                               primaryAnalysisParticle.KineticEnergy());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromRange,                  primaryAnalysisParticle.KineticEnergyFromRangeFraction());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromCorrectedTrackCharge,   primaryAnalysisParticle.KineticEnergyFromCorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromUncorrectedTrackCharge, primaryAnalysisParticle.KineticEnergyFromUncorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromShowerCharge,           primaryAnalysisParticle.KineticEnergyFromShowerChargeFraction());
    PUSH_TREE_RECORD(m_primary_VertexX,                                     primaryAnalysisParticle.VertexPosition().GetX());
    PUSH_TREE_RECORD(m_primary_VertexY,                                     primaryAnalysisParticle.VertexPosition().GetY());
    PUSH_TREE_RECORD(m_primary_VertexZ,                                     primaryAnalysisParticle.VertexPosition().GetZ());
    PUSH_TREE_RECORD(m_primary_DirectionCosineX,                            primaryAnalysisParticle.DirectionCosines().GetX());
    PUSH_TREE_RECORD(m_primary_DirectionCosineY,                            primaryAnalysisParticle.DirectionCosines().GetY());
    PUSH_TREE_RECORD(m_primary_DirectionCosineZ,                            primaryAnalysisParticle.DirectionCosines().GetZ());
    PUSH_TREE_RECORD(m_primary_TypeTree,                                    LArAnalysisParticleHelper::TypeTreeAsString(primaryAnalysisParticle.GetTypeTree()));
    PUSH_TREE_RECORD(m_primary_IsShower,                                    primaryAnalysisParticle.IsShower());
    PUSH_TREE_RECORD(m_primary_IsTrack,                                    !primaryAnalysisParticle.IsShower());
    PUSH_TREE_RECORD(m_primary_IsProton,                                   (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::PROTON));
    PUSH_TREE_RECORD(m_primary_IsPionOrMuon,                               (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::PION_MUON) || (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_primary_NumberOf3dHits,                              primaryAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_primary_NumberOfCollectionPlaneHits,                 primaryAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_primary_NumberOfDownstreamParticles,                 primaryAnalysisParticle.NumberOfDownstreamParticles());

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

        PUSH_TREE_RECORD(m_primary_mc_IsParticleSplitByReco,          particleSplitByReco);
        PUSH_TREE_RECORD(m_primary_mc_Energy,                         primaryAnalysisParticle.McEnergy());
        PUSH_TREE_RECORD(m_primary_mc_KineticEnergy,                  primaryAnalysisParticle.McKineticEnergy());
        PUSH_TREE_RECORD(m_primary_mc_Mass,                           primaryAnalysisParticle.McMass());
        PUSH_TREE_RECORD(m_primary_mc_VertexX,                        primaryAnalysisParticle.McVertexPosition().GetX());
        PUSH_TREE_RECORD(m_primary_mc_VertexY,                        primaryAnalysisParticle.McVertexPosition().GetY());
        PUSH_TREE_RECORD(m_primary_mc_VertexZ,                        primaryAnalysisParticle.McVertexPosition().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,               primaryAnalysisParticle.McDirectionCosines().GetX());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,               primaryAnalysisParticle.McDirectionCosines().GetY());
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,               primaryAnalysisParticle.McDirectionCosines().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_Momentum,                       primaryAnalysisParticle.McMomentum().GetMagnitude());
        PUSH_TREE_RECORD(m_primary_mc_MomentumX,                      primaryAnalysisParticle.McMomentum().GetX());
        PUSH_TREE_RECORD(m_primary_mc_MomentumY,                      primaryAnalysisParticle.McMomentum().GetY());
        PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                      primaryAnalysisParticle.McMomentum().GetZ());
        PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,               primaryAnalysisParticle.McIsVertexFiducial());
        PUSH_TREE_RECORD(m_primary_mc_IsContained,                   (primaryAnalysisParticle.McContainmentFraction() >= m_mcContainmentFractionLowerBound));
        PUSH_TREE_RECORD(m_primary_mc_ContainmentFraction,            primaryAnalysisParticle.McContainmentFraction());
        PUSH_TREE_RECORD(m_primary_mc_TypeTree,                       LArAnalysisParticleHelper::TypeTreeAsString(primaryAnalysisParticle.GetMcTypeTree()));
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                       primaryAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_primary_mc_IsTrack,                       !primaryAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_primary_mc_IsProton,                      (primaryAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PROTON));
        PUSH_TREE_RECORD(m_primary_mc_IsPionOrMuon,                  (primaryAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PION_MUON));
        PUSH_TREE_RECORD(m_primary_mc_IsCosmicRay,                   (primaryAnalysisParticle.McType() == LArAnalysisParticle::TYPE::COSMIC_RAY));
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                        primaryAnalysisParticle.McPdgCode());
        PUSH_TREE_RECORD(m_primary_mc_HitPurity,                      primaryAnalysisParticle.McHitPurity());
        PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,                primaryAnalysisParticle.McHitCompleteness());
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitPurity,       primaryAnalysisParticle.McCollectionPlaneHitPurity());
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitCompleteness, primaryAnalysisParticle.McCollectionPlaneHitCompleteness());
    }

    else
    {
        PUSH_TREE_RECORD(m_primary_mc_IsParticleSplitByReco,          false);
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid,                  0ULL);
        PUSH_TREE_RECORD(m_primary_mc_Energy,                         0.f);
        PUSH_TREE_RECORD(m_primary_mc_KineticEnergy,                  0.f);
        PUSH_TREE_RECORD(m_primary_mc_Mass,                           0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexX,                        0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexY,                        0.f);
        PUSH_TREE_RECORD(m_primary_mc_VertexZ,                        0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,               0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,               0.f);
        PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,               0.f);
        PUSH_TREE_RECORD(m_primary_mc_Momentum,                       0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumX,                      0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumY,                      0.f);
        PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                      0.f);
        PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,               false);
        PUSH_TREE_RECORD(m_primary_mc_IsContained,                    false);
        PUSH_TREE_RECORD(m_primary_mc_ContainmentFraction,            0.f);
        PUSH_TREE_RECORD(m_primary_mc_TypeTree,                       "-");
        PUSH_TREE_RECORD(m_primary_mc_IsShower,                       false);
        PUSH_TREE_RECORD(m_primary_mc_IsTrack,                        false);
        PUSH_TREE_RECORD(m_primary_mc_IsProton,                       false);
        PUSH_TREE_RECORD(m_primary_mc_IsPionOrMuon,                   false);
        PUSH_TREE_RECORD(m_primary_mc_IsCosmicRay,                    false);
        PUSH_TREE_RECORD(m_primary_mc_PdgCode,                        0);
        PUSH_TREE_RECORD(m_primary_mc_HitPurity,                      0.f);
        PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,                0.f);
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitPurity,       0.f);
        PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitCompleteness, 0.f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyPrimaryDaughterRecord(const MCParticle *const pMainMcParticle, const float mcEnergy,
    const float mcKineticEnergy, const float mcMass, const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines,
    const CartesianVector &mcMomentum, const bool mcIsVertexFiducial, const float mcContainmentFraction,
    const LArAnalysisParticle::TYPE mcType, const bool mcIsShower, const int mcPdgCode,
    const LArAnalysisParticle::TypeTree mcTypeTree) const
{
    PUSH_TREE_RECORD(m_primary_WasReconstructed,                            false);
    PUSH_TREE_RECORD(m_primary_IsVertexFiducial,                            false);
    PUSH_TREE_RECORD(m_primary_IsContained,                                 false);
    PUSH_TREE_RECORD(m_primary_FiducialHitFraction,                         0.f);
    PUSH_TREE_RECORD(m_primary_HasMcInfo,                                   true);
    PUSH_TREE_RECORD(m_primary_KineticEnergy,                               0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromRange,                  0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromCorrectedTrackCharge,   0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromUncorrectedTrackCharge, 0.f);
    PUSH_TREE_RECORD(m_primary_KineticEnergyFracFromShowerCharge,           0.f);
    PUSH_TREE_RECORD(m_primary_VertexX,                                     0.f);
    PUSH_TREE_RECORD(m_primary_VertexY,                                     0.f);
    PUSH_TREE_RECORD(m_primary_VertexZ,                                     0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineX,                            0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineY,                            0.f);
    PUSH_TREE_RECORD(m_primary_DirectionCosineZ,                            0.f);
    PUSH_TREE_RECORD(m_primary_TypeTree,                                    "-");
    PUSH_TREE_RECORD(m_primary_IsShower,                                    false);
    PUSH_TREE_RECORD(m_primary_IsTrack,                                     false);
    PUSH_TREE_RECORD(m_primary_IsProton,                                    false);
    PUSH_TREE_RECORD(m_primary_IsPionOrMuon,                                false);
    PUSH_TREE_RECORD(m_primary_NumberOf3dHits,                              0U);
    PUSH_TREE_RECORD(m_primary_NumberOfCollectionPlaneHits,                 0U);
    PUSH_TREE_RECORD(m_primary_NumberOfDownstreamParticles,                 0U);

    if (pMainMcParticle)
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, reinterpret_cast<std::uint64_t>(pMainMcParticle->GetUid()));

    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: primary neutrino daughter had MC info but no associated main MC particle" << std::endl;
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, 0ULL);
    }

    PUSH_TREE_RECORD(m_primary_mc_IsParticleSplitByReco,          false); // not covered by any reconstructed particle
    PUSH_TREE_RECORD(m_primary_mc_Energy,                         mcEnergy);
    PUSH_TREE_RECORD(m_primary_mc_KineticEnergy,                  mcKineticEnergy);
    PUSH_TREE_RECORD(m_primary_mc_Mass,                           mcMass);
    PUSH_TREE_RECORD(m_primary_mc_VertexX,                        mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_primary_mc_VertexY,                        mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_primary_mc_VertexZ,                        mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,               mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,               mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,               mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_Momentum,                       mcMomentum.GetMagnitude());
    PUSH_TREE_RECORD(m_primary_mc_MomentumX,                      mcMomentum.GetX());
    PUSH_TREE_RECORD(m_primary_mc_MomentumY,                      mcMomentum.GetY());
    PUSH_TREE_RECORD(m_primary_mc_MomentumZ,                      mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,               mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_primary_mc_IsContained,                   (mcContainmentFraction >= m_mcContainmentFractionLowerBound));
    PUSH_TREE_RECORD(m_primary_mc_ContainmentFraction,            mcContainmentFraction);
    PUSH_TREE_RECORD(m_primary_mc_TypeTree,                       LArAnalysisParticleHelper::TypeTreeAsString(mcTypeTree));
    PUSH_TREE_RECORD(m_primary_mc_IsShower,                       mcIsShower);
    PUSH_TREE_RECORD(m_primary_mc_IsTrack,                       !mcIsShower);
    PUSH_TREE_RECORD(m_primary_mc_IsProton,                      (mcType == LArAnalysisParticle::TYPE::PROTON));
    PUSH_TREE_RECORD(m_primary_mc_IsPionOrMuon,                  (mcType == LArAnalysisParticle::TYPE::PION_MUON));
    PUSH_TREE_RECORD(m_primary_mc_IsCosmicRay,                   (mcType == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_primary_mc_PdgCode,                        mcPdgCode);
    PUSH_TREE_RECORD(m_primary_mc_HitPurity,                      0.f); // undefined
    PUSH_TREE_RECORD(m_primary_mc_HitCompleteness,                0.f); // undefined
    PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitPurity,       0.f); // undefined
    PUSH_TREE_RECORD(m_primary_mc_CollectionPlaneHitCompleteness, 0.f); // undefined
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddCosmicRayRecord(const LArAnalysisParticle &cosmicRayAnalysisParticle,
    const MCPrimaryMap &coveredMCPrimaries) const
{
    // m_cr_Number is dealt with by the calling method.

    PUSH_TREE_RECORD(m_cr_WasReconstructed,                            true);
    PUSH_TREE_RECORD(m_cr_IsVertexFiducial,                            cosmicRayAnalysisParticle.IsVertexFiducial());
    PUSH_TREE_RECORD(m_cr_IsContained,                                (cosmicRayAnalysisParticle.FiducialHitFraction() >= m_fiducialHitFractionLowerBound));
    PUSH_TREE_RECORD(m_cr_FiducialHitFraction,                         cosmicRayAnalysisParticle.FiducialHitFraction());
    PUSH_TREE_RECORD(m_cr_HasMcInfo,                                   cosmicRayAnalysisParticle.HasMcInfo());
    PUSH_TREE_RECORD(m_cr_KineticEnergy,                               cosmicRayAnalysisParticle.KineticEnergy());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromRange,                  cosmicRayAnalysisParticle.KineticEnergyFromRangeFraction());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromCorrectedTrackCharge,   cosmicRayAnalysisParticle.KineticEnergyFromCorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromUncorrectedTrackCharge, cosmicRayAnalysisParticle.KineticEnergyFromUncorrectedTrackChargeFraction());
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromShowerCharge,           cosmicRayAnalysisParticle.KineticEnergyFromShowerChargeFraction());
    PUSH_TREE_RECORD(m_cr_VertexX,                                     cosmicRayAnalysisParticle.VertexPosition().GetX());
    PUSH_TREE_RECORD(m_cr_VertexY,                                     cosmicRayAnalysisParticle.VertexPosition().GetY());
    PUSH_TREE_RECORD(m_cr_VertexZ,                                     cosmicRayAnalysisParticle.VertexPosition().GetZ());
    PUSH_TREE_RECORD(m_cr_DirectionCosineX,                            cosmicRayAnalysisParticle.DirectionCosines().GetX());
    PUSH_TREE_RECORD(m_cr_DirectionCosineY,                            cosmicRayAnalysisParticle.DirectionCosines().GetY());
    PUSH_TREE_RECORD(m_cr_DirectionCosineZ,                            cosmicRayAnalysisParticle.DirectionCosines().GetZ());
    PUSH_TREE_RECORD(m_cr_TypeTree,                                    LArAnalysisParticleHelper::TypeTreeAsString(cosmicRayAnalysisParticle.GetTypeTree()));
    PUSH_TREE_RECORD(m_cr_NumberOf3dHits,                              cosmicRayAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_cr_NumberOfCollectionPlaneHits,                 cosmicRayAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_cr_NumberOfDownstreamParticles,                 cosmicRayAnalysisParticle.NumberOfDownstreamParticles());

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

        PUSH_TREE_RECORD(m_cr_mc_IsParticleSplitByReco,          particleSplitByReco);
        PUSH_TREE_RECORD(m_cr_mc_Energy,                         cosmicRayAnalysisParticle.McEnergy());
        PUSH_TREE_RECORD(m_cr_mc_KineticEnergy,                  cosmicRayAnalysisParticle.McKineticEnergy());
        PUSH_TREE_RECORD(m_cr_mc_Mass,                           cosmicRayAnalysisParticle.McMass());
        PUSH_TREE_RECORD(m_cr_mc_VertexX,                        cosmicRayAnalysisParticle.McVertexPosition().GetX());
        PUSH_TREE_RECORD(m_cr_mc_VertexY,                        cosmicRayAnalysisParticle.McVertexPosition().GetY());
        PUSH_TREE_RECORD(m_cr_mc_VertexZ,                        cosmicRayAnalysisParticle.McVertexPosition().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,               cosmicRayAnalysisParticle.McDirectionCosines().GetX());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,               cosmicRayAnalysisParticle.McDirectionCosines().GetY());
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,               cosmicRayAnalysisParticle.McDirectionCosines().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_Momentum,                       cosmicRayAnalysisParticle.McMomentum().GetMagnitude());
        PUSH_TREE_RECORD(m_cr_mc_MomentumX,                      cosmicRayAnalysisParticle.McMomentum().GetX());
        PUSH_TREE_RECORD(m_cr_mc_MomentumY,                      cosmicRayAnalysisParticle.McMomentum().GetY());
        PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                      cosmicRayAnalysisParticle.McMomentum().GetZ());
        PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,               cosmicRayAnalysisParticle.McIsVertexFiducial());
        PUSH_TREE_RECORD(m_cr_mc_IsContained,                   (cosmicRayAnalysisParticle.McContainmentFraction() >= m_mcContainmentFractionLowerBound));
        PUSH_TREE_RECORD(m_cr_mc_ContainmentFraction,            cosmicRayAnalysisParticle.McContainmentFraction());
        PUSH_TREE_RECORD(m_cr_mc_TypeTree,                       LArAnalysisParticleHelper::TypeTreeAsString(cosmicRayAnalysisParticle.GetMcTypeTree()));
        PUSH_TREE_RECORD(m_cr_mc_IsShower,                       cosmicRayAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_cr_mc_IsTrack,                       !cosmicRayAnalysisParticle.McIsShower());
        PUSH_TREE_RECORD(m_cr_mc_IsProton,                      (cosmicRayAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PROTON));
        PUSH_TREE_RECORD(m_cr_mc_IsPionOrMuon,                  (cosmicRayAnalysisParticle.McType() == LArAnalysisParticle::TYPE::PION_MUON));
        PUSH_TREE_RECORD(m_cr_mc_IsCosmicRay,                   (cosmicRayAnalysisParticle.McType() == LArAnalysisParticle::TYPE::COSMIC_RAY));
        PUSH_TREE_RECORD(m_cr_mc_PdgCode,                        cosmicRayAnalysisParticle.McPdgCode());
        PUSH_TREE_RECORD(m_cr_mc_HitPurity,                      cosmicRayAnalysisParticle.McHitPurity());
        PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,                cosmicRayAnalysisParticle.McHitCompleteness());
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitPurity,       cosmicRayAnalysisParticle.McCollectionPlaneHitPurity());
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitCompleteness, cosmicRayAnalysisParticle.McCollectionPlaneHitCompleteness());
    }

    else
    {
        PUSH_TREE_RECORD(m_cr_mc_IsParticleSplitByReco,          false);
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid,                  0ULL);
        PUSH_TREE_RECORD(m_cr_mc_Energy,                         0.f);
        PUSH_TREE_RECORD(m_cr_mc_KineticEnergy,                  0.f);
        PUSH_TREE_RECORD(m_cr_mc_Mass,                           0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexX,                        0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexY,                        0.f);
        PUSH_TREE_RECORD(m_cr_mc_VertexZ,                        0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,               0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,               0.f);
        PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,               0.f);
        PUSH_TREE_RECORD(m_cr_mc_Momentum,                       0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumX,                      0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumY,                      0.f);
        PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                      0.f);
        PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,               false);
        PUSH_TREE_RECORD(m_cr_mc_IsContained,                    false);
        PUSH_TREE_RECORD(m_cr_mc_ContainmentFraction,            0.f);
        PUSH_TREE_RECORD(m_cr_mc_TypeTree,                       "-");
        PUSH_TREE_RECORD(m_cr_mc_IsShower,                       false);
        PUSH_TREE_RECORD(m_cr_mc_IsTrack,                        false);
        PUSH_TREE_RECORD(m_cr_mc_IsProton,                       false);
        PUSH_TREE_RECORD(m_cr_mc_IsPionOrMuon,                   false);
        PUSH_TREE_RECORD(m_cr_mc_IsCosmicRay,                    false);
        PUSH_TREE_RECORD(m_cr_mc_PdgCode,                        0);
        PUSH_TREE_RECORD(m_cr_mc_HitPurity,                      0.f);
        PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,                0.f);
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitPurity,       0.f);
        PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitCompleteness, 0.f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyCosmicRayRecord(const MCParticle *const pMainMcParticle, const float mcEnergy,
    const float mcKineticEnergy, const float mcMass, const CartesianVector &mcVertexPosition, const CartesianVector &mcDirectionCosines,
    const CartesianVector &mcMomentum, const bool mcIsVertexFiducial, const float mcContainmentFraction,
    const LArAnalysisParticle::TYPE mcType, const bool mcIsShower, const int mcPdgCode,
    const LArAnalysisParticle::TypeTree mcTypeTree) const
{
    // m_cr_Number is dealt with by the calling method.

    PUSH_TREE_RECORD(m_cr_WasReconstructed,                            false);
    PUSH_TREE_RECORD(m_cr_IsVertexFiducial,                            false);
    PUSH_TREE_RECORD(m_cr_IsContained,                                 false);
    PUSH_TREE_RECORD(m_cr_FiducialHitFraction,                         0.f);
    PUSH_TREE_RECORD(m_cr_HasMcInfo,                                   true);
    PUSH_TREE_RECORD(m_cr_KineticEnergy,                               0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromRange,                  0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromCorrectedTrackCharge,   0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromUncorrectedTrackCharge, 0.f);
    PUSH_TREE_RECORD(m_cr_KineticEnergyFracFromShowerCharge,           0.f);
    PUSH_TREE_RECORD(m_cr_VertexX,                                     0.f);
    PUSH_TREE_RECORD(m_cr_VertexY,                                     0.f);
    PUSH_TREE_RECORD(m_cr_VertexZ,                                     0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineX,                            0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineY,                            0.f);
    PUSH_TREE_RECORD(m_cr_DirectionCosineZ,                            0.f);
    PUSH_TREE_RECORD(m_cr_TypeTree,                                    "-");
    PUSH_TREE_RECORD(m_cr_NumberOf3dHits,                              0U);
    PUSH_TREE_RECORD(m_cr_NumberOfCollectionPlaneHits,                 0U);
    PUSH_TREE_RECORD(m_cr_NumberOfDownstreamParticles,                 0U);

    if (pMainMcParticle)
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid, reinterpret_cast<std::uint64_t>(pMainMcParticle->GetUid()));

    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: cosmic ray had MC info but no associated main MC particle" << std::endl;
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid, 0ULL);
    }

    PUSH_TREE_RECORD(m_cr_mc_IsParticleSplitByReco,          false); // not covered by any reconstructed particle
    PUSH_TREE_RECORD(m_cr_mc_Energy,                         mcEnergy);
    PUSH_TREE_RECORD(m_cr_mc_KineticEnergy,                  mcKineticEnergy);
    PUSH_TREE_RECORD(m_cr_mc_Mass,                           mcMass);
    PUSH_TREE_RECORD(m_cr_mc_VertexX,                        mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_cr_mc_VertexY,                        mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_cr_mc_VertexZ,                        mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,               mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,               mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,               mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_Momentum,                       mcMomentum.GetMagnitude());
    PUSH_TREE_RECORD(m_cr_mc_MomentumX,                      mcMomentum.GetX());
    PUSH_TREE_RECORD(m_cr_mc_MomentumY,                      mcMomentum.GetY());
    PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                      mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,               mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_cr_mc_IsContained,                   (mcContainmentFraction >= m_mcContainmentFractionLowerBound));
    PUSH_TREE_RECORD(m_cr_mc_ContainmentFraction,            mcContainmentFraction);
    PUSH_TREE_RECORD(m_cr_mc_TypeTree,                       LArAnalysisParticleHelper::TypeTreeAsString(mcTypeTree));
    PUSH_TREE_RECORD(m_cr_mc_IsShower,                       mcIsShower);
    PUSH_TREE_RECORD(m_cr_mc_IsTrack,                       !mcIsShower);
    PUSH_TREE_RECORD(m_cr_mc_IsProton,                      (mcType == LArAnalysisParticle::TYPE::PROTON));
    PUSH_TREE_RECORD(m_cr_mc_IsPionOrMuon,                  (mcType == LArAnalysisParticle::TYPE::PION_MUON) || (mcType == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_cr_mc_IsCosmicRay,                   (mcType == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_cr_mc_PdgCode,                        mcPdgCode);
    PUSH_TREE_RECORD(m_cr_mc_HitPurity,                      0.f); // undefined
    PUSH_TREE_RECORD(m_cr_mc_HitCompleteness,                0.f); // undefined
    PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitPurity,       0.f); // undefined
    PUSH_TREE_RECORD(m_cr_mc_CollectionPlaneHitCompleteness, 0.f); // undefined
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::DumpTree() const
{
    const std::string nuLabel = " - [nu]        ";

    std::cout << "Pandora Tree dump:\n";
    std::cout << nuLabel << "Was reconstructed:                " << std::boolalpha << m_treeParameters.m_nu_WasReconstructed << std::noboolalpha << '\n';
    std::cout << nuLabel << "Is vertex fiducial:               " << std::boolalpha << m_treeParameters.m_nu_IsVertexFiducial << std::noboolalpha << '\n';
    std::cout << nuLabel << "Is contained:                     " << std::boolalpha << m_treeParameters.m_nu_IsContained << std::noboolalpha << '\n';
    std::cout << nuLabel << "Fiducial hit fraction:            " << 100.f * m_treeParameters.m_nu_FiducialHitFraction << "%\n";
    std::cout << nuLabel << "Has MC info:                      " << std::boolalpha << m_treeParameters.m_nu_HasMcInfo << std::noboolalpha << '\n';
    std::cout << nuLabel << "Visible energy:                   " << 1000.f * m_treeParameters.m_nu_VisibleEnergy << " MeV\n";
    std::cout << nuLabel << "Longitudinal visible energy:      " << 1000.f * m_treeParameters.m_nu_VisibleLongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "Transverse visible energy:        " << 1000.f * m_treeParameters.m_nu_VisibleTransverseEnergy << " MeV\n";
    std::cout << nuLabel << "E frac from range:                " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromRange << "%\n";
    std::cout << nuLabel << "E frac from cor. track charge:    " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromCorrectedTrackCharge << "%\n";
    std::cout << nuLabel << "E frac from uncor. track charge:  " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromUncorrectedTrackCharge << "%\n";
    std::cout << nuLabel << "E frac from shower charge:        " << 100.f * m_treeParameters.m_nu_VisibleEnergyFracFromShowerCharge << "%\n";
    std::cout << nuLabel << "Vertex x:                         " << m_treeParameters.m_nu_VertexX << " cm\n";
    std::cout << nuLabel << "Vertex y:                         " << m_treeParameters.m_nu_VertexY << " cm\n";
    std::cout << nuLabel << "Vertex z:                         " << m_treeParameters.m_nu_VertexZ << " cm\n";
    std::cout << nuLabel << "Dir cosine x:                     " << m_treeParameters.m_nu_DirectionCosineX << '\n';
    std::cout << nuLabel << "Dir cosine y:                     " << m_treeParameters.m_nu_DirectionCosineY << '\n';
    std::cout << nuLabel << "Dir cosine z:                     " << m_treeParameters.m_nu_DirectionCosineZ << '\n';
    std::cout << nuLabel << "Type tree:                        " << m_treeParameters.m_nu_TypeTree << '\n';
    std::cout << nuLabel << "Number of 3D hits:                " << m_treeParameters.m_nu_NumberOf3dHits << '\n';
    std::cout << nuLabel << "Number of coll plane hits:        " << m_treeParameters.m_nu_NumberOfCollectionPlaneHits << '\n';
    std::cout << nuLabel << "Number of downstream pfos:        " << m_treeParameters.m_nu_NumberOfDownstreamParticles << '\n';
    std::cout << nuLabel << "MC particle UID:                  " << m_treeParameters.m_nu_mc_McParticleUid << '\n';
    std::cout << nuLabel << "MC energy:                        " << 1000.f * m_treeParameters.m_nu_mc_Energy << " MeV\n";
    std::cout << nuLabel << "MC longitudinal energy:           " << 1000.f * m_treeParameters.m_nu_mc_LongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "MC transverse energy:             " << 1000.f * m_treeParameters.m_nu_mc_TransverseEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible energy:                " << 1000.f * m_treeParameters.m_nu_mc_VisibleEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible long energy:           " << 1000.f * m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy << " MeV\n";
    std::cout << nuLabel << "MC visible trans energy:          " << 1000.f * m_treeParameters.m_nu_mc_VisibleTransverseEnergy << " MeV\n";
    std::cout << nuLabel << "MC vertex x:                      " << m_treeParameters.m_nu_mc_VertexX << " cm\n";
    std::cout << nuLabel << "MC vertex y:                      " << m_treeParameters.m_nu_mc_VertexY << " cm\n";
    std::cout << nuLabel << "MC vertex z:                      " << m_treeParameters.m_nu_mc_VertexZ << " cm\n";
    std::cout << nuLabel << "MC dir cosine x:                  " << m_treeParameters.m_nu_mc_DirectionCosineX << '\n';
    std::cout << nuLabel << "MC dir cosine y:                  " << m_treeParameters.m_nu_mc_DirectionCosineY << '\n';
    std::cout << nuLabel << "MC dir cosine z:                  " << m_treeParameters.m_nu_mc_DirectionCosineZ << '\n';
    std::cout << nuLabel << "MC momentum:                      " << 1000.f * m_treeParameters.m_nu_mc_Momentum << " MeV/c\n";
    std::cout << nuLabel << "MC momentum x:                    " << 1000.f * m_treeParameters.m_nu_mc_MomentumX << " MeV/c\n";
    std::cout << nuLabel << "MC momentum y:                    " << 1000.f * m_treeParameters.m_nu_mc_MomentumY << " MeV/c\n";
    std::cout << nuLabel << "MC momentum z:                    " << 1000.f * m_treeParameters.m_nu_mc_MomentumZ << " MeV/c\n";
    std::cout << nuLabel << "MC is vertex fiducial:            " << std::boolalpha << m_treeParameters.m_nu_mc_IsVertexFiducial << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC is contained:                  " << std::boolalpha << m_treeParameters.m_nu_mc_IsContained << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC containment fraction:          " << 100.f * m_treeParameters.m_nu_mc_ContainmentFraction << "%\n";
    std::cout << nuLabel << "MC type tree:                     " << m_treeParameters.m_nu_mc_TypeTree << '\n';
    std::cout << nuLabel << "MC interaction type:              " << m_treeParameters.m_nu_mc_InteractionType << '\n';
    std::cout << nuLabel << "MC is charged-current:            " << std::boolalpha << m_treeParameters.m_nu_mc_IsChargedCurrent << std::noboolalpha << '\n';
    std::cout << nuLabel << "MC visible energy frac:           " << 100.f * m_treeParameters.m_nu_mc_VisibleEnergyFraction << "%\n";
    std::cout << nuLabel << "MC PDG code:                      " << m_treeParameters.m_nu_mc_PdgCode << '\n';
    std::cout << nuLabel << "MC hit purity:                    " << 100.f * m_treeParameters.m_nu_mc_HitPurity << "%\n";
    std::cout << nuLabel << "MC hit completeness:              " << 100.f * m_treeParameters.m_nu_mc_HitCompleteness << "%\n";
    std::cout << nuLabel << "MC hit purity (coll plane):       " << 100.f * m_treeParameters.m_nu_mc_CollectionPlaneHitPurity << "%\n";
    std::cout << nuLabel << "MC hit completeness (coll plane): " << 100.f * m_treeParameters.m_nu_mc_CollectionPlaneHitCompleteness << "%\n";

    std::cout << " - [primary]   Number of primaries:              " << m_treeParameters.m_primary_Number << '\n';

    for (int i = 0; i < m_treeParameters.m_primary_Number; ++i)
    {
        const std::string label = " - [primary " + std::to_string(i) + "] ";

        std::cout << label << "Was reconstructed:                " << std::boolalpha << m_treeParameters.m_primary_WasReconstructed.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is vertex fiducial:               " << std::boolalpha << m_treeParameters.m_primary_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is contained:                     " << std::boolalpha << m_treeParameters.m_primary_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Fiducial hit fraction:            " << 100.f * m_treeParameters.m_primary_FiducialHitFraction.at(i) << "%\n";
        std::cout << label << "Has MC info:                      " << std::boolalpha << m_treeParameters.m_primary_HasMcInfo.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Kinetic energy:                   " << 1000.f * m_treeParameters.m_primary_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "KE frac from range:               " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromRange.at(i) << "%\n";
        std::cout << label << "KE frac from cor. track charge:   " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromCorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from uncor. track charge: " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromUncorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from shower charge:       " << 100.f * m_treeParameters.m_primary_KineticEnergyFracFromShowerCharge.at(i) << "%\n";
        std::cout << label << "Vertex x:                         " << m_treeParameters.m_primary_VertexX.at(i) << " cm\n";
        std::cout << label << "Vertex y:                         " << m_treeParameters.m_primary_VertexY.at(i) << " cm\n";
        std::cout << label << "Vertex z:                         " << m_treeParameters.m_primary_VertexZ.at(i) << " cm\n";
        std::cout << label << "Dir cosine x:                     " << m_treeParameters.m_primary_DirectionCosineX.at(i) << '\n';
        std::cout << label << "Dir cosine y:                     " << m_treeParameters.m_primary_DirectionCosineY.at(i) << '\n';
        std::cout << label << "Dir cosine z:                     " << m_treeParameters.m_primary_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "Type tree:                        " << m_treeParameters.m_primary_TypeTree.at(i) << '\n';
        std::cout << label << "Is shower:                        " << std::boolalpha << m_treeParameters.m_primary_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is track:                         " << std::boolalpha << m_treeParameters.m_primary_IsTrack.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is proton:                        " << std::boolalpha << m_treeParameters.m_primary_IsProton.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is pion or muon:                  " << std::boolalpha << m_treeParameters.m_primary_IsPionOrMuon.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Number of 3D hits:                " << m_treeParameters.m_primary_NumberOf3dHits.at(i) << '\n';
        std::cout << label << "Number of coll plane hits:        " << m_treeParameters.m_primary_NumberOfCollectionPlaneHits.at(i) << '\n';
        std::cout << label << "Number of downstream pfos:        " << m_treeParameters.m_primary_NumberOfDownstreamParticles.at(i) << '\n';
        std::cout << label << "MC particle UID:                  " << m_treeParameters.m_primary_mc_McParticleUid.at(i) << '\n';
        std::cout << label << "MC is particle split by reco:     " << std::boolalpha << m_treeParameters.m_primary_mc_IsParticleSplitByReco.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC energy:                        " << 1000.f * m_treeParameters.m_primary_mc_Energy.at(i) << " MeV\n";
        std::cout << label << "MC kinetic energy:                " << 1000.f * m_treeParameters.m_primary_mc_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "MC mass:                          " << 1000.f * m_treeParameters.m_primary_mc_Mass.at(i) << " MeV/c^2\n";
        std::cout << label << "MC vertex x:                      " << m_treeParameters.m_primary_mc_VertexX.at(i) << " cm\n";
        std::cout << label << "MC vertex y:                      " << m_treeParameters.m_primary_mc_VertexY.at(i) << " cm\n";
        std::cout << label << "MC vertex z:                      " << m_treeParameters.m_primary_mc_VertexZ.at(i) << " cm\n";
        std::cout << label << "MC dir cosine x:                  " << m_treeParameters.m_primary_mc_DirectionCosineX.at(i) << '\n';
        std::cout << label << "MC dir cosine y:                  " << m_treeParameters.m_primary_mc_DirectionCosineY.at(i) << '\n';
        std::cout << label << "MC dir cosine z:                  " << m_treeParameters.m_primary_mc_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "MC momentum:                      " << 1000.f * m_treeParameters.m_primary_mc_Momentum.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum x:                    " << 1000.f * m_treeParameters.m_primary_mc_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum y:                    " << 1000.f * m_treeParameters.m_primary_mc_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum z:                    " << 1000.f * m_treeParameters.m_primary_mc_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "MC is vertex fiducial:            " << std::boolalpha << m_treeParameters.m_primary_mc_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is contained:                  " << std::boolalpha << m_treeParameters.m_primary_mc_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC containment fraction:          " << 100.f * m_treeParameters.m_primary_mc_ContainmentFraction.at(i) << "%\n";
        std::cout << label << "MC type tree:                     " << m_treeParameters.m_primary_mc_TypeTree.at(i) << '\n';
        std::cout << label << "MC is shower:                     " << std::boolalpha << m_treeParameters.m_primary_mc_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is track:                      " << std::boolalpha << m_treeParameters.m_primary_mc_IsTrack.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is proton:                     " << std::boolalpha << m_treeParameters.m_primary_mc_IsProton.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is pion or muon:               " << std::boolalpha << m_treeParameters.m_primary_mc_IsPionOrMuon.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is cosmic ray:                 " << std::boolalpha << m_treeParameters.m_primary_mc_IsCosmicRay.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC PDG code:                      " << m_treeParameters.m_primary_mc_PdgCode.at(i) << '\n';
        std::cout << label << "MC hit purity:                    " << 100.f * m_treeParameters.m_primary_mc_HitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness:              " << 100.f * m_treeParameters.m_primary_mc_HitCompleteness.at(i) << "%\n";
        std::cout << label << "MC hit purity (coll plane):       " << 100.f * m_treeParameters.m_primary_mc_CollectionPlaneHitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness (coll plane): " << 100.f * m_treeParameters.m_primary_mc_CollectionPlaneHitCompleteness.at(i) << "%\n";
    }

    std::cout << " - [cr]        Number of cosmic rays:            " << m_treeParameters.m_cr_Number << '\n';

    for (int i = 0; i < m_treeParameters.m_cr_Number; ++i)
    {
        const std::string label = " - [cr " + std::to_string(i) + "]      ";

        std::cout << label << "Was reconstructed:                " << std::boolalpha << m_treeParameters.m_cr_WasReconstructed.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is vertex fiducial:               " << std::boolalpha << m_treeParameters.m_cr_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Is contained:                     " << std::boolalpha << m_treeParameters.m_cr_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Fiducial hit fraction:            " << 100.f * m_treeParameters.m_cr_FiducialHitFraction.at(i) << "%\n";
        std::cout << label << "Has MC info:                      " << std::boolalpha << m_treeParameters.m_cr_HasMcInfo.at(i) << std::noboolalpha << '\n';
        std::cout << label << "Kinetic energy:                   " << 1000.f * m_treeParameters.m_cr_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "KE frac from range:               " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromRange.at(i) << "%\n";
        std::cout << label << "KE frac from cor. track charge:   " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromCorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from uncor. track charge: " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromUncorrectedTrackCharge.at(i) << "%\n";
        std::cout << label << "KE frac from shower charge:       " << 100.f * m_treeParameters.m_cr_KineticEnergyFracFromShowerCharge.at(i) << "%\n";
        std::cout << label << "Vertex x:                         " << m_treeParameters.m_cr_VertexX.at(i) << " cm\n";
        std::cout << label << "Vertex y:                         " << m_treeParameters.m_cr_VertexY.at(i) << " cm\n";
        std::cout << label << "Vertex z:                         " << m_treeParameters.m_cr_VertexZ.at(i) << " cm\n";
        std::cout << label << "Dir cosine x:                     " << m_treeParameters.m_cr_DirectionCosineX.at(i) << '\n';
        std::cout << label << "Dir cosine y:                     " << m_treeParameters.m_cr_DirectionCosineY.at(i) << '\n';
        std::cout << label << "Dir cosine z:                     " << m_treeParameters.m_cr_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "Type tree:                        " << m_treeParameters.m_cr_TypeTree.at(i) << '\n';
        std::cout << label << "Number of 3D hits:                " << m_treeParameters.m_cr_NumberOf3dHits.at(i) << '\n';
        std::cout << label << "Number of coll plane hits:        " << m_treeParameters.m_cr_NumberOfCollectionPlaneHits.at(i) << '\n';
        std::cout << label << "Number of downstream pfos:        " << m_treeParameters.m_cr_NumberOfDownstreamParticles.at(i) << '\n';
        std::cout << label << "MC particle UID:                  " << m_treeParameters.m_cr_mc_McParticleUid.at(i) << '\n';
        std::cout << label << "MC is particle split by reco:     " << std::boolalpha << m_treeParameters.m_cr_mc_IsParticleSplitByReco.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC energy:                        " << 1000.f * m_treeParameters.m_cr_mc_Energy.at(i) << " MeV\n";
        std::cout << label << "MC kinetic energy:                " << 1000.f * m_treeParameters.m_cr_mc_KineticEnergy.at(i) << " MeV\n";
        std::cout << label << "MC mass:                          " << 1000.f * m_treeParameters.m_cr_mc_Mass.at(i) << " MeV/c^2\n";
        std::cout << label << "MC vertex x:                      " << m_treeParameters.m_cr_mc_VertexX.at(i) << " cm\n";
        std::cout << label << "MC vertex y:                      " << m_treeParameters.m_cr_mc_VertexY.at(i) << " cm\n";
        std::cout << label << "MC vertex z:                      " << m_treeParameters.m_cr_mc_VertexZ.at(i) << " cm\n";
        std::cout << label << "MC dir cosine x:                  " << m_treeParameters.m_cr_mc_DirectionCosineX.at(i) << '\n';
        std::cout << label << "MC dir cosine y:                  " << m_treeParameters.m_cr_mc_DirectionCosineY.at(i) << '\n';
        std::cout << label << "MC dir cosine z:                  " << m_treeParameters.m_cr_mc_DirectionCosineZ.at(i) << '\n';
        std::cout << label << "MC momentum:                      " << 1000.f * m_treeParameters.m_cr_mc_Momentum.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum x:                    " << 1000.f * m_treeParameters.m_cr_mc_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum y:                    " << 1000.f * m_treeParameters.m_cr_mc_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "MC momentum z:                    " << 1000.f * m_treeParameters.m_cr_mc_MomentumZ.at(i) << " MeV/c\n";
        std::cout << label << "MC is vertex fiducial:            " << std::boolalpha << m_treeParameters.m_cr_mc_IsVertexFiducial.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is contained:                  " << std::boolalpha << m_treeParameters.m_cr_mc_IsContained.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC containment fraction:          " << 100.f * m_treeParameters.m_cr_mc_ContainmentFraction.at(i) << "%\n";
        std::cout << label << "MC type tree:                     " << m_treeParameters.m_cr_mc_TypeTree.at(i) << '\n';
        std::cout << label << "MC is shower:                     " << std::boolalpha << m_treeParameters.m_cr_mc_IsShower.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is track:                      " << std::boolalpha << m_treeParameters.m_cr_mc_IsTrack.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is proton:                     " << std::boolalpha << m_treeParameters.m_cr_mc_IsProton.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is pion or muon:               " << std::boolalpha << m_treeParameters.m_cr_mc_IsPionOrMuon.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC is cosmic ray:                 " << std::boolalpha << m_treeParameters.m_cr_mc_IsCosmicRay.at(i) << std::noboolalpha << '\n';
        std::cout << label << "MC PDG code:                      " << m_treeParameters.m_cr_mc_PdgCode.at(i) << '\n';
        std::cout << label << "MC hit purity:                    " << 100.f * m_treeParameters.m_cr_mc_HitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness:              " << 100.f * m_treeParameters.m_cr_mc_HitCompleteness.at(i) << "%\n";
        std::cout << label << "MC hit purity (coll plane):       " << 100.f * m_treeParameters.m_cr_mc_CollectionPlaneHitPurity.at(i) << "%\n";
        std::cout << label << "MC hit completeness (coll plane): " << 100.f * m_treeParameters.m_cr_mc_CollectionPlaneHitCompleteness.at(i) << "%\n";
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_outputFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Verbose", m_verbose));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowXMargin", m_fiducialCutLowXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighXMargin", m_fiducialCutHighXMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowYMargin", m_fiducialCutLowYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighYMargin", m_fiducialCutHighYMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowZMargin", m_fiducialCutLowZMargin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighZMargin", m_fiducialCutHighZMargin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McContainmentFractionLowerBound", m_mcContainmentFractionLowerBound));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialHitFractionLowerBound", m_fiducialHitFractionLowerBound));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McOnlyParticleContainmentCut", m_mcOnlyParticleContainmentCut));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McOnlyParticleEnergyCut", m_mcOnlyParticleEnergyCut));

    // Use the detector geometry and the margins to get the maximum and minimum fiducial volume coordinates.
    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), m_fiducialCutLowXMargin, m_fiducialCutHighXMargin,
        m_fiducialCutLowYMargin, m_fiducialCutHighYMargin, m_fiducialCutLowZMargin, m_fiducialCutHighZMargin, m_minCoordinates, m_maxCoordinates);

    m_pOutputTFile = new TFile(m_outputFile.c_str(), "UPDATE");
    m_pOutputTree = new TTree("PandoraTree", "PandoraTree");

    // Neutrino parameters.
    m_pOutputTree->Branch("nu_WasReconstructed",                                 &m_treeParameters.m_nu_WasReconstructed);
    m_pOutputTree->Branch("nu_IsVertexFiducial",                                 &m_treeParameters.m_nu_IsVertexFiducial);
    m_pOutputTree->Branch("nu_IsContained",                                      &m_treeParameters.m_nu_IsContained);
    m_pOutputTree->Branch("nu_FiducialHitFraction",                              &m_treeParameters.m_nu_FiducialHitFraction);
    m_pOutputTree->Branch("nu_HasMcInfo",                                        &m_treeParameters.m_nu_HasMcInfo);
    m_pOutputTree->Branch("nu_VisibleEnergy",                                    &m_treeParameters.m_nu_VisibleEnergy);
    m_pOutputTree->Branch("nu_VisibleLongitudinalEnergy",                        &m_treeParameters.m_nu_VisibleLongitudinalEnergy);
    m_pOutputTree->Branch("nu_VisibleTransverseEnergy",                          &m_treeParameters.m_nu_VisibleTransverseEnergy);
    m_pOutputTree->Branch("nu_VisibleEnergyFracFromRange",                       &m_treeParameters.m_nu_VisibleEnergyFracFromRange);
    m_pOutputTree->Branch("nu_VisibleEnergyFracFromCorrectedTrackCharge",        &m_treeParameters.m_nu_VisibleEnergyFracFromCorrectedTrackCharge);
    m_pOutputTree->Branch("nu_VisibleEnergyFracFromUncorrectedTrackCharge",      &m_treeParameters.m_nu_VisibleEnergyFracFromUncorrectedTrackCharge);
    m_pOutputTree->Branch("nu_VisibleEnergyFracFromShowerCharge",                &m_treeParameters.m_nu_VisibleEnergyFracFromShowerCharge);
    m_pOutputTree->Branch("nu_VertexX",                                          &m_treeParameters.m_nu_VertexX);
    m_pOutputTree->Branch("nu_VertexY",                                          &m_treeParameters.m_nu_VertexY);
    m_pOutputTree->Branch("nu_VertexZ",                                          &m_treeParameters.m_nu_VertexZ);
    m_pOutputTree->Branch("nu_DirectionCosineX",                                 &m_treeParameters.m_nu_DirectionCosineX);
    m_pOutputTree->Branch("nu_DirectionCosineY",                                 &m_treeParameters.m_nu_DirectionCosineY);
    m_pOutputTree->Branch("nu_DirectionCosineZ",                                 &m_treeParameters.m_nu_DirectionCosineZ);
    m_pOutputTree->Branch("nu_TypeTree",                                         &m_treeParameters.m_nu_TypeTree);
    m_pOutputTree->Branch("nu_NumberOf3dHits",                                   &m_treeParameters.m_nu_NumberOf3dHits);
    m_pOutputTree->Branch("nu_NumberOfCollectionPlaneHits",                      &m_treeParameters.m_nu_NumberOfCollectionPlaneHits);
    m_pOutputTree->Branch("nu_NumberOfDownstreamParticles",                      &m_treeParameters.m_nu_NumberOfDownstreamParticles);
    m_pOutputTree->Branch("nu_mc_McParticleUid",                                 &m_treeParameters.m_nu_mc_McParticleUid);
    m_pOutputTree->Branch("nu_mc_Energy",                                        &m_treeParameters.m_nu_mc_Energy);
    m_pOutputTree->Branch("nu_mc_LongitudinalEnergy",                            &m_treeParameters.m_nu_mc_LongitudinalEnergy);
    m_pOutputTree->Branch("nu_mc_TransverseEnergy",                              &m_treeParameters.m_nu_mc_TransverseEnergy);
    m_pOutputTree->Branch("nu_mc_VisibleEnergy",                                 &m_treeParameters.m_nu_mc_VisibleEnergy);
    m_pOutputTree->Branch("nu_mc_VisibleLongitudinalEnergy",                     &m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy);
    m_pOutputTree->Branch("nu_mc_VisibleTransverseEnergy",                       &m_treeParameters.m_nu_mc_VisibleTransverseEnergy);
    m_pOutputTree->Branch("nu_mc_VertexX",                                       &m_treeParameters.m_nu_mc_VertexX);
    m_pOutputTree->Branch("nu_mc_VertexY",                                       &m_treeParameters.m_nu_mc_VertexY);
    m_pOutputTree->Branch("nu_mc_VertexZ",                                       &m_treeParameters.m_nu_mc_VertexZ);
    m_pOutputTree->Branch("nu_mc_DirectionCosineX",                              &m_treeParameters.m_nu_mc_DirectionCosineX);
    m_pOutputTree->Branch("nu_mc_DirectionCosineY",                              &m_treeParameters.m_nu_mc_DirectionCosineY);
    m_pOutputTree->Branch("nu_mc_DirectionCosineZ",                              &m_treeParameters.m_nu_mc_DirectionCosineZ);
    m_pOutputTree->Branch("nu_mc_Momentum",                                      &m_treeParameters.m_nu_mc_Momentum);
    m_pOutputTree->Branch("nu_mc_MomentumX",                                     &m_treeParameters.m_nu_mc_MomentumX);
    m_pOutputTree->Branch("nu_mc_MomentumY",                                     &m_treeParameters.m_nu_mc_MomentumY);
    m_pOutputTree->Branch("nu_mc_MomentumZ",                                     &m_treeParameters.m_nu_mc_MomentumZ);
    m_pOutputTree->Branch("nu_mc_IsVertexFiducial",                              &m_treeParameters.m_nu_mc_IsVertexFiducial);
    m_pOutputTree->Branch("nu_mc_IsContained",                                   &m_treeParameters.m_nu_mc_IsContained);
    m_pOutputTree->Branch("nu_mc_ContainmentFraction",                           &m_treeParameters.m_nu_mc_ContainmentFraction);
    m_pOutputTree->Branch("nu_mc_TypeTree",                                      &m_treeParameters.m_nu_mc_TypeTree);
    m_pOutputTree->Branch("nu_mc_InteractionType",                               &m_treeParameters.m_nu_mc_InteractionType);
    m_pOutputTree->Branch("nu_mc_IsChargedCurrent",                              &m_treeParameters.m_nu_mc_IsChargedCurrent);
    m_pOutputTree->Branch("nu_mc_VisibleEnergyFraction",                         &m_treeParameters.m_nu_mc_VisibleEnergyFraction);
    m_pOutputTree->Branch("nu_mc_PdgCode",                                       &m_treeParameters.m_nu_mc_PdgCode);
    m_pOutputTree->Branch("nu_mc_HitPurity",                                     &m_treeParameters.m_nu_mc_HitPurity);
    m_pOutputTree->Branch("nu_mc_HitCompleteness",                               &m_treeParameters.m_nu_mc_HitCompleteness);
    m_pOutputTree->Branch("nu_mc_CollectionPlaneHitPurity",                      &m_treeParameters.m_nu_mc_CollectionPlaneHitPurity);
    m_pOutputTree->Branch("nu_mc_CollectionPlaneHitCompleteness",                &m_treeParameters.m_nu_mc_CollectionPlaneHitCompleteness);

    // Primary neutrino daughter parameters.
    m_pOutputTree->Branch("primary_Number",                                      &m_treeParameters.m_primary_Number);
    m_pOutputTree->Branch("primary_WasReconstructed",                            &m_treeParameters.m_primary_WasReconstructed);
    m_pOutputTree->Branch("primary_IsVertexFiducial",                            &m_treeParameters.m_primary_IsVertexFiducial);
    m_pOutputTree->Branch("primary_IsContained",                                 &m_treeParameters.m_primary_IsContained);
    m_pOutputTree->Branch("primary_FiducialHitFraction",                         &m_treeParameters.m_primary_FiducialHitFraction);
    m_pOutputTree->Branch("primary_HasMcInfo",                                   &m_treeParameters.m_primary_HasMcInfo);
    m_pOutputTree->Branch("primary_KineticEnergy",                               &m_treeParameters.m_primary_KineticEnergy);
    m_pOutputTree->Branch("primary_KineticEnergyFracFromRange",                  &m_treeParameters.m_primary_KineticEnergyFracFromRange);
    m_pOutputTree->Branch("primary_KineticEnergyFracFromCorrectedTrackCharge",   &m_treeParameters.m_primary_KineticEnergyFracFromCorrectedTrackCharge);
    m_pOutputTree->Branch("primary_KineticEnergyFracFromUncorrectedTrackCharge", &m_treeParameters.m_primary_KineticEnergyFracFromUncorrectedTrackCharge);
    m_pOutputTree->Branch("primary_KineticEnergyFracFromShowerCharge",           &m_treeParameters.m_primary_KineticEnergyFracFromShowerCharge);
    m_pOutputTree->Branch("primary_VertexX",                                     &m_treeParameters.m_primary_VertexX);
    m_pOutputTree->Branch("primary_VertexY",                                     &m_treeParameters.m_primary_VertexY);
    m_pOutputTree->Branch("primary_VertexZ",                                     &m_treeParameters.m_primary_VertexZ);
    m_pOutputTree->Branch("primary_DirectionCosineX",                            &m_treeParameters.m_primary_DirectionCosineX);
    m_pOutputTree->Branch("primary_DirectionCosineY",                            &m_treeParameters.m_primary_DirectionCosineY);
    m_pOutputTree->Branch("primary_DirectionCosineZ",                            &m_treeParameters.m_primary_DirectionCosineZ);
    m_pOutputTree->Branch("primary_TypeTree",                                    &m_treeParameters.m_primary_TypeTree);
    m_pOutputTree->Branch("primary_IsShower",                                    &m_treeParameters.m_primary_IsShower);
    m_pOutputTree->Branch("primary_IsTrack",                                     &m_treeParameters.m_primary_IsTrack);
    m_pOutputTree->Branch("primary_IsProton",                                    &m_treeParameters.m_primary_IsProton);
    m_pOutputTree->Branch("primary_IsPionOrMuon",                                &m_treeParameters.m_primary_IsPionOrMuon);
    m_pOutputTree->Branch("primary_NumberOf3dHits",                              &m_treeParameters.m_primary_NumberOf3dHits);
    m_pOutputTree->Branch("primary_NumberOfCollectionPlaneHits",                 &m_treeParameters.m_primary_NumberOfCollectionPlaneHits);
    m_pOutputTree->Branch("primary_NumberOfDownstreamParticles",                 &m_treeParameters.m_primary_NumberOfDownstreamParticles);
    m_pOutputTree->Branch("primary_mc_McParticleUid",                            &m_treeParameters.m_primary_mc_McParticleUid);
    m_pOutputTree->Branch("primary_mc_IsParticleSplitByReco",                    &m_treeParameters.m_primary_mc_IsParticleSplitByReco);
    m_pOutputTree->Branch("primary_mc_Energy",                                   &m_treeParameters.m_primary_mc_Energy);
    m_pOutputTree->Branch("primary_mc_KineticEnergy",                            &m_treeParameters.m_primary_mc_KineticEnergy);
    m_pOutputTree->Branch("primary_mc_Mass",                                     &m_treeParameters.m_primary_mc_Mass);
    m_pOutputTree->Branch("primary_mc_VertexX",                                  &m_treeParameters.m_primary_mc_VertexX);
    m_pOutputTree->Branch("primary_mc_VertexY",                                  &m_treeParameters.m_primary_mc_VertexY);
    m_pOutputTree->Branch("primary_mc_VertexZ",                                  &m_treeParameters.m_primary_mc_VertexZ);
    m_pOutputTree->Branch("primary_mc_DirectionCosineX",                         &m_treeParameters.m_primary_mc_DirectionCosineX);
    m_pOutputTree->Branch("primary_mc_DirectionCosineY",                         &m_treeParameters.m_primary_mc_DirectionCosineY);
    m_pOutputTree->Branch("primary_mc_DirectionCosineZ",                         &m_treeParameters.m_primary_mc_DirectionCosineZ);
    m_pOutputTree->Branch("primary_mc_Momentum",                                 &m_treeParameters.m_primary_mc_Momentum);
    m_pOutputTree->Branch("primary_mc_MomentumX",                                &m_treeParameters.m_primary_mc_MomentumX);
    m_pOutputTree->Branch("primary_mc_MomentumY",                                &m_treeParameters.m_primary_mc_MomentumY);
    m_pOutputTree->Branch("primary_mc_MomentumZ",                                &m_treeParameters.m_primary_mc_MomentumZ);
    m_pOutputTree->Branch("primary_mc_IsVertexFiducial",                         &m_treeParameters.m_primary_mc_IsVertexFiducial);
    m_pOutputTree->Branch("primary_mc_IsContained",                              &m_treeParameters.m_primary_mc_IsContained);
    m_pOutputTree->Branch("primary_mc_ContainmentFraction",                      &m_treeParameters.m_primary_mc_ContainmentFraction);
    m_pOutputTree->Branch("primary_mc_TypeTree",                                 &m_treeParameters.m_primary_mc_TypeTree);
    m_pOutputTree->Branch("primary_mc_IsShower",                                 &m_treeParameters.m_primary_mc_IsShower);
    m_pOutputTree->Branch("primary_mc_IsTrack",                                  &m_treeParameters.m_primary_mc_IsTrack);
    m_pOutputTree->Branch("primary_mc_IsProton",                                 &m_treeParameters.m_primary_mc_IsProton);
    m_pOutputTree->Branch("primary_mc_IsPionOrMuon",                             &m_treeParameters.m_primary_mc_IsPionOrMuon);
    m_pOutputTree->Branch("primary_mc_IsCosmicRay",                              &m_treeParameters.m_primary_mc_IsCosmicRay);
    m_pOutputTree->Branch("primary_mc_PdgCode",                                  &m_treeParameters.m_primary_mc_PdgCode);
    m_pOutputTree->Branch("primary_mc_HitPurity",                                &m_treeParameters.m_primary_mc_HitPurity);
    m_pOutputTree->Branch("primary_mc_HitCompleteness",                          &m_treeParameters.m_primary_mc_HitCompleteness);
    m_pOutputTree->Branch("primary_mc_CollectionPlaneHitPurity",                 &m_treeParameters.m_primary_mc_CollectionPlaneHitPurity);
    m_pOutputTree->Branch("primary_mc_CollectionPlaneHitCompleteness",           &m_treeParameters.m_primary_mc_CollectionPlaneHitCompleteness);

    // Cosmic ray parameters.
    m_pOutputTree->Branch("cr_Number",                                           &m_treeParameters.m_cr_Number);
    m_pOutputTree->Branch("cr_WasReconstructed",                                 &m_treeParameters.m_cr_WasReconstructed);
    m_pOutputTree->Branch("cr_IsVertexFiducial",                                 &m_treeParameters.m_cr_IsVertexFiducial);
    m_pOutputTree->Branch("cr_IsContained",                                      &m_treeParameters.m_cr_IsContained);
    m_pOutputTree->Branch("cr_FiducialHitFraction",                              &m_treeParameters.m_cr_FiducialHitFraction);
    m_pOutputTree->Branch("cr_HasMcInfo",                                        &m_treeParameters.m_cr_HasMcInfo);
    m_pOutputTree->Branch("cr_KineticEnergy",                                    &m_treeParameters.m_cr_KineticEnergy);
    m_pOutputTree->Branch("cr_KineticEnergyFracFromRange",                       &m_treeParameters.m_cr_KineticEnergyFracFromRange);
    m_pOutputTree->Branch("cr_KineticEnergyFracFromCorrectedTrackCharge",        &m_treeParameters.m_cr_KineticEnergyFracFromCorrectedTrackCharge);
    m_pOutputTree->Branch("cr_KineticEnergyFracFromUncorrectedTrackCharge",      &m_treeParameters.m_cr_KineticEnergyFracFromUncorrectedTrackCharge);
    m_pOutputTree->Branch("cr_KineticEnergyFracFromShowerCharge",                &m_treeParameters.m_cr_KineticEnergyFracFromShowerCharge);
    m_pOutputTree->Branch("cr_VertexX",                                          &m_treeParameters.m_cr_VertexX);
    m_pOutputTree->Branch("cr_VertexY",                                          &m_treeParameters.m_cr_VertexY);
    m_pOutputTree->Branch("cr_VertexZ",                                          &m_treeParameters.m_cr_VertexZ);
    m_pOutputTree->Branch("cr_DirectionCosineX",                                 &m_treeParameters.m_cr_DirectionCosineX);
    m_pOutputTree->Branch("cr_DirectionCosineY",                                 &m_treeParameters.m_cr_DirectionCosineY);
    m_pOutputTree->Branch("cr_DirectionCosineZ",                                 &m_treeParameters.m_cr_DirectionCosineZ);
    m_pOutputTree->Branch("cr_TypeTree",                                         &m_treeParameters.m_cr_TypeTree);
    m_pOutputTree->Branch("cr_NumberOf3dHits",                                   &m_treeParameters.m_cr_NumberOf3dHits);
    m_pOutputTree->Branch("cr_NumberOfCollectionPlaneHits",                      &m_treeParameters.m_cr_NumberOfCollectionPlaneHits);
    m_pOutputTree->Branch("cr_NumberOfDownstreamParticles",                      &m_treeParameters.m_cr_NumberOfDownstreamParticles);
    m_pOutputTree->Branch("cr_mc_McParticleUid",                                 &m_treeParameters.m_cr_mc_McParticleUid);
    m_pOutputTree->Branch("cr_mc_IsParticleSplitByReco",                         &m_treeParameters.m_cr_mc_IsParticleSplitByReco);
    m_pOutputTree->Branch("cr_mc_Energy",                                        &m_treeParameters.m_cr_mc_Energy);
    m_pOutputTree->Branch("cr_mc_KineticEnergy",                                 &m_treeParameters.m_cr_mc_KineticEnergy);
    m_pOutputTree->Branch("cr_mc_Mass",                                          &m_treeParameters.m_cr_mc_Mass);
    m_pOutputTree->Branch("cr_mc_VertexX",                                       &m_treeParameters.m_cr_mc_VertexX);
    m_pOutputTree->Branch("cr_mc_VertexY",                                       &m_treeParameters.m_cr_mc_VertexY);
    m_pOutputTree->Branch("cr_mc_VertexZ",                                       &m_treeParameters.m_cr_mc_VertexZ);
    m_pOutputTree->Branch("cr_mc_DirectionCosineX",                              &m_treeParameters.m_cr_mc_DirectionCosineX);
    m_pOutputTree->Branch("cr_mc_DirectionCosineY",                              &m_treeParameters.m_cr_mc_DirectionCosineY);
    m_pOutputTree->Branch("cr_mc_DirectionCosineZ",                              &m_treeParameters.m_cr_mc_DirectionCosineZ);
    m_pOutputTree->Branch("cr_mc_Momentum",                                      &m_treeParameters.m_cr_mc_Momentum);
    m_pOutputTree->Branch("cr_mc_MomentumX",                                     &m_treeParameters.m_cr_mc_MomentumX);
    m_pOutputTree->Branch("cr_mc_MomentumY",                                     &m_treeParameters.m_cr_mc_MomentumY);
    m_pOutputTree->Branch("cr_mc_MomentumZ",                                     &m_treeParameters.m_cr_mc_MomentumZ);
    m_pOutputTree->Branch("cr_mc_IsVertexFiducial",                              &m_treeParameters.m_cr_mc_IsVertexFiducial);
    m_pOutputTree->Branch("cr_mc_IsContained",                                   &m_treeParameters.m_cr_mc_IsContained);
    m_pOutputTree->Branch("cr_mc_ContainmentFraction",                           &m_treeParameters.m_cr_mc_ContainmentFraction);
    m_pOutputTree->Branch("cr_mc_TypeTree",                                      &m_treeParameters.m_cr_mc_TypeTree);
    m_pOutputTree->Branch("cr_mc_IsShower",                                      &m_treeParameters.m_cr_mc_IsShower);
    m_pOutputTree->Branch("cr_mc_IsTrack",                                       &m_treeParameters.m_cr_mc_IsTrack);
    m_pOutputTree->Branch("cr_mc_IsProton",                                      &m_treeParameters.m_cr_mc_IsProton);
    m_pOutputTree->Branch("cr_mc_IsPionOrMuon",                                  &m_treeParameters.m_cr_mc_IsPionOrMuon);
    m_pOutputTree->Branch("cr_mc_IsCosmicRay",                                   &m_treeParameters.m_cr_mc_IsCosmicRay);
    m_pOutputTree->Branch("cr_mc_PdgCode",                                       &m_treeParameters.m_cr_mc_PdgCode);
    m_pOutputTree->Branch("cr_mc_HitPurity",                                     &m_treeParameters.m_cr_mc_HitPurity);
    m_pOutputTree->Branch("cr_mc_HitCompleteness",                               &m_treeParameters.m_cr_mc_HitCompleteness);
    m_pOutputTree->Branch("cr_mc_CollectionPlaneHitPurity",                      &m_treeParameters.m_cr_mc_CollectionPlaneHitPurity);
    m_pOutputTree->Branch("cr_mc_CollectionPlaneHitCompleteness",                &m_treeParameters.m_cr_mc_CollectionPlaneHitCompleteness);

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
