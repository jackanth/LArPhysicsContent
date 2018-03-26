/**
 *  @file   larphysicscontent/WriteAnalysisParticlesAlgorithm.cc
 *
 *  @brief  Implementation of the write AnalysisParticles algorithm class.
 *
 *  $Log: $
 */
#include "larphysicscontent/WriteAnalysisParticlesAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

WriteAnalysisParticlesAlgorithm::WriteAnalysisParticlesAlgorithm() :
    m_pfoListName(),
    m_outputFile(),
    m_treeName("PandoraTree"),
    m_treeTitle("PandoraTree"),
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
    m_mcOnlyParticleContainmentCut(0.f),
    m_mcOnlyParticleEnergyCut(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

WriteAnalysisParticlesAlgorithm::~WriteAnalysisParticlesAlgorithm()
{
    // Write any remaining data in the TTree.
    if (m_pOutputTree)
        m_pOutputTree->Write();

    // If we have a TFile, close and delete it. The TTree is then deleted by its owning TFile.
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

    // Construct a multimap from MCParticles to AnalysisParticles and use this to fill the is-split boolean later.
    const MCPrimaryMap mainMcParticleMap = this->GetMainMcParticleMap(*pPfoList);

    for (const ParticleFlowObject * const pPfo : *pPfoList)
    {
        if (const LArAnalysisParticle *const pAnalysisParticle = dynamic_cast<const LArAnalysisParticle *>(pPfo))
            this->ProcessAnalysisParticle(pAnalysisParticle, mainMcParticleMap, pMCParticleList);
    }

    // Use the multimap entries to add records for unreconstructed particles.
    if (pMCParticleList)
        this->RecordUnreconstructedParticles(pMCParticleList, mainMcParticleMap);

    this->CheckTreeVectorSizes();

    if (m_verbose)
        this->PrintTree();

    m_pOutputTree->Fill();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

WriteAnalysisParticlesAlgorithm::MCPrimaryMap WriteAnalysisParticlesAlgorithm::GetMainMcParticleMap(const PfoList &pfoList) const
{
    MCPrimaryMap mainMcParticleMap;

    // Find all the analysis particles with associated MC particles and add them to the map.
    for (const ParticleFlowObject * const pPfo : pfoList)
    {
        if (const LArAnalysisParticle *const pAnalysisParticle = dynamic_cast<const LArAnalysisParticle *const>(pPfo))
        {
            if (pAnalysisParticle->HasMcInfo() && pAnalysisParticle->McMainMCParticle())
                mainMcParticleMap.emplace(pAnalysisParticle->McMainMCParticle(), pAnalysisParticle);
        }
    }

    return mainMcParticleMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool WriteAnalysisParticlesAlgorithm::ProcessAnalysisParticle(const LArAnalysisParticle *const pAnalysisParticle,
    const MCPrimaryMap &mainMcParticleMap, const MCParticleList *const pMCParticleList) const
{
    if (LArAnalysisParticleHelper::IsNeutrino(pAnalysisParticle))
    {
        // We should only have up to one neutrino.
        if (m_treeParameters.m_nu_WasReconstructed)
        {
            std::cout << "WriteAnalysisParticlesAlgorithm: multiple neutrinos found - only recording one" << std::endl;
            return false;
        }

        this->PopulateNeutrinoParameters(*pAnalysisParticle, pMCParticleList);
        m_treeParameters.m_nu_WasReconstructed = true;
    }

    else if (LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(pAnalysisParticle))
    {
        this->AddPrimaryDaughterRecord(*pAnalysisParticle, mainMcParticleMap);
        ++m_treeParameters.m_primary_Number;
    }

    else if (LArAnalysisParticleHelper::IsCosmicRay(pAnalysisParticle))
    {
        this->AddCosmicRayRecord(*pAnalysisParticle, mainMcParticleMap);
        ++m_treeParameters.m_cr_Number;
    }

    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: analysis particle was not a cosmic ray, neutrino or primary daughter" << std::endl;
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::RecordUnreconstructedParticles(const MCParticleList *const pMCParticleList,
    const MCPrimaryMap &mainMcParticleMap) const
{
    for (const MCParticle *const pMCPrimary : this->GetAllMcPrimaries(pMCParticleList))
    {
        // Check if it's been recorded by a reco particle.
        if (mainMcParticleMap.find(pMCPrimary) != mainMcParticleMap.end())
            continue;

        // Check if it's a neutrino, primary daughter or CR.
        if (!LArMCParticleHelper::IsNeutrino(pMCPrimary) && !LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) &&
            !LArMCParticleHelper::IsCosmicRay(pMCPrimary))
        {
            continue;
        }

        this->RecordUnreconstructedParticle(pMCParticleList, pMCPrimary);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::RecordUnreconstructedParticle(const MCParticleList *const pMCParticleList, const MCParticle *const pMCPrimary) const
{
    // Perform any quality cuts before recording.
    const LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo = LArAnalysisParticleHelper::GetMcInformation(pMCPrimary, m_minCoordinates, m_maxCoordinates,
        m_mcContainmentFractionLowerBound);

    if ((pfoMcInfo.m_mcEnergy < m_mcOnlyParticleEnergyCut) || (pfoMcInfo.m_mcContainmentFraction < m_mcOnlyParticleContainmentCut))
        return;

    if (LArMCParticleHelper::IsNeutrino(pMCPrimary))
    {
        // We should only have up to one neutrino.
        if (m_treeParameters.m_nu_HasMcInfo)
        {
            std::cout << "WriteAnalysisParticlesAlgorithm: found neutrino that was not covered by MC particles but neutrino MC properties already exist" << std::endl;
            return;
        }

        m_treeParameters.m_nu_HasMcInfo = true;
        this->PopulateNeutrinoMcParameters(pfoMcInfo, pMCParticleList);
    }

    else if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary))
    {
        ++m_treeParameters.m_primary_Number;
        this->AddMcOnlyPrimaryDaughterRecord(pfoMcInfo);
    }

    else if (pMCPrimary->GetParticleId() == MU_MINUS)
    {
        ++m_treeParameters.m_cr_Number;
        this->AddMcOnlyCosmicRayRecord(pfoMcInfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleSet WriteAnalysisParticlesAlgorithm::GetAllMcPrimaries(const MCParticleList *const pMCParticleList) const
{
    MCParticleSet mcPrimarySet;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        // LArMCParticleHelper::GetPrimaryMCParticle will throw if there's no visible particle up the chain.
        try
        {
            if (pMCParticle && (pMCParticle->GetParentList().size() == 0) && LArMCParticleHelper::IsNeutrino(pMCParticle))
                mcPrimarySet.insert(pMCParticle);

            else
                mcPrimarySet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle));
        }

        catch (...)
        {
            continue;
        }
    }

    return mcPrimarySet;
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
    m_treeParameters.m_nu_Momentum                                    = neutrinoAnalysisParticle.AnalysisMomentum().GetMagnitude();
    m_treeParameters.m_nu_MomentumX                                   = neutrinoAnalysisParticle.AnalysisMomentum().GetX();
    m_treeParameters.m_nu_MomentumY                                   = neutrinoAnalysisParticle.AnalysisMomentum().GetY();
    m_treeParameters.m_nu_MomentumZ                                   = neutrinoAnalysisParticle.AnalysisMomentum().GetZ();
    m_treeParameters.m_nu_TypeTree                                    = LArAnalysisParticleHelper::TypeTreeAsString(neutrinoAnalysisParticle.GetTypeTree());
    m_treeParameters.m_nu_NumberOf3dHits                              = neutrinoAnalysisParticle.NumberOf3dHits();
    m_treeParameters.m_nu_NumberOfCollectionPlaneHits                 = neutrinoAnalysisParticle.NumberOfCollectionPlaneHits();
    m_treeParameters.m_nu_NumberOfRecoParticles                       = neutrinoAnalysisParticle.NumberOfDownstreamParticles();

    // Recurse through the analysis particle hierarchy, counting the numbers of tracks and showers.
    unsigned int numberOfRecoTracks(0UL), numberOfRecoShowers(0UL);
    this->CountRecoTracksAndShowers(neutrinoAnalysisParticle, numberOfRecoTracks, numberOfRecoShowers);

    m_treeParameters.m_nu_NumberOfRecoTracks  = numberOfRecoTracks;
    m_treeParameters.m_nu_NumberOfRecoShowers = numberOfRecoShowers;

    // Use the 'momentum' to get the transverse and longitudinal visible energy components.
    const CartesianVector zAxis(0.f, 0.f, 1.f);
    m_treeParameters.m_nu_VisibleLongitudinalEnergy = neutrinoAnalysisParticle.AnalysisMomentum().GetDotProduct(zAxis);
    m_treeParameters.m_nu_VisibleTransverseEnergy   = neutrinoAnalysisParticle.AnalysisMomentum().GetCrossProduct(zAxis).GetMagnitude();

    // Populate the MC parameters if the info exists - if not, they can be left at default values.
    if (neutrinoAnalysisParticle.HasMcInfo() && neutrinoAnalysisParticle.McMainMCParticle())
    {
        const LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo = LArAnalysisParticleHelper::GetMcInformation(neutrinoAnalysisParticle.McMainMCParticle(), m_minCoordinates,
            m_maxCoordinates, m_mcContainmentFractionLowerBound);

        this->PopulateNeutrinoMcParameters(pfoMcInfo, pMCParticleList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::CountRecoTracksAndShowers(const LArAnalysisParticle &currentAnalysisParticle,
    unsigned int &numberOfRecoTracks, unsigned int &numberOfRecoShowers) const
{
    for (const ParticleFlowObject *const pDaughterPfo : currentAnalysisParticle.GetDaughterPfoList())
    {
        if (const LArAnalysisParticle *const pAnalysisDaughter = dynamic_cast<const LArAnalysisParticle *const>(pDaughterPfo))
        {
            if (pAnalysisDaughter->IsShower())
                ++numberOfRecoShowers;

            else
                ++numberOfRecoTracks;

            // Recurse through the analysis particle hierarchy.
            this->CountRecoTracksAndShowers(*pAnalysisDaughter, numberOfRecoTracks, numberOfRecoShowers);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PopulateNeutrinoMcParameters(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo,
    const MCParticleList *const pMCParticleList) const
{
    if (pfoMcInfo.m_pMCParticle)
        m_treeParameters.m_nu_mc_McParticleUid = reinterpret_cast<std::uint64_t>(pfoMcInfo.m_pMCParticle->GetUid());

    else // it can stay at its default value
        std::cout << "WriteAnalysisParticlesAlgorithm: could not populate MC particle UID" << std::endl;

    m_treeParameters.m_nu_mc_Energy              = pfoMcInfo.m_mcEnergy;
    m_treeParameters.m_nu_mc_VertexX             = pfoMcInfo.m_mcVertexPosition.GetX();
    m_treeParameters.m_nu_mc_VertexY             = pfoMcInfo.m_mcVertexPosition.GetY();
    m_treeParameters.m_nu_mc_VertexZ             = pfoMcInfo.m_mcVertexPosition.GetZ();
    m_treeParameters.m_nu_mc_DirectionCosineX    = pfoMcInfo.m_mcDirectionCosines.GetX();
    m_treeParameters.m_nu_mc_DirectionCosineY    = pfoMcInfo.m_mcDirectionCosines.GetY();
    m_treeParameters.m_nu_mc_DirectionCosineZ    = pfoMcInfo.m_mcDirectionCosines.GetZ();
    m_treeParameters.m_nu_mc_Momentum            = pfoMcInfo.m_mcMomentum.GetMagnitude();
    m_treeParameters.m_nu_mc_MomentumX           = pfoMcInfo.m_mcMomentum.GetX();
    m_treeParameters.m_nu_mc_MomentumY           = pfoMcInfo.m_mcMomentum.GetY();
    m_treeParameters.m_nu_mc_MomentumZ           = pfoMcInfo.m_mcMomentum.GetZ();
    m_treeParameters.m_nu_mc_IsVertexFiducial    = pfoMcInfo.m_mcIsVertexFiducial;
    m_treeParameters.m_nu_mc_IsContained         = pfoMcInfo.m_mcIsContained;
    m_treeParameters.m_nu_mc_ContainmentFraction = pfoMcInfo.m_mcContainmentFraction;
    m_treeParameters.m_nu_mc_TypeTree            = LArAnalysisParticleHelper::TypeTreeAsString(pfoMcInfo.m_mcTypeTree);
    m_treeParameters.m_nu_mc_PdgCode             = pfoMcInfo.m_mcPdgCode;

    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: no valid MC particle list name provided but neutrino analysis particle has MC info" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    // Use the MC neutrino's visible primary daughters to construct the MC analogues of the energy and momentum measurements.
    CartesianVector visibleMomentum(0.f, 0.f, 0.f);
    float visibleEnergy(0.f);
    this->CalculateNeutrinoMcVisibleMomentum(pMCParticleList, visibleEnergy, visibleMomentum);

    const LArInteractionTypeHelper::InteractionType interactionType = this->GetInteractionType(pMCParticleList);

    // To align with the non-MC definition, we define the longitudinal/transverse components with respect to the z-direction.
    const CartesianVector zDirectionVector(0.f, 0.f, 1.f);

    m_treeParameters.m_nu_mc_LongitudinalEnergy        = pfoMcInfo.m_mcMomentum.GetDotProduct(zDirectionVector);
    m_treeParameters.m_nu_mc_TransverseEnergy          = pfoMcInfo.m_mcMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    m_treeParameters.m_nu_mc_VisibleEnergy             = visibleEnergy;
    m_treeParameters.m_nu_mc_VisibleLongitudinalEnergy = visibleMomentum.GetDotProduct(zDirectionVector);
    m_treeParameters.m_nu_mc_VisibleTransverseEnergy   = visibleMomentum.GetCrossProduct(zDirectionVector).GetMagnitude();
    m_treeParameters.m_nu_mc_InteractionType           = LArInteractionTypeHelper::ToString(interactionType);
    m_treeParameters.m_nu_mc_IsChargedCurrent          = this->IsChargedCurrent(interactionType);

    float visibleEnergyFraction(0.f);

    if (pfoMcInfo.m_mcEnergy > std::numeric_limits<float>::epsilon())
        visibleEnergyFraction = visibleEnergy / pfoMcInfo.m_mcEnergy;

    m_treeParameters.m_nu_mc_VisibleEnergyFraction = visibleEnergyFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::CalculateNeutrinoMcVisibleMomentum(const MCParticleList *const pMCParticleList, float &visibleEnergy,
    CartesianVector &visibleMomentum) const
{
    visibleMomentum = CartesianVector(0.f, 0.f, 0.f);
    visibleEnergy   = 0.f;

    // Sum over the visible primary daughters of the MC neutrino.
    for (const MCParticle *const pMCPrimary : this->GetAllMcPrimaryDaughters(pMCParticleList))
    {
        const float primaryVisibleEnergy = pMCPrimary->GetEnergy() - PdgTable::GetParticleMass(pMCPrimary->GetParticleId());
        visibleEnergy += primaryVisibleEnergy;

        if (pMCPrimary->GetMomentum().GetMagnitude() > std::numeric_limits<float>::epsilon())
            visibleMomentum += pMCPrimary->GetMomentum().GetUnitVector() * primaryVisibleEnergy;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleSet WriteAnalysisParticlesAlgorithm::GetAllMcPrimaryDaughters(const MCParticleList *const pMCParticleList) const
{
    MCParticleSet mcPrimarySet;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        // LArMCParticleHelper::GetPrimaryMCParticle will throw if there's no visible particle up the chain.
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

    return mcPrimarySet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionType WriteAnalysisParticlesAlgorithm::GetInteractionType(const MCParticleList *const pMCParticleList) const
{
    if (!pMCParticleList)
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: could not get interaction type because there was no MC particle list" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

    MCParticleList mcPrimaryList;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))
            mcPrimaryList.push_back(pMCParticle);
    }

    if (mcPrimaryList.empty())
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: could not get interaction type because the MC particle list was empty" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

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
    PUSH_TREE_RECORD(m_primary_Momentum,                                    primaryAnalysisParticle.AnalysisMomentum().GetMagnitude());
    PUSH_TREE_RECORD(m_primary_MomentumX,                                   primaryAnalysisParticle.AnalysisMomentum().GetX());
    PUSH_TREE_RECORD(m_primary_MomentumY,                                   primaryAnalysisParticle.AnalysisMomentum().GetY());
    PUSH_TREE_RECORD(m_primary_MomentumZ,                                   primaryAnalysisParticle.AnalysisMomentum().GetZ());
    PUSH_TREE_RECORD(m_primary_TypeTree,                                    LArAnalysisParticleHelper::TypeTreeAsString(primaryAnalysisParticle.GetTypeTree()));
    PUSH_TREE_RECORD(m_primary_IsShower,                                    primaryAnalysisParticle.IsShower());
    PUSH_TREE_RECORD(m_primary_IsTrack,                                    !primaryAnalysisParticle.IsShower());
    PUSH_TREE_RECORD(m_primary_IsProton,                                   (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::PROTON));
    PUSH_TREE_RECORD(m_primary_IsPionOrMuon,                               (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::PION_MUON) || (primaryAnalysisParticle.Type() == LArAnalysisParticle::TYPE::COSMIC_RAY));
    PUSH_TREE_RECORD(m_primary_NumberOf3dHits,                              primaryAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_primary_NumberOfCollectionPlaneHits,                 primaryAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_primary_NumberOfDownstreamParticles,                 primaryAnalysisParticle.NumberOfDownstreamParticles());

    // Add the MC record if the info is there, otherwise populate a blank record.
    if (primaryAnalysisParticle.HasMcInfo() && primaryAnalysisParticle.McMainMCParticle())
    {
        const LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo = LArAnalysisParticleHelper::GetMcInformation(primaryAnalysisParticle.McMainMCParticle(),
            m_minCoordinates, m_maxCoordinates, m_mcContainmentFractionLowerBound);

        const bool particleSplitByReco = (coveredMCPrimaries.count(primaryAnalysisParticle.McMainMCParticle()) > 1U);
        this->AddPrimaryDaughterMcRecord(pfoMcInfo, particleSplitByReco);
    }

    else
        this->AddBlankPrimaryDaughterMcRecord();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddPrimaryDaughterMcRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo, const bool particleSplitByReco) const
{
    if (pfoMcInfo.m_pMCParticle)
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, reinterpret_cast<std::uint64_t>(pfoMcInfo.m_pMCParticle->GetUid()));

    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: primary neutrino daughter had MC info but no associated main MC particle" << std::endl;
        PUSH_TREE_RECORD(m_primary_mc_McParticleUid, 0ULL);
    }

    PUSH_TREE_RECORD(m_primary_mc_IsParticleSplitByReco, particleSplitByReco);
    PUSH_TREE_RECORD(m_primary_mc_Energy,                pfoMcInfo.m_mcEnergy);
    PUSH_TREE_RECORD(m_primary_mc_KineticEnergy,         pfoMcInfo.m_mcKineticEnergy);
    PUSH_TREE_RECORD(m_primary_mc_Mass,                  pfoMcInfo.m_mcMass);
    PUSH_TREE_RECORD(m_primary_mc_VertexX,               pfoMcInfo.m_mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_primary_mc_VertexY,               pfoMcInfo.m_mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_primary_mc_VertexZ,               pfoMcInfo.m_mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineX,      pfoMcInfo.m_mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineY,      pfoMcInfo.m_mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_primary_mc_DirectionCosineZ,      pfoMcInfo.m_mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_Momentum,              pfoMcInfo.m_mcMomentum.GetMagnitude());
    PUSH_TREE_RECORD(m_primary_mc_MomentumX,             pfoMcInfo.m_mcMomentum.GetX());
    PUSH_TREE_RECORD(m_primary_mc_MomentumY,             pfoMcInfo.m_mcMomentum.GetY());
    PUSH_TREE_RECORD(m_primary_mc_MomentumZ,             pfoMcInfo.m_mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_primary_mc_IsVertexFiducial,      pfoMcInfo.m_mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_primary_mc_IsContained,           pfoMcInfo.m_mcIsContained);
    PUSH_TREE_RECORD(m_primary_mc_ContainmentFraction,   pfoMcInfo.m_mcContainmentFraction);
    PUSH_TREE_RECORD(m_primary_mc_TypeTree,              LArAnalysisParticleHelper::TypeTreeAsString(pfoMcInfo.m_mcTypeTree));
    PUSH_TREE_RECORD(m_primary_mc_IsShower,              pfoMcInfo.m_mcIsShower);
    PUSH_TREE_RECORD(m_primary_mc_IsTrack,              !pfoMcInfo.m_mcIsShower);
    PUSH_TREE_RECORD(m_primary_mc_IsProton,              pfoMcInfo.m_mcIsProton);
    PUSH_TREE_RECORD(m_primary_mc_IsPionOrMuon,          pfoMcInfo.m_mcIsPionOrMuon);
    PUSH_TREE_RECORD(m_primary_mc_IsCosmicRay,           pfoMcInfo.m_mcIsCosmicRay);
    PUSH_TREE_RECORD(m_primary_mc_PdgCode,               pfoMcInfo.m_mcPdgCode);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddBlankPrimaryDaughterMcRecord() const
{
    TREE_VECTOR_MEMBERS_PRIMARY_MC(VECTOR_MEMBER_PUSH_DEFAULT)
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyPrimaryDaughterRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo) const
{
    // m_primary_Number is dealt with by the calling method.

    TREE_VECTOR_MEMBERS_PRIMARY(VECTOR_MEMBER_PUSH_DEFAULT)

    this->AddPrimaryDaughterMcRecord(pfoMcInfo, false); // cannot refer to a split particle by construction
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
    PUSH_TREE_RECORD(m_cr_Momentum,                                    cosmicRayAnalysisParticle.AnalysisMomentum().GetMagnitude());
    PUSH_TREE_RECORD(m_cr_MomentumX,                                   cosmicRayAnalysisParticle.AnalysisMomentum().GetX());
    PUSH_TREE_RECORD(m_cr_MomentumY,                                   cosmicRayAnalysisParticle.AnalysisMomentum().GetY());
    PUSH_TREE_RECORD(m_cr_MomentumZ,                                   cosmicRayAnalysisParticle.AnalysisMomentum().GetZ());
    PUSH_TREE_RECORD(m_cr_TypeTree,                                    LArAnalysisParticleHelper::TypeTreeAsString(cosmicRayAnalysisParticle.GetTypeTree()));
    PUSH_TREE_RECORD(m_cr_NumberOf3dHits,                              cosmicRayAnalysisParticle.NumberOf3dHits());
    PUSH_TREE_RECORD(m_cr_NumberOfCollectionPlaneHits,                 cosmicRayAnalysisParticle.NumberOfCollectionPlaneHits());
    PUSH_TREE_RECORD(m_cr_NumberOfDownstreamParticles,                 cosmicRayAnalysisParticle.NumberOfDownstreamParticles());

    // Add the MC record if the info is there, otherwise populate a blank record.
    if (cosmicRayAnalysisParticle.HasMcInfo())
    {
        const LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo = LArAnalysisParticleHelper::GetMcInformation(cosmicRayAnalysisParticle.McMainMCParticle(),
            m_minCoordinates, m_maxCoordinates, m_mcContainmentFractionLowerBound);

        const bool particleSplitByReco = (coveredMCPrimaries.count(cosmicRayAnalysisParticle.McMainMCParticle()) > 1U);
        this->AddCosmicRayMcRecord(pfoMcInfo, particleSplitByReco);
    }

    else
        this->AddBlankCosmicRayMcRecord();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddCosmicRayMcRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo,
    const bool particleSplitByReco) const
{
    if (pfoMcInfo.m_pMCParticle)
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid, reinterpret_cast<std::uint64_t>(pfoMcInfo.m_pMCParticle->GetUid()));

    else
    {
        std::cout << "WriteAnalysisParticlesAlgorithm: cosmic ray had MC info but no associated main MC particle" << std::endl;
        PUSH_TREE_RECORD(m_cr_mc_McParticleUid, 0ULL);
    }

    PUSH_TREE_RECORD(m_cr_mc_IsParticleSplitByReco,          particleSplitByReco);
    PUSH_TREE_RECORD(m_cr_mc_Energy,                         pfoMcInfo.m_mcEnergy);
    PUSH_TREE_RECORD(m_cr_mc_KineticEnergy,                  pfoMcInfo.m_mcKineticEnergy);
    PUSH_TREE_RECORD(m_cr_mc_Mass,                           pfoMcInfo.m_mcMass);
    PUSH_TREE_RECORD(m_cr_mc_VertexX,                        pfoMcInfo.m_mcVertexPosition.GetX());
    PUSH_TREE_RECORD(m_cr_mc_VertexY,                        pfoMcInfo.m_mcVertexPosition.GetY());
    PUSH_TREE_RECORD(m_cr_mc_VertexZ,                        pfoMcInfo.m_mcVertexPosition.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineX,               pfoMcInfo.m_mcDirectionCosines.GetX());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineY,               pfoMcInfo.m_mcDirectionCosines.GetY());
    PUSH_TREE_RECORD(m_cr_mc_DirectionCosineZ,               pfoMcInfo.m_mcDirectionCosines.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_Momentum,                       pfoMcInfo.m_mcMomentum.GetMagnitude());
    PUSH_TREE_RECORD(m_cr_mc_MomentumX,                      pfoMcInfo.m_mcMomentum.GetX());
    PUSH_TREE_RECORD(m_cr_mc_MomentumY,                      pfoMcInfo.m_mcMomentum.GetY());
    PUSH_TREE_RECORD(m_cr_mc_MomentumZ,                      pfoMcInfo.m_mcMomentum.GetZ());
    PUSH_TREE_RECORD(m_cr_mc_IsVertexFiducial,               pfoMcInfo.m_mcIsVertexFiducial);
    PUSH_TREE_RECORD(m_cr_mc_IsContained,                    pfoMcInfo.m_mcIsContained);
    PUSH_TREE_RECORD(m_cr_mc_ContainmentFraction,            pfoMcInfo.m_mcContainmentFraction);
    PUSH_TREE_RECORD(m_cr_mc_TypeTree,                       LArAnalysisParticleHelper::TypeTreeAsString(pfoMcInfo.m_mcTypeTree));
    PUSH_TREE_RECORD(m_cr_mc_IsShower,                       pfoMcInfo.m_mcIsShower);
    PUSH_TREE_RECORD(m_cr_mc_IsTrack,                       !pfoMcInfo.m_mcIsShower);
    PUSH_TREE_RECORD(m_cr_mc_IsProton,                       pfoMcInfo.m_mcIsProton);
    PUSH_TREE_RECORD(m_cr_mc_IsPionOrMuon,                   pfoMcInfo.m_mcIsPionOrMuon);
    PUSH_TREE_RECORD(m_cr_mc_IsCosmicRay,                    pfoMcInfo.m_mcIsCosmicRay);
    PUSH_TREE_RECORD(m_cr_mc_PdgCode,                        pfoMcInfo.m_mcPdgCode);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddBlankCosmicRayMcRecord() const
{
    TREE_VECTOR_MEMBERS_CR_MC(VECTOR_MEMBER_PUSH_DEFAULT)
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::AddMcOnlyCosmicRayRecord(const LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo) const
{
    // m_cr_Number is dealt with by the calling method.

    TREE_VECTOR_MEMBERS_CR(VECTOR_MEMBER_PUSH_DEFAULT)

    this->AddCosmicRayMcRecord(pfoMcInfo, false); // cannot refer to a split particle by construction
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::CheckTreeVectorSizes() const
{
    TREE_VECTOR_MEMBERS_PRIMARY(CHECK_VECTOR_MEMBER_SIZE)
    TREE_VECTOR_MEMBERS_PRIMARY_MC(CHECK_VECTOR_MEMBER_SIZE)
    TREE_VECTOR_MEMBERS_CR(CHECK_VECTOR_MEMBER_SIZE)
    TREE_VECTOR_MEMBERS_CR_MC(CHECK_VECTOR_MEMBER_SIZE)
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::PrintTree() const
{
    const std::string nuLabel = " - [nu]        ";

    std::cout << "Pandora Tree:\n";
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
    std::cout << nuLabel << "Momentum:                         " << 1000.f * m_treeParameters.m_nu_Momentum << " MeV/c\n";
    std::cout << nuLabel << "Momentum x:                       " << 1000.f * m_treeParameters.m_nu_MomentumX << " MeV/c\n";
    std::cout << nuLabel << "Momentum y:                       " << 1000.f * m_treeParameters.m_nu_MomentumY << " MeV/c\n";
    std::cout << nuLabel << "Momentum z:                       " << 1000.f * m_treeParameters.m_nu_MomentumZ << " MeV/c\n";
    std::cout << nuLabel << "Type tree:                        " << m_treeParameters.m_nu_TypeTree << '\n';
    std::cout << nuLabel << "Number of 3D hits:                " << m_treeParameters.m_nu_NumberOf3dHits << '\n';
    std::cout << nuLabel << "Number of coll plane hits:        " << m_treeParameters.m_nu_NumberOfCollectionPlaneHits << '\n';
    std::cout << nuLabel << "Number of reco particles:         " << m_treeParameters.m_nu_NumberOfRecoParticles << '\n';
    std::cout << nuLabel << "Number of reco tracks:            " << m_treeParameters.m_nu_NumberOfRecoTracks << '\n';
    std::cout << nuLabel << "Number of reco showers:           " << m_treeParameters.m_nu_NumberOfRecoShowers << '\n';
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
        std::cout << label << "Momentum:                         " << 1000.f * m_treeParameters.m_primary_Momentum.at(i) << " MeV/c\n";
        std::cout << label << "Momentum x:                       " << 1000.f * m_treeParameters.m_primary_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "Momentum y:                       " << 1000.f * m_treeParameters.m_primary_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "Momentum z:                       " << 1000.f * m_treeParameters.m_primary_MomentumZ.at(i) << " MeV/c\n";
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
        std::cout << label << "Momentum:                         " << m_treeParameters.m_cr_Momentum.at(i) << " MeV/c\n";
        std::cout << label << "Momentum x:                       " << m_treeParameters.m_cr_MomentumX.at(i) << " MeV/c\n";
        std::cout << label << "Momentum y:                       " << m_treeParameters.m_cr_MomentumY.at(i) << " MeV/c\n";
        std::cout << label << "Momentum z:                       " << m_treeParameters.m_cr_MomentumZ.at(i) << " MeV/c\n";
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TreeTitle", m_treeTitle));
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

    this->InitializeTree();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteAnalysisParticlesAlgorithm::InitializeTree()
{
    m_pOutputTFile = new TFile(m_outputFile.c_str(), "UPDATE");
    m_pOutputTree  = new TTree(m_treeName.c_str(), m_treeTitle.c_str());

    TREE_SCALAR_MEMBERS(SET_SCALAR_MEMBER_BRANCH)
    TREE_VECTOR_MEMBERS_PRIMARY(SET_VECTOR_MEMBER_BRANCH)
    TREE_VECTOR_MEMBERS_PRIMARY_MC(SET_VECTOR_MEMBER_BRANCH)
    TREE_VECTOR_MEMBERS_CR(SET_VECTOR_MEMBER_BRANCH)
    TREE_VECTOR_MEMBERS_CR_MC(SET_VECTOR_MEMBER_BRANCH)
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TreeParameters::TreeParameters() noexcept :
    TREE_SCALAR_MEMBERS(INITIALIZE_SCALAR_MEMBER)
    TREE_VECTOR_MEMBERS_PRIMARY(INITIALIZE_VECTOR_MEMBER)
    TREE_VECTOR_MEMBERS_PRIMARY_MC(INITIALIZE_VECTOR_MEMBER)
    TREE_VECTOR_MEMBERS_CR(INITIALIZE_VECTOR_MEMBER)
    TREE_VECTOR_MEMBERS_CR_MC(INITIALIZE_VECTOR_MEMBER)
    m_dummy(false)
{
}

} // namespace lar_physics_content
