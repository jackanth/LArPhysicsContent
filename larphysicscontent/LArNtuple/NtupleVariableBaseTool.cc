/**
 *  @file   larphysicscontent/LArNtuple/NtupleVariableBaseTool.cc
 *
 *  @brief  Implementation of the ntuple variable base tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

#include "larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"
#include "larphysicscontent/LArNtuple/LArNtuple.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

NtupleVariableBaseTool::NtupleVariableBaseTool() noexcept :
    AlgorithmTool(),
    m_spNtuple(),
    m_eventPrefix("evt_"),
    m_neutrinoPrefix("nu_"),
    m_primaryPrefix("primary_"),
    m_particlePrefix("pfo_"),
    m_cosmicPrefix("cr_"),
    m_minFiducialCoordinates(0.f, 0.f, 0.f),
    m_maxFiducialCoordinates(0.f, 0.f, 0.f),
    m_pAlgorithm(nullptr),
    m_spPlotsRegistry(nullptr),
    m_spTmpRegistry(nullptr),
    m_isSetup(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &NtupleVariableBaseTool::GetAllDownstreamThreeDHits(const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetAllDownstreamThreeDHits(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList NtupleVariableBaseTool::GetAllDownstreamTwoDHits(const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetAllDownstreamTwoDHits(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &NtupleVariableBaseTool::GetAllDownstreamUHits(const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetAllDownstreamUHits(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &NtupleVariableBaseTool::GetAllDownstreamVHits(const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetAllDownstreamVHits(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &NtupleVariableBaseTool::GetAllDownstreamWHits(const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetAllDownstreamWHits(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList &NtupleVariableBaseTool::GetAllDownstreamPfos(const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetAllDownstreamPfos(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArNtupleHelper::TrackFitSharedPtr &NtupleVariableBaseTool::GetTrackFit(const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetTrackFit(this->GetPandora(), pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NtupleVariableBaseTool::PrepareEventWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const PfoList &pfoList,
    const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo)
{
    if (!m_isSetup)
    {
        std::cerr << "NtupleVariableBaseTool: Tool was not setup" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    this->PrepareEvent(pfoList, eventValidationInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessEventWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::vector<LArNtupleRecord> records = this->ProcessEvent(pfoList, eventValidationInfo);

    for (LArNtupleRecord &record : records)
        record.AddBranchNamePrefix(m_eventPrefix);

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessNeutrinoWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo)
{
    return this->ProcessImpl(pAlgorithm, m_neutrinoPrefix, pPfo, spInteractionInfo->GetMcNeutrino(),
        [&]() { return this->ProcessNeutrino(pPfo, pfoList, spInteractionInfo); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessPrimaryWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    return this->ProcessImpl(
        pAlgorithm, m_primaryPrefix, pPfo, spMcTarget->GetMCParticle(), [&]() { return this->ProcessPrimary(pPfo, pfoList, spMcTarget); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessCosmicRayWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    return this->ProcessImpl(
        pAlgorithm, m_cosmicPrefix, pPfo, spMcTarget->GetMCParticle(), [&]() { return this->ProcessCosmicRay(pPfo, pfoList, spMcTarget); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessImpl(const AnalysisNtupleAlgorithm *const pAlgorithm, const std::string &prefix,
    const ParticleFlowObject *const pPfo, const MCParticle *const pMcParticle, const Processor &processor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Add the prefix and particles to the records
    std::vector<LArNtupleRecord> records = processor();

    for (LArNtupleRecord &record : records)
    {
        record.AddBranchNamePrefix(prefix);

        if (pPfo)
            record.SetPfo(pPfo);

        if (pMcParticle)
            record.SetMCParticle(pMcParticle);
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NtupleVariableBaseTool::Setup(std::shared_ptr<LArNtuple> spNtuple, const pandora::Algorithm *const pAlgorithm,
    pandora::CartesianVector minFiducialCoordinates, pandora::CartesianVector maxFiducialCoordinates,
    std::shared_ptr<LArRootRegistry> spPlotsRegistry, std::shared_ptr<LArRootRegistry> spTmpRegistry)
{
    if (!pAlgorithm || !spNtuple || !spPlotsRegistry || !spTmpRegistry)
    {
        std::cerr << "NtupleVariableBaseTool: Failed to perform setup" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    m_spNtuple               = std::move(spNtuple);
    m_pAlgorithm             = pAlgorithm;
    m_minFiducialCoordinates = std::move(minFiducialCoordinates);
    m_maxFiducialCoordinates = std::move(maxFiducialCoordinates);
    m_spPlotsRegistry        = std::move(spPlotsRegistry);
    m_spTmpRegistry          = std::move(spTmpRegistry);
    m_isSetup                = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArBranchPlaceholder::NtupleRecordSPtr NtupleVariableBaseTool::GetScalarRecord(const std::string &branchName) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetScalarRecord(branchName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArBranchPlaceholder::NtupleRecordSPtr NtupleVariableBaseTool::GetVectorRecordElement(
    LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetVectorRecordElement(type, branchName, pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArBranchPlaceholder::NtupleRecordSPtr NtupleVariableBaseTool::GetVectorRecordElement(
    LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const MCParticle *const pMCParticle) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetVectorRecordElement(type, branchName, pMCParticle);
}

} // namespace lar_physics_content
