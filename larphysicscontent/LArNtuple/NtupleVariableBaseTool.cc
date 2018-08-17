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

NtupleVariableBaseTool::NtupleVariableBaseTool() noexcept : AlgorithmTool(), m_spNtuple()
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

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessEventWrapper(
    const AnalysisNtupleAlgorithm *const pAlgorithm, const PfoList &pfoList, const MCParticleList *const pMCParticleList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::vector<LArNtupleRecord> records = ProcessEvent(pfoList, pMCParticleList);

    for (LArNtupleRecord &record : records)
        record.AddBranchNamePrefix("evt");

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessParticleWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    
    return this->ProcessImpl(pAlgorithm, "pfo", [&]() { return this->ProcessParticle(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessNeutrinoWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    return this->ProcessImpl(pAlgorithm, "nu", [&]() { return this->ProcessNeutrino(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessPrimaryWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    return this->ProcessImpl(pAlgorithm, "primary", [&]() { return this->ProcessPrimary(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessCosmicRayWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    return this->ProcessImpl(pAlgorithm, "cr", [&]() { return this->ProcessCosmicRay(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessImpl(
    const AnalysisNtupleAlgorithm *const pAlgorithm, const std::string &prefix, const Processor &processor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Add the prefix to the records and return them
    std::vector<LArNtupleRecord> records = processor();

    for (LArNtupleRecord &record : records)
        record.AddBranchNamePrefix(prefix);

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::MCParticle *NtupleVariableBaseTool::GetMCParticle(
    const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return m_spNtuple->GetMCParticleWrapper(pPfo, pMCParticleList);
}

} // namespace lar_physics_content
