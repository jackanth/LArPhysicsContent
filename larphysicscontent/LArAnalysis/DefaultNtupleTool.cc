/**
 *  @file   larphysicscontent/LArAnalysis/DefaultNtupleTool.cc
 *
 *  @brief  Implementation of the default ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/DefaultNtupleTool.h"
#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{
DefaultNtupleTool::DefaultNtupleTool() :
    NtupleVariableBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DefaultNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessEvent(const pandora::PfoList &, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const, const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessParticle(
    const ParticleFlowObject *const, const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const, const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

} // namespace lar_physics_content
