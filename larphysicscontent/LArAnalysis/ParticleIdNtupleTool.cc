/**
 *  @file   larphysicscontent/LArAnalysis/ParticleIdNtupleTool.cc
 *
 *  @brief  Implementation of the particle ID ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/ParticleIdNtupleTool.h"
#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{
ParticleIdNtupleTool::ParticleIdNtupleTool() : NtupleVariableBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleIdNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> ParticleIdNtupleTool::ProcessEvent(const PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> ParticleIdNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> ParticleIdNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> ParticleIdNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

} // namespace lar_physics_content
