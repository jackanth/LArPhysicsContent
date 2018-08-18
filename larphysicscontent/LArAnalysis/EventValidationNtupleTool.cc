/**
 *  @file   larphysicscontent/LArAnalysis/EventValidationNtupleTool.cc
 *
 *  @brief  Implementation of the event validation ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/EventValidationNtupleTool.h"
#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{
EventValidationNtupleTool::EventValidationNtupleTool() : NtupleVariableBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessEvent(const PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

} // namespace lar_physics_content
