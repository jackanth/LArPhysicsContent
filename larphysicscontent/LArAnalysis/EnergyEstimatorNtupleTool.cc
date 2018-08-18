/**
 *  @file   larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.cc
 *
 *  @brief  Implementation of the energy estimator ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.h"
#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{
EnergyEstimatorNtupleTool::EnergyEstimatorNtupleTool() : NtupleVariableBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyEstimatorNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessEvent(const PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

} // namespace lar_physics_content
