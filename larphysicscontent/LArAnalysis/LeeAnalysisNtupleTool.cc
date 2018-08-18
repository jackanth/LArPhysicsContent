/**
 *  @file   larphysicscontent/LArAnalysis/LeeAnalysisNtupleTool.cc
 *
 *  @brief  Implementation of the LEE analysis ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/LeeAnalysisNtupleTool.h"
#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{
LeeAnalysisNtupleTool::LeeAnalysisNtupleTool() : NtupleVariableBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LeeAnalysisNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> LeeAnalysisNtupleTool::ProcessEvent(const PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> LeeAnalysisNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> LeeAnalysisNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> LeeAnalysisNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

} // namespace lar_physics_content
