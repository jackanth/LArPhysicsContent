/**
 *  @file   larphysicscontent/LArAnalysis/EventValidationNtupleTool.cc
 *
 *  @brief  Implementation of the event validation ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/EventValidationNtupleTool.h"

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;
using namespace lar_content;

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

void EventValidationNtupleTool::PrepareEvent(const PfoList &, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessEvent(const PfoList &, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const std::shared_ptr<LArInteractionValidationInfo> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const, const PfoList &, const std::shared_ptr<LArMCTargetValidationInfo> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const, const PfoList &, const std::shared_ptr<LArMCTargetValidationInfo> &)
{
    return {};
}

} // namespace lar_physics_content
