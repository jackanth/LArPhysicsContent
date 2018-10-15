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

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessEvent(
    const PfoList &, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo)
{
    std::vector<LArNtupleRecord> records;

    if (!eventValidationInfo.empty())
    {
        std::size_t numCosmicInteractions(0UL), numNuInteractions(0UL);

        for (const std::shared_ptr<LArInteractionValidationInfo> &spMcInteraction : eventValidationInfo)
        {
            if (spMcInteraction->IsCosmicRay())
                ++numCosmicInteractions;

            else
                ++numNuInteractions;
        }

        records.emplace_back("mc_NumInteractions", static_cast<LArNtupleRecord::RUInt>(eventValidationInfo.size()));
        records.emplace_back("mc_NumNeutrinoInteractions", static_cast<LArNtupleRecord::RUInt>(numNuInteractions));
        records.emplace_back("mc_NumCosmicRayInteractions", static_cast<LArNtupleRecord::RUInt>(numCosmicInteractions));
    }

    else
    {
        records.emplace_back("mc_NumInteractions", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("mc_NumNeutrinoInteractions", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("mc_NumCosmicRayInteractions", static_cast<LArNtupleRecord::RUInt>(0U));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const std::shared_ptr<LArInteractionValidationInfo> &spMcInteraction)
{
    return this->WriteInteractionRecords(spMcInteraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const pPfo, const PfoList &, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    return this->WriteMatchRecords(pPfo, spMcTarget);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const pPfo, const PfoList &, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    std::vector<LArNtupleRecord> records = this->WriteInteractionRecords(spMcTarget ? spMcTarget->GetParentInteractionInfo() : nullptr);

    std::vector<LArNtupleRecord> matchRecords = this->WriteMatchRecords(pPfo, spMcTarget);
    records.insert(records.end(), std::make_move_iterator(matchRecords.begin()), std::make_move_iterator(matchRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::WriteInteractionRecords(const std::shared_ptr<LArInteractionValidationInfo> &spMcInteraction) const
{
    std::vector<LArNtupleRecord> records;

    if (spMcInteraction)
    {
        records.emplace_back("mc_NuanceCode", static_cast<LArNtupleRecord::RUInt>(spMcInteraction->GetNuanceCode()));
        records.emplace_back("mc_IsCorrect", static_cast<LArNtupleRecord::RBool>(spMcInteraction->IsCorrect()));
        records.emplace_back("mc_IsFake", static_cast<LArNtupleRecord::RBool>(spMcInteraction->IsFake()));
        records.emplace_back("mc_IsSplit", static_cast<LArNtupleRecord::RBool>(spMcInteraction->IsSplit()));
        records.emplace_back("mc_IsLost", static_cast<LArNtupleRecord::RBool>(spMcInteraction->IsLost()));
        records.emplace_back(
            "mc_InteractionType", LArNtupleRecord::RTString(LArInteractionTypeHelper::ToString(spMcInteraction->GetInteractionType())));
    }

    else
    {
        records.emplace_back("mc_NuanceCode", static_cast<LArNtupleRecord::RUInt>(0UL));
        records.emplace_back("mc_IsCorrect", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_IsFake", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_IsSplit", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_IsLost", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_InteractionType", LArNtupleRecord::RTString(""));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EventValidationNtupleTool::WriteMatchRecords(
    const ParticleFlowObject *const pPfo, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) const
{
    std::vector<LArNtupleRecord> records;

    if (pPfo && spMcTarget)
    {
        const std::shared_ptr<LArMCMatchValidationInfo> &spMcMatch = spMcTarget->GetMatch(pPfo);

        if (!spMcMatch)
        {
            std::cerr << "EventValidationNtupleTool: Could not find PFO match in parent target" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        records.emplace_back("mc_MatchPurity", static_cast<LArNtupleRecord::RFloat>(spMcMatch->GetPurity()));
        records.emplace_back("mc_MatchCompleteness", static_cast<LArNtupleRecord::RFloat>(spMcMatch->GetCompleteness()));
        records.emplace_back("mc_IsGoodMatch", static_cast<LArNtupleRecord::RBool>(spMcMatch->IsGoodMatch()));
    }

    else
    {
        records.emplace_back("mc_MatchPurity", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_MatchCompleteness", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_IsGoodMatch", static_cast<LArNtupleRecord::RBool>(false));
    }

    return records;
}

} // namespace lar_physics_content
