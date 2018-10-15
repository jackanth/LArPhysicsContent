/**
 *  @file   larphysicscontent/LArNtuple/LArMCTargetValidationInfo.cc
 *
 *  @brief  Implementation of the lar MC target validation info class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArObjects/LArMCTargetValidationInfo.h"

using namespace pandora;

namespace lar_physics_content
{
LArMCTargetValidationInfo::LArMCTargetValidationInfo(std::shared_ptr<LArInteractionValidationInfo> spInteractionValidationInfo,
    const MCParticle *const pMCPrimary, const CaloHitList &mcPrimaryHitList, const bool isTargetMCPrimary) noexcept :
    m_spParentInteractionInfo(std::move_if_noexcept(spInteractionValidationInfo)),
    m_daughterMatches(),
    m_pMCParticle(pMCPrimary),
    m_isTargetMCPrimary(isTargetMCPrimary),
    m_mcHits(std::move_if_noexcept(mcPrimaryHitList))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<LArMCMatchValidationInfo> LArMCTargetValidationInfo::GetMatch(const pandora::ParticleFlowObject *const pPfo) const noexcept
{
    for (const std::shared_ptr<LArMCMatchValidationInfo> &spMatch : m_daughterMatches)
    {
        if (spMatch->GetPfo() == pPfo)
            return spMatch;
    }

    return nullptr;
}

} // namespace lar_physics_content
