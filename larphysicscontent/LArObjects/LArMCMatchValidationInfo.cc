/**
 *  @file   larphysicscontent/LArNtuple/LArMCMatchValidationInfo.cc
 *
 *  @brief  Implementation of the lar MC match validation info class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArObjects/LArMCMatchValidationInfo.h"

using namespace pandora;

namespace lar_physics_content
{
LArMCMatchValidationInfo::LArMCMatchValidationInfo(std::shared_ptr<LArMCTargetValidationInfo> spParentTarget,
    const ParticleFlowObject *const pPfo, const bool isRecoCosmicRay, CaloHitList sharedHits, CaloHitList pfoHits, const float purity,
    const float completeness, const bool isGoodMatch, const bool isBestMatch) noexcept :
    m_spParentTarget(std::move_if_noexcept(spParentTarget)),
    m_pPfo(pPfo),
    m_isRecoCosmicRay(isRecoCosmicRay),
    m_sharedHits(std::move_if_noexcept(sharedHits)),
    m_pfoHits(std::move_if_noexcept(pfoHits)),
    m_purity(purity),
    m_completeness(completeness),
    m_isGoodMatch(isGoodMatch),
    m_isBestMatch(isBestMatch)
{
}

} // namespace lar_physics_content
