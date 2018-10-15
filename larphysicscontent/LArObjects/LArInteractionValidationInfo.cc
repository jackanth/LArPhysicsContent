/**
 *  @file   larphysicscontent/LArNtuple/LArInteractionValidationInfo.cc
 *
 *  @brief  Implementation of the lar interaction validation info class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArObjects/LArInteractionValidationInfo.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

void LArInteractionValidationInfo::SetParameters(const lar_content::LArInteractionTypeHelper::InteractionType interactionType, const bool isCorrect,
    const bool isFake, const bool isSplit, const bool isLost, const ParticleFlowObject *const pRecoNeutrino, const MCParticle *const pMcNeutrino)
{
    if (m_parametersSet)
    {
        std::cerr << "LArInteractionValidationInfo: Could not set parameters because they had already been set" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    m_interactionType = interactionType;
    m_isCorrect       = isCorrect;
    m_isFake          = isFake;
    m_isSplit         = isSplit;
    m_isLost          = isLost;
    m_pRecoNeutrino   = pRecoNeutrino;
    m_pMcNeutrino     = pMcNeutrino;
    m_parametersSet   = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::shared_ptr<LArMCTargetValidationInfo> &LArInteractionValidationInfo::AddDaughterTarget(
    const MCParticle *const pMCPrimary, const CaloHitList &mcPrimaryHitList, const bool isTargetPrimary)
{
    auto spDaughterTarget = std::shared_ptr<LArMCTargetValidationInfo>(
        new LArMCTargetValidationInfo(this->shared_from_this(), pMCPrimary, mcPrimaryHitList, isTargetPrimary));
    return *m_daughterTargets.emplace(std::move(spDaughterTarget)).first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

lar_content::LArInteractionTypeHelper::InteractionType LArInteractionValidationInfo::GetInteractionType() const
{
    if (m_parametersSet)
        return m_interactionType;

    std::cerr << "LArInteractionValidationInfo: Could not get parameter because it has not been set" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArInteractionValidationInfo::IsCorrect() const
{
    if (m_parametersSet)
        return m_isCorrect;

    std::cerr << "LArInteractionValidationInfo: Could not get parameter because it has not been set" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArInteractionValidationInfo::IsFake() const
{
    if (m_parametersSet)
        return m_isFake;

    std::cerr << "LArInteractionValidationInfo: Could not get parameter because it has not been set" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArInteractionValidationInfo::IsSplit() const
{
    if (m_parametersSet)
        return m_isSplit;

    std::cerr << "LArInteractionValidationInfo: Could not get parameter because it has not been set" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArInteractionValidationInfo::IsLost() const
{
    if (m_parametersSet)
        return m_isLost;

    std::cerr << "LArInteractionValidationInfo: Could not get parameter because it has not been set" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *LArInteractionValidationInfo::GetRecoNeutrino() const
{
    if (m_parametersSet)
        return m_pRecoNeutrino;

    std::cerr << "LArInteractionValidationInfo: Could not get parameter because it has not been set" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArInteractionValidationInfo::GetMcNeutrino() const
{
    if (m_parametersSet)
        return m_pMcNeutrino;

    std::cerr << "LArInteractionValidationInfo: Could not get parameter because it has not been set" << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionValidationInfo::LArInteractionValidationInfo(const int mcNuanceCode, const bool isCosmicRay) noexcept :
    m_daughterTargets(),
    m_mcNuanceCode(mcNuanceCode),
    m_isCosmicRay(isCosmicRay),
    m_interactionType(LArInteractionTypeHelper::InteractionType::ALL_INTERACTIONS),
    m_isCorrect(false),
    m_isFake(false),
    m_isSplit(false),
    m_isLost(false),
    m_pRecoNeutrino(nullptr),
    m_pMcNeutrino(nullptr),
    m_parametersSet(false)
{
}

} // namespace lar_physics_content
