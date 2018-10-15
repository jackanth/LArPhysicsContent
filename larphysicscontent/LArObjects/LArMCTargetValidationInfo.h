/**
 *  @file   larphysicscontent/LArObjects/LArMCTargetValidationInfo.h
 *
 *  @brief  Header file for the lar MC target validation info class.
 *
 *  $Log: $
 */
#ifndef LAR_MC_TARGET_VALIDATION_INFO_H
#define LAR_MC_TARGET_VALIDATION_INFO_H 1

#include "larphysicscontent/LArObjects/LArMCMatchValidationInfo.h"

#include <memory>

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the LArInteractionValidationInfo class
 */
class LArInteractionValidationInfo;

/**
 *  @brief  LArMCTargetValidationInfo class
 */
class LArMCTargetValidationInfo : public std::enable_shared_from_this<LArMCTargetValidationInfo>
{
public:
    using MatchSet = std::set<std::shared_ptr<LArMCMatchValidationInfo>, std::owner_less<std::shared_ptr<LArMCMatchValidationInfo>>>; ///< Alias for a set of shared pointers to matches

    /**
     * @brief  Default copy constructor
     */
    LArMCTargetValidationInfo(const LArMCTargetValidationInfo &) = default;

    /**
     * @brief  Default move constructor
     */
    LArMCTargetValidationInfo(LArMCTargetValidationInfo &&) = default;

    /**
     * @brief  Default copy assignment operator
     */
    LArMCTargetValidationInfo &operator=(const LArMCTargetValidationInfo &) = default;

    /**
     * @brief  Default move assignment operator
     */
    LArMCTargetValidationInfo &operator=(LArMCTargetValidationInfo &&) = default;

    /**
     * @brief  Default destructor
     */
    ~LArMCTargetValidationInfo() = default;

    /**
     *  @brief  Add a daughter match
     *
     *  @param  pPfo address of the PFO
     *  @param  isRecoCosmicRay whether it is a reco cosmic ray
     *  @param  sharedHits the shared hits
     *  @param  pfoHits the PFO hits
     *  @param  purity the match purity
     *  @param  completeness the match completeness
     *  @param  isGoodMatch whether it is a good match
     *  @param  isBestMatch whether it is the best match
     *
     *  @return shared pointer to the match info object
     */
    const std::shared_ptr<LArMCMatchValidationInfo> &AddDaughterMatch(const pandora::ParticleFlowObject *const pPfo,
        const bool isRecoCosmicRay, pandora::CaloHitList sharedHits, pandora::CaloHitList pfoHits, const float purity,
        const float completeness, const bool isGoodMatch, const bool isBestMatch);

    /**
     *  @brief  Get the parent interaction info
     *
     *  @return shared pointer to the parent interaction info
     */
    const std::shared_ptr<LArInteractionValidationInfo> &GetParentInteractionInfo() const noexcept;

    /**
     *  @brief  Get the daughter matches
     *
     *  @return the daughter matches
     */
    const MatchSet &GetDaughterMatches() const noexcept;

    /**
     *  @brief  Get the MCParticle
     *
     *  @return address of the MCParticle
     */
    const pandora::MCParticle *GetMCParticle() const noexcept;

    /**
     *  @brief  Get whether this is a target MC primary
     *
     *  @return whether this is a target MC primary
     */
    bool IsTargetMCPrimary() const noexcept;

    /**
     *  @brief  Get the MC hits
     *
     *  @return the MC hits
     */
    const pandora::CaloHitList &GetMCHits() const noexcept;

    /**
     *  @brief  Get a match for a specified  PFO
     * 
     *  @param  pPfo address of the PFO
     *
     *  @return shared pointer to the match, or nullptr if no match
     */
    std::shared_ptr<LArMCMatchValidationInfo> GetMatch(const pandora::ParticleFlowObject *const pPfo) const noexcept;

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  spInteractionValidationInfo shared pointer to the interaction validation info
     *  @param  pMCPrimary address of the MC primary
     *  @param  mcPrimaryHitList the hit list
     *  @param  isTargetMCPrimary whether this is a target MC primary
     */
    LArMCTargetValidationInfo(std::shared_ptr<LArInteractionValidationInfo> spInteractionValidationInfo,
        const pandora::MCParticle *const pMCPrimary, const pandora::CaloHitList &mcPrimaryHitList, const bool isTargetMCPrimary) noexcept;

    friend class LArInteractionValidationInfo;

private:
    std::shared_ptr<LArInteractionValidationInfo> m_spParentInteractionInfo; ///< Shared pointer to the parent evinteractionent info
    MatchSet                                      m_daughterMatches;         ///< The daughter matches
    const pandora::MCParticle *                   m_pMCParticle;             ///< Address of the MCParticle
    bool                                          m_isTargetMCPrimary;       ///< Whether this is a target MC primary
    pandora::CaloHitList                          m_mcHits;                  ///< The set of MC hits
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::shared_ptr<LArMCMatchValidationInfo> &LArMCTargetValidationInfo::AddDaughterMatch(const pandora::ParticleFlowObject *const pPfo,
    const bool isRecoCosmicRay, pandora::CaloHitList sharedHits, pandora::CaloHitList pfoHits, const float purity, const float completeness,
    const bool isGoodMatch, const bool isBestMatch)
{
    auto spDaughterMatch = std::shared_ptr<LArMCMatchValidationInfo>(new LArMCMatchValidationInfo(this->shared_from_this(), pPfo,
        isRecoCosmicRay, std::move(sharedHits), std::move(pfoHits), purity, completeness, isGoodMatch, isBestMatch));

    return *m_daughterMatches.emplace(std::move(spDaughterMatch)).first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::shared_ptr<LArInteractionValidationInfo> &LArMCTargetValidationInfo::GetParentInteractionInfo() const noexcept
{
    return m_spParentInteractionInfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArMCTargetValidationInfo::MatchSet &LArMCTargetValidationInfo::GetDaughterMatches() const noexcept
{
    return m_daughterMatches;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArMCTargetValidationInfo::GetMCParticle() const noexcept
{
    return m_pMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArMCTargetValidationInfo::IsTargetMCPrimary() const noexcept
{
    return m_isTargetMCPrimary;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArMCTargetValidationInfo::GetMCHits() const noexcept
{
    return m_mcHits;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_MC_TARGET_VALIDATION_INFO_H