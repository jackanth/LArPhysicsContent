/**
 *  @file   larphysicscontent/LArObjects/LArMCMatchValidationInfo.h
 *
 *  @brief  Header file for the lar MC match validation info class.
 *
 *  $Log: $
 */
#ifndef LAR_MC_MATCH_VALIDATION_INFO_H
#define LAR_MC_MATCH_VALIDATION_INFO_H 1

#include "Objects/MCParticle.h"
#include "Objects/ParticleFlowObject.h"

#include <memory>

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the LArMCTargetValidationInfo class
 */
class LArMCTargetValidationInfo;

/**
 *  @brief  LArMCMatchValidationInfo class
 */
class LArMCMatchValidationInfo : public std::enable_shared_from_this<LArMCMatchValidationInfo>
{
public:
    /**
     *  @brief  Default copy constructor
     */
    LArMCMatchValidationInfo(const LArMCMatchValidationInfo &) = default;

    /**
     *  @brief  Default move constructor
     */
    LArMCMatchValidationInfo(LArMCMatchValidationInfo &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    LArMCMatchValidationInfo &operator=(const LArMCMatchValidationInfo &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    LArMCMatchValidationInfo &operator=(LArMCMatchValidationInfo &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~LArMCMatchValidationInfo() = default;

    /**
     *  @brief  Get the parent target
     *
     *  @return shared pointer to the parent target
     */
    const std::shared_ptr<LArMCTargetValidationInfo> &ParentTarget() const noexcept;

    /**
     *  @brief  Get the PFO
     *
     *  @return address of the PFO
     */
    const pandora::Pfo *GetPfo() const noexcept;

    /**
     *  @brief  Get whether this is a reco cosmic ray
     *
     *  @return whether this is a reco cosmic ray
     */
    bool IsRecoCosmicRay() const noexcept;

    /**
     *  @brief  Get the shared hits
     *
     *  @return the shared hits
     */
    const pandora::CaloHitList &GetSharedHits() const noexcept;

    /**
     *  @brief  Get the PFO hits
     *
     *  @return the PFO hits
     */
    const pandora::CaloHitList &GetPfoHits() const noexcept;

    /**
     *  @brief  Get the purity
     *
     *  @return the purity
     */
    float GetPurity() const noexcept;

    /**
     *  @brief  Get the completeness
     *
     *  @return the completeness
     */
    float GetCompleteness() const noexcept;

    /**
     *  @brief  Get whether this is a good match
     *
     *  @return whether this is a good match
     */
    bool IsGoodMatch() const noexcept;

    /**
     *  @brief  Get whether this is the best match
     *
     *  @return whether this is the best match
     */
    bool IsBestMatch() const noexcept;

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  spParentTarget shared pointer to the parent target
     *  @param  pPfo address of the PFO
     *  @param  isRecoCosmicRay whether this is a reco cosmic ray
     *  @param  sharedHits the set of shared hits
     *  @param  pfoHits the set of PFO hits
     *  @param  purity the match purity
     *  @param  completeness the match completeness
     *  @param  isGoodMatch whether this is a good match
     *  @param  isBestMatch whether it is the best match
     */
    LArMCMatchValidationInfo(std::shared_ptr<LArMCTargetValidationInfo> spParentTarget, const pandora::ParticleFlowObject *const pPfo,
        const bool isRecoCosmicRay, pandora::CaloHitList sharedHits, pandora::CaloHitList pfoHits, const float purity,
        const float completeness, const bool isGoodMatch, const bool isBestMatch) noexcept;

    friend class LArMCTargetValidationInfo;

private:
    std::shared_ptr<LArMCTargetValidationInfo> m_spParentTarget;  ///< Shared pointer to the parent target
    const pandora::Pfo *                       m_pPfo;            ///< Address of the PFO
    bool                                       m_isRecoCosmicRay; ///< Whether this is a reco cosmic ray
    pandora::CaloHitList                       m_sharedHits;      ///< The set of shared hits
    pandora::CaloHitList                       m_pfoHits;         ///< The set of PFO hits
    float                                      m_purity;          ///< The match purity
    float                                      m_completeness;    ///< The match completeness
    bool                                       m_isGoodMatch;     ///< Whether it is a good match
    bool                                       m_isBestMatch;     ///< Whether it is the best match
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::shared_ptr<LArMCTargetValidationInfo> &LArMCMatchValidationInfo::ParentTarget() const noexcept
{
    return m_spParentTarget;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Pfo *LArMCMatchValidationInfo::GetPfo() const noexcept
{
    return m_pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArMCMatchValidationInfo::IsRecoCosmicRay() const noexcept
{
    return m_isRecoCosmicRay;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArMCMatchValidationInfo::GetSharedHits() const noexcept
{
    return m_sharedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArMCMatchValidationInfo::GetPfoHits() const noexcept
{
    return m_pfoHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArMCMatchValidationInfo::GetPurity() const noexcept
{
    return m_purity;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArMCMatchValidationInfo::GetCompleteness() const noexcept
{
    return m_completeness;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArMCMatchValidationInfo::IsGoodMatch() const noexcept
{
    return m_isGoodMatch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArMCMatchValidationInfo::IsBestMatch() const noexcept
{
    return m_isBestMatch;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_MC_MATCH_VALIDATION_INFO_H