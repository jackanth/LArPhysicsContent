/**
 *  @file   larphysicscontent/LArObjects/LArInteractionValidationInfo.h
 *
 *  @brief  Header file for the lar interaction validation info class.
 *
 *  $Log: $
 */

#ifndef LAR_INTERACTION_VALIDATION_INFO_H
#define LAR_INTERACTION_VALIDATION_INFO_H 1

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larphysicscontent/LArObjects/LArMCTargetValidationInfo.h"

#include <memory>
#include <set>

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the EventValidationTool class
 */
class EventValidationTool;

/**
 *  @brief  LArInteractionValidationInfo class
 */
class LArInteractionValidationInfo : public std::enable_shared_from_this<LArInteractionValidationInfo>
{
public:
    using TargetSet = std::set<std::shared_ptr<LArMCTargetValidationInfo>, std::owner_less<std::shared_ptr<LArMCTargetValidationInfo>>>; ///< Alias for a set of shared pointers to targets

    /**
     *  @brief  Default copy constructor
     */
    LArInteractionValidationInfo(const LArInteractionValidationInfo &) = default;

    /**
     *  @brief  Default move constructor
     */
    LArInteractionValidationInfo(LArInteractionValidationInfo &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    LArInteractionValidationInfo &operator=(const LArInteractionValidationInfo &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    LArInteractionValidationInfo &operator=(LArInteractionValidationInfo &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~LArInteractionValidationInfo() = default;

    /**
     *  @brief  Add a daughter target
     *
     *  @param  pMCPrimary address of the MC primary
     *  @param  mcNuanceCode the nuance code
     *  @param  mcPrimaryHitList the hit list
     *  @param  isTargetMCPrimary whether this is a target MC primary
     *
     *  @return shared pointer to the daughter target
     */
    const std::shared_ptr<LArMCTargetValidationInfo> &AddDaughterTarget(
        const pandora::MCParticle *const pMCPrimary, const pandora::CaloHitList &mcPrimaryHitList, const bool isTargetMCPrimary);

    /**
     *  @brief  Get the daughter targets
     *
     *  @return the daughter targets
     */
    const TargetSet &GetDaughterTargets() const noexcept;

    /**
     *  @brief  Get the nuance code
     *
     *  @return the nuance code
     */
    int GetNuanceCode() const noexcept;

    /**
     *  @brief  Get whether this is a cosmic ray
     *
     *  @return whether this is a cosmic ray
     */
    bool IsCosmicRay() const noexcept;

    /**
     *  @brief  Get the interaction type
     *
     *  @return the daughter targets
     */
    lar_content::LArInteractionTypeHelper::InteractionType GetInteractionType() const;

    /**
     *  @brief  Get whether this is correct
     *
     *  @return whether this is correct
     */
    bool IsCorrect() const;

    /**
     *  @brief  Get whether this is fake
     *
     *  @return whether this is fake
     */
    bool IsFake() const;

    /**
     *  @brief  Get whether this is split
     *
     *  @return whether this is split
     */
    bool IsSplit() const;

    /**
     *  @brief  Get whether this is lost
     *
     *  @return whether this is lost
     */
    bool IsLost() const;

    /**
     *  @brief  Get the address of the reco neutrino, if appropriate
     *
     *  @return address of the reco neutrinos
     */
    const pandora::ParticleFlowObject *GetRecoNeutrino() const;

    /**
     *  @brief  Get the address of the MC neutrino, if appropriate
     *
     *  @return address of the MC neutrinos
     */
    const pandora::MCParticle *GetMcNeutrino() const;

    /**
     *  @brief  Get whether the parameters are set
     *
     *  @return whether the parameters are set
     */
    bool AreParametersSet() const noexcept;

    /**
     *  @brief  Set the parameters
     *
     *  @param  interactionType the interaction type
     *  @param  isCorrect whether it is correct
     *  @param  isFake whether it is fake
     *  @param  isSplit whether it is split
     *  @param  isLost whether it is lost
     *  @param  pRecoNeutrino address of the reco neutrino
     *  @param  pMcNeutrino address of the MC neutrino
     */
    void SetParameters(const lar_content::LArInteractionTypeHelper::InteractionType interactionType, const bool isCorrect, const bool isFake,
        const bool isSplit, const bool isLost, const pandora::ParticleFlowObject *const pRecoNeutrino, const pandora::MCParticle *const pMcNeutrino);

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  mcNuanceCode the nuance code
     *  @param  isCosmicRay whether it is a cosmic ray
     */
    LArInteractionValidationInfo(const int mcNuanceCode, const bool isCosmicRay) noexcept;

    friend class EventValidationTool;

private:
    TargetSet                                              m_daughterTargets; ///< The daughter targets
    int                                                    m_mcNuanceCode;    ///< The nuance code
    bool                                                   m_isCosmicRay;     ///< Whether this is a cosmic ray
    lar_content::LArInteractionTypeHelper::InteractionType m_interactionType; ///< The interaction type
    bool                                                   m_isCorrect;       ///< Whether this is correct
    bool                                                   m_isFake;          ///< Whether this is fake
    bool                                                   m_isSplit;         ///< Whether this is split
    bool                                                   m_isLost;          ///< Whether this is lost
    const pandora::ParticleFlowObject *                    m_pRecoNeutrino;   ///< Address of the reco neutrino, if appropriate
    const pandora::MCParticle *                            m_pMcNeutrino;     ///< Address of the reco neutrino, if appropriate
    bool                                                   m_parametersSet;   ///< Whether the parameters have been set
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArInteractionValidationInfo::TargetSet &LArInteractionValidationInfo::GetDaughterTargets() const noexcept
{
    return m_daughterTargets;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArInteractionValidationInfo::GetNuanceCode() const noexcept
{
    return m_mcNuanceCode;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArInteractionValidationInfo::IsCosmicRay() const noexcept
{
    return m_isCosmicRay;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArInteractionValidationInfo::AreParametersSet() const noexcept
{
    return m_parametersSet;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_INTERACTION_VALIDATION_INFO_H