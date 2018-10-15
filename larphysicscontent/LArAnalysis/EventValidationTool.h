/**
 *  @file   larphysicscontent/LArAnalysis/EventValidationTool.h
 *
 *  @brief  Header file for the event validation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_VALIDATION_TOOL_H
#define LAR_EVENT_VALIDATION_TOOL_H 1

#include "larphysicscontent/LArObjects/LArInteractionValidationInfo.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include <memory>

namespace lar_physics_content
{
/**
 *  @brief  EventValidationTool class
 */
class EventValidationTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Constructor
     */
    EventValidationTool();

    /**
     *  @brief  Default copy constructor
     */
    EventValidationTool(const EventValidationTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    EventValidationTool(EventValidationTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    EventValidationTool &operator=(const EventValidationTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    EventValidationTool &operator=(EventValidationTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~EventValidationTool() = default;

    /**
     *  @brief  Run the validation
     *
     *  @param  pfoList the PFO list
     *  @param  caloHitList the CaloHit list
     *  @param  mcParticleList the MCParticle list
     *
     *  @return vector of shared pointers to the interaction validation info
     */
    std::vector<std::shared_ptr<LArInteractionValidationInfo>> RunValidation(
        const pandora::PfoList &pfoList, const pandora::CaloHitList &caloHitList, const pandora::MCParticleList &mcParticleList) const;

private:
    using PfoToIdMap = std::unordered_map<const pandora::ParticleFlowObject *, unsigned int>; ///< Alias for a map from PFOs to IDs

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  ValidationInfo class
     */
    class ValidationInfo
    {
    public:
        /**
         *  @brief  Get the all mc particle to hits map
         *
         *  @return the all mc particle to hits map
         */
        const lar_content::LArMCParticleHelper::MCContributionMap &GetAllMCParticleToHitsMap() const;

        /**
         *  @brief  Get the target mc particle to hits map
         *
         *  @return the target mc particle to hits map
         */
        const lar_content::LArMCParticleHelper::MCContributionMap &GetTargetMCParticleToHitsMap() const;

        /**
         *  @brief  Get the pfo to hits map
         *
         *  @return the pfo to hits map
         */
        const lar_content::LArMCParticleHelper::PfoContributionMap &GetPfoToHitsMap() const;

        /**
         *  @brief  Get the mc to pfo hit sharing map
         *
         *  @return the mc to pfo hit sharing map
         */
        const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &GetMCToPfoHitSharingMap() const;

        /**
         *  @brief  Get the interpreted mc to pfo hit sharing map
         *
         *  @return the interpreted mc to pfo hit sharing map
         */
        const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &GetInterpretedMCToPfoHitSharingMap() const;

        /**
         *  @brief  Set the all mc particle to hits map
         *
         *  @param  allMCParticleToHitsMap the all mc particle to hits map
         */
        void SetAllMCParticleToHitsMap(const lar_content::LArMCParticleHelper::MCContributionMap &allMCParticleToHitsMap);

        /**
         *  @brief  Set the target mc particle to hits map
         *
         *  @param  targetMCParticleToHitsMap the target mc particle to hits map
         */
        void SetTargetMCParticleToHitsMap(const lar_content::LArMCParticleHelper::MCContributionMap &targetMCParticleToHitsMap);

        /**
         *  @brief  Set the pfo to hits map
         *
         *  @param  pfoToHitsMap the pfo to hits map
         */
        void SetPfoToHitsMap(const lar_content::LArMCParticleHelper::PfoContributionMap &pfoToHitsMap);

        /**
         *  @brief  Set the mc to pfo hit sharing map
         *
         *  @param  mcToPfoHitSharingMap the mc to pfo hit sharing map
         */
        void SetMCToPfoHitSharingMap(const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap);

        /**
         *  @brief  Set the interpreted mc to pfo hit sharing map
         *
         *  @param  interpretedMCToPfoHitSharingMap the interpreted mc to pfo hit sharing map
         */
        void SetInterpretedMCToPfoHitSharingMap(const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap);

    private:
        lar_content::LArMCParticleHelper::MCContributionMap            m_allMCParticleToHitsMap;    ///< The all mc particle to hits map
        lar_content::LArMCParticleHelper::MCContributionMap            m_targetMCParticleToHitsMap; ///< The target mc particle to hits map
        lar_content::LArMCParticleHelper::PfoContributionMap           m_pfoToHitsMap;              ///< The pfo to hits map
        lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap m_mcToPfoHitSharingMap;      ///< The mc to pfo hit sharing map
        lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap m_interpretedMCToPfoHitSharingMap; ///< The interpreted mc to pfo hit sharing map
    };

    /**
     *  @brief  Fill the validation info containers
     *
     *  @param  mcParticleList the mc particle list
     *  @param  caloHitList the calo hit list
     *  @param  pfoList the pfo list
     *  @param  validationInfo to receive the validation info
     */
    void FillValidationInfo(const pandora::MCParticleList &mcParticleList, const pandora::CaloHitList &caloHitList,
        const pandora::PfoList &pfoList, ValidationInfo &validationInfo) const;

    /**
     *  @brief  Apply an interpretative matching procedure to the comprehensive matches in the provided validation info object
     *
     *  @param  validationInfo the validation info
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     */
    void InterpretMatching(const ValidationInfo &                       validationInfo,
        lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Get the strongest pfo match (most matched hits) between an available mc primary and an available pfo
     *
     *  @param  validationInfo the validation info
     *  @param  mcPrimaryVector the mc primary vector
     *  @param  usedPfos the set of previously used pfos
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     *
     *  @return whether a strong match was identified
     */
    bool GetStrongestPfoMatch(const ValidationInfo &validationInfo, const pandora::MCParticleVector &mcPrimaryVector,
        pandora::PfoSet &usedPfos, lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Get the best matches for any pfos left-over after the strong matching procedure
     *
     *  @param  validationInfo the validation info
     *  @param  mcPrimaryVector the mc primary vector
     *  @param  usedPfos the set of previously used pfos
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     */
    void GetRemainingPfoMatches(const ValidationInfo &validationInfo, const pandora::MCParticleVector &mcPrimaryVector,
        const pandora::PfoSet &usedPfos, lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Create the event validation info object
     *
     *  @param  validationInfo the validation info
     *
     *  @return the vector of shared pointers to the interaction validation objects
     */
    std::vector<std::shared_ptr<LArInteractionValidationInfo>> GetEventValidationInfo(const ValidationInfo &validationInfo) const;

    /**
     *  @brief  Test if a match is good
     *
     *  @param  trueHits the true hits
     *  @param  recoHits the reco hits
     *  @param  sharedHits the shared hits
     *
     *  @return whether it is good
     */
    bool IsGoodMatch(const pandora::CaloHitList &trueHits, const pandora::CaloHitList &recoHits, const pandora::CaloHitList &sharedHits) const;

    /**
     *  @brief  Get the purity
     *
     *  @param  recoHits the reco hits
     *  @param  sharedHits the shared hits
     *
     *  @return the purity
     */
    float GetPurity(const pandora::CaloHitList &recoHits, const pandora::CaloHitList &sharedHits) const;

    /**
     *  @brief  Get the completeness
     *
     *  @param  trueHits the true hits
     *  @param  sharedHits the shared hits
     *
     *  @return the completeness
     */
    float GetCompleteness(const pandora::CaloHitList &trueHits, const pandora::CaloHitList &sharedHits) const;

    bool         m_useTrueNeutrinosOnly;    ///< Whether to consider only mc particles that were neutrino induced
    bool         m_selectInputHits;         ///< Whether to use only hits passing mc-based quality (is "reconstructable") checks
    float        m_minHitSharingFraction;   ///< Minimum fraction of energy deposited by selected primary in a single "good" hit
    float        m_maxPhotonPropagation;    ///< Maximum distance travelled by photon, downstream of a track, in mc particle hierarchy
    bool         m_useSmallPrimaries;       ///< Whether to consider matches to mc primaries with fewer than m_matchingMinPrimaryHits
    unsigned int m_matchingMinSharedHits;   ///< The minimum number of shared hits used in matching scheme
    float        m_matchingMinCompleteness; ///< The minimum particle completeness to declare a match
    float        m_matchingMinPurity;       ///< The minimum particle purity to declare a match
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const lar_content::LArMCParticleHelper::MCContributionMap &EventValidationTool::ValidationInfo::GetAllMCParticleToHitsMap() const
{
    return m_allMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const lar_content::LArMCParticleHelper::MCContributionMap &EventValidationTool::ValidationInfo::GetTargetMCParticleToHitsMap() const
{
    return m_targetMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const lar_content::LArMCParticleHelper::PfoContributionMap &EventValidationTool::ValidationInfo::GetPfoToHitsMap() const
{
    return m_pfoToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &EventValidationTool::ValidationInfo::GetMCToPfoHitSharingMap() const
{
    return m_mcToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &EventValidationTool::ValidationInfo::GetInterpretedMCToPfoHitSharingMap() const
{
    return m_interpretedMCToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationTool::ValidationInfo::SetAllMCParticleToHitsMap(const lar_content::LArMCParticleHelper::MCContributionMap &allMCParticleToHitsMap)
{
    m_allMCParticleToHitsMap = allMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationTool::ValidationInfo::SetTargetMCParticleToHitsMap(const lar_content::LArMCParticleHelper::MCContributionMap &targetMCParticleToHitsMap)
{
    m_targetMCParticleToHitsMap = targetMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationTool::ValidationInfo::SetPfoToHitsMap(const lar_content::LArMCParticleHelper::PfoContributionMap &pfoToHitsMap)
{
    m_pfoToHitsMap = pfoToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationTool::ValidationInfo::SetMCToPfoHitSharingMap(const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap)
{
    m_mcToPfoHitSharingMap = mcToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationTool::ValidationInfo::SetInterpretedMCToPfoHitSharingMap(
    const lar_content::LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap)
{
    m_interpretedMCToPfoHitSharingMap = interpretedMCToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float EventValidationTool::GetPurity(const pandora::CaloHitList &recoHits, const pandora::CaloHitList &sharedHits) const
{
    return (recoHits.size() > 0UL) ? static_cast<float>(sharedHits.size()) / static_cast<float>(recoHits.size()) : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float EventValidationTool::GetCompleteness(const pandora::CaloHitList &trueHits, const pandora::CaloHitList &sharedHits) const
{
    return (trueHits.size() > 0UL) ? static_cast<float>(sharedHits.size()) / static_cast<float>(trueHits.size()) : 0.f;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_EVENT_VALIDATION_TOOL_H
