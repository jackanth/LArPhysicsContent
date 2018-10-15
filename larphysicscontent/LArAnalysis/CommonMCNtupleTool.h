/**
 *  @file   larphysicscontent/LArAnalysis/CommonMCNtupleTool.h
 *
 *  @brief  Header file for the common MC ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_COMMON_MC_NTUPLE_TOOL_H
#define LAR_COMMON_MC_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

namespace lar_physics_content
{
/**
 *  @brief  CommonMCNtupleTool class
 */
class CommonMCNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CommonMCNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    CommonMCNtupleTool(const CommonMCNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    CommonMCNtupleTool(CommonMCNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    CommonMCNtupleTool &operator=(const CommonMCNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    CommonMCNtupleTool &operator=(CommonMCNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~CommonMCNtupleTool() = default;

protected:
    void PrepareEvent(const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) override;

    std::vector<LArNtupleRecord> ProcessEvent(
        const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

private:
    using HitSelector = std::function<bool(const pandora::CaloHit *const)>; ///< Alias for a hit selector function
    using HitGetter   = std::function<pandora::CaloHitList()>;              ///< Alias for a hit getter function

    std::string                 m_twoDCaloHitListName; ///< The 2D CaloHit list name
    const pandora::CaloHitList *m_pTwoDCaloHitList;    ///< Address of the 2D CaloHit list

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Produce MC records generic to every considered particle class
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional address of the MC particle
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceGenericPfoMCRecords(
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Produce MC records generic to every considered particle class except neutrinos
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional address of the MC particle
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceNonNeutrinoPfoMCRecords(
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Get the kinetic-energy-weighted contained PFO fraction
     *
     *  @param  pMCParticle optional address of the MC particle
     *
     *  @return the fraction
     */
    float CalculateContainmentFraction(const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Recurse through the MC particle hierarchy to get the amount of contained and total visible kinetic energy
     *
     *  @param  pMCParticle optional address of the MC particle
     *  @param  containedEnergy the contained visible kinetic energy (to be populated)
     *  @param  totalEnergy the total visible kinetic energy (to be populated)
     */
    void RecursivelyGetContainedEnergy(const pandora::MCParticle *const pMCParticle, float &containedEnergy, float &totalEnergy) const;

    /**
     *  @brief  Recurse through the MC particle hierarchy, summing the visible momentum
     *
     *  @param  visibleMomentum the visible momentum (to be populated)
     */
    void RecursivelySumVisibleMomentum(const pandora::MCParticle *const pMCParticle, pandora::CartesianVector &visibleMomentum) const;

    /**
     *  @brief  Get the reco/MC match quality parameters
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the match purity and the match completeness
     */
    std::tuple<float, float> GetRecoMcMatchQuality(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Get the reco/MC match quality parameters for collection-plane hits only
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the match purity and the match completeness
     */
    std::tuple<float, float> GetCollectionPlaneRecoMcMatchQuality(
        const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Get the reco/MC match quality parameters (implementation method)
     *
     *  @param  pMCParticle address of the MC particle
     *  @param  pfoHitGetter the PFO-associated hit-getter function
     *  @parma  hitSelector the hit-selector function
     *
     *  @return the match purity and the match completeness
     */
    std::tuple<float, float> GetRecoMcMatchQualityImpl(
        const pandora::MCParticle *const pMCParticle, const HitGetter &pfoHitGetter, const HitSelector &hitSelector) const;

    /**
     *  @brief  Get the hit weight associated with a given MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *  @parma  hitSelector the hit-selector function
     *
     *  @return the hit weight associated with the MC particle
     */
    float GetHitWeightAssociatedWithMcParticle(const pandora::MCParticle *const pMCParticle, const HitSelector &hitSelector) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float CommonMCNtupleTool::CalculateContainmentFraction(const pandora::MCParticle *const pMCParticle) const
{
    float containedEnergy(0.f), totalEnergy(0.f);
    this->RecursivelyGetContainedEnergy(pMCParticle, containedEnergy, totalEnergy);
    return totalEnergy > std::numeric_limits<float>::epsilon() ? containedEnergy / totalEnergy : 0.f;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_COMMON_MC_NTUPLE_TOOL_H
