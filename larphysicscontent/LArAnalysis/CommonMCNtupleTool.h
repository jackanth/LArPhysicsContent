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
    std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Produce MC records generic to every considered particle class
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional address of the MC particle
     *  @param  pMCParticleList optional list of all MC particles
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceGenericPfoMCRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Produce MC records generic to every considered particle class except neutrinos
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional address of the MC particle
     *  @param  pMCParticleList optional list of all MC particles
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceNonNeutrinoPfoMCRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) const;

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
