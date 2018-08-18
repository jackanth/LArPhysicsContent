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
};

} // namespace lar_physics_content

#endif // #ifndef LAR_COMMON_MC_NTUPLE_TOOL_H
