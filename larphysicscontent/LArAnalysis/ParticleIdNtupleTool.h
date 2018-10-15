/**
 *  @file   larphysicscontent/LArAnalysis/ParticleIdNtupleTool.h
 *
 *  @brief  Header file for the particle ID ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_PARTICLE_ID_NTUPLE_TOOL_H
#define LAR_PARTICLE_ID_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

namespace lar_physics_content
{
/**
 *  @brief  ParticleIdNtupleTool class
 */
class ParticleIdNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ParticleIdNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    ParticleIdNtupleTool(const ParticleIdNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    ParticleIdNtupleTool(ParticleIdNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    ParticleIdNtupleTool &operator=(const ParticleIdNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    ParticleIdNtupleTool &operator=(ParticleIdNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~ParticleIdNtupleTool() = default;

protected:
    std::vector<LArNtupleRecord> ProcessEvent(
        const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_PARTICLE_ID_NTUPLE_TOOL_H
