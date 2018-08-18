/**
 *  @file   larphysicscontent/LArAnalysis/EventValidationNtupleTool.h
 *
 *  @brief  Header file for the event validation ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_VALIDATION_NTUPLE_TOOL_H
#define LAR_EVENT_VALIDATION_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

namespace lar_physics_content
{
/**
 *  @brief  EventValidationNtupleTool class
 */
class EventValidationNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    EventValidationNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    EventValidationNtupleTool(const EventValidationNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    EventValidationNtupleTool(EventValidationNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    EventValidationNtupleTool &operator=(const EventValidationNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    EventValidationNtupleTool &operator=(EventValidationNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~EventValidationNtupleTool() = default;

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
};

} // namespace lar_physics_content

#endif // #ifndef LAR_EVENT_VALIDATION_NTUPLE_TOOL_H
