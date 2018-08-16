/**
 *  @file   larphysicscontent/LArAnalysis/DefaultNtupleTool.h
 *
 *  @brief  Header file for the default ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_DEFAULT_NTUPLE_TOOL_H
#define LAR_DEFAULT_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief  DefaultNtupleTool class
 */
class DefaultNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DefaultNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    DefaultNtupleTool(const DefaultNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    DefaultNtupleTool(DefaultNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    DefaultNtupleTool &operator=(const DefaultNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    DefaultNtupleTool &operator=(DefaultNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~DefaultNtupleTool() = default;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessParticle(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_DEFAULT_NTUPLE_TOOL_H
