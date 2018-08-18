/**
 *  @file   larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.h
 *
 *  @brief  Header file for the energy estimator ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
#define LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

namespace lar_physics_content
{
/**
 *  @brief  EnergyEstimatorNtupleTool class
 */
class EnergyEstimatorNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    EnergyEstimatorNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    EnergyEstimatorNtupleTool(const EnergyEstimatorNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    EnergyEstimatorNtupleTool(EnergyEstimatorNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    EnergyEstimatorNtupleTool &operator=(const EnergyEstimatorNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    EnergyEstimatorNtupleTool &operator=(EnergyEstimatorNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~EnergyEstimatorNtupleTool() = default;

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

#endif // #ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
