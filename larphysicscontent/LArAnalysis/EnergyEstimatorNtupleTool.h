/**
 *  @file   larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.h
 *
 *  @brief  Header file for the energy estimator ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
#define LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H 1

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
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
    bool m_writeEnergiesToNtuple; ///< Whether to write the energies to the ntuple
    bool m_useParticleId;         ///< Whether to use particle ID
    bool m_trainingSetMode;       ///< Whether to run in training set mode

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get a 3D distance from a cell
     *
     *  @param  hitWidth the hit width
     *  @param  wirePitch the wire pitch
     *  @param  fitDirection the fit direction at the hit
     *
     *  @return the 3D distance
     */
    float CellToThreeDDistance(const float hitWidth, const float wirePitch, const pandora::CartesianVector &fitDirection) const;

    /**
     *  @brief  Get the 3D distance from a CaloHit
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pCaloHit address of the CaloHit
     *  @param  trackFit the 3D track fit to which the hit belongs
     *
     *  @return the 3D distance
     */
    float CaloHitToThreeDDistance(const pandora::CaloHit *const pCaloHit, const lar_content::ThreeDSlidingFitResult &trackFit) const;

    /**
     *  @brief  Turn a direction in polar and azimuthal angles
     *
     *  @param  direction the Cartesian direction
     *
     *  @return the polar and azimuthal angles
     */
    std::tuple<float, float> GetPolarAnglesFromDirection(const pandora::CartesianVector &direction) const;

    /**
     *  @brief  Produce training records for a cosmic or primary PFO
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional pointer to the corresponding MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the training records
     */
    std::vector<LArNtupleRecord> ProduceTrainingRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &,
        const pandora::MCParticle *const pMcParticle, const pandora::MCParticleList *const);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
