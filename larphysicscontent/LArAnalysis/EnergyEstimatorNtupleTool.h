/**
 *  @file   larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.h
 *
 *  @brief  Header file for the energy estimator ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
#define LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArHelpers/LArRootHelper.h"
#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "bethe-faster/BetheFaster.h"

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
    std::vector<LArNtupleRecord> ProcessEvent(
        const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

private:
    /**
     *  @brief  Struct containing hit calorimetry info
     */
    struct HitCalorimetryInfo
    {
        /**
         *  @brief  Constructor
         */
        HitCalorimetryInfo();

        bool                     m_projectionSuccessful; ///< Whether the projection was successful
        pandora::CartesianVector m_threeDPosition;       ///< If projection successful, the 3D position
        float                    m_projectionError;      ///< If projection successful, the projection error
        double                   m_coordinate;           ///< If projection successful, the 3D track coordinate
        double                   m_dQ;                   ///< The hit charge
        double                   m_dX;                   ///< The 3D dx
        bool                     m_isNoise;              ///< Whether this hit is tagged as noisy
        bool                     m_isFake;               ///< Whether this is a fake hit
    };

    using DoubleVector          = std::vector<double>;                 ///< Alias for a vector of doubles
    using HitCalorimetryInfoPtr = std::shared_ptr<HitCalorimetryInfo>; ///< Alias for a shared to a HitCalorimetryInfo object
    using HitCalorimetryInfoMap =
        std::unordered_map<const pandora::CaloHit *, HitCalorimetryInfoPtr>; ///< Alias for a map from CaloHits to HitCalorimetryInfo shared pointers
    using CaloHitMap = std::unordered_map<const pandora::CaloHit *, const pandora::CaloHit *>; ///< Alias for a map between CaloHits

    template <typename T>
    using HitDataGetter = std::function<std::optional<std::decay_t<T>>(const pandora::CaloHit *const, const HitCalorimetryInfo &hitInfo)>; ///< Alias for a CaloHit info getter function

    bool  m_braggGradientTrainingMode; ///< Whether to run in Bragg gradient training mode
    bool  m_makePlots;                 ///< Make plots
    float m_modboxRho;                 ///< The ModBox rho parameter
    float m_modboxA;                   ///< The ModBox A parameter
    float m_modboxB;                   ///< The ModBox B parameter
    float m_modboxEpsilon;             ///< The ModBox epsilon parameter
    float m_modboxWion;                ///< The ModBox W_ion parameter
    float m_modboxC;                   ///< The ModBox C parameter
    float m_modboxFactor;              ///< The ModBox (rho * epsilon / B) value

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get energy estimator records for a PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticle optional address of the MC particle
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> GetEnergyEstimatorRecords(
        const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Get CaloHit calorimetry info
     *
     *  @param  pPfo address of the PFO
     *  @param  pMcParticle optional address of the MCParticle
     *
     *  @return the dQ/dx vector, the dx vector, the shower charge, and the number of hits lost to errors
     */
    std::tuple<LArNtupleRecord::RFloatVector, LArNtupleRecord::RFloatVector, LArNtupleRecord::RFloat, LArNtupleRecord::RUInt> GetHitCalorimetryInfo(
        const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMcParticle) const;

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
     *  @brief  Get the 3D distance from CaloHit info
     *
     *  @param  hitWidth the instance of Pandora
     *  @param  trackFit the 3D track fit to which the hit belongs
     *  @param  threeDPosition the 3D position
     *  @param  threeDDistance the 3D distance (to populate)
     *
     *  @return the status code
     */
    pandora::StatusCode CaloHitToThreeDDistance(const float hitWidth, const lar_content::ThreeDSlidingFitResult &trackFit,
        const pandora::CartesianVector &threeDPosition, float &threeDDistance) const;

    /**
     *  @brief  Turn a direction in polar and azimuthal angles
     *
     *  @param  direction the Cartesian direction
     *
     *  @return the polar and azimuthal angles
     */
    std::tuple<float, float> GetPolarAnglesFromDirection(const pandora::CartesianVector &direction) const;

    /**
     *  @brief  Produce Bragg gradient training records for a PFO
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional pointer to the corresponding MCParticle
     *
     *  @return the training records
     */
    std::vector<LArNtupleRecord> ProduceBraggGradientTrainingRecords(
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const pandora::MCParticle *const pMcParticle);

    /**
     *  @brief  Get the Bragg gradient parameters
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pMcParticle optional pointer to the corresponding MCParticle
     *  @param  firstOrderGradient the first-order gradient parameter (to populate)
     *  @param  firstOrderIntercept the first-order intercept parameter (to populate)
     *  @param  secondOrderGradient the second-order gradient parameter (to populate)
     *  @param  secondOrderIntercept the second-order intercept parameter (to populate)
     *  @param  averageDetectorThickness the average detector thickness (to populate)
     *  @param  pida the PIDA value
     *
     *  @return success
     */
    bool GetBraggGradientParameters(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMcParticle,
        float &firstOrderGradient, float &firstOrderIntercept, float &secondOrderGradient, float &secondOrderIntercept, float &averageDetectorThickness, float &pida, float &medianFilteredEnergyLossRate, float &medianUnfilteredEnergyLossRate) const;

    /**
     *  @brief  Get the dE/dx distribution
     *
     *  @param  trackFit the 3D track fit
     *  @param  caloHitList the CaloHit list
     *  @param  caloHitMap the map from 2D to 3D CaloHits
     *  @param  isBackwards whether the particle is going backwards
     *  @param  pMcParticle address of the MC particle
     *
     *  @return the vector of hit charge objects
     */
    std::vector<bf::HitCharge> GetdEdxDistribution(const lar_content::ThreeDSlidingFitResult &trackFit, pandora::CaloHitList caloHitList,
        const CaloHitMap &caloHitMap, const bool isBackwards, const pandora::MCParticle * const pMcParticle) const;

    /**
     *  @brief  Generate the hit calorimetry info map
     *
     *  @param  trackFit the track fit object
     *  @param  caloHitList the CaloHit list
     *  @param  caloHitMap the map from 2D to 3D CaloHits
     *  @param  isBackwards whether the particle is going backwards
     *
     *  @return the hit calorimetry info vector
     */
    std::vector<HitCalorimetryInfoPtr> CalculateHitCalorimetryInfo(const lar_content::ThreeDSlidingFitResult &trackFit,
        const pandora::CaloHitList &caloHitList, const CaloHitMap &caloHitMap, const bool isBackwards) const;

    /**
     *  @brief  Calculate the hit calorimetry info
     *
     *  @param  trackFit the track fit object
     *  @param  twoDPositionVector the 2D position
     *  @param  dQ the hit charge
     *  @param  hitWidth the hit width
     *  @param  caloHitMap the map from 2D to 3D CaloHits
     *  @param  pCaloHit address of the CaloHit
     *
     *  @return shared pointer to the hit calorimetry info object
     */
    HitCalorimetryInfoPtr CalculateHitCalorimetryInfo(const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CartesianVector &twoDPositionVector,
        const float dQ, const float hitWidth, const CaloHitMap &caloHitMap, const pandora::CaloHit *const pCaloHit) const;

    /**
     *  @brief  Esimate the energy of a track hit
     *
     *  @param  dQdX the dQ/dx value
     *  @param  dX the dx value
     *
     *  @return the estimated energy
     */
    float EstimateTrackHitEnergy(const float dQdX, const float dX) const;

    /**
     *  @brief  Estimate the energy of some shower charge
     *
     *  @param  showerCharge the shower charge
     *
     *  @return the estimated energy
     */
    float EstimateShowerEnergy(const float showerCharge) const;

    /**
     *  @brief  Apply the ModBox correction to map dQ/dx to dE/dx
     *
     *  @param  dQdx the dQ/dx value
     *
     *  @return the dE/dx value
     */
    float ApplyModBoxCorrection(const float dQdx) const;

    /**
     *  @brief  Apply a scaling from charge to energy in the absence of recombination
     *
     *  @param  dQdx the dQ/dx value
     *
     *  @return the dE/dx value
     */
    float ApplyChargeScaling(const float dQdx) const;

    float CorrectChargeDeposition(const float dQdxUncorrected, const pandora::CartesianVector &threeDPosition) const;

    float GetDriftCoordinateCorrection(const float xPosition) const;

    float GetYZCoordinateCorrection(const float yPosition, const float zPosition) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float EnergyEstimatorNtupleTool::ApplyModBoxCorrection(const float dQdx) const
{
    return m_modboxFactor * (std::exp(this->ApplyChargeScaling(dQdx) / m_modboxFactor) - m_modboxA);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float EnergyEstimatorNtupleTool::ApplyChargeScaling(const float chargeValue) const
{
    return chargeValue * m_modboxWion / m_modboxC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float EnergyEstimatorNtupleTool::EstimateTrackHitEnergy(const float dQdX, const float dX) const
{
    return this->ApplyModBoxCorrection(dQdX) * dX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float EnergyEstimatorNtupleTool::EstimateShowerEnergy(const float showerCharge) const
{
    return this->ApplyChargeScaling(showerCharge);
}

} // namespace lar_physics_content

#endif // #ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
