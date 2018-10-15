/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisHelper.h
 *
 *  @brief  Header file for the lar analysis helper class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_HELPER_H
#define LAR_ANALYSIS_HELPER_H 1

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larphysicscontent/LArObjects/LArInteractionValidationInfo.h"
#include "larphysicscontent/LArObjects/LArMCTargetValidationInfo.h"

#include "Objects/CartesianVector.h"
#include "Objects/MCParticle.h"
#include "Pandora/Pandora.h"
#include "Pandora/PdgTable.h"

#include <tuple>

namespace lar_physics_content
{
/**
 *  @brief  LArAnalysisHelper class
 */
class LArAnalysisHelper
{
public:
    using PfoToTargetMap =
        std::unordered_map<const pandora::ParticleFlowObject *, std::shared_ptr<LArMCTargetValidationInfo>>; ///< Alias for a map from PFOs to MC targets
    using PfoToInteractionMap = std::unordered_map<const pandora::ParticleFlowObject *, std::shared_ptr<LArInteractionValidationInfo>>; ///< Alias for a map from PFOs to MC interactions

    /**
     *  @brief  Deleted copy constructor
     */
    LArAnalysisHelper(const LArAnalysisHelper &) = delete;

    /**
     *  @brief  Deleted move constructor
     */
    LArAnalysisHelper(LArAnalysisHelper &&) = delete;

    /**
     *  @brief  Deleted copy assignment operator
     */
    LArAnalysisHelper &operator=(const LArAnalysisHelper &) = delete;

    /**
     *  @brief  Deleted move assignment operator
     */
    LArAnalysisHelper &operator=(LArAnalysisHelper &&) = delete;

    /**
     *  @brief  Deleted destructor
     */
    ~LArAnalysisHelper() = delete;

    /**
     *  @brief  Get the minimum and maximum fiducial cut coordinates
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  fiducialCutLowMargins the fiducial cut low margins
     *  @param  fiducialCutHighMargins the fiducial cut high margins
     *
     *  @return the minimum and maximum fiducial cut coordinates
     */
    static std::tuple<pandora::CartesianVector, pandora::CartesianVector> GetFiducialCutCoordinates(const pandora::Pandora &pandoraInstance,
        const pandora::CartesianVector &fiducialCutLowMargins, const pandora::CartesianVector &fiducialCutHighMargins);

    /**
     *  @brief  Find out whether a given point lies in the fiducial region of the detector
     *
     *  @param  point the point
     *  @param  minCoordinates the minimum fiducial volume coordinates
     *  @param  maxCoordinates the maximum fiducial volume coordinates
     *
     *  @return whether the point is fiducial
     */
    static bool IsPointFiducial(const pandora::CartesianVector &point, const pandora::CartesianVector &minCoordinates,
        const pandora::CartesianVector &maxCoordinates);

    /**
     *  @brief  Get the true kinetic energy for an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the true kinetic energy (GeV)
     */
    static float GetTrueKineticEnergy(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the true mass for an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the true mass (GeV)
     */
    static float GetTrueMass(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Find out whether an MC particle is a true shower
     *
     *  @param  pMCParticle address of the MC particle
     *
     *  @return whether it is a true shower
     */
    static bool IsTrueShower(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Project a 2D position to a 3D position using a track fit
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  trackFit the track fit object
     *  @param  twoDPosition the 2D position
     *  @param  hitType the hit type
     *  @param  linearlyProjectEnds whether to linearly project the ends of the fit
     *  @param  threeDPosition the 3D position (to populate)
     *  @param  projectionError the projection error (to populate)
     *
     *  @return the status code
     */
    static pandora::StatusCode ProjectTwoDPositionOntoTrackFit(const pandora::Pandora &pandoraInstance,
        const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CartesianVector &twoDPosition, const pandora::HitType hitType,
        const bool linearlyProjectEnds, pandora::CartesianVector &threeDPosition, float &projectionError);

    /**
     *  @brief  Get the fitted track direction at a given 2D position
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  trackFit the track fit object
     *  @param  twoDPosition the 2D position
     *  @param  hitType the hit type
     *  @param  linearlyProjectEnds whether to linearly project the ends of the fit
     *  @param  direction the direction (to populate)
     *  @param  projectionError the projection error (to populate)
     *
     *  @return the status code
     */
    static pandora::StatusCode GetFittedDirectionAtTwoDPosition(const pandora::Pandora &pandoraInstance,
        const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CartesianVector &twoDPosition, const pandora::HitType hitType,
        const bool linearlyProjectEnds, pandora::CartesianVector &direction, float &projectionError);

    /**
     *  @brief  Get the fitted track direction at a given 3D position
     *
     *  @param  trackFit the track fit object
     *  @param  threeDPosition the 3D position
     *  @param  snapToEnds whether to snap to the ends of the fit if the point is out of bounds
     *  @param  direction the direction (to populate)
     *
     *  @return the status code
     */
    static pandora::StatusCode GetFittedDirectionAtThreeDPosition(const lar_content::ThreeDSlidingFitResult &trackFit,
        const pandora::CartesianVector &threeDPosition, const bool snapToEnds, pandora::CartesianVector &direction);

private:
    static pandora::CartesianVector AlignVectorWithFit(
        const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CartesianVector &threeDPosition, const bool antiAlign);

    static pandora::StatusCode LinearlyExtrapolateVectorFromFitEnd(const lar_content::ThreeDSlidingFitResult &trackFit,
        const float longCoord, const bool maxLayer, pandora::CartesianVector &extrapolatedVector);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArAnalysisHelper::GetTrueKineticEnergy(const pandora::MCParticle *const pMCParticle)
{
    return pMCParticle->GetEnergy() - LArAnalysisHelper::GetTrueMass(pMCParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArAnalysisHelper::GetTrueMass(const pandora::MCParticle *const pMCParticle)
{
    return pandora::PdgTable::GetParticleMass(pMCParticle->GetParticleId());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArAnalysisHelper::GetFittedDirectionAtTwoDPosition(const pandora::Pandora &pandoraInstance,
    const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CartesianVector &twoDPosition, const pandora::HitType hitType,
    const bool linearlyProjectEnds, pandora::CartesianVector &fittedDirection, float &projectionError)
{
    pandora::CartesianVector threeDPosition(0.f, 0.f, 0.f);

    const pandora::StatusCode status = LArAnalysisHelper::ProjectTwoDPositionOntoTrackFit(
        pandoraInstance, trackFit, twoDPosition, hitType, linearlyProjectEnds, threeDPosition, projectionError);

    if (status == pandora::STATUS_CODE_SUCCESS)
        return LArAnalysisHelper::GetFittedDirectionAtThreeDPosition(trackFit, threeDPosition, linearlyProjectEnds, fittedDirection);

    return status;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_HELPER_H
