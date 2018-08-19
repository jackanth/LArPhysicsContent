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
     *  @brief  Get the fitted track direction at a given position
     *
     *  @param  trackFit the track fit object
     *  @param  position the position
     *
     *  @return the track fit direction
     */
    static pandora::CartesianVector GetFittedDirectionAtPosition(
        const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CartesianVector &position);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArAnalysisHelper::GetTrueKineticEnergy(const pandora::MCParticle *const pMCParticle)
{
    return pMCParticle->GetEnergy() - GetTrueMass(pMCParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArAnalysisHelper::GetTrueMass(const pandora::MCParticle *const pMCParticle)
{
    return pandora::PdgTable::GetParticleMass(pMCParticle->GetParticleId());
}

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_HELPER_H
