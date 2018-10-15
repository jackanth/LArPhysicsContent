/**
 *  @file   larphysicscontent/LArHelpers/LArNtupleHelper.h
 *
 *  @brief  Header file for the lar ntuple helper class.
 *
 *  $Log: $
 */
#ifndef LAR_NTUPLE_HELPER_H
#define LAR_NTUPLE_HELPER_H 1

#include "Objects/ParticleFlowObject.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <memory>

namespace lar_physics_content
{
/**
 *  @brief  LArNtupleHelper class
 */
class LArNtupleHelper
{
public:
    using TrackFitSharedPtr = std::shared_ptr<lar_content::ThreeDSlidingFitResult>; ///< Alias for a shared pointer to a track fit object

    /**
     *  @brief  The particle class
     */
    enum class PARTICLE_CLASS : unsigned
    {
        NEUTRINO   = 0U, ///< A beam neutrino
        PRIMARY    = 1U, ///< A primary daughter of a beam neutrino
        COSMIC_RAY = 2U, ///< A primary cosmic ray
        OTHER      = 3U  ///< Everything else
    };

    /**
     *  @brief  The vector branch type
     */
    enum class VECTOR_BRANCH_TYPE : unsigned
    {
        NEUTRINO   = 0U, ///< A neutrino-type vector
        PRIMARY    = 1U, ///< A primary-type vector
        COSMIC_RAY = 2U, ///< A cosmic-type vector
        PARTICLE   = 3U  ///< A particle-type vector
    };

    /**
     *  @brief  Deleted copy constructor
     */
    LArNtupleHelper(const LArNtupleHelper &) = delete;

    /**
     *  @brief  Deleted move constructor
     */
    LArNtupleHelper(LArNtupleHelper &&) = delete;

    /**
     *  @brief  Deleted copy assignment operator
     */
    LArNtupleHelper &operator=(const LArNtupleHelper &) = delete;

    /**
     *  @brief  Deleted move assignment operator
     */
    LArNtupleHelper &operator=(LArNtupleHelper &&) = delete;

    /**
     *  @brief  Deleted destructor
     */
    ~LArNtupleHelper() = delete;

    /**
     *  @brief  Get the class of a PFO
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the PFO class
     */
    static PARTICLE_CLASS GetParticleClass(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Translate a particle class to a string
     *
     *  @param  particleClass the particle class
     *
     *  @return the class as a string
     */
    static std::string ToString(const PARTICLE_CLASS particleClass);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_HELPER_H
