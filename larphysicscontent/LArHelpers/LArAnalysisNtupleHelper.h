/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisNtupleHelper.h
 *
 *  @brief  Header file for the lar analysis ntuple helper class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_NTUPLE_HELPER_H
#define LAR_ANALYSIS_NTUPLE_HELPER_H 1

#include "Objects/ParticleFlowObject.h"

namespace lar_physics_content
{
/**
 *  @brief  LArAnalysisNtupleHelper class
 */
class LArAnalysisNtupleHelper
{
public:
    /**
     *  @brief  The particle class
     */
    enum class PARTICLE_CLASS : unsigned
    {
        NEUTRINO   = 0, ///< A neutrino
        PRIMARY    = 1, ///< A primary daughter of the neutrino
        COSMIC_RAY = 2, ///< A cosmic ray
        OTHER      = 3  ///< Another class
    };

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

#endif // #ifndef LAR_ANALYSIS_NTUPLE_HELPER_H
