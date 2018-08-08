/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisNtupleHelper.cxx
 *
 *  @brief  Implementation of the lar analysis ntuple helper class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArHelpers/LArAnalysisNtupleHelper.h"

using namespace pandora;

namespace lar_physics_content
{

LArAnalysisNtupleHelper::PARTICLE_CLASS LArAnalysisNtupleHelper::GetParticleClass(const ParticleFlowObject *const pPfo)
{
    (void)pPfo;
    return PARTICLE_CLASS::OTHER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArAnalysisNtupleHelper::ToString(const PARTICLE_CLASS particleClass)
{
    switch (particleClass)
    {
        case PARTICLE_CLASS::NEUTRINO:
            return "NEUTRINO";

        case PARTICLE_CLASS::PRIMARY:
            return "PRIMARY";

        case PARTICLE_CLASS::COSMIC_RAY:
            return "COSMIC_RAY";

        case PARTICLE_CLASS::OTHER:
            return "OTHER";

        default:
            break;
    }

    std::cout << "LArAnalysisNtupleHelper: Unknown particle class" << std::endl;
    throw STATUS_CODE_INVALID_PARAMETER;
}

} // namespace lar_physics_content
