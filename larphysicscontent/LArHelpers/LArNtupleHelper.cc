/**
 *  @file   larphysicscontent/LArHelpers/LArNtupleHelper.cxx
 *
 *  @brief  Implementation of the lar ntuple helper class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

LArNtupleHelper::PARTICLE_CLASS LArNtupleHelper::GetParticleClass(const ParticleFlowObject *const pPfo)
{
    if (pPfo->GetParentPfoList().empty() && !LArPfoHelper::IsNeutrino(pPfo))
        return PARTICLE_CLASS::COSMIC_RAY;

    if (LArPfoHelper::IsNeutrinoFinalState(pPfo))
        return PARTICLE_CLASS::PRIMARY;

    if (LArPfoHelper::IsNeutrino(pPfo) && pPfo->GetParentPfoList().empty())
        return PARTICLE_CLASS::NEUTRINO;

    return PARTICLE_CLASS::OTHER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArNtupleHelper::ToString(const PARTICLE_CLASS particleClass)
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

    std::cout << "LArNtupleHelper: Unknown particle class" << std::endl;
    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_physics_content
