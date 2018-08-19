/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisHelper.cxx
 *
 *  @brief  Implementation of the lar analysis helper class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"

#include "Geometry/LArTPC.h"
#include "Managers/GeometryManager.h"

using namespace pandora;

namespace lar_physics_content
{

std::tuple<CartesianVector, CartesianVector> LArAnalysisHelper::GetFiducialCutCoordinates(
    const Pandora &pandoraInstance, const CartesianVector &fiducialCutLowMargins, const CartesianVector &fiducialCutHighMargins)
{
    const LArTPCMap &larTPCMap(pandoraInstance.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.size() != 1UL)
    {
        std::cout << "LArAnalysisHelper: the number of LArTPCs was not equal to 1" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

    const LArTPC *const pLArTPC(larTPCMap.begin()->second);

    // The extremal coordinates are at the centre +- half of the widths
    const float xMin(pLArTPC->GetCenterX() - (0.5f * pLArTPC->GetWidthX()) + fiducialCutLowMargins.GetX());
    const float yMin(pLArTPC->GetCenterY() - (0.5f * pLArTPC->GetWidthY()) + fiducialCutLowMargins.GetY());
    const float zMin(pLArTPC->GetCenterZ() - (0.5f * pLArTPC->GetWidthZ()) + fiducialCutLowMargins.GetZ());

    const float xMax(pLArTPC->GetCenterX() + (0.5f * pLArTPC->GetWidthX()) - fiducialCutHighMargins.GetX());
    const float yMax(pLArTPC->GetCenterY() + (0.5f * pLArTPC->GetWidthY()) - fiducialCutHighMargins.GetY());
    const float zMax(pLArTPC->GetCenterZ() + (0.5f * pLArTPC->GetWidthZ()) - fiducialCutHighMargins.GetZ());

    return {CartesianVector(xMin, yMin, zMin), CartesianVector(xMax, yMax, zMax)};
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisHelper::IsPointFiducial(const CartesianVector &point, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates)
{
    const float xPosition = point.GetX();
    const float yPosition = point.GetY();
    const float zPosition = point.GetZ();

    if (xPosition < minCoordinates.GetX() || xPosition > maxCoordinates.GetX())
        return false;

    if (yPosition < minCoordinates.GetY() || yPosition > maxCoordinates.GetY())
        return false;

    if (zPosition < minCoordinates.GetZ() || zPosition > maxCoordinates.GetZ())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisHelper::IsTrueShower(const MCParticle *const pMCParticle)
{
    switch (pMCParticle->GetParticleId())
    {
        case PHOTON:
        case E_MINUS:
        case E_PLUS:
            return true;

        default:
            break;
    }

    return false;
}

} // namespace lar_physics_content
