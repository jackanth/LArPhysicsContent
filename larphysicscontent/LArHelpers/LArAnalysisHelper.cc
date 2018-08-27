/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisHelper.cxx
 *
 *  @brief  Implementation of the lar analysis helper class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "Geometry/LArTPC.h"
#include "Managers/GeometryManager.h"

using namespace pandora;
using namespace lar_content;

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

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArAnalysisHelper::ProjectTwoDPositionOntoTrackFit(const Pandora &pandoraInstance, const ThreeDSlidingFitResult &trackFit,
    const CartesianVector &twoDPosition, const HitType hitType, const bool linearlyProjectEnds, CartesianVector &threeDPosition, float &projectionError)
{
    // Calculate the effective longitudinal coordinate along the 3D fit
    const float longCoordScalingFactor = LArGeometryHelper::ProjectDirection(pandoraInstance, trackFit.GetAxisDirection(), hitType).GetDotProduct(trackFit.GetAxisDirection());

    if (longCoordScalingFactor <= std::numeric_limits<float>::epsilon())
        return STATUS_CODE_FAILURE;

    const CartesianVector &projectedIntercept = LArGeometryHelper::ProjectPosition(pandoraInstance, trackFit.GetAxisIntercept(), hitType);
    const CartesianVector &projectedDirection = LArGeometryHelper::ProjectDirection(pandoraInstance, trackFit.GetAxisDirection(), hitType);

    const float projectedLongCoord = (twoDPosition - projectedIntercept).GetDotProduct(projectedDirection);
    const float longCoord          = projectedLongCoord / longCoordScalingFactor;

    // If this is within bounds, then get the fit position
    if (STATUS_CODE_SUCCESS == trackFit.GetGlobalFitPosition(longCoord, threeDPosition))
    {
        projectionError = (LArGeometryHelper::ProjectPosition(pandoraInstance, threeDPosition, hitType) - twoDPosition).GetMagnitude();
        return STATUS_CODE_SUCCESS;
    }

    // If out of bounds and we're not willing to try to project the ends, this is a failure
    if (!linearlyProjectEnds)
        return STATUS_CODE_FAILURE;

    const CartesianVector &minPosition = trackFit.GetGlobalMinLayerPosition();
    const CartesianVector &maxPosition = trackFit.GetGlobalMaxLayerPosition();

    const CartesianVector projectedMinPosition = LArGeometryHelper::ProjectPosition(pandoraInstance, minPosition, hitType);
    const CartesianVector projectedMaxPosition = LArGeometryHelper::ProjectPosition(pandoraInstance, maxPosition, hitType);

    const bool closerToMax = ((twoDPosition - projectedMinPosition).GetMagnitude() > (twoDPosition - projectedMaxPosition).GetMagnitude());

    // Check the quality
    if (STATUS_CODE_SUCCESS == LArAnalysisHelper::LinearlyExtrapolateVectorFromFitEnd(trackFit, longCoord, closerToMax, threeDPosition))
    {
        projectionError = (LArGeometryHelper::ProjectPosition(pandoraInstance, threeDPosition, hitType) - twoDPosition).GetMagnitude();
        return STATUS_CODE_SUCCESS;
    }

    return STATUS_CODE_FAILURE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArAnalysisHelper::GetFittedDirectionAtThreeDPosition(
    const ThreeDSlidingFitResult &trackFit, const CartesianVector &threeDPosition, const bool snapToEnds, CartesianVector &direction)
{
    // Get the extremal fit parameters
    const CartesianVector &minPosition = trackFit.GetGlobalMinLayerPosition();
    const CartesianVector &maxPosition = trackFit.GetGlobalMaxLayerPosition();

    const float longCoord   = trackFit.GetLongitudinalDisplacement(threeDPosition);
    const bool  closerToMin = ((threeDPosition - minPosition).GetMagnitude() < (threeDPosition - maxPosition).GetMagnitude());

    if (trackFit.GetGlobalFitDirection(longCoord, direction) != STATUS_CODE_SUCCESS)
    {
        if (!snapToEnds)
        {
            std::cerr << "LArAnalysisHelper: Could not get fitted direction because position was not within fit" << std::endl;
            return STATUS_CODE_FAILURE;
        }

        // Get the fit direction at the fit position closest to the point
        // ATTN snapping to ends isn't ideal for multivalued end values: would rather find the point on the fit closest to the desired point
        direction = closerToMin ? trackFit.GetGlobalMinLayerDirection() : trackFit.GetGlobalMaxLayerDirection();
    }

    // Make sure the fit direction points more towards the fit centre
    direction = LArAnalysisHelper::AlignVectorWithFit(trackFit, direction, closerToMin);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArAnalysisHelper::AlignVectorWithFit(const ThreeDSlidingFitResult &trackFit, const CartesianVector &threeDPosition, const bool antiAlign)
{
    // Get the extremal fit parameters
    const CartesianVector &minPosition  = trackFit.GetGlobalMinLayerPosition();
    const CartesianVector &maxPosition  = trackFit.GetGlobalMaxLayerPosition();
    const CartesianVector  fitDirection = (maxPosition - minPosition).GetUnitVector();

    if (antiAlign)
    {
        if (threeDPosition.GetDotProduct(fitDirection) < 0.f)
            return threeDPosition;

        return threeDPosition * -1.f;
    }

    if (threeDPosition.GetDotProduct(fitDirection) > 0.f)
        return threeDPosition;

    return threeDPosition * -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArAnalysisHelper::LinearlyExtrapolateVectorFromFitEnd(
    const ThreeDSlidingFitResult &trackFit, const float longCoord, const bool maxLayer, CartesianVector &extrapolatedVector)
{
    const CartesianVector &endPosition      = maxLayer ? trackFit.GetGlobalMaxLayerPosition() : trackFit.GetGlobalMinLayerPosition();
    const CartesianVector &endDirection     = maxLayer ? trackFit.GetGlobalMaxLayerDirection() : trackFit.GetGlobalMinLayerDirection();
    const float            coordinateAtEnd  = trackFit.GetLongitudinalDisplacement(endDirection);
    const float            excessCoordinate = maxLayer ? longCoord - coordinateAtEnd : coordinateAtEnd - longCoord;

    if ((excessCoordinate < 0.f) || (endDirection.GetDotProduct(trackFit.GetAxisDirection()) < 0.f))
        return STATUS_CODE_FAILURE;

    // Check if the fit end direction points the wrong way
    if (maxLayer && (endDirection.GetDotProduct(trackFit.GetAxisDirection()) < 0.f))
        return STATUS_CODE_FAILURE;

    if (!maxLayer && (endDirection.GetDotProduct(trackFit.GetAxisDirection()) > 0.f))
        return STATUS_CODE_FAILURE;

    const float excessCoordinateScalingFactor = std::fabs(endDirection.GetDotProduct(trackFit.GetAxisDirection()));

    if (excessCoordinateScalingFactor <= std::numeric_limits<float>::epsilon())
        return STATUS_CODE_FAILURE;

    extrapolatedVector = endPosition + endDirection * (excessCoordinate / excessCoordinateScalingFactor);
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
