/**
 *  @file   larphysicscontent/TrackHitEnergyTool.cc
 *
 *  @brief  Implementation of the track hit energy tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/TrackHitEnergyTool.h"

#include "larphysicscontent/LArFittedTrackInfo.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_physics_content
{

TrackHitEnergyTool::TrackHitEnergyTool() :
    m_trackSlidingFitWindow(25U)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackHitEnergyTool::Run(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pPfo,
    LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap, float &excessCharge,
    const HitPurityToolCallback &hitPurityToolCallback)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    this->RecursivelyAppendMap(pPfo, excessCharge, hitPurityToolCallback, fittedTrackInfoMap);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackHitEnergyTool::RecursivelyAppendMap(const ParticleFlowObject *const pPfo, float &excessCharge,
    HitPurityToolCallback hitPurityToolCallback, LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap) const
{
    if (LArPfoHelper::IsTrack(pPfo))
    {
        try
        {
            const ThreeDSlidingFitResult fitResult = this->PerformSlidingTrackFit(pPfo);
            const LArFittedTrackInfo::TrackHitValueVector hitChargeVector = this->AppendLArTrackHitEnergyMap(pPfo, fitResult, excessCharge, hitPurityToolCallback);
            const float range = this->GetParticleRange(pPfo, fitResult);

            fittedTrackInfoMap.emplace(pPfo, LArFittedTrackInfo(pPfo, hitChargeVector, fitResult, range));
        }

        catch (...)
        {
            std::cout << "TrackHitEnergyTool: failed to perform track fit" << std::endl;
        }
    }

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        this->RecursivelyAppendMap(pDaughterPfo, excessCharge, hitPurityToolCallback, fittedTrackInfoMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDSlidingFitResult TrackHitEnergyTool::PerformSlidingTrackFit(const ParticleFlowObject *const pPfo) const
{
    // Get the 3D cluster and make sure there's at least one.
    ClusterList threeDClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, threeDClusterList);

    if (threeDClusterList.empty())
    {
        std::cout << "TrackHitEnergyTool: could not get track length because there were no 3D clusters" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

    // Get the complete coordinate vector and make sure it's not empty.
    CartesianPointVector coordinateVector;
    for (const Cluster *const pCluster : threeDClusterList)
    {
        CartesianPointVector clusterCoordinateVector;
        LArClusterHelper::GetCoordinateVector(pCluster, clusterCoordinateVector);
        coordinateVector.insert(coordinateVector.end(), std::make_move_iterator(clusterCoordinateVector.begin()),
                                std::make_move_iterator(clusterCoordinateVector.end()));
    }

    if (coordinateVector.empty())
    {
        std::cout << "TrackHitEnergyTool: Could not get track length because there were no points in the 3D clusters" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

    // Make a 3D sliding linear fit.
    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    return ThreeDSlidingFitResult{&coordinateVector, m_trackSlidingFitWindow, layerPitch};
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArFittedTrackInfo::TrackHitValueVector TrackHitEnergyTool::AppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo,
    const ThreeDSlidingFitResult &trackFit, float &excessCharge, HitPurityToolCallback hitPurityToolCallback) const
{
    // Get all hits and order them by projection along the track fit.
    const CaloHitList collectionPlaneHits = LArAnalysisParticleHelper::GetHitsOfType(pPfo, TPC_VIEW_W, false);

    const HitProjectionVector orderedHitProjections = this->OrderHitsByProjectionOnToTrackFit(collectionPlaneHits, trackFit);
    LArFittedTrackInfo::TrackHitValueVector trackHitChargeVector;

    for (const HitProjectionPair &projectionPair : orderedHitProjections)
    {
        const CaloHit *const pCaloHit = projectionPair.first;
        const float coordinate        = projectionPair.second;

        const float threeDDistance    = this->CaloHitToThreeDDistance(pCaloHit, trackFit);
        trackHitChargeVector.emplace_back(pCaloHit, coordinate, threeDDistance, pCaloHit->GetInputEnergy());
    }

    hitPurityToolCallback(trackHitChargeVector, excessCharge);
    return trackHitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackHitEnergyTool::CaloHitToThreeDDistance(const CaloHit *const pCaloHit, const ThreeDSlidingFitResult &trackFit) const
{
    const CartesianVector fitDirection = LArAnalysisParticleHelper::GetFittedDirectionAtPosition(trackFit, pCaloHit->GetPositionVector(), true);
    const float wirePitch = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W);
    return this->CellToThreeDDistance(pCaloHit->GetCellSize1(), wirePitch, fitDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackHitEnergyTool::CellToThreeDDistance(const float hitWidth, const float wirePitch, const CartesianVector &fitDirection) const
{
    float polarAngle = 0.f, azimuthalAngle = 0.f;
    std::tie(polarAngle, azimuthalAngle) = this->GetPolarAnglesFromDirection(fitDirection);

    const float cosPhi_sinTheta = std::fabs(std::cos(azimuthalAngle) * std::sin(polarAngle));
    const float sinPhi_sinTheta = std::fabs(std::sin(azimuthalAngle) * std::sin(polarAngle));

    float dx_p = std::numeric_limits<float>::max();
    float dx_w = std::numeric_limits<float>::max();

    if (cosPhi_sinTheta > std::numeric_limits<float>::epsilon())
        dx_p = wirePitch / cosPhi_sinTheta;

    if (sinPhi_sinTheta > std::numeric_limits<float>::epsilon())
        dx_w = hitWidth / sinPhi_sinTheta;

    return std::min(dx_p, dx_w);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<float, float> TrackHitEnergyTool::GetPolarAnglesFromDirection(const CartesianVector &direction) const
{
    const float polarAngle     = std::acos(std::fabs(direction.GetY()));
    const float azimuthalAngle = std::asin(std::fabs(direction.GetX() / std::sin(polarAngle)));

    return {polarAngle, azimuthalAngle};
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackHitEnergyTool::GetParticleRange(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit) const
{
    // Get the 3D CaloHits and order them by projection along the track fit.
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList);

    const HitProjectionVector orderedHitProjectionVector = this->OrderHitsByProjectionOnToTrackFit(caloHitList, trackFit);

    // Use the projections to create vector positions along the track.
    CartesianPointVector pointVector;

    for (const auto &projectionPair : orderedHitProjectionVector)
    {
        CartesianVector position{0.f, 0.f, 0.f};

        if (trackFit.GetGlobalFitPosition(projectionPair.second, position) != STATUS_CODE_SUCCESS)
        {
            std::cout << "LArAnalysisParticleHelper: could not add CaloHit position because it was not within the fit" << std::endl;
            continue;
        }

        pointVector.push_back(position);
    }

    if (pointVector.empty())
    {
        std::cout << "LArAnalysisParticleHelper: could not get track length because point vector was empty" << std::endl;
        return 0.f;
    }

    // Add up all the inter-hit distances to get a track length.
    CartesianVector currentPosition(pointVector.front());
    float totalLength(0.f);

    for (const CartesianVector &position : pointVector)
    {
        const float increment = (position - currentPosition).GetMagnitude();

        if (increment >= 0.f && increment < 1000.f)
            totalLength += increment;

        else
            std::cout << "LArAnalysisParticleHelper: Range increment was " << increment << std::endl;

        currentPosition = position;
    }

    return totalLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackHitEnergyTool::HitProjectionVector TrackHitEnergyTool::OrderHitsByProjectionOnToTrackFit(const CaloHitList &caloHitList,
    const ThreeDSlidingFitResult &trackFit) const
{
    HitProjectionVector orderedHitProjectionVector;

    for (const CaloHit *const pCaloHit : caloHitList)
        orderedHitProjectionVector.emplace_back(pCaloHit, trackFit.GetLongitudinalDisplacement(pCaloHit->GetPositionVector()));

    std::sort(orderedHitProjectionVector.begin(), orderedHitProjectionVector.end(),
        [](const HitProjectionPair &lhs, const HitProjectionPair &rhs)
        {
            return lhs.second < rhs.second;
        });

    return orderedHitProjectionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackHitEnergyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackSlidingFitWindow", m_trackSlidingFitWindow));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
