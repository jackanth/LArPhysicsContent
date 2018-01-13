/**
 *  @file   LArPhysicsContent/src/LArAnalysisParticleHelper.cxx
 *
 *  @brief  Implementation of the LEE analysis helper class.
 *
 *  $Log: $
 */

#include "LArAnalysisParticleHelper.h"

#include "DebugDefinitions.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Managers/GeometryManager.h"
#include "Geometry/LArTPC.h"
#include "Api/PandoraContentApi.h"
#include "PandoraMonitoringApi.h"

#include "TFile.h"
#include "TCanvas.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
void LArAnalysisParticleHelper::GetFiducialCutParameters(const Pandora &pandoraInstance, const float fiducialCutXMargin, const float fiducialCutYMargin,
    const float fiducialCutZMargin, CartesianVector &minCoordinates, CartesianVector &maxCoordinates)
{
    const LArTPCMap &larTPCMap(pandoraInstance.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.size() != 1)
    {
        CERR("The number of LArTPCs was not equal to 1 (" << larTPCMap.size() << ")");
        throw STATUS_CODE_NOT_FOUND;
    }

    const LArTPC *const pLArTPC(larTPCMap.begin()->second);

    const float xMin = pLArTPC->GetCenterX() - (0.5f * pLArTPC->GetWidthX()) + fiducialCutXMargin;
    const float yMin = pLArTPC->GetCenterY() - (0.5f * pLArTPC->GetWidthY()) + fiducialCutYMargin;
    const float zMin = pLArTPC->GetCenterZ() - (0.5f * pLArTPC->GetWidthZ()) + fiducialCutZMargin;

    const float xMax = pLArTPC->GetCenterX() + (0.5f * pLArTPC->GetWidthX()) - fiducialCutXMargin;
    const float yMax = pLArTPC->GetCenterY() + (0.5f * pLArTPC->GetWidthY()) - fiducialCutYMargin;
    const float zMax = pLArTPC->GetCenterZ() + (0.5f * pLArTPC->GetWidthZ()) - fiducialCutZMargin;
    
    minCoordinates = CartesianVector(xMin, yMin, zMin);
    maxCoordinates = CartesianVector(xMax, yMax, zMax);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticleHelper::RecursivelyAppendTrackFitMap(const Pandora &pandoraInstance, const ParticleFlowObject *const pPfo,
                                                     TrackFitMap &trackFitMap, const unsigned int slidingFitWindow)
{
    if (LArPfoHelper::IsTrack(pPfo))
    {
        try
        {
            trackFitMap.emplace(pPfo, LArAnalysisParticleHelper::PerformSlidingTrackFit(pandoraInstance, pPfo, slidingFitWindow));
        }
        
        catch (...)
        {
            CERR("Failed to perform track fit");
        }
    }
    
    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        LArAnalysisParticleHelper::RecursivelyAppendTrackFitMap(pandoraInstance, pDaughterPfo, trackFitMap, slidingFitWindow);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDSlidingFitResult LArAnalysisParticleHelper::PerformSlidingTrackFit(const Pandora &pandoraInstance, const ParticleFlowObject *const pPfo,
                                                                 const unsigned int slidingFitWindow)
{
    // Get the 3D cluster and make sure there's at least one.
    ClusterList threeDClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, threeDClusterList);

    if (threeDClusterList.empty())
    {
        CERR("Could not get track length because there were no 3D clusters");
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
        CERR("Could not get track length because there were no points in the 3D clusters");
        throw STATUS_CODE_NOT_FOUND;
    }

    // Make a 3D sliding linear fit.
    const float layerPitch(LArGeometryHelper::GetWireZPitch(pandoraInstance));
    return ThreeDSlidingFitResult{&coordinateVector, slidingFitWindow, layerPitch};
}

//------------------------------------------------------------------------------------------------------------------------------------------

// ATTN this method is copied from elsewhere.
PfoList LArAnalysisParticleHelper::GetRecoNeutrinoList(const Algorithm &algorithm, const std::string &pfoListName)
{
    const PfoList *pPfoList = nullptr;
    PandoraContentApi::GetList(algorithm, pfoListName, pPfoList);

    // Obtain all reco particles.
    PfoList allRecoParticleList(pPfoList ? *pPfoList : PfoList());

    PfoList allRecoNeutrinoList;
    LArPfoHelper::GetRecoNeutrinos(&allRecoParticleList, allRecoNeutrinoList);
    allRecoNeutrinoList.sort(LArAnalysisParticleHelper::SortRecoNeutrinos);
    
    return allRecoNeutrinoList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArAnalysisParticleHelper::CaloHitToThreeDDistance(const Pandora &pandoraInstance, const CaloHit *const pCaloHit,
                                                 const ThreeDSlidingFitResult &trackFit)
{
    const CartesianVector fitDirection = LArAnalysisParticleHelper::GetFittedDirectionAtPosition(trackFit, pCaloHit->GetPositionVector());
    return LArAnalysisParticleHelper::CellToThreeDDistance(pCaloHit->GetCellSize1(), LArGeometryHelper::GetWirePitch(pandoraInstance, TPC_VIEW_W),
                                                   fitDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float> LArAnalysisParticleHelper::GetPolarAnglesFromDirection(const CartesianVector &direction)
{
    const float polarAngle     = std::acos(std::fabs(direction.GetY()));
    const float azimuthalAngle = std::asin(std::fabs(direction.GetX() / std::sin(polarAngle)));
    
    return std::make_tuple(polarAngle, azimuthalAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArAnalysisParticleHelper::GetHitsOfType(const ParticleFlowObject *const pPfo, const HitType hitType, const bool recurseOverDaughters)
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, hitType, caloHitList);
    
    if (recurseOverDaughters)
    {
        for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
        {
            CaloHitList daughterCaloHitList = GetHitsOfType(pDaughterPfo, hitType, recurseOverDaughters);
            caloHitList.insert(caloHitList.end(), std::make_move_iterator(daughterCaloHitList.begin()),
                               std::make_move_iterator(daughterCaloHitList.end()));
        }
    }
    
    return caloHitList;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArAnalysisParticleHelper::GetFittedDirectionAtPosition(const ThreeDSlidingFitResult &trackFit, const CartesianVector &position)
{
    const CartesianVector &minPosition = trackFit.GetGlobalMinLayerPosition();
    const CartesianVector &maxPosition = trackFit.GetGlobalMaxLayerPosition();
    
    const float minCoordinate = trackFit.GetLongitudinalDisplacement(minPosition);
    const float maxCoordinate = trackFit.GetLongitudinalDisplacement(maxPosition);
    
    const CartesianVector &minDirection = trackFit.GetGlobalMinLayerDirection();
    const CartesianVector &maxDirection = trackFit.GetGlobalMaxLayerDirection();
 
    CartesianVector fitDirection{0.f, 0.f, 0.f};
    const float displacementAlongFittedTrack = trackFit.GetLongitudinalDisplacement(position);
        
    if (trackFit.GetGlobalFitDirection(displacementAlongFittedTrack, fitDirection) != STATUS_CODE_SUCCESS)
    {
        if (displacementAlongFittedTrack <= minCoordinate)
            fitDirection = minDirection;
            
        else if (displacementAlongFittedTrack >= maxCoordinate)
            fitDirection = maxDirection;
            
        else
        {
            const float distanceToMinPosition = (position - minPosition).GetMagnitude();
            const float distanceToMaxPosition = (position - maxPosition).GetMagnitude();
            
            if (distanceToMinPosition < distanceToMaxPosition)
                fitDirection = minDirection;
                
            else
                fitDirection = maxDirection;
        }
    }
    
    return fitDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisParticleHelper::RecursivelyCheckFiducialCut(const ParticleFlowObject *const pPfo, const CartesianVector &minCoordinates,
    const CartesianVector &maxCoordinates)
{
    if (!LArAnalysisParticleHelper::CheckFiducialCut(pPfo, minCoordinates, maxCoordinates))
        return false;
    
    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        if (!LArAnalysisParticleHelper::RecursivelyCheckFiducialCut(pDaughterPfo, minCoordinates, maxCoordinates))
            return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisParticleHelper::CheckFiducialCut(const ParticleFlowObject *const pPfo, const CartesianVector &minCoordinates,
    const CartesianVector &maxCoordinates)
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList);
    
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float hitXPosition = pCaloHit->GetPositionVector().GetX();
        const float hitYPosition = pCaloHit->GetPositionVector().GetY();
        const float hitZPosition = pCaloHit->GetPositionVector().GetZ();
        
        if (hitXPosition < minCoordinates.GetX() || hitXPosition > maxCoordinates.GetX())
            return false;
            
        if (hitYPosition < minCoordinates.GetY() || hitYPosition > maxCoordinates.GetY())
            return false;
            
        if (hitZPosition < minCoordinates.GetZ() || hitZPosition > maxCoordinates.GetZ())
            return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArAnalysisParticleHelper::GetParticleRange(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit)
{
    // Get the 3D CaloHits and order them by projection along the track fit.
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList);
    
    const HitProjectionVector orderedHitProjectionVector = LArAnalysisParticleHelper::OrderHitsByProjectionOnToTrackFit(caloHitList, trackFit);

    // Use the projections to create vector positions along the track.
    CartesianPointVector pointVector;
    
    for (const auto &projectionPair : orderedHitProjectionVector)
    {
        CartesianVector position{0.f, 0.f, 0.f};

        if (trackFit.GetGlobalFitPosition(projectionPair.second, position) != STATUS_CODE_SUCCESS)
        {
            CERR("Could not add CaloHit position because it was not within the fit");
            continue;
        }

        pointVector.push_back(position);
    }

    if (pointVector.empty())
    {
        CERR("Could not get track length because point vector was empty");
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
            CERR("Range increment was " << increment);
            
        currentPosition = position;
    }
    
    return totalLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticleHelper::HitProjectionVector LArAnalysisParticleHelper::OrderHitsByProjectionOnToTrackFit(const CaloHitList &caloHitList,
                                                                                            const ThreeDSlidingFitResult &trackFit)
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

void LArAnalysisParticleHelper::WriteNTuple(TNtuple *const pNtuple, const std::string &fileName, const bool verboseMode)
{
    if (verboseMode)
        pNtuple->Print();
    
    TFile *pFile = new TFile(fileName.c_str(), "NEW");
    pNtuple->Write();
    
    pFile->Close();
    delete pFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TNtuple * LArAnalysisParticleHelper::LoadNTupleFromFile(const std::string &filePath, const std::string &nTupleName)
{
    TFile *pFile = new TFile(filePath.c_str(), "READ");
    
    if (!pFile->IsOpen())
    {
        CERR("Failed to open file at " << filePath);
        return NULL;
    }
    
    if (!pFile->GetListOfKeys()->Contains(nTupleName.c_str()))
    {
        CERR("Data file at " << filePath << " did not contain key '" << nTupleName << "'");
        return NULL;
    }
    
    return (TNtuple *)pFile->Get(nTupleName.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

// ATTN this method is copied from elsewhere.
bool LArAnalysisParticleHelper::SortRecoNeutrinos(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    if (!LArPfoHelper::IsNeutrino(pLhs) || !LArPfoHelper::IsNeutrino(pRhs))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PfoList downstreamPfosLhs, downstreamPfosRhs;
    LArPfoHelper::GetAllDownstreamPfos(pLhs, downstreamPfosLhs);
    LArPfoHelper::GetAllDownstreamPfos(pRhs, downstreamPfosRhs);

    // If just left with the neutrino pfos themselves
    if ((1 == downstreamPfosLhs.size()) && (1 == downstreamPfosRhs.size()))
    {
        // ATTN Not a good pair of tie-breakers, but this should be rare (ideally shouldn't have any neutrinos without daughter pfos)
        if (!pLhs->GetVertexList().empty() && !pRhs->GetVertexList().empty())
            return ((*(pLhs->GetVertexList().begin()))->GetPosition().GetZ() < (*(pRhs->GetVertexList().begin()))->GetPosition().GetZ());

        return (pLhs->GetParticleId() < pRhs->GetParticleId());
    }

    PfoVector pfoVectorLhs(downstreamPfosLhs.begin(), downstreamPfosLhs.end());
    PfoVector pfoVectorRhs(downstreamPfosRhs.begin(), downstreamPfosRhs.end());
    std::sort(pfoVectorLhs.begin(), pfoVectorLhs.end(), LArPfoHelper::SortByNHits);
    std::sort(pfoVectorRhs.begin(), pfoVectorRhs.end(), LArPfoHelper::SortByNHits);

    if (pfoVectorLhs.empty() || pfoVectorRhs.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return (LArPfoHelper::SortByNHits(pfoVectorLhs.front(), pfoVectorRhs.front()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArAnalysisParticleHelper::CellToThreeDDistance(const float hitWidth, const float wirePitch, const CartesianVector &fitDirection)
{
    float polarAngle = 0.f, azimuthalAngle = 0.f;
    std::tie(polarAngle, azimuthalAngle) = LArAnalysisParticleHelper::GetPolarAnglesFromDirection(fitDirection);
    
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
} // namespace lar_physics_content
