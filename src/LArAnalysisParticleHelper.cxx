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
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "Pandora/PdgTable.h"
#include "Managers/GeometryManager.h"
#include "Geometry/LArTPC.h"
#include "Api/PandoraContentApi.h"
#include "PandoraMonitoringApi.h"
#include "Helpers/MCParticleHelper.h"

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
        if (!LArAnalysisParticleHelper::IsPointFiducial(pCaloHit->GetPositionVector(), minCoordinates, maxCoordinates))
            return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisParticleHelper::IsPointFiducial(const CartesianVector &point, const CartesianVector &minCoordinates,
    const CartesianVector &maxCoordinates)
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

const MCParticle *LArAnalysisParticleHelper::GetMainMCParticle(const ParticleFlowObject *const pPfo)
{
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
    
    McParticleVotingMap mcParticleVotingMap;

    for (const Cluster *const pCluster : clusterList)
    {
        try
        {
            if (const MCParticle *const pThisMainMCParticle = MCParticleHelper::GetMainMCParticle(pCluster))
            {
                const MCParticle *const pThisMCPrimary = LArMCParticleHelper::GetPrimaryMCParticle(pThisMainMCParticle);
                const auto findIter = mcParticleVotingMap.find(pThisMCPrimary);
                
                if (findIter == mcParticleVotingMap.end())
                    mcParticleVotingMap.emplace(pThisMCPrimary, 1U);
                    
                else
                    ++findIter->second;
            }
        }
        
        catch (...)
        {
            continue;
        }
    }
    
    if (mcParticleVotingMap.empty())
        return nullptr;
        
    if (mcParticleVotingMap.size() == 1)
        return mcParticleVotingMap.begin()->first;
    
    // There are at least two different candidates, so pick the one with the most votes.
    McParticleVotingVector mcParticleVotingVector;
    mcParticleVotingVector.insert(mcParticleVotingVector.end(), std::make_move_iterator(mcParticleVotingMap.begin()), 
                                  std::make_move_iterator(mcParticleVotingMap.end()));
    
    std::sort(mcParticleVotingVector.begin(), mcParticleVotingVector.end(),
        [](const McParticleVotingPair &lhs, const McParticleVotingPair &rhs)
        {
            return lhs.second > rhs.second;
        }
    );
    
    return mcParticleVotingVector.front().first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisParticleHelper::IsNeutrino(const ParticleFlowObject *const pPfo)
{
    return (LArPfoHelper::IsNeutrino(pPfo) && pPfo->GetParentPfoList().empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisParticleHelper::IsCosmicRay(const ParticleFlowObject *const pPfo)
{
    return (pPfo->GetParentPfoList().empty() && !LArPfoHelper::IsNeutrino(pPfo) && LArPfoHelper::IsTrack(pPfo));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisParticleHelper::IsPrimaryNeutrinoDaughter(const ParticleFlowObject *const pPfo)
{
    return LArPfoHelper::IsNeutrinoFinalState(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArAnalysisParticleHelper::GetMcInformation(const MCParticle *const pMCParticle, float &mcEnergy, 
    LArAnalysisParticle::TypeTree &typeTree, LArAnalysisParticle::TYPE &mcType, CartesianVector &mcVertexPosition,
    CartesianVector &mcMomentum, int &mcPdgCode, float &mcContainmentFraction, const CartesianVector &minCoordinates,
    const CartesianVector &maxCoordinates)
{
    if (!pMCParticle)
        return false;
    
    typeTree         = LArAnalysisParticleHelper::CreateMcTypeTree(pMCParticle);
    mcEnergy         = pMCParticle->GetEnergy() - PdgTable::GetParticleMass(pMCParticle->GetParticleId());
    mcType           = typeTree.Type();
    mcVertexPosition = pMCParticle->GetVertex();
    mcMomentum       = pMCParticle->GetMomentum();
    mcPdgCode        = pMCParticle->GetParticleId();
    
    float escapedEnergy(0.f), totalEnergy(0.f);
    LArAnalysisParticleHelper::RecursivelyAddEscapedEnergy(pMCParticle, escapedEnergy, totalEnergy, minCoordinates, maxCoordinates);
    mcContainmentFraction = (totalEnergy - escapedEnergy) / totalEnergy;
    
    return true;
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

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticleHelper::RecursivelyAddEscapedEnergy(const MCParticle *const pCurrentMCParticle, float &escapedEnergy,
    float &totalEnergy, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates)
{
    bool allEnergyContained(true);
    
    switch (pCurrentMCParticle->GetParticleId())
    {
        case E_MINUS:
        case E_PLUS:
        case MU_MINUS:
        case MU_PLUS:
        case PI_PLUS:
        case PI_MINUS:
        case PROTON:
        {
            const float mcParticleEnergy = pCurrentMCParticle->GetEnergy() - PdgTable::GetParticleMass(pCurrentMCParticle->GetParticleId());
            const CartesianVector displacementVector = pCurrentMCParticle->GetEndpoint() - pCurrentMCParticle->GetVertex();
            
            if (displacementVector.GetMagnitude() < std::numeric_limits<float>::epsilon())
                break;
            
            float muMin(0.f), muMax(1.f);
            bool forceZeroContainment(false);
            
            LArAnalysisParticleHelper::AdjustMusForContainmentFraction(CartesianVector(minCoordinates.GetX(), 0.f, 0.f), 
                CartesianVector(-1.f, 0.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);
                
            LArAnalysisParticleHelper::AdjustMusForContainmentFraction(CartesianVector(maxCoordinates.GetX(), 0.f, 0.f), 
                CartesianVector(1.f, 0.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);
                
            LArAnalysisParticleHelper::AdjustMusForContainmentFraction(CartesianVector(0.f, minCoordinates.GetY(), 0.f), 
                CartesianVector(0.f, -1.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);
                
            LArAnalysisParticleHelper::AdjustMusForContainmentFraction(CartesianVector(0.f, maxCoordinates.GetY(), 0.f),
                CartesianVector(0.f, 1.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);
                
            LArAnalysisParticleHelper::AdjustMusForContainmentFraction(CartesianVector(0.f, 0.f, minCoordinates.GetZ()),
                CartesianVector(0.f, 0.f, -1.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);
                
            LArAnalysisParticleHelper::AdjustMusForContainmentFraction(CartesianVector(0.f, 0.f, maxCoordinates.GetZ()),
                CartesianVector(0.f, 0.f, 1.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);
                
            if (forceZeroContainment)
            {
                escapedEnergy += mcParticleEnergy;
                allEnergyContained = false;
            }
                
            else
            {
                const float containmentFraction = std::max(0.f, muMax - muMin);
                
                if (containmentFraction < 1.f)
                {
                    escapedEnergy += (1.f - containmentFraction) * mcParticleEnergy;
                    allEnergyContained = false;
                }
            }
                
            totalEnergy += mcParticleEnergy;
        }
        
        default: break;
    }
    
    if (allEnergyContained)
    {
        for (const MCParticle *const pDaughterParticle : pCurrentMCParticle->GetDaughterList())
        {
            LArAnalysisParticleHelper::RecursivelyAddEscapedEnergy(pDaughterParticle, escapedEnergy, totalEnergy, minCoordinates,
                maxCoordinates);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticleHelper::AdjustMusForContainmentFraction(const CartesianVector &planePoint, const CartesianVector &planeNormal, 
    const CartesianVector &vertexPosition, const CartesianVector &originalDisplacementVector, float &muMin, float &muMax,
    bool &forceZeroContainment)
{
    const float projectedDisplacement = originalDisplacementVector.GetDotProduct(planeNormal);
    
    if (std::fabs(projectedDisplacement) < std::numeric_limits<float>::epsilon())
        return;
        
    const float muIntercept = (planePoint - vertexPosition).GetDotProduct(planeNormal) / projectedDisplacement;
    const bool isAligned = (projectedDisplacement > 0.f);
    
    if (isAligned)
    {
        if (muIntercept > std::numeric_limits<float>::min() && muIntercept < std::numeric_limits<float>::max())
            muMax = std::min(muMax, muIntercept);
            
        else if (muIntercept < 0.f)
            forceZeroContainment = true;
    }
    
    else
    {
        if (muIntercept > std::numeric_limits<float>::min() && muIntercept < std::numeric_limits<float>::max())
            muMin = std::max(muMin, muIntercept);
            
        else if (muIntercept > 0.f)
            forceZeroContainment = true;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
LArAnalysisParticle::TypeTree LArAnalysisParticleHelper::CreateMcTypeTree(const MCParticle *const pMCParticle)
{
    LArAnalysisParticle::TypeTree::List daughterTypeTrees;
    
    const LArAnalysisParticle::TYPE type = LArAnalysisParticleHelper::GetMcParticleType(pMCParticle);
    
    if (type != LArAnalysisParticle::TYPE::SHOWER && type != LArAnalysisParticle::TYPE::UNKNOWN)
    {
        for (const MCParticle *const pDaughterParticle : pMCParticle->GetDaughterList())
        {
            if (pDaughterParticle->GetEnergy() > 0.05f)
                daughterTypeTrees.push_back(LArAnalysisParticleHelper::CreateMcTypeTree(pDaughterParticle));
        }
    }
    
    return {type, daughterTypeTrees};
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticle::TYPE LArAnalysisParticleHelper::GetMcParticleType(const MCParticle *const pMCParticle)
{
    switch (pMCParticle->GetParticleId())
    {
        case PROTON:   return LArAnalysisParticle::TYPE::PROTON;
        case MU_MINUS: 
        case PI_MINUS:
        case PI_PLUS:  return LArAnalysisParticle::TYPE::PION_MUON;
        case PHOTON:
        case E_MINUS:
        case E_PLUS:   return LArAnalysisParticle::TYPE::SHOWER;
        case NEUTRON:  return LArAnalysisParticle::TYPE::UNKNOWN;
        default: break;
    }
    
    return LArAnalysisParticle::TYPE::UNKNOWN;
}

} // namespace lar_physics_content
