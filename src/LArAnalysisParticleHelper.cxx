/**
 *  @file   LArPhysicsContent/src/LArAnalysisParticleHelper.cxx
 *
 *  @brief  Implementation of the LEE analysis helper class.
 *
 *  $Log: $
 */

#include "LArAnalysisParticleHelper.h"

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

std::string LArAnalysisParticleHelper::TypeAsString(const LArAnalysisParticle::TYPE type)
{
    switch (type)
    {
        case LArAnalysisParticle::TYPE::PION_MUON:  return "PION_MUON";
        case LArAnalysisParticle::TYPE::PROTON:     return "PROTON";
        case LArAnalysisParticle::TYPE::SHOWER:     return "SHOWER";
        case LArAnalysisParticle::TYPE::TRACK:      return "TRACK";
        case LArAnalysisParticle::TYPE::NEUTRINO:   return "NEUTRINO";
        case LArAnalysisParticle::TYPE::COSMIC_RAY: return "COSMIC_RAY";
        default:               break;
    }

    return "UNKNOWN";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticleHelper::GetFiducialCutParameters(const Pandora &pandoraInstance, const float fiducialCutLowXMargin,
    const float fiducialCutHighXMargin, const float fiducialCutLowYMargin, const float fiducialCutHighYMargin,
    const float fiducialCutLowZMargin, const float fiducialCutHighZMargin, CartesianVector &minCoordinates, CartesianVector &maxCoordinates)
{
    const LArTPCMap &larTPCMap(pandoraInstance.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.size() != 1)
    {
        std::cout << "LArAnalysisParticleHelper: the number of LArTPCs was not equal to 1" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

    const LArTPC *const pLArTPC(larTPCMap.begin()->second);

    const float xMin = pLArTPC->GetCenterX() - (0.5f * pLArTPC->GetWidthX()) + fiducialCutLowXMargin;
    const float yMin = pLArTPC->GetCenterY() - (0.5f * pLArTPC->GetWidthY()) + fiducialCutLowYMargin;
    const float zMin = pLArTPC->GetCenterZ() - (0.5f * pLArTPC->GetWidthZ()) + fiducialCutLowZMargin;

    const float xMax = pLArTPC->GetCenterX() + (0.5f * pLArTPC->GetWidthX()) - fiducialCutHighXMargin;
    const float yMax = pLArTPC->GetCenterY() + (0.5f * pLArTPC->GetWidthY()) - fiducialCutHighYMargin;
    const float zMax = pLArTPC->GetCenterZ() + (0.5f * pLArTPC->GetWidthZ()) - fiducialCutHighZMargin;

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
            std::cout << "LArAnalysisParticleHelper: failed to perform track fit" << std::endl;
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
        std::cout << "LArAnalysisParticleHelper: could not get track length because there were no 3D clusters" << std::endl;
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
        std::cout << "LArAnalysisParticleHelper: Could not get track length because there were no points in the 3D clusters" << std::endl;
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
    const CartesianVector fitDirection = LArAnalysisParticleHelper::GetFittedDirectionAtPosition(trackFit, pCaloHit->GetPositionVector(), true);
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

CartesianVector LArAnalysisParticleHelper::GetFittedDirectionAtPosition(const ThreeDSlidingFitResult &trackFit, const CartesianVector &position,
    const bool pointTowardsMiddle)
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

    if (pointTowardsMiddle)
    {
        if ((maxCoordinate - displacementAlongFittedTrack) < (displacementAlongFittedTrack - minCoordinate))
            fitDirection *= -1.f;
    }

    else
    {
        if (fitDirection.GetDotProduct(CartesianVector(0.f, -1.f, 0.f)) < 0.f)
            fitDirection *= -1.f;
    }

    return fitDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArAnalysisParticleHelper::GetFractionOfFiducialHits(const ParticleFlowObject *const pPfo, const CartesianVector &minCoordinates,
    const CartesianVector &maxCoordinates)
{
    PfoList downstreamPfos;
    LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfos);

    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(downstreamPfos, TPC_3D, caloHitList);

    unsigned fiducialHits(0U);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (LArAnalysisParticleHelper::IsPointFiducial(pCaloHit->GetPositionVector(), minCoordinates, maxCoordinates))
            ++fiducialHits;
    }

    return static_cast<float>(fiducialHits) / static_cast<float>(caloHitList.size());
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
        std::cout << "LArAnalysisParticleHelper: failed to open file at " << filePath << std::endl;
        return NULL;
    }

    if (!pFile->GetListOfKeys()->Contains(nTupleName.c_str()))
    {
        std::cout << "LArAnalysisParticleHelper: data file at " << filePath << " did not contain key '" << nTupleName << "'" << std::endl;
        return NULL;
    }

    return (TNtuple *)pFile->Get(nTupleName.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArAnalysisParticleHelper::GetMainMCParticle(const ParticleFlowObject *const pPfo)
{
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);

    CaloHitList caloHitList;

    for (const Cluster *const pCluster : clusterList)
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        MCParticleVector mcParticleVector;

        for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap)
            mcParticleVector.push_back(mapEntry.first);

        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            try
            {
                const MCParticle *const pMCPrimary = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
                mcParticleWeightMap[pMCPrimary] += hitMCParticleWeightMap.at(pMCParticle);
            }

            catch (...)
            {
                continue;
            }
        }
    }

    float bestWeight(0.f);
    const MCParticle *pBestMCParticle(nullptr);

    MCParticleVector mcParticleVector;

    for (const MCParticleWeightMap::value_type &mapEntry : mcParticleWeightMap)
        mcParticleVector.push_back(mapEntry.first);

    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

    for (const MCParticle *const pCurrentMCParticle : mcParticleVector)
    {
        const float currentWeight(mcParticleWeightMap.at(pCurrentMCParticle));

        if (currentWeight > bestWeight)
        {
            pBestMCParticle = pCurrentMCParticle;
            bestWeight = currentWeight;
        }
    }

    if (!pBestMCParticle)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pBestMCParticle;
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

bool LArAnalysisParticleHelper::GetMcInformation(const MCParticle *const pMCParticle, float &mcEnergy, float &mcKineticEnergy, float &mcMass,
    LArAnalysisParticle::TypeTree &typeTree, LArAnalysisParticle::TYPE &mcType, CartesianVector &mcVertexPosition,
    CartesianVector &mcMomentum, int &mcPdgCode, float &mcContainmentFraction, const CartesianVector &minCoordinates,
    const CartesianVector &maxCoordinates)
{
    if (!pMCParticle)
        return false;

    LArAnalysisParticleHelper::CreateMcTypeTree(pMCParticle, typeTree);

    mcEnergy         = pMCParticle->GetEnergy();
    mcKineticEnergy  = pMCParticle->GetEnergy() - PdgTable::GetParticleMass(pMCParticle->GetParticleId());
    mcMass           = PdgTable::GetParticleMass(pMCParticle->GetParticleId());
    mcType           = typeTree.Type();
    mcVertexPosition = pMCParticle->GetVertex();
    mcMomentum       = pMCParticle->GetMomentum();
    mcPdgCode        = pMCParticle->GetParticleId();

    float escapedEnergy(0.f), totalEnergy(0.f);
    LArAnalysisParticleHelper::RecursivelyAddEscapedEnergy(pMCParticle, escapedEnergy, totalEnergy, minCoordinates, maxCoordinates);

    if (totalEnergy > std::numeric_limits<float>::epsilon())
        mcContainmentFraction = (totalEnergy - escapedEnergy) / totalEnergy;

    else
        mcContainmentFraction = 0.f;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticleHelper::CalculateHitPurityAndCompleteness(const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle,
    const CaloHitList *const pCaloHitList, const bool isNeutrino, float &hitPurity, float &hitCompleteness, float &collectionPlaneHitPurity,
    float &collectionPlaneHitCompleteness)
{
    PfoList downstreamPfos;
    LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfos);

    CaloHitList pfoAssociatedCaloHits;
    LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_U, pfoAssociatedCaloHits);
    LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_V, pfoAssociatedCaloHits);

    CaloHitList pfoAssociatedWCaloHits;
    LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_W, pfoAssociatedWCaloHits);
    pfoAssociatedCaloHits.insert(pfoAssociatedCaloHits.end(), pfoAssociatedWCaloHits.begin(), pfoAssociatedWCaloHits.end());

    LArAnalysisParticleHelper::CalculateHitPurityAndCompleteness(pfoAssociatedCaloHits, pMCParticle, pCaloHitList, isNeutrino, hitPurity, hitCompleteness, false);
    LArAnalysisParticleHelper::CalculateHitPurityAndCompleteness(pfoAssociatedWCaloHits, pMCParticle, pCaloHitList, isNeutrino, collectionPlaneHitPurity,
        collectionPlaneHitCompleteness, true);
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

bool LArAnalysisParticleHelper::CreateMcTypeTree(const MCParticle *const pMCParticle, LArAnalysisParticle::TypeTree &typeTree)
{
    LArAnalysisParticle::TypeTree::List daughterTypeTrees;

    const LArAnalysisParticle::TYPE type = LArAnalysisParticleHelper::GetMcParticleType(pMCParticle);

    if (type == LArAnalysisParticle::TYPE::UNKNOWN)
        return false;

    if (type != LArAnalysisParticle::TYPE::SHOWER)
    {
        for (const MCParticle *const pDaughterParticle : pMCParticle->GetDaughterList())
        {
            if (pDaughterParticle->GetEnergy() > 0.05f)
            {
                LArAnalysisParticle::TypeTree daughterTypeTree;

                if (LArAnalysisParticleHelper::CreateMcTypeTree(pDaughterParticle, daughterTypeTree))
                    daughterTypeTrees.push_back(daughterTypeTree);
            }
        }
    }

    typeTree = LArAnalysisParticle::TypeTree(type, daughterTypeTrees);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticle::TYPE LArAnalysisParticleHelper::GetMcParticleType(const MCParticle *const pMCParticle)
{
    if (LArMCParticleHelper::IsNeutrino(pMCParticle))
        return LArAnalysisParticle::TYPE::NEUTRINO;

    if (LArMCParticleHelper::IsCosmicRay(pMCParticle))
        return LArAnalysisParticle::TYPE::COSMIC_RAY;

    switch (pMCParticle->GetParticleId())
    {
        case PROTON:   return LArAnalysisParticle::TYPE::PROTON;
        case MU_MINUS:
        case MU_PLUS:
        case PI_MINUS:
        case PI_PLUS:  return LArAnalysisParticle::TYPE::PION_MUON;
        case PHOTON:
        case E_MINUS:
        case E_PLUS:   return LArAnalysisParticle::TYPE::SHOWER;
        case NEUTRON:
        default: break;
    }

    return LArAnalysisParticle::TYPE::UNKNOWN;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArAnalysisParticleHelper::CalculateHitPurityAndCompleteness(const CaloHitList  &pfoAssociatedCaloHits, const MCParticle *const pMCParticle,
    const CaloHitList *const pCaloHitList, const bool isNeutrino, float &hitPurity, float &hitCompleteness,
    const bool useCollectionPlaneOnly)
{
    // Purity       = (num 2D hits assoc with PFO or its descendents \cap assoc with MC particle or its descendents) /
    //                (num 2D hits assoc with PFO or its descendents)

    // Completeness = (num 2D hits assoc with PFO or its descendents \cap assoc with MC particle or its descendents) /
    //                (num 2D hits assoc with MC particle or its descendents)

    std::unordered_map<const CaloHit *, float> mcAssociatedCaloHits;
    float totalMcHitWeight(0.f);

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (useCollectionPlaneOnly && (pCaloHit->GetHitType() != TPC_VIEW_W))
            continue;

        for (const MCParticleWeightMap::value_type &mapPair : pCaloHit->GetMCParticleWeightMap())
        {
            try
            {
                if ((isNeutrino && LArMCParticleHelper::IsBeamNeutrinoFinalState(mapPair.first)) ||
                    (!isNeutrino && LArMCParticleHelper::GetPrimaryMCParticle(mapPair.first) == pMCParticle))
                {
                    mcAssociatedCaloHits.emplace(pCaloHit, mapPair.second);
                    totalMcHitWeight += mapPair.second;
                }
            }

            catch (...)
            {
                continue;
            }
        }
    }

    float numerator(0.f);

    for (const CaloHit *const pPfoAssocCaloHit : pfoAssociatedCaloHits)
    {
        const auto findIter = mcAssociatedCaloHits.find(pPfoAssocCaloHit);

        if (findIter != mcAssociatedCaloHits.end())
            numerator += findIter->second;
    }

    if (pfoAssociatedCaloHits.empty() || (totalMcHitWeight < std::numeric_limits<float>::epsilon()))
    {
        hitPurity       = 0.f;
        hitCompleteness = 0.f;
        return;
    }

    hitPurity       = numerator / static_cast<float>(pfoAssociatedCaloHits.size());
    hitCompleteness = numerator / totalMcHitWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArAnalysisParticleHelper::TypeTreeAsStringImpl(const LArAnalysisParticle::TypeTree &typeTree,
    const bool printTrailingDelimiter)
{
    const std::string delimiter = " - ";
    std::string typeTreeString = TypeAsString(typeTree.Type());

    if (!typeTree.Daughters().empty())
    {
        typeTreeString += delimiter;
        typeTreeString += "[ ";

        for (auto iter = typeTree.Daughters().begin(); iter != typeTree.Daughters().end(); ++iter)
        {
            const bool isLast = (std::next(iter, 1) == typeTree.Daughters().end());
            typeTreeString += TypeTreeAsStringImpl(*iter, !isLast);
        }

        typeTreeString += " ]";
    }

    if (printTrailingDelimiter)
        typeTreeString += delimiter;

    return typeTreeString;
}

} // namespace lar_physics_content
