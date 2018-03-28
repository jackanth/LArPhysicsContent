/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisParticleHelper.cxx
 *
 *  @brief  Implementation of the lar analysis particle helper class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArHelpers/LArAnalysisParticleHelper.h"

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
void LArAnalysisParticleHelper::GetFiducialCutParameters(const Pandora &pandoraInstance, const CartesianVector &fiducialCutLowMargins,
    const CartesianVector &fiducialCutHighMargins, CartesianVector &minCoordinates, CartesianVector &maxCoordinates)
{
    const LArTPCMap &larTPCMap(pandoraInstance.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.size() != 1UL)
    {
        std::cout << "LArAnalysisParticleHelper: the number of LArTPCs was not equal to 1" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

    const LArTPC *const pLArTPC(larTPCMap.begin()->second);

    // The extremal coordinates are at the centre +- half of the widths.
    const float xMin(pLArTPC->GetCenterX() - (0.5f * pLArTPC->GetWidthX()) + fiducialCutLowMargins.GetX());
    const float yMin(pLArTPC->GetCenterY() - (0.5f * pLArTPC->GetWidthY()) + fiducialCutLowMargins.GetY());
    const float zMin(pLArTPC->GetCenterZ() - (0.5f * pLArTPC->GetWidthZ()) + fiducialCutLowMargins.GetZ());

    const float xMax(pLArTPC->GetCenterX() + (0.5f * pLArTPC->GetWidthX()) - fiducialCutHighMargins.GetX());
    const float yMax(pLArTPC->GetCenterY() + (0.5f * pLArTPC->GetWidthY()) - fiducialCutHighMargins.GetY());
    const float zMax(pLArTPC->GetCenterZ() + (0.5f * pLArTPC->GetWidthZ()) - fiducialCutHighMargins.GetZ());

    minCoordinates = CartesianVector(xMin, yMin, zMin);
    maxCoordinates = CartesianVector(xMax, yMax, zMax);
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
        const CartesianVector &minToMaxVector = maxPosition - minPosition;

        if (std::fabs(maxCoordinate - displacementAlongFittedTrack) < std::fabs(displacementAlongFittedTrack - minCoordinate))
        {
            // If closer to the maximum coordinate, we want the fit direction and the min-to-max vector to be more anti-aligned.
            if (minToMaxVector.GetDotProduct(fitDirection) > 0.f)
                fitDirection *= -1.f;

        }

        else
        {
            // If closer to the minimum coordinate, we want the fit direction and the min-to-max vector to be more aligned.
            if (minToMaxVector.GetDotProduct(fitDirection) < 0.f)
                fitDirection *= -1.f;
        }
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

std::string LArAnalysisParticleHelper::TypeTreeAsStringImpl(const LArAnalysisParticle::TypeTree &typeTree,
    const bool printTrailingDelimiter)
{
    const std::string delimiter = " - ";
    std::string typeTreeString = LArAnalysisParticle::TypeAsString(typeTree.Type());

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
