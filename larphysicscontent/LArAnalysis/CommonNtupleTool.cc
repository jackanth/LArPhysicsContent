/**
 *  @file   larphysicscontent/LArAnalysis/CommonNtupleTool.cc
 *
 *  @brief  Implementation of the common ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/CommonNtupleTool.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"

#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
CommonNtupleTool::CommonNtupleTool() : NtupleVariableBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CommonNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProcessEvent(const PfoList &pfoList, const MCParticleList *const)
{
    std::vector<LArNtupleRecord> records;
    std::size_t                  numPrimaryTracks(0UL), numPrimaryShowers(0UL);

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        if (LArNtupleHelper::GetParticleClass(pPfo) != LArNtupleHelper::PARTICLE_CLASS::PRIMARY)
            continue;

        if (LArPfoHelper::IsShower(pPfo))
            ++numPrimaryShowers;

        else
            ++numPrimaryTracks;
    }

    records.emplace_back("NumberOfRecoPrimaryTracks", static_cast<LArNtupleRecord::RUInt>(numPrimaryTracks));
    records.emplace_back("NumberOfRecoPrimaryShowers", static_cast<LArNtupleRecord::RUInt>(numPrimaryShowers));
    records.emplace_back("NumOfRecoPfos", static_cast<LArNtupleRecord::RUInt>(pfoList.size()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const, const MCParticleList *const)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoRecords = this->ProduceGenericPfoRecords(pPfo, pfoList);
    records.insert(records.end(), std::make_move_iterator(genericPfoRecords.begin()), std::make_move_iterator(genericPfoRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const, const MCParticleList *const)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoRecords = this->ProduceGenericPfoRecords(pPfo, pfoList);
    records.insert(records.end(), std::make_move_iterator(genericPfoRecords.begin()), std::make_move_iterator(genericPfoRecords.end()));

    std::vector<LArNtupleRecord> nonNuPfoRecords = this->ProduceNonNeutrinoPfoRecords(pPfo, pfoList);
    records.insert(records.end(), std::make_move_iterator(nonNuPfoRecords.begin()), std::make_move_iterator(nonNuPfoRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const, const MCParticleList *const)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoRecords = this->ProduceGenericPfoRecords(pPfo, pfoList);
    records.insert(records.end(), std::make_move_iterator(genericPfoRecords.begin()), std::make_move_iterator(genericPfoRecords.end()));

    std::vector<LArNtupleRecord> nonNuPfoRecords = this->ProduceNonNeutrinoPfoRecords(pPfo, pfoList);
    records.insert(records.end(), std::make_move_iterator(nonNuPfoRecords.begin()), std::make_move_iterator(nonNuPfoRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProduceGenericPfoRecords(const ParticleFlowObject *const pPfo, const PfoList &) const
{
    std::vector<LArNtupleRecord> records;
    records.emplace_back("WasReconstructed", static_cast<LArNtupleRecord::RBool>(pPfo));

    if (pPfo)
    {
        const CartesianVector &vertexPosition = this->GetVertexPosition(pPfo);

        records.emplace_back("IsVertexFiducial", static_cast<LArNtupleRecord::RBool>(this->IsPointFiducial(vertexPosition)));
        records.emplace_back("VertexX", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetX()));
        records.emplace_back("VertexY", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetY()));
        records.emplace_back("VertexZ", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetZ()));
        records.emplace_back("FiducialThreeDHitFraction", static_cast<LArNtupleRecord::RFloat>(this->GetFractionOfFiducialThreeDHits(pPfo)));
        records.emplace_back("NumberOfThreeDHits", static_cast<LArNtupleRecord::RUInt>(this->GetAllDownstreamThreeDHits(pPfo).size()));
        records.emplace_back("NumberOfTwoDHits", static_cast<LArNtupleRecord::RUInt>(this->GetAllDownstreamTwoDHits(pPfo).size()));
        records.emplace_back("NumberOfCollectionPlaneHits", static_cast<LArNtupleRecord::RUInt>(this->GetAllDownstreamWHits(pPfo).size()));
        records.emplace_back("NumberOfDownstreamPfos", static_cast<LArNtupleRecord::RUInt>(this->GetAllDownstreamPfos(pPfo).size()));
    }

    else // null values for size consistency
    {
        records.emplace_back("IsVertexFiducial", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("VertexX", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("VertexY", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("VertexZ", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("FiducialThreeDHitFraction", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("NumberOfThreeDHits", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("NumberOfTwoDHits", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("NumberOfCollectionPlaneHits", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("NumberOfDownstreamPfos", static_cast<LArNtupleRecord::RUInt>(0U));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProduceNonNeutrinoPfoRecords(const ParticleFlowObject *const pPfo, const PfoList &) const
{
    std::vector<LArNtupleRecord> records;

    if (pPfo)
    {
        const bool            isShower         = LArPfoHelper::IsShower(pPfo);
        const CartesianVector directionCosines = isShower ? this->GetShowerDirectionAtVertex(pPfo) : this->GetTrackDirectionAtVertex(pPfo);

        records.emplace_back("IsShower", static_cast<LArNtupleRecord::RBool>(isShower));
        records.emplace_back("IsTrack", static_cast<LArNtupleRecord::RBool>(!isShower));
        records.emplace_back("DirectionCosineX", static_cast<LArNtupleRecord::RFloat>(directionCosines.GetX()));
        records.emplace_back("DirectionCosineY", static_cast<LArNtupleRecord::RFloat>(directionCosines.GetY()));
        records.emplace_back("DirectionCosineZ", static_cast<LArNtupleRecord::RFloat>(directionCosines.GetZ()));
    }

    else // null values for size consistency
    {
        records.emplace_back("IsShower", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("IsTrack", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("DirectionCosineX", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("DirectionCosineY", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("DirectionCosineZ", static_cast<LArNtupleRecord::RFloat>(0.f));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CommonNtupleTool::GetFractionOfFiducialThreeDHits(const ParticleFlowObject *const pPfo) const
{
    const auto &caloHitList = this->GetAllDownstreamThreeDHits(pPfo);
    std::size_t fiducialHits(0UL);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (this->IsPointFiducial(pCaloHit->GetPositionVector()))
            ++fiducialHits;
    }

    return static_cast<float>(fiducialHits) / static_cast<float>(caloHitList.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector CommonNtupleTool::GetShowerDirectionAtVertex(const ParticleFlowObject *const pPfo) const
{
    const CaloHitList &downstreamThreeDHits = this->GetAllDownstreamThreeDHits(pPfo);

    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::EigenValues  eigenValues(0.f, 0.f, 0.f);
    CartesianVector            centroid(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(downstreamThreeDHits, centroid, eigenValues, eigenVectors);

    if (eigenVectors.empty())
    {
        std::cout << "AnalysisAlgorithm: PCA eigenvectors were empty" << std::endl;
        return CartesianVector(0.f, 0.f, 0.f);
    }

    CartesianVector       fitDirection     = eigenVectors.at(0UL);
    const CartesianVector vertexToCentroid = centroid - this->GetVertexPosition(pPfo);

    // We want the fit direction to be mostly aligned with the vertex-to-centroid direction.
    if (vertexToCentroid.GetDotProduct(fitDirection) < 0.f)
        fitDirection *= -1.f;

    return fitDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &CommonNtupleTool::GetVertexPosition(const ParticleFlowObject *const pPfo) const
{
    const VertexList &vertexList = pPfo->GetVertexList();

    if (vertexList.empty())
    {
        std::cerr << "CommonNtupleTool: Could not get vertex because list was empty" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (vertexList.size() > 1UL)
    {
        std::cerr << "CommonNtupleTool: Could not get vertex because there were more than one" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    return vertexList.front()->GetPosition();
}

} // namespace lar_physics_content
