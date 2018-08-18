/**
 *  @file   larphysicscontent/LArAnalysis/CommonNtupleTool.cc
 *
 *  @brief  Implementation of the common ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/CommonNtupleTool.h"

#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
CommonNtupleTool::CommonNtupleTool() :
    NtupleVariableBaseTool(),
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minFiducialCoordinates(0.f, 0.f, 0.f),
    m_maxFiducialCoordinates(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CommonNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowMargins", m_fiducialCutLowMargins));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighMargins", m_fiducialCutHighMargins));

    std::tie(m_minFiducialCoordinates, m_maxFiducialCoordinates) =
        LArAnalysisHelper::GetFiducialCutCoordinates(this->GetPandora(), m_fiducialCutLowMargins, m_fiducialCutHighMargins);

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

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const, const MCParticleList *const)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoRecords = this->ProduceGenericPfoRecords(pPfo, pfoList);
    records.insert(records.end(), std::make_move_iterator(genericPfoRecords.begin()), std::make_move_iterator(genericPfoRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonNtupleTool::ProduceGenericPfoRecords(const ParticleFlowObject *const pPfo, const PfoList &) const
{
    std::vector<LArNtupleRecord> records;
    records.emplace_back("WasReconstructed", static_cast<LArNtupleRecord::RBool>(pPfo));

    if (pPfo)
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

        const CartesianVector &vertexPosition(vertexList.front()->GetPosition());

        records.emplace_back("IsVertexFiducial", static_cast<LArNtupleRecord::RBool>(LArAnalysisHelper::IsPointFiducial(
                                                     vertexPosition, m_minFiducialCoordinates, m_maxFiducialCoordinates)));
        records.emplace_back("VertexX", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetX()));
        records.emplace_back("VertexY", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetY()));
        records.emplace_back("VertexZ", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetZ()));

        records.emplace_back("FiducialThreeDHitFraction", static_cast<LArNtupleRecord::RFloat>(this->GetFractionOfFiducialThreeDHits(pPfo)));
        records.emplace_back("NumberOfThreeDHits", static_cast<LArNtupleRecord::RUInt>(this->GetAllDownstreamThreeDHits(pPfo).size()));
        records.emplace_back("NumberOfTwoDHits", static_cast<LArNtupleRecord::RUInt>(this->GetAllDownstreamTwoDHits(pPfo).size()));
        records.emplace_back("NumberOfCollectionPlaneHits", static_cast<LArNtupleRecord::RUInt>(this->GetAllDownstreamWHits(pPfo).size()));
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
        if (LArAnalysisHelper::IsPointFiducial(pCaloHit->GetPositionVector(), m_minFiducialCoordinates, m_maxFiducialCoordinates))
            ++fiducialHits;
    }

    return static_cast<float>(fiducialHits) / static_cast<float>(caloHitList.size());
}

} // namespace lar_physics_content
