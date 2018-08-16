/**
 *  @file   larphysicscontent/LArAnalysis/DefaultNtupleTool.cc
 *
 *  @brief  Implementation of the default ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/DefaultNtupleTool.h"
#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
DefaultNtupleTool::DefaultNtupleTool() :
    NtupleVariableBaseTool(),
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minFiducialCoordinates(0.f, 0.f, 0.f),
    m_maxFiducialCoordinates(0.f, 0.f, 0.f),
    m_cacheDownstreamThreeDHits()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DefaultNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
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

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessEvent(const PfoList &pfoList, const MCParticleList *const)
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

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessNeutrino(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoRecords   = this->ProduceGenericPfoRecords(pPfo, pfoList);
    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);

    records.insert(records.end(), std::make_move_iterator(genericPfoRecords.begin()), std::make_move_iterator(genericPfoRecords.end()));
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessPrimary(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoRecords   = this->ProduceGenericPfoRecords(pPfo, pfoList);
    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);

    records.insert(records.end(), std::make_move_iterator(genericPfoRecords.begin()), std::make_move_iterator(genericPfoRecords.end()));
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProcessCosmicRay(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoRecords   = this->ProduceGenericPfoRecords(pPfo, pfoList);
    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);

    records.insert(records.end(), std::make_move_iterator(genericPfoRecords.begin()), std::make_move_iterator(genericPfoRecords.end()));
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> DefaultNtupleTool::ProduceGenericPfoRecords(const ParticleFlowObject *const pPfo, const PfoList &) const
{
    std::vector<LArNtupleRecord> records;
    records.emplace_back("WasReconstructed", static_cast<LArNtupleRecord::RBool>(pPfo));

    if (pPfo)
    {
        const VertexList &vertexList = pPfo->GetVertexList();

        if (vertexList.empty())
        {
            std::cerr << "DefaultNtupleTool: Could not get vertex because list was empty" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        if (vertexList.size() > 1UL)
        {
            std::cerr << "DefaultNtupleTool: Could not get vertex because there were more than one" << std::endl;
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

    else
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

std::vector<LArNtupleRecord> DefaultNtupleTool::ProduceGenericPfoMCRecords(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const pMCParticle, const MCParticleList *const) const
{
    std::vector<LArNtupleRecord> records;
    records.emplace_back("HasMCInfo", static_cast<LArNtupleRecord::RBool>(pMCParticle));

    if (pMCParticle)
    {
        records.emplace_back("mc_McParticleUid", reinterpret_cast<LArNtupleRecord::RULong64>(pMCParticle->GetUid()));
    }

    else
    {
        records.emplace_back("mc_McParticleUid", reinterpret_cast<LArNtupleRecord::RULong64>(0ULL));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList DefaultNtupleTool::GetAllDownstreamTwoDHits(const ParticleFlowObject *const pPfo) const
{
    // Get all the (possibly-cached) U-, V-, and W-hits and add them all up.
    CaloHitList hitsU(this->GetAllDownstreamUHits(pPfo));
    CaloHitList hitsV(this->GetAllDownstreamVHits(pPfo));
    CaloHitList hitsW(this->GetAllDownstreamWHits(pPfo));

    CaloHitList hitsTwoD;
    hitsTwoD.insert(hitsTwoD.end(), std::make_move_iterator(hitsU.begin()), std::make_move_iterator(hitsU.end()));
    hitsTwoD.insert(hitsTwoD.end(), std::make_move_iterator(hitsV.begin()), std::make_move_iterator(hitsV.end()));
    hitsTwoD.insert(hitsTwoD.end(), std::make_move_iterator(hitsW.begin()), std::make_move_iterator(hitsW.end()));

    return hitsTwoD;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &DefaultNtupleTool::GetAllDownstreamHitsImpl(
    const ParticleFlowObject *const pPfo, const HitType hitType, PfoCache<CaloHitList> &cache) const
{
    // Check the cache first
    const auto findIter = cache.find(pPfo);

    if (findIter != cache.end())
        return findIter->second;

    // Not in cache, so work it out and cache it
    PfoList downstreamPfos;
    LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfos);

    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(downstreamPfos, hitType, caloHitList);

    return cache.emplace(pPfo, std::move(caloHitList)).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DefaultNtupleTool::GetFractionOfFiducialThreeDHits(const ParticleFlowObject *const pPfo) const
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
