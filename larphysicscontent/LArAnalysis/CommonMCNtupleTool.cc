/**
 *  @file   larphysicscontent/LArAnalysis/CommonMCNtupleTool.cc
 *
 *  @brief  Implementation of the common MC ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/CommonMCNtupleTool.h"
#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{
CommonMCNtupleTool::CommonMCNtupleTool() : NtupleVariableBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CommonMCNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProcessEvent(const PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProcessNeutrino(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProcessPrimary(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;

    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProcessCosmicRay(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;
    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProduceGenericPfoMCRecords(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const pMCParticle, const MCParticleList *const) const
{
    std::vector<LArNtupleRecord> records;
    records.emplace_back("HasMCInfo", static_cast<LArNtupleRecord::RBool>(pMCParticle));

    if (pMCParticle)
    {
        records.emplace_back("mc_McParticleUid", reinterpret_cast<LArNtupleRecord::RULong64>(pMCParticle->GetUid()));
    }

    else // null values for size consistency
    {
        records.emplace_back("mc_McParticleUid", reinterpret_cast<LArNtupleRecord::RULong64>(0ULL));
    }

    return records;
}

} // namespace lar_physics_content
