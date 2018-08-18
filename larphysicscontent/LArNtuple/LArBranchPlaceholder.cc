/**
 *  @file   larphysicscontent/LArNtuple/LArBranchPlaceholder.cc
 *
 *  @brief  Implementation of the lar branch placeholder class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArNtuple/LArBranchPlaceholder.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

using namespace pandora;

namespace lar_physics_content
{

void LArBranchPlaceholder::SetNtupleScalarRecord(const LArNtupleRecord &record)
{
    if (m_valueType != record.ValueType())
    {
        std::cerr << "LArBranchPlaceholder: Could not replace value of branch '" << record.BranchName()
                  << "' because the types did not match" << std::endl;
        throw pandora::STATUS_CODE_NOT_ALLOWED;
    }

    m_spNtupleScalarRecord = NtupleRecordSPtr(new LArNtupleRecord(record));
    this->PopulateMaps();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArBranchPlaceholder::PushNtupleScalarRecord()
{
    if (!m_spNtupleScalarRecord)
    {
        std::cerr << "LArBranchPlaceholder: Failed to push ntuple scalar record because it was null" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    m_ntupleVectorRecord.emplace_back(std::move(m_spNtupleScalarRecord));
    m_spNtupleScalarRecord = nullptr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArBranchPlaceholder::LArBranchPlaceholder(const LArNtupleRecord &record) noexcept :
    m_spNtupleScalarRecord(new LArNtupleRecord(record)),
    m_ntupleVectorRecord(),
    m_valueType(record.ValueType()),
    m_pCacheElement(nullptr),
    m_pfoRecordMap(),
    m_mcParticleRecordMap()
{
    this->PopulateMaps();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArBranchPlaceholder::PopulateMaps()
{
    if (const ParticleFlowObject *const pPfo = m_spNtupleScalarRecord->GetPfo())
        m_pfoRecordMap.emplace(pPfo, m_spNtupleScalarRecord);

    if (const MCParticle *const pMCParticle = m_spNtupleScalarRecord->GetMCParticle())
        m_mcParticleRecordMap.emplace(pMCParticle, m_spNtupleScalarRecord);
}

} // namespace lar_physics_content
