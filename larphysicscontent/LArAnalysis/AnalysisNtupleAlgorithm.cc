/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.cc
 *
 *  @brief  Implementation of the analysis ntuple algorithm class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"
#include "larphysicscontent/LArNtuple/LArNtuple.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{

AnalysisNtupleAlgorithm::AnalysisNtupleAlgorithm() :
    m_pfoListName(),
    m_ntupleVariableTools(),
    m_ntupleOutputFile(),
    m_ntupleTreeName("PandoraNtuple"),
    m_ntupleTreeTitle("Pandora Ntuple"),
    m_spNtuple(nullptr),
    m_fileIdentifier(0),
    m_eventNumber(0),
    m_appendNtuple(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisNtupleAlgorithm::Run()
{
    ++m_eventNumber;

    // Get input PFO list
    const PfoList *pPfoList(nullptr);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_pfoListName, pPfoList) || !pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "AnalysisNtupleAlgorithm: cannot find pfo list " << m_pfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // Try to get the MC particle list
    const MCParticleList *pMCParticleList(nullptr);
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);

    // Prepare the ntuple state in case previous instance encountered an exception
    m_spNtuple->Reset();

    const PfoList &pfoList                        = *pPfoList;
    const auto [neutrinos, cosmicRays, primaries] = this->GetParticleLists(pfoList);

    this->RegisterNtupleRecords(neutrinos, cosmicRays, primaries, pfoList, pMCParticleList);
    m_spNtuple->Fill();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<PfoList, PfoList, PfoList> AnalysisNtupleAlgorithm::GetParticleLists(const PfoList &pfoList) const
{
    PfoList particles, neutrinos, cosmicRays, primaries;

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        const LArNtupleHelper::PARTICLE_CLASS particleClass = LArNtupleHelper::GetParticleClass(pPfo);

        switch (particleClass)
        {
            case LArNtupleHelper::PARTICLE_CLASS::NEUTRINO:
                neutrinos.push_back(pPfo);
                break;

            case LArNtupleHelper::PARTICLE_CLASS::PRIMARY:
                primaries.push_back(pPfo);
                break;

            case LArNtupleHelper::PARTICLE_CLASS::COSMIC_RAY:
                cosmicRays.push_back(pPfo);
                break;

            default:
                break;
        }
    }

    return {neutrinos, cosmicRays, primaries};
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::RegisterVectorRecords(const PfoList &particles, const VectorRecordProcessor &processor) const
{
    for (const ParticleFlowObject *const pPfo : particles)
    {
        for (const LArNtupleRecord &record : processor(pPfo))
            m_spNtuple->AddVectorRecordElement(record);

        m_spNtuple->FillVectors();
    }

    m_spNtuple->PushVectors();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::RegisterNtupleRecords(const PfoList &neutrinos, const PfoList &cosmicRays, const PfoList &primaries,
    const PfoList &pfoList, const MCParticleList *const pMCParticleList) const
{
    // Register the reserved records.
    m_spNtuple->AddScalarRecord(LArNtupleRecord("fileId", static_cast<LArNtupleRecord::RInt>(m_fileIdentifier)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("eventNum", static_cast<LArNtupleRecord::RInt>(m_eventNumber - 1)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numNeutrinos", static_cast<LArNtupleRecord::RUInt>(neutrinos.size())));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numCosmicRays", static_cast<LArNtupleRecord::RUInt>(cosmicRays.size())));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numPrimaries", static_cast<LArNtupleRecord::RUInt>(primaries.size())));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numPfos", static_cast<LArNtupleRecord::RUInt>(pfoList.size())));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("hasMcInfo", static_cast<LArNtupleRecord::RBool>(pMCParticleList)));

    // Register the per-event records.
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        for (const LArNtupleRecord &record : pNtupleTool->ProcessEventWrapper(this, pfoList, pMCParticleList))
            m_spNtuple->AddScalarRecord(record);
    }

    // Register the vector records for the neutrinos.
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        this->RegisterVectorRecords(neutrinos,
            [&](const ParticleFlowObject *const pPfo) { return pNtupleTool->ProcessNeutrinoWrapper(this, pPfo, pfoList, pMCParticleList); });
    }

    // Register the vector records for the primaries.
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        this->RegisterVectorRecords(primaries,
            [&](const ParticleFlowObject *const pPfo) { return pNtupleTool->ProcessPrimaryWrapper(this, pPfo, pfoList, pMCParticleList); });
    }

    // Register the vector records for the cosmics.
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        this->RegisterVectorRecords(cosmicRays,
            [&](const ParticleFlowObject *const pPfo) { return pNtupleTool->ProcessCosmicRayWrapper(this, pPfo, pfoList, pMCParticleList); });
    }

    // Register the vector records for all the particles.
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        this->RegisterVectorRecords(pfoList,
            [&](const ParticleFlowObject *const pPfo) { return pNtupleTool->ProcessParticleWrapper(this, pPfo, pfoList, pMCParticleList); });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisNtupleAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NtupleOutputFile", m_ntupleOutputFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NtupleTreeName", m_ntupleTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NtupleTreeTitle", m_ntupleTreeTitle));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileIdentifier", m_fileIdentifier));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AppendNtuple", m_appendNtuple));

    m_spNtuple = std::shared_ptr<LArNtuple>(new LArNtuple(m_ntupleOutputFile, m_ntupleTreeName, m_ntupleTreeTitle, m_appendNtuple));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "NtupleTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
    {
        if (NtupleVariableBaseTool *const pNtupleTool = dynamic_cast<NtupleVariableBaseTool *const>(pAlgorithmTool))
            m_ntupleVariableTools.push_back(pNtupleTool);

        else
        {
            std::cerr << "AnalysisNtupleAlgorithm: Failed to cast algorithm tool as ntuple tool" << std::endl;
            throw STATUS_CODE_FAILURE;
        }
    }

    return STATUS_CODE_SUCCESS;
}
} // namespace lar_physics_content
