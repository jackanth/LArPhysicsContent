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
    m_appendNtuple(false),
    m_minUnmatchedMcParticleEnergy(0.f)
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

    const PfoList &pfoList                              = *pPfoList;
    const auto [neutrinos, cosmicRays, primaries]       = this->GetParticleLists(pfoList);
    const auto [mcNeutrinos, mcCosmicRays, mcPrimaries] = this->GetMCParticleLists(pMCParticleList);

    this->RegisterNtupleRecords(neutrinos, cosmicRays, primaries, pfoList, mcNeutrinos, mcCosmicRays, mcPrimaries, pMCParticleList);
    m_spNtuple->Fill();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<PfoList, PfoList, PfoList> AnalysisNtupleAlgorithm::GetParticleLists(const PfoList &pfoList) const
{
    PfoList neutrinos, cosmicRays, primaries;

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

            default: // everything that isn't a beam neutrino, primary cosmic ray, or primary neutrino daughter
                break;
        }
    }

    return {neutrinos, cosmicRays, primaries};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<MCParticleList, MCParticleList, MCParticleList> AnalysisNtupleAlgorithm::GetMCParticleLists(const MCParticleList *const pMCParticleList) const
{
    MCParticleList mcNeutrinos, mcCosmicRays, mcPrimaries;

    if (!pMCParticleList) // no MC info for the event - return empty lists
        return {mcNeutrinos, mcCosmicRays, mcPrimaries};

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (!pMCParticle)
            continue;

        const LArNtupleHelper::PARTICLE_CLASS particleClass = LArNtupleHelper::GetParticleClass(pMCParticle);

        switch (particleClass)
        {
            case LArNtupleHelper::PARTICLE_CLASS::NEUTRINO:
                mcNeutrinos.push_back(pMCParticle);
                break;

            case LArNtupleHelper::PARTICLE_CLASS::PRIMARY:
                mcPrimaries.push_back(pMCParticle);
                break;

            case LArNtupleHelper::PARTICLE_CLASS::COSMIC_RAY:
                mcCosmicRays.push_back(pMCParticle);
                break;

            default: // everything that isn't a beam neutrino, primary cosmic ray, or primary neutrino daughter
                break;
        }
    }

    return {mcNeutrinos, mcCosmicRays, mcPrimaries};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::size_t AnalysisNtupleAlgorithm::RegisterVectorRecords(NtupleVariableBaseTool *const pNtupleTool, const PfoList &particles,
    const MCParticleList &mcParticleList, const pandora::MCParticleList *const pMCParticleList, const VectorRecordProcessor &processor,
    const MCParticleRetriever &mcParticleRetriever) const
{
    // Run over the reco PFOs, matching to MC particles where possible
    MCParticleSet encounteredMCParticles;

    for (const ParticleFlowObject *const pPfo : particles)
    {
        const MCParticle *const pMCParticle = pMCParticleList ? mcParticleRetriever(pPfo) : nullptr;

        for (const LArNtupleRecord &record : processor(pNtupleTool, pPfo, pMCParticle))
            m_spNtuple->AddVectorRecordElement(record);

        if (pMCParticle)
            encounteredMCParticles.insert(pMCParticle);

        m_spNtuple->FillVectors();
    }

    // Find all the MC particles that are in our main classes but not matched to a PFO
    std::size_t extraMCRecords(0UL);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        if (encounteredMCParticles.find(pMCParticle) != encounteredMCParticles.end())
            continue;

        if (!this->TestMCParticleQuality(pMCParticle))
            continue;

        for (const LArNtupleRecord &record : processor(pNtupleTool, nullptr, pMCParticle))
            m_spNtuple->AddVectorRecordElement(record);

        ++extraMCRecords;
        m_spNtuple->FillVectors();
    }

    m_spNtuple->PushVectors();
    return particles.size() + extraMCRecords;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::RegisterNtupleRecords(const PfoList &neutrinos, const PfoList &cosmicRays, const PfoList &primaries,
    const PfoList &pfoList, const MCParticleList &mcNeutrinos, const MCParticleList &mcCosmicRays, const MCParticleList &mcPrimaries,
    const MCParticleList *const pMCParticleList) const
{
    // Register the vector records for all the particles
    const std::size_t numPfoEntries = this->CheckAndRegisterVectorRecords(pfoList, pMCParticleList ? *pMCParticleList : MCParticleList(), pMCParticleList,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
            return pNtupleTool->ProcessParticleWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
        },
        [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCParticleWrapper(pPfo, pMCParticleList); });

    // Register the vector records for all the cosmics
    const std::size_t numCosmicRayEntries = this->CheckAndRegisterVectorRecords(cosmicRays, mcCosmicRays, pMCParticleList,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
            return pNtupleTool->ProcessCosmicRayWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
        },
        [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCCosmicWrapper(pPfo, pMCParticleList); });

    // Register the vector records for all the primaries
    const std::size_t numPrimaryEntries = this->CheckAndRegisterVectorRecords(primaries, mcPrimaries, pMCParticleList,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
            return pNtupleTool->ProcessPrimaryWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
        },
        [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCPrimaryWrapper(pPfo, pMCParticleList); });

    // Register the vector records for all the neutrinos
    const std::size_t numNeutrinoEntries = this->CheckAndRegisterVectorRecords(neutrinos, mcNeutrinos, pMCParticleList,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
            return pNtupleTool->ProcessNeutrinoWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
        },
        [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCNeutrinoWrapper(pPfo, pMCParticleList); });

    // Register the per-event records
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        for (const LArNtupleRecord &record : pNtupleTool->ProcessEventWrapper(this, pfoList, pMCParticleList))
            m_spNtuple->AddScalarRecord(record);
    }

    // Register the standard per-event records (no prefix)
    m_spNtuple->AddScalarRecord(LArNtupleRecord("fileId", static_cast<LArNtupleRecord::RInt>(m_fileIdentifier)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("eventNum", static_cast<LArNtupleRecord::RInt>(m_eventNumber - 1)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numNeutrinoEntries", static_cast<LArNtupleRecord::RUInt>(numNeutrinoEntries)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numCosmicRayEntries", static_cast<LArNtupleRecord::RUInt>(numCosmicRayEntries)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numPrimaryEntries", static_cast<LArNtupleRecord::RUInt>(numPrimaryEntries)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numPfoEntries", static_cast<LArNtupleRecord::RUInt>(numPfoEntries)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("hasMcInfo", static_cast<LArNtupleRecord::RBool>(pMCParticleList)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::size_t AnalysisNtupleAlgorithm::CheckAndRegisterVectorRecords(const PfoList &particleList, const MCParticleList &mcParticleList,
    const MCParticleList *const pMCParticleList, const VectorRecordProcessor &processor, const MCParticleRetriever &mcParticleRetriever) const
{
    std::size_t numEntries(0UL);

    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        const std::size_t newNumEntries =
            this->RegisterVectorRecords(pNtupleTool, particleList, mcParticleList, pMCParticleList, processor, mcParticleRetriever);

        if ((numEntries > 0UL) && (numEntries != newNumEntries))
        {
            std::cerr << "NtupleVariableBaseTool: There was an inconsistent vector record size due to an internal error" << std::endl;
            throw STATUS_CODE_FAILURE;
        }

        numEntries = newNumEntries;
    }

    return numEntries;
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinUnmatchedMcParticleEnergy", m_minUnmatchedMcParticleEnergy));

    m_spNtuple = std::shared_ptr<LArNtuple>(new LArNtuple(m_ntupleOutputFile, m_ntupleTreeName, m_ntupleTreeTitle, m_appendNtuple));

    // Downcast and store the algorithm tools
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "NtupleTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
    {
        if (NtupleVariableBaseTool *const pNtupleTool = dynamic_cast<NtupleVariableBaseTool *const>(pAlgorithmTool))
        {
            pNtupleTool->SetNtuple(m_spNtuple);
            m_ntupleVariableTools.push_back(pNtupleTool);
        }

        else
        {
            std::cerr << "AnalysisNtupleAlgorithm: Failed to cast algorithm tool as ntuple tool" << std::endl;
            throw STATUS_CODE_FAILURE;
        }
    }

    return STATUS_CODE_SUCCESS;
}
} // namespace lar_physics_content
