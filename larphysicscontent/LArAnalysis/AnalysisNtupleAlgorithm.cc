/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.cc
 *
 *  @brief  Implementation of the analysis ntuple algorithm class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h"

#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"
#include "larphysicscontent/LArNtuple/LArNtuple.h"

#include "Pandora/AlgorithmHeaders.h"

#include "TROOT.h"
#include "TF1.h"

using namespace pandora;

namespace lar_physics_content
{

AnalysisNtupleAlgorithm::AnalysisNtupleAlgorithm() :
    m_pfoListName(),
    m_ntupleVariableTools(),
    m_ntupleOutputFile(),
    m_ntupleTreeName("PandoraNtuple"),
    m_ntupleTreeTitle("Pandora Ntuple"),
    m_plotsOutputFile(),
    m_tmpOutputFile(),
    m_spNtuple(nullptr),
    m_fileIdentifier(0),
    m_eventNumber(0),
    m_appendNtuple(false),
    m_minUnmatchedMcParticleEnergy(0.f),
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minFiducialCoordinates(0.f, 0.f, 0.f),
    m_maxFiducialCoordinates(0.f, 0.f, 0.f),
    m_spTmpRegistry(nullptr),
    m_spPlotsRegistry(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisNtupleAlgorithm::Run()
{
    ++m_eventNumber;
    gROOT->Reset();
    
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

    gSystem->ProcessEvents();
    m_spNtuple->Fill();
    m_spTmpRegistry->Clear();
    m_spPlotsRegistry->Write();
    m_spPlotsRegistry->ClearMemory();

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

std::size_t AnalysisNtupleAlgorithm::RegisterVectorRecords(const PfoList &particles, const MCParticleList &mcParticleList,
    const pandora::MCParticleList *const pMCParticleList, const LArNtupleHelper::VECTOR_BRANCH_TYPE type,
    const VectorRecordProcessor &processor, const MCParticleRetriever &mcParticleRetriever) const
{
    // Run over the reco PFOs, matching to MC particles where possible
    MCParticleSet encounteredMCParticles;

    for (const ParticleFlowObject *const pPfo : particles)
    {
        const MCParticle *const pMCParticle = pMCParticleList ? mcParticleRetriever(pPfo) : nullptr;

        for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
        {
            for (const LArNtupleRecord &record : processor(pNtupleTool, pPfo, pMCParticle))
                m_spNtuple->AddVectorRecordElement(record, type);

            if (pMCParticle)
                encounteredMCParticles.insert(pMCParticle);
        }

        m_spNtuple->FillVectors(type);
    }

    // Find all the MC particles that are in our main classes but not matched to a PFO
    std::size_t extraMCRecords(0UL);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        if (encounteredMCParticles.find(pMCParticle) != encounteredMCParticles.end())
            continue;

        if (!this->TestMCParticleQuality(pMCParticle))
            continue;

        for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
        {
            for (const LArNtupleRecord &record : processor(pNtupleTool, nullptr, pMCParticle))
                m_spNtuple->AddVectorRecordElement(record, type);
        }

        ++extraMCRecords;
        m_spNtuple->FillVectors(type);
    }

    m_spNtuple->PushVectors(type);
    return particles.size() + extraMCRecords;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::RegisterNtupleRecords(const PfoList &neutrinos, const PfoList &cosmicRays, const PfoList &primaries,
    const PfoList &pfoList, const MCParticleList &mcNeutrinos, const MCParticleList &mcCosmicRays, const MCParticleList &mcPrimaries,
    const MCParticleList *const pMCParticleList) const
{
    std::cout << "AnalysisNtupleAlgorithm: Preparing ntuple tools for new event" << std::endl;

    // Prepare the tools
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
        pNtupleTool->PrepareEventWrapper(this, pfoList, pMCParticleList);

    std::cout << "AnalysisNtupleAlgorithm: Registering PFO records" << std::endl;

    // Register the vector records for all the particles
    const std::size_t numPfoEntries = this->RegisterVectorRecords(pfoList, pMCParticleList ? *pMCParticleList : MCParticleList(),
        pMCParticleList, LArNtupleHelper::VECTOR_BRANCH_TYPE::PARTICLE,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
            return pNtupleTool->ProcessParticleWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
        },
        [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCParticleWrapper(pPfo, pMCParticleList); });

    std::cout << "AnalysisNtupleAlgorithm: Registering cosmic records" << std::endl;

    // Register the vector records for all the cosmics
    const std::size_t numCosmicRayEntries =
        this->RegisterVectorRecords(cosmicRays, mcCosmicRays, pMCParticleList, LArNtupleHelper::VECTOR_BRANCH_TYPE::COSMIC_RAY,
            [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
                return pNtupleTool->ProcessCosmicRayWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
            },
            [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCCosmicWrapper(pPfo, pMCParticleList); });

    std::cout << "AnalysisNtupleAlgorithm: Registering primary records" << std::endl;

    // Register the vector records for all the primaries
    const std::size_t numPrimaryEntries =
        this->RegisterVectorRecords(primaries, mcPrimaries, pMCParticleList, LArNtupleHelper::VECTOR_BRANCH_TYPE::PRIMARY,
            [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
                return pNtupleTool->ProcessPrimaryWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
            },
            [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCPrimaryWrapper(pPfo, pMCParticleList); });

    std::cout << "AnalysisNtupleAlgorithm: Registering neutrino records" << std::endl;

    // Register the vector records for all the neutrinos
    const std::size_t numNeutrinoEntries =
        this->RegisterVectorRecords(neutrinos, mcNeutrinos, pMCParticleList, LArNtupleHelper::VECTOR_BRANCH_TYPE::NEUTRINO,
            [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) {
                return pNtupleTool->ProcessNeutrinoWrapper(this, pPfo, pfoList, pMCParticle, pMCParticleList);
            },
            [&](const ParticleFlowObject *const pPfo) { return m_spNtuple->GetMCNeutrinoWrapper(pPfo, pMCParticleList); });

    std::cout << "AnalysisNtupleAlgorithm: Registering event records" << std::endl;

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

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PlotsOutputFile", m_plotsOutputFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TmpOutputFile", m_tmpOutputFile));

    gROOT->SetBatch(kTRUE);
    TF1::DefaultAddToGlobalList(false);

    m_spTmpRegistry  = std::shared_ptr<LArRootRegistry>(new LArRootRegistry(m_tmpOutputFile, LArRootRegistry::FILE_MODE::OVERWRITE));
    m_spPlotsRegistry = std::shared_ptr<LArRootRegistry>(new LArRootRegistry(m_plotsOutputFile, LArRootRegistry::FILE_MODE::APPEND));

    // Get the minimum and maximum fiducial coordinates
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowMargins", m_fiducialCutLowMargins));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighMargins", m_fiducialCutHighMargins));

    std::tie(m_minFiducialCoordinates, m_maxFiducialCoordinates) =
        LArAnalysisHelper::GetFiducialCutCoordinates(this->GetPandora(), m_fiducialCutLowMargins, m_fiducialCutHighMargins);

    m_spNtuple = std::shared_ptr<LArNtuple>(new LArNtuple(m_ntupleOutputFile, m_ntupleTreeName, m_ntupleTreeTitle, m_appendNtuple));

    // Downcast and store the algorithm tools
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "NtupleTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
    {
        if (NtupleVariableBaseTool *const pNtupleTool = dynamic_cast<NtupleVariableBaseTool *const>(pAlgorithmTool))
        {
            pNtupleTool->Setup(m_spNtuple, this, m_minFiducialCoordinates, m_maxFiducialCoordinates, m_spPlotsRegistry, m_spTmpRegistry);
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
