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

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "TF1.h"
#include "TROOT.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

AnalysisNtupleAlgorithm::AnalysisNtupleAlgorithm() :
    m_eventNumber(0UL),
    m_pEventValidationTool(nullptr),
    m_caloHitListName(),
    m_mcParticleListName(),
    m_printValidation(false),
    m_produceAllOutcomes(true),
    m_pfoListName(),
    m_ntupleOutputFile(),
    m_ntupleTreeName("PandoraNtuple"),
    m_ntupleTreeTitle("Pandora Ntuple"),
    m_plotsOutputFile(),
    m_tmpOutputFile(),
    m_spNtuple(nullptr),
    m_fileIdentifier(0),
    m_appendNtuple(false),
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minFiducialCoordinates(0.f, 0.f, 0.f),
    m_maxFiducialCoordinates(0.f, 0.f, 0.f),
    m_spTmpRegistry(nullptr),
    m_spPlotsRegistry(nullptr),
    m_ntupleVariableTools(),
    m_batchMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisNtupleAlgorithm::Run()
{
    ++m_eventNumber;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const MCParticleList *pMCParticleList(nullptr);

    if (!m_mcParticleListName.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    if (m_produceAllOutcomes) // treat each hypothesis as its own event
    {
        PfoVector        clearCosmics;
        PfoHypothesisMap pfoHypotheses;
        this->CollectAllPfoOutcomes(clearCosmics, pfoHypotheses);

        for (unsigned int hypothesisId = 0UL, numHypotheses = pfoHypotheses.size(); hypothesisId < numHypotheses; ++hypothesisId)
        {
            const auto findIter = pfoHypotheses.find(hypothesisId);

            if (findIter == pfoHypotheses.end())
            {
                std::cerr << "AnalysisNtupleAlgorithm: Could not find hypothesis with index " << hypothesisId << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            const PfoVector &hypothesisPfos = findIter->second;
            PfoList          allPfos(hypothesisPfos.begin(), hypothesisPfos.end());
            allPfos.insert(allPfos.end(), clearCosmics.begin(), clearCosmics.end());

            this->ProcessEventHypothesis(hypothesisId, allPfos, *pCaloHitList, pMCParticleList);
        }
    }

    else // use the PFO list name
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

        PfoVector allConnectedPfosVector;
        this->CollectPfos(*pPfoList, allConnectedPfosVector);
        PfoList allConnectedPfos(allConnectedPfosVector.begin(), allConnectedPfosVector.end());

        this->ProcessEventHypothesis(-1, allConnectedPfos, *pCaloHitList, pMCParticleList);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::CollectAllPfoOutcomes(PfoVector &clearCosmics, PfoHypothesisMap &pfoHypotheses) const
{
    clearCosmics.clear();
    pfoHypotheses.clear();

    const PfoList *pParentPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pParentPfoList));

    const auto [pSlicingWorker, pSliceNuWorker, pSliceCRWorker] = this->GetPandoraWorkers();
    this->BuildHypotheses(pfoHypotheses, pSlicingWorker, pSliceNuWorker, pSliceCRWorker);

    // Collect clear cosmic-rays
    PfoList clearCosmicsList;

    for (const ParticleFlowObject *const pPfo : *pParentPfoList)
    {
        bool        isClearCosmic(false);
        const auto &properties(pPfo->GetPropertiesMap());
        const auto  it(properties.find("IsClearCosmic"));

        if (it != properties.end())
            isClearCosmic = static_cast<bool>(std::round(pPfo->GetPropertiesMap().at("IsClearCosmic")));

        if (!isClearCosmic)
            continue;

        clearCosmicsList.push_back(pPfo);
    }

    this->CollectPfos(clearCosmicsList, clearCosmics);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<const Pandora *, const Pandora *, const Pandora *> AnalysisNtupleAlgorithm::GetPandoraWorkers() const
{
    // Identify the pandora worker instances by their name
    const Pandora *pSlicingWorker(nullptr);
    const Pandora *pSliceNuWorker(nullptr);
    const Pandora *pSliceCRWorker(nullptr);

    for (const Pandora *const pPandora : MultiPandoraApi::GetDaughterPandoraInstanceList(&this->GetPandora()))
    {
        const std::string &name(pPandora->GetName());

        if (name == "SlicingWorker")
        {
            if (pSlicingWorker)
            {
                std::cerr << "AnalysisNtupleAlgorithm: Multiple slice worker instances!" << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            pSlicingWorker = pPandora;
        }

        else if (name == "SliceNuWorker")
        {
            if (pSliceNuWorker)
            {
                std::cerr << "AnalysisNtupleAlgorithm: Multiple neutrino slice worker instances!" << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            pSliceNuWorker = pPandora;
        }

        else if (name == "SliceCRWorker")
        {
            if (pSliceCRWorker)
            {
                std::cerr << "AnalysisNtupleAlgorithm: Multiple cosmic ray slice worker instances!" << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            pSliceCRWorker = pPandora;
        }
    }

    if (!pSlicingWorker || !pSliceNuWorker || !pSliceCRWorker)
    {
        std::cerr << "AnalysisNtupleAlgorithm: Can't produce all outcomes for a non-consolidated pandora producer" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    return {pSlicingWorker, pSliceNuWorker, pSliceCRWorker};
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::BuildHypotheses(
    PfoHypothesisMap &pfoHypotheses, const Pandora *pSlicingWorker, const Pandora *pSliceNuWorker, const Pandora *pSliceCRWorker) const
{
    // Collect slices under both reconstruction outcomes
    const PfoList *pSlicePfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pSlicingWorker, pSlicePfoList));

    for (unsigned int nuSliceIndex = 0UL; nuSliceIndex < pSlicePfoList->size(); ++nuSliceIndex)
    {
        PfoList        pfoList;
        const PfoList *pNuPfoList(nullptr);

        if (STATUS_CODE_SUCCESS == PandoraApi::GetPfoList(*pSliceNuWorker, "NeutrinoParticles3D" + std::to_string(nuSliceIndex), pNuPfoList))
            pfoList.insert(pfoList.end(), pNuPfoList->begin(), pNuPfoList->end());

        for (unsigned int crSliceIndex = 0UL; crSliceIndex < pSlicePfoList->size(); ++crSliceIndex)
        {
            if (crSliceIndex == nuSliceIndex)
                continue;

            const PfoList *pCRPfoList(nullptr);

            if (STATUS_CODE_SUCCESS == PandoraApi::GetPfoList(*pSliceCRWorker, "MuonParticles3D" + std::to_string(crSliceIndex), pCRPfoList))
                pfoList.insert(pfoList.end(), pCRPfoList->begin(), pCRPfoList->end());
        }

        PfoVector pfoVector;
        this->CollectPfos(pfoList, pfoVector);

        pfoHypotheses.emplace(nuSliceIndex + 1U, std::move(pfoVector));
    }

    // There's one more hypothesis: that everything is a cosmic ray - call this hypothesis 0
    PfoList allCosmicsPfoList;

    for (unsigned int crSliceIndex = 0UL; crSliceIndex < pSlicePfoList->size(); ++crSliceIndex)
    {
        const PfoList *pCRPfoList(nullptr);

        if (STATUS_CODE_SUCCESS == PandoraApi::GetPfoList(*pSliceCRWorker, "MuonParticles3D" + std::to_string(crSliceIndex), pCRPfoList))
            allCosmicsPfoList.insert(allCosmicsPfoList.end(), pCRPfoList->begin(), pCRPfoList->end());
    }

    PfoVector allCosmicsPfoVector;
    this->CollectPfos(allCosmicsPfoList, allCosmicsPfoVector);

    pfoHypotheses.emplace(0U, std::move(allCosmicsPfoVector));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::CollectPfos(const PfoList &parentPfoList, PfoVector &pfoVector) const
{
    if (!pfoVector.empty())
    {
        std::cerr << "AnalysisNtupleAlgorithm: Trying to collect pfos into a non-empty list" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    PfoList pfoList;
    LArPfoHelper::GetAllConnectedPfos(parentPfoList, pfoList);

    PfoSet pfoSet(pfoList.begin(), pfoList.end());

    pfoVector.insert(pfoVector.end(), pfoSet.begin(), pfoSet.end());
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::ProcessEventHypothesis(
    const int hypothesisId, const PfoList &allPfos, const CaloHitList &caloHitList, const MCParticleList *const pMCParticleList) const
{
    // Prepare the ntuple state in case previous instance encountered an exception
    gROOT->Reset();
    m_spNtuple->Reset();

    std::vector<std::shared_ptr<LArInteractionValidationInfo>> eventValidationInfo;

    if (pMCParticleList && m_pEventValidationTool)
    {
        eventValidationInfo = m_pEventValidationTool->RunValidation(allPfos, caloHitList, *pMCParticleList);

        if (m_printValidation)
            this->PrintValidation(eventValidationInfo);
    }

    const auto [neutrinos, cosmicRays, primaries]                     = this->GetParticleLists(allPfos);
    const auto [mcNeutrinoInts, mcCosmicRayTargets, mcPrimaryTargets] = this->GetMCParticleLists(eventValidationInfo);
    const auto [pfoToTargetMap, pfoToInteractionMap]                  = this->GetPfoToMcObjectMaps(eventValidationInfo);

    try
    {
        this->RegisterNtupleRecords(hypothesisId, neutrinos, cosmicRays, primaries, allPfos, mcNeutrinoInts, mcCosmicRayTargets,
            mcPrimaryTargets, pfoToTargetMap, pfoToInteractionMap, eventValidationInfo);
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "AnalysisNtupleAlgorithm: Error: " << err.what() << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
    catch (...)
    {
        std::cerr << "AnalysisNtupleAlgorithm: Unknown error" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    gSystem->ProcessEvents();
    m_spNtuple->Fill();
    m_spTmpRegistry->Clear();
    m_spPlotsRegistry->Write();
    m_spPlotsRegistry->ClearMemory();

    for (TObject *pFunction : *gROOT->GetListOfFunctions())
    {
        if (pFunction)
        {
            delete pFunction;
            pFunction = nullptr;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::PrintValidation(const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const
{
    std::cout << "---NTUPLE-VALIDATION-OUTPUT---------------------------------------------------------------------" << std::endl;

    for (const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo : eventValidationInfo)
    {
        std::cout << LArInteractionTypeHelper::ToString(spInteractionInfo->GetInteractionType()) << " (Nuance "
                  << spInteractionInfo->GetNuanceCode() << ", Nu " << !spInteractionInfo->IsCosmicRay() << ", CR "
                  << spInteractionInfo->IsCosmicRay() << ")" << std::endl;

        for (const std::shared_ptr<LArMCTargetValidationInfo> &spTargetInfo : spInteractionInfo->GetDaughterTargets())
        {
            if (spTargetInfo->GetDaughterMatches().empty() && !spTargetInfo->IsTargetMCPrimary())
                continue;

            const MCParticle *const pMCPrimary = spTargetInfo->GetMCParticle();

            // clang-format off
            std::cout << (!spTargetInfo->IsTargetMCPrimary() ? "(Non target) " : "")
                    << "PrimaryId ?"
                    << ", Nu " << !spInteractionInfo->IsCosmicRay()
                    << ", TB 0"
                    << ", CR " << spInteractionInfo->IsCosmicRay()
                    << ", MCPDG " << pMCPrimary->GetParticleId()
                    << ", Energy " << pMCPrimary->GetEnergy()
                    << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude()
                    << ", nMCHits " << spTargetInfo->GetMCHits().size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, spTargetInfo->GetMCHits())
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, spTargetInfo->GetMCHits())
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, spTargetInfo->GetMCHits()) << ")" << std::endl;
            // clang-format on

            for (const std::shared_ptr<LArMCMatchValidationInfo> &spMatchInfo : spTargetInfo->GetDaughterMatches())
            {
                std::cout << "-" << (!spMatchInfo->IsGoodMatch() ? "(Below threshold) " : "") << "MatchedPfoId ?, Nu "
                          << !spMatchInfo->IsRecoCosmicRay();
                const ParticleFlowObject *const pPfo = spMatchInfo->GetPfo();

                if (!spMatchInfo->IsRecoCosmicRay())
                    std::cout << " [NuId: ?]";

                std::cout << ", CR " << spMatchInfo->IsRecoCosmicRay() << ", PDG " << pPfo->GetParticleId() << ", nMatchedHits "
                          << spMatchInfo->GetSharedHits().size() << " ("
                          << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, spMatchInfo->GetSharedHits()) << ", "
                          << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, spMatchInfo->GetSharedHits()) << ", "
                          << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, spMatchInfo->GetSharedHits()) << ")"
                          << ", nPfoHits " << spMatchInfo->GetPfoHits().size() << " ("
                          << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, spMatchInfo->GetPfoHits()) << ", "
                          << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, spMatchInfo->GetPfoHits()) << ", "
                          << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, spMatchInfo->GetPfoHits()) << ")" << std::endl;
            }
        }

        std::cout << std::endl;
    }
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

std::tuple<AnalysisNtupleAlgorithm::McInteractionVector, AnalysisNtupleAlgorithm::McTargetVector, AnalysisNtupleAlgorithm::McTargetVector>
AnalysisNtupleAlgorithm::GetMCParticleLists(const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const
{
    McInteractionVector mcNeutrinoInt;
    McTargetVector      mcCosmicRayTargets, mcPrimaryTargets;

    // Collect all the neutrino interactions
    for (const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo : eventValidationInfo)
    {
        if (spInteractionInfo->IsCosmicRay())
            continue;

        mcNeutrinoInt.push_back(spInteractionInfo);
    }

    // Collect all the primary- and CR-targets
    for (const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo : eventValidationInfo)
    {
        for (const std::shared_ptr<LArMCTargetValidationInfo> &spTargetInfo : spInteractionInfo->GetDaughterTargets())
        {
            if (spInteractionInfo->IsCosmicRay())
                mcCosmicRayTargets.push_back(spTargetInfo);

            else
                mcPrimaryTargets.push_back(spTargetInfo);
        }
    }

    return {mcNeutrinoInt, mcCosmicRayTargets, mcPrimaryTargets};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<LArAnalysisHelper::PfoToTargetMap, LArAnalysisHelper::PfoToInteractionMap> AnalysisNtupleAlgorithm::GetPfoToMcObjectMaps(
    const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const
{
    LArAnalysisHelper::PfoToTargetMap      pfoToTargetMap;
    LArAnalysisHelper::PfoToInteractionMap pfoToInteractionMap;

    // Collect the primary/CR matches
    for (const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo : eventValidationInfo)
    {
        for (const std::shared_ptr<LArMCTargetValidationInfo> &spTargetInfo : spInteractionInfo->GetDaughterTargets())
        {
            for (const std::shared_ptr<LArMCMatchValidationInfo> &spMatchInfo : spTargetInfo->GetDaughterMatches())
            {
                // Add best matches to the map - and check they are unique
                if (spMatchInfo->IsBestMatch())
                {
                    if (!pfoToTargetMap.emplace(spMatchInfo->GetPfo(), spTargetInfo).second)
                    {
                        std::cerr << "AnalysisNtupleAlgorithm: Best-matched PFO features in map more than once" << std::endl;
                        throw StatusCodeException(STATUS_CODE_FAILURE);
                    }
                }
            }
        }
    }

    // Collect the neutrino interactions
    for (const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo : eventValidationInfo)
    {
        if (spInteractionInfo->IsCosmicRay())
            continue;

        if (const ParticleFlowObject *const pRecoNeutrino = spInteractionInfo->GetRecoNeutrino())
        {
            if (!pfoToInteractionMap.emplace(pRecoNeutrino, spInteractionInfo).second)
            {
                std::cerr << "AnalysisNtupleAlgorithm: Reco neutrino features in map more than once" << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }
    }

    return {pfoToTargetMap, pfoToInteractionMap};
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::RegisterNtupleRecords(const int hypothesisId, const PfoList &neutrinos, const PfoList &cosmicRays,
    const PfoList &primaries, const PfoList &pfoList, const McInteractionVector &mcNeutrinoInts, const McTargetVector &mcCosmicRayTargets,
    const McTargetVector &mcPrimaryTargets, const LArAnalysisHelper::PfoToTargetMap &pfoToTargetMap,
    const LArAnalysisHelper::PfoToInteractionMap &                    pfoToInteractionMap,
    const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const
{
    std::cout << "AnalysisNtupleAlgorithm: Preparing ntuple tools for new event" << std::endl;

    // Prepare the tools
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
        pNtupleTool->PrepareEventWrapper(this, pfoList, eventValidationInfo);

    std::cout << "AnalysisNtupleAlgorithm: Registering cosmic records" << std::endl;

    // Register the vector records for all the cosmics
    const std::size_t numCosmicRayEntries = this->RegisterVectorRecords<LArMCTargetValidationInfo>(cosmicRays, pfoToTargetMap,
        mcCosmicRayTargets, LArNtupleHelper::VECTOR_BRANCH_TYPE::COSMIC_RAY,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) {
            return pNtupleTool->ProcessCosmicRayWrapper(this, pPfo, pfoList, spMcTarget);
        });

    std::cout << "AnalysisNtupleAlgorithm: Registering primary records" << std::endl;

    // Register the vector records for all the primaries
    const std::size_t numPrimaryEntries = this->RegisterVectorRecords<LArMCTargetValidationInfo>(primaries, pfoToTargetMap,
        mcPrimaryTargets, LArNtupleHelper::VECTOR_BRANCH_TYPE::PRIMARY,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) {
            return pNtupleTool->ProcessPrimaryWrapper(this, pPfo, pfoList, spMcTarget);
        });

    std::cout << "AnalysisNtupleAlgorithm: Registering neutrino records" << std::endl;

    // Register the vector records for all the neutrinos
    const std::size_t numNeutrinoEntries = this->RegisterVectorRecords<LArInteractionValidationInfo>(neutrinos, pfoToInteractionMap,
        mcNeutrinoInts, LArNtupleHelper::VECTOR_BRANCH_TYPE::NEUTRINO,
        [&](NtupleVariableBaseTool *const pNtupleTool, const ParticleFlowObject *const pPfo,
            const std::shared_ptr<LArInteractionValidationInfo> &spMcInteraction) {
            return pNtupleTool->ProcessNeutrinoWrapper(this, pPfo, pfoList, spMcInteraction);
        });

    std::cout << "AnalysisNtupleAlgorithm: Registering event records" << std::endl;

    // Register the per-event records
    for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
    {
        for (const LArNtupleRecord &record : pNtupleTool->ProcessEventWrapper(this, pfoList, eventValidationInfo))
            m_spNtuple->AddScalarRecord(record);
    }

    // Register the standard per-event records (no prefix)
    m_spNtuple->AddScalarRecord(LArNtupleRecord("fileId", static_cast<LArNtupleRecord::RInt>(m_fileIdentifier)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("eventNum", static_cast<LArNtupleRecord::RInt>(m_eventNumber - 1)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("hypothesisId", static_cast<LArNtupleRecord::RInt>(hypothesisId)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numNeutrinoEntries", static_cast<LArNtupleRecord::RUInt>(numNeutrinoEntries)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numCosmicRayEntries", static_cast<LArNtupleRecord::RUInt>(numCosmicRayEntries)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("numPrimaryEntries", static_cast<LArNtupleRecord::RUInt>(numPrimaryEntries)));
    m_spNtuple->AddScalarRecord(LArNtupleRecord("hasMcInfo", static_cast<LArNtupleRecord::RBool>(!eventValidationInfo.empty())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisNtupleAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrintValidation", m_printValidation));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProduceAllOutcomes", m_produceAllOutcomes));

    if (!m_produceAllOutcomes)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NtupleOutputFile", m_ntupleOutputFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NtupleTreeName", m_ntupleTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NtupleTreeTitle", m_ntupleTreeTitle));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileIdentifier", m_fileIdentifier));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AppendNtuple", m_appendNtuple));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PlotsOutputFile", m_plotsOutputFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TmpOutputFile", m_tmpOutputFile));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "BatchMode", m_batchMode));
    gROOT->SetBatch(m_batchMode);

    m_spTmpRegistry   = std::shared_ptr<LArRootRegistry>(new LArRootRegistry(m_tmpOutputFile, LArRootRegistry::FILE_MODE::OVERWRITE));
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
    AlgorithmToolVector validationToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "EventValidationTools", validationToolVector));

    if (validationToolVector.size() != 1UL)
    {
        std::cerr << "AnalysisNtupleAlgorithm: Must specify exactly one event validation tool" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    if (EventValidationTool *const pEventValidationTool = dynamic_cast<EventValidationTool *const>(validationToolVector.front()))
        m_pEventValidationTool = pEventValidationTool;

    else
    {
        std::cerr << "AnalysisNtupleAlgorithm: Failed to cast event validation tool" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    // Downcast and store the algorithm tools
    AlgorithmToolVector ntupleToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "NtupleTools", ntupleToolVector));

    for (AlgorithmTool *const pAlgorithmTool : ntupleToolVector)
    {
        if (NtupleVariableBaseTool *const pNtupleTool = dynamic_cast<NtupleVariableBaseTool *const>(pAlgorithmTool))
        {
            pNtupleTool->Setup(m_spNtuple, this, m_minFiducialCoordinates, m_maxFiducialCoordinates, m_spPlotsRegistry, m_spTmpRegistry);
            m_ntupleVariableTools.push_back(pNtupleTool);
        }

        else
        {
            std::cerr << "AnalysisNtupleAlgorithm: Failed to cast algorithm tool as ntuple tool" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
