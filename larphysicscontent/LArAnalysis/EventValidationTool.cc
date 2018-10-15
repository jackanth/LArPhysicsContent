/**
 *  @file   larphysicscontent/LArAnalysis/EventValidationTool.cc
 *
 *  @brief  Implementation of the event validation tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/EventValidationTool.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

EventValidationTool::EventValidationTool() :
    AlgorithmTool(),
    m_useTrueNeutrinosOnly(false),
    m_selectInputHits(true),
    m_minHitSharingFraction(0.9f),
    m_maxPhotonPropagation(2.5f),
    m_useSmallPrimaries(true),
    m_matchingMinSharedHits(5),
    m_matchingMinCompleteness(0.1f),
    m_matchingMinPurity(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::shared_ptr<LArInteractionValidationInfo>> EventValidationTool::RunValidation(
    const PfoList &pfoList, const CaloHitList &caloHitList, const MCParticleList &mcParticleList) const
{
    ValidationInfo validationInfo;
    this->FillValidationInfo(mcParticleList, caloHitList, pfoList, validationInfo);
    return this->GetEventValidationInfo(validationInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationTool::FillValidationInfo(
    const MCParticleList &mcParticleList, const CaloHitList &caloHitList, const PfoList &pfoList, ValidationInfo &validationInfo) const
{
    LArMCParticleHelper::PrimaryParameters parameters;

    parameters.m_selectInputHits       = m_selectInputHits;
    parameters.m_minHitSharingFraction = m_minHitSharingFraction;
    parameters.m_maxPhotonPropagation  = m_maxPhotonPropagation;

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(
        &mcParticleList, &caloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    if (!m_useTrueNeutrinosOnly)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(
            &mcParticleList, &caloHitList, parameters, LArMCParticleHelper::IsBeamParticle, targetMCParticleToHitsMap);
    }

    if (!m_useTrueNeutrinosOnly)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(
            &mcParticleList, &caloHitList, parameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);
    }

    parameters.m_minPrimaryGoodHits    = 0;
    parameters.m_minHitsForGoodView    = 0;
    parameters.m_minHitSharingFraction = 0.f;

    LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(
        &mcParticleList, &caloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, allMCParticleToHitsMap);

    if (!m_useTrueNeutrinosOnly)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(
            &mcParticleList, &caloHitList, parameters, LArMCParticleHelper::IsBeamParticle, allMCParticleToHitsMap);
    }

    if (!m_useTrueNeutrinosOnly)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(
            &mcParticleList, &caloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);
    }

    validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
    validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);

    PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(pfoList, allConnectedPfos);

    PfoList finalStatePfos;
    for (const ParticleFlowObject *const pPfo : allConnectedPfos)
    {
        if (LArPfoHelper::IsFinalState(pPfo))
            finalStatePfos.push_back(pPfo);
    }

    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap);
    validationInfo.SetPfoToHitsMap(pfoToHitsMap);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;

    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(
        validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);

    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationTool::InterpretMatching(
    const ValidationInfo &validationInfo, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetAllMCParticleToHitsMap()}, mcPrimaryVector);

    PfoSet usedPfos;
    while (this->GetStrongestPfoMatch(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap))
    {
    }
    this->GetRemainingPfoMatches(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap);

    // Ensure all primaries have an entry, and sorting is as desired
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        LArMCParticleHelper::PfoToSharedHitsVector &pfoHitPairs(interpretedMCToPfoHitSharingMap[pMCPrimary]);
        std::sort(pfoHitPairs.begin(), pfoHitPairs.end(),
            [](const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool {
                return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArPfoHelper::SortByNHits(a.first, b.first));
            });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationTool::GetStrongestPfoMatch(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    const MCParticle *                      pBestMCParticle(nullptr);
    LArMCParticleHelper::PfoCaloHitListPair bestPfoHitPair(nullptr, CaloHitList());

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (interpretedMCToPfoHitSharingMap.count(pMCPrimary))
            continue;

        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            if (!this->IsGoodMatch(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary),
                    validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first), pfoToSharedHits.second))
                continue;

            if (pfoToSharedHits.second.size() > bestPfoHitPair.second.size())
            {
                pBestMCParticle = pMCPrimary;
                bestPfoHitPair  = pfoToSharedHits;
            }
        }
    }

    if (!pBestMCParticle || !bestPfoHitPair.first)
        return false;

    interpretedMCToPfoHitSharingMap[pBestMCParticle].push_back(bestPfoHitPair);
    usedPfos.insert(bestPfoHitPair.first);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationTool::GetRemainingPfoMatches(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    const PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            const LArMCParticleHelper::MCParticleCaloHitListPair        mcParticleToHits(pMCPrimary, pfoToSharedHits.second);
            LArMCParticleHelper::PfoToMCParticleHitSharingMap::iterator iter(pfoToMCParticleHitSharingMap.find(pfoToSharedHits.first));

            if (pfoToMCParticleHitSharingMap.end() == iter)
            {
                pfoToMCParticleHitSharingMap[pfoToSharedHits.first].push_back(mcParticleToHits);
            }
            else
            {
                if (1 != iter->second.size())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                LArMCParticleHelper::MCParticleCaloHitListPair &originalMCParticleToHits(iter->second.at(0));

                if (mcParticleToHits.second.size() > originalMCParticleToHits.second.size())
                    originalMCParticleToHits = mcParticleToHits;
            }
        }
    }

    for (const auto &mapEntry : pfoToMCParticleHitSharingMap)
    {
        const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleToHits(mapEntry.second.at(0));
        interpretedMCToPfoHitSharingMap[mcParticleToHits.first].push_back(
            LArMCParticleHelper::PfoCaloHitListPair(mapEntry.first, mcParticleToHits.second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::shared_ptr<LArInteractionValidationInfo>> EventValidationTool::GetEventValidationInfo(const ValidationInfo &validationInfo) const
{
    std::vector<std::shared_ptr<LArInteractionValidationInfo>> interactionValidationVector;
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &  mcToPfoHitSharingMap(validationInfo.GetInterpretedMCToPfoHitSharingMap());

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    int nNeutrinoPrimaries(0);
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            ++nNeutrinoPrimaries;
    }

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    int        pfoIndex(0), neutrinoPfoIndex(0);
    PfoToIdMap pfoToIdMap, neutrinoPfoToIdMap;
    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
    {
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, ++pfoIndex));
        const Pfo *const pRecoNeutrino(LArPfoHelper::IsNeutrinoFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentNeutrino(pPrimaryPfo) : nullptr);

        if (pRecoNeutrino && !neutrinoPfoToIdMap.count(pRecoNeutrino))
            neutrinoPfoToIdMap.insert(PfoToIdMap::value_type(pRecoNeutrino, ++neutrinoPfoIndex));
    }

    PfoSet         recoNeutrinos;
    MCParticleList associatedMCPrimaries;

    int nCorrectNu(0), nTotalNu(0), nCorrectCR(0), nTotalCR(0), nFakeNu(0), nFakeCR(0), nSplitNu(0), nSplitCR(0), nLost(0);
    int mcPrimaryIndex(0), nTargetMatches(0), nTargetNuMatches(0), nTargetCRMatches(0), nTargetGoodNuMatches(0), nTargetNuSplits(0),
        nTargetNuLosses(0);

    std::shared_ptr<LArInteractionValidationInfo> spNeutrinoInteraction = nullptr;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary));

        if (!hasMatch && !isTargetPrimary)
            continue;

        associatedMCPrimaries.push_back(pMCPrimary);
        const int  nTargetPrimaries(associatedMCPrimaries.size());
        const bool isLastNeutrinoPrimary(++mcPrimaryIndex == nNeutrinoPrimaries);

        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        const int mcNuanceCode(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        const int isBeamNeutrinoFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary));
        const int isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));

        // Print the MC targets
        std::shared_ptr<LArInteractionValidationInfo> spInteractionInfo = nullptr;

        if (mcPrimaryIndex <= nNeutrinoPrimaries)
        {
            if (!spNeutrinoInteraction)
            {
                spNeutrinoInteraction = std::shared_ptr<LArInteractionValidationInfo>(new LArInteractionValidationInfo(mcNuanceCode, isCosmicRay));
                interactionValidationVector.push_back(spNeutrinoInteraction);
            }

            spInteractionInfo = spNeutrinoInteraction;
        }

        else
            spInteractionInfo = std::shared_ptr<LArInteractionValidationInfo>(new LArInteractionValidationInfo(mcNuanceCode, isCosmicRay));

        if (!spInteractionInfo)
        {
            std::cerr << "EventValidationTool: Could not find/create interaction" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        if (mcPrimaryIndex > nNeutrinoPrimaries)
            interactionValidationVector.push_back(spInteractionInfo);

        const std::shared_ptr<LArMCTargetValidationInfo> &spTargetInfo =
            spInteractionInfo->AddDaughterTarget(pMCPrimary, mcPrimaryHitList, isTargetPrimary);
        int nPrimaryMatches(0), nPrimaryNuMatches(0), nPrimaryCRMatches(0), nPrimaryGoodNuMatches(0), nPrimaryNuSplits(0);

        bool isBestMatch(true);

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool  isRecoNeutrinoFinalState(LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first));
            const float purity       = this->GetPurity(pfoHitList, sharedHitList);
            const float completeness = this->GetCompleteness(mcPrimaryHitList, sharedHitList);
            const bool  isGoodMatch  = this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList);

            if (isGoodMatch)
                ++nPrimaryMatches;

            if (isRecoNeutrinoFinalState)
            {
                const Pfo *const pRecoNeutrino(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first));
                const bool       isSplitRecoNeutrino(!recoNeutrinos.empty() && !recoNeutrinos.count(pRecoNeutrino));
                if (!isSplitRecoNeutrino && isGoodMatch)
                    ++nPrimaryGoodNuMatches;
                if (isSplitRecoNeutrino && isBeamNeutrinoFinalState && isGoodMatch)
                    ++nPrimaryNuSplits;
                recoNeutrinos.insert(pRecoNeutrino);
            }

            if (isRecoNeutrinoFinalState && isGoodMatch)
                ++nPrimaryNuMatches;
            if (!isRecoNeutrinoFinalState && isGoodMatch)
                ++nPrimaryCRMatches;

            // Print the matched PFO
            spTargetInfo->AddDaughterMatch(
                pfoToSharedHits.first, !isRecoNeutrinoFinalState, sharedHitList, pfoHitList, purity, completeness, isGoodMatch, isBestMatch);
            isBestMatch = false;
        }

        nTargetMatches += nPrimaryMatches;
        nTargetNuMatches += nPrimaryNuMatches;
        nTargetCRMatches += nPrimaryCRMatches;
        nTargetGoodNuMatches += nPrimaryGoodNuMatches;
        nTargetNuSplits += nPrimaryNuSplits;

        if (0 == nPrimaryMatches)
            ++nTargetNuLosses;

        if (isLastNeutrinoPrimary || isCosmicRay)
        {
            const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(associatedMCPrimaries));

            // ATTN Some redundancy introduced to contributing variables
            const int isCorrectNu(isBeamNeutrinoFinalState && (nTargetGoodNuMatches == nTargetNuMatches) && (nTargetGoodNuMatches == nTargetPrimaries) &&
                                  (nTargetCRMatches == 0) && (nTargetNuSplits == 0) && (nTargetNuLosses == 0));
            const int isCorrectCR(isCosmicRay && (nTargetNuMatches == 0) && (nTargetCRMatches == 1));
            const int isFakeNu(isCosmicRay && (nTargetNuMatches > 0));
            const int isFakeCR(!isCosmicRay && (nTargetCRMatches > 0));
            const int isSplitNu(!isCosmicRay && ((nTargetNuMatches > nTargetPrimaries) || (nTargetNuSplits > 0)));
            const int isSplitCR(isCosmicRay && (nTargetCRMatches > 1));
            const int isLost(nTargetMatches == 0);

            const ParticleFlowObject *pRecoNeutrino(nullptr);
            const MCParticle *pMcNeutrino(nullptr);

            if (isLastNeutrinoPrimary)
            {
                if (recoNeutrinos.size() > 1UL)
                {
                    std::cerr << "EventValidationTool: More than one reco neutrino" << std::endl;
                    throw StatusCodeException(STATUS_CODE_FAILURE);
                }

                pRecoNeutrino = recoNeutrinos.empty() ? nullptr : *recoNeutrinos.begin();

                pMcNeutrino = LArMCParticleHelper::GetParentMCParticle(pMCPrimary);

                if (!LArMCParticleHelper::IsNeutrino(pMcNeutrino))
                {
                    std::cerr << "EventValidationTool: Neutrino primary parent was not neutrino" << std::endl;
                    throw StatusCodeException(STATUS_CODE_FAILURE);
                }
            }

            spInteractionInfo->SetParameters(interactionType, (isCorrectNu || isCorrectCR), (isFakeNu || isFakeCR),
                (isSplitNu || isSplitCR), isLost, pRecoNeutrino, pMcNeutrino);

            if (isLastNeutrinoPrimary)
                ++nTotalNu;
            if (isCosmicRay)
                ++nTotalCR;
            if (isCorrectNu)
                ++nCorrectNu;
            if (isCorrectCR)
                ++nCorrectCR;
            if (isFakeNu)
                ++nFakeNu;
            if (isFakeCR)
                ++nFakeCR;
            if (isSplitNu)
                ++nSplitNu;
            if (isSplitCR)
                ++nSplitCR;
            if (isLost)
                ++nLost;

            recoNeutrinos.clear();
            associatedMCPrimaries.clear();
            nTargetMatches       = 0;
            nTargetNuMatches     = 0;
            nTargetCRMatches     = 0;
            nTargetGoodNuMatches = 0;
            nTargetNuSplits      = 0;
            nTargetNuLosses      = 0;
        }
    }

    for (const std::shared_ptr<LArInteractionValidationInfo> &spInteractionValidationInfo : interactionValidationVector)
    {
        if (!spInteractionValidationInfo->AreParametersSet())
        {
            std::cerr << "EventValidationTool: Not all interactions had their parameters set" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }

    return interactionValidationVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationTool::IsGoodMatch(const CaloHitList &trueHits, const CaloHitList &recoHits, const CaloHitList &sharedHits) const
{
    const float purity       = this->GetPurity(recoHits, sharedHits);
    const float completeness = this->GetCompleteness(trueHits, sharedHits);
    return ((sharedHits.size() >= m_matchingMinSharedHits) && (purity >= m_matchingMinPurity) && (completeness >= m_matchingMinCompleteness));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectInputHits", m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHitSharingFraction", m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseSmallPrimaries", m_useSmallPrimaries));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchingMinSharedHits", m_matchingMinSharedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MatchingMinCompleteness", m_matchingMinCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchingMinPurity", m_matchingMinPurity));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
