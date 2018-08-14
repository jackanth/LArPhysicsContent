/**
 *  @file   larphysicscontent/LArNtuple/NtupleVariableBaseTool.cc
 *
 *  @brief  Implementation of the ntuple variable base tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"
#include "larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

NtupleVariableBaseTool::NtupleVariableBaseTool() noexcept : AlgorithmTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle * NtupleVariableBaseTool::GetMCParticle(const pandora::ParticleFlowObject *const pPfo,
    const MCParticleList *const)
{
    // Get all the PFO's 2D hits
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);

    CaloHitList caloHitList;

    for (const Cluster *const pCluster : clusterList)
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    // Find the MC particle with the highest weight
    const MCParticleWeightMap mcParticleWeightMap = this->GetMCParticleWeightMap(caloHitList);

    float bestWeight(0.f);
    const MCParticle *pBestMCParticle(nullptr);

    MCParticleVector mcParticleVector;

    for (const MCParticleWeightMap::value_type &mapEntry : mcParticleWeightMap)
        mcParticleVector.push_back(mapEntry.first);

    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

    for (const MCParticle *const pCurrentMCParticle : mcParticleVector)
    {
        const float currentWeight(mcParticleWeightMap.at(pCurrentMCParticle));

        if (currentWeight > bestWeight)
        {
            pBestMCParticle = pCurrentMCParticle;
            bestWeight      = currentWeight;
        }
    }

    return pBestMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessEventWrapper(
    const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::vector<LArNtupleRecord> records = ProcessEvent(pfoList, pMCParticleList);

    for (LArNtupleRecord &record : records)
        record.AddBranchNamePrefix("evt");

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessParticleWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList)
{
    return this->ProcessImpl(pAlgorithm, pPfo, pMCParticleList, "pfo",
        [&](const MCParticle *const pMCParticle) { return this->ProcessParticle(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessNeutrinoWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList)
{
    return this->ProcessImpl(pAlgorithm, pPfo, pMCParticleList, "nu",
        [&](const MCParticle *const pMCParticle) { return this->ProcessNeutrino(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessPrimaryWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList)
{
    return this->ProcessImpl(pAlgorithm, pPfo, pMCParticleList, "primary",
        [&](const MCParticle *const pMCParticle) { return this->ProcessPrimary(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessCosmicRayWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList)
{
    return this->ProcessImpl(pAlgorithm, pPfo, pMCParticleList, "cr",
        [&](const MCParticle *const pMCParticle) { return this->ProcessCosmicRay(pPfo, pfoList, pMCParticle, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessImpl(const AnalysisNtupleAlgorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pPfo, const MCParticleList *const pMCParticleList,
    const std::string &prefix, const Processor &processor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    const MCParticle *const pMCParticle  = pMCParticleList ? this->GetMCParticle(pPfo, pMCParticleList) : nullptr;
    std::vector<LArNtupleRecord> records = processor(pMCParticle);

    for (LArNtupleRecord &record : records)
        record.AddBranchNamePrefix(prefix);

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleWeightMap NtupleVariableBaseTool::GetMCParticleWeightMap(const CaloHitList &caloHitList) const
{
    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        // Synthesize the weights from every hit into the main map
        const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        MCParticleVector mcParticleVector;

        for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap)
            mcParticleVector.push_back(mapEntry.first);

        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            try
            {
                const MCParticle *const pMCPrimary = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
                mcParticleWeightMap[pMCPrimary] += hitMCParticleWeightMap.at(pMCParticle);
            }

            catch (...)
            {
                continue;
            }
        }
    }

    return mcParticleWeightMap;
}
} // namespace lar_physics_content
