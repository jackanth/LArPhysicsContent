/**
 *  @file   larphysicscontent/LArAnalysis/CommonMCNtupleTool.cc
 *
 *  @brief  Implementation of the common MC ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/CommonMCNtupleTool.h"
#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
CommonMCNtupleTool::CommonMCNtupleTool() : NtupleVariableBaseTool(), m_twoDCaloHitListName(), m_pTwoDCaloHitList(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CommonMCNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TwoDCaloHitListName", m_twoDCaloHitListName));
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CommonMCNtupleTool::PrepareEvent(const PfoList &, const MCParticleList *const)
{
    if ((PandoraContentApi::GetList(*this->GetAlgorithm(), m_twoDCaloHitListName, m_pTwoDCaloHitList) != STATUS_CODE_SUCCESS) || !m_pTwoDCaloHitList)
    {
        std::cerr << "CommonMCNtupleTool: Failed to get list of 2D calo hits '" << m_twoDCaloHitListName << "'" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }
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

    if (pMCParticle)
    {
        // Calculate the true energy components
        const CartesianVector zDirection(0.f, 0.f, 1.f);
        const float           longitudinalEnergy = pMCParticle->GetMomentum().GetDotProduct(zDirection);
        const float           tranverseEnergy    = pMCParticle->GetMomentum().GetCrossProduct(zDirection).GetMagnitude();

        // Calculate the visible analogues of the energy/direction parameters
        CartesianVector visibleMomentum(0.f, 0.f, 0.f);
        this->RecursivelySumVisibleMomentum(pMCParticle, visibleMomentum);

        const float visibleEnergy             = visibleMomentum.GetMagnitude();
        const float visibleLongitudinalEnergy = visibleMomentum.GetDotProduct(zDirection);
        const float visibleTranverseEnergy    = visibleMomentum.GetCrossProduct(zDirection).GetMagnitude();

        const CartesianVector visibleInitialDirection =
            (visibleEnergy > std::numeric_limits<float>::epsilon()) ? visibleMomentum.GetUnitVector() : CartesianVector(0.f, 0.f, 0.f);

        records.emplace_back("mc_LongitudinalEnergy", static_cast<LArNtupleRecord::RFloat>(longitudinalEnergy));
        records.emplace_back("mc_TransverseEnergy", static_cast<LArNtupleRecord::RFloat>(tranverseEnergy));
        records.emplace_back("mc_VisibleEnergy", static_cast<LArNtupleRecord::RFloat>(visibleEnergy));
        records.emplace_back("mc_VisibleLongitudinalEnergy", static_cast<LArNtupleRecord::RFloat>(visibleLongitudinalEnergy));
        records.emplace_back("mc_VisibleTranverseEnergy", static_cast<LArNtupleRecord::RFloat>(visibleTranverseEnergy));
        records.emplace_back("mc_VisibleInitialDirectionX", static_cast<LArNtupleRecord::RFloat>(visibleInitialDirection.GetX()));
        records.emplace_back("mc_VisibleInitialDirectionY", static_cast<LArNtupleRecord::RFloat>(visibleInitialDirection.GetY()));
        records.emplace_back("mc_VisibleInitialDirectionZ", static_cast<LArNtupleRecord::RFloat>(visibleInitialDirection.GetZ()));
    }

    else // null values for size consistency
    {
        records.emplace_back("mc_LongitudinalEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_TransverseEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VisibleEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VisibleLongitudinalEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VisibleTranverseEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VisibleInitialDirectionX", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VisibleInitialDirectionY", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VisibleInitialDirectionZ", static_cast<LArNtupleRecord::RFloat>(0.f));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProcessPrimary(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;
    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    std::vector<LArNtupleRecord> nonNuPfoMCRecords = this->ProduceNonNeutrinoPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);
    records.insert(records.end(), std::make_move_iterator(nonNuPfoMCRecords.begin()), std::make_move_iterator(nonNuPfoMCRecords.end()));

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProcessCosmicRay(const ParticleFlowObject *const pPfo, const PfoList &pfoList,
    const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList)
{
    std::vector<LArNtupleRecord> records;
    std::vector<LArNtupleRecord> genericPfoMCRecords = this->ProduceGenericPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);
    records.insert(records.end(), std::make_move_iterator(genericPfoMCRecords.begin()), std::make_move_iterator(genericPfoMCRecords.end()));

    std::vector<LArNtupleRecord> nonNuPfoMCRecords = this->ProduceNonNeutrinoPfoMCRecords(pPfo, pfoList, pMCParticle, pMCParticleList);
    records.insert(records.end(), std::make_move_iterator(nonNuPfoMCRecords.begin()), std::make_move_iterator(nonNuPfoMCRecords.end()));

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
        const CartesianVector &vertexPosition    = pMCParticle->GetVertex();
        const CartesianVector &momentum          = pMCParticle->GetMomentum();
        const float            momentumMagnitude = momentum.GetMagnitude();
        const CartesianVector  initialDirection =
            (momentumMagnitude > std::numeric_limits<float>::epsilon()) ? momentum.GetUnitVector() : CartesianVector(0.f, 0.f, 0.f);

        records.emplace_back("mc_McParticleUid", reinterpret_cast<LArNtupleRecord::RULong64>(pMCParticle->GetUid()));
        records.emplace_back("mc_Energy", static_cast<LArNtupleRecord::RFloat>(pMCParticle->GetEnergy()));
        records.emplace_back("mc_VertexX", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetX()));
        records.emplace_back("mc_VertexY", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetY()));
        records.emplace_back("mc_VertexZ", static_cast<LArNtupleRecord::RFloat>(vertexPosition.GetZ()));
        records.emplace_back("mc_IsVertexFiducial", static_cast<LArNtupleRecord::RBool>(this->IsPointFiducial(vertexPosition)));
        records.emplace_back("mc_PdgCode", static_cast<LArNtupleRecord::RInt>(pMCParticle->GetParticleId()));
        records.emplace_back("mc_NuanceCode", static_cast<LArNtupleRecord::RUInt>(LArMCParticleHelper::GetNuanceCode(pMCParticle)));
        records.emplace_back("mc_MomentumX", static_cast<LArNtupleRecord::RFloat>(momentum.GetX()));
        records.emplace_back("mc_MomentumY", static_cast<LArNtupleRecord::RFloat>(momentum.GetY()));
        records.emplace_back("mc_MomentumZ", static_cast<LArNtupleRecord::RFloat>(momentum.GetZ()));
        records.emplace_back("mc_DirectionCosineX", static_cast<LArNtupleRecord::RFloat>(initialDirection.GetX()));
        records.emplace_back("mc_DirectionCosineY", static_cast<LArNtupleRecord::RFloat>(initialDirection.GetY()));
        records.emplace_back("mc_DirectionCosineZ", static_cast<LArNtupleRecord::RFloat>(initialDirection.GetZ()));
        records.emplace_back(
            "mc_EnergyWeightedContainedPfoFraction", static_cast<LArNtupleRecord::RFloat>(this->CalculateContainmentFraction(pMCParticle)));
    }

    else // null values for size consistency
    {
        records.emplace_back("mc_McParticleUid", reinterpret_cast<LArNtupleRecord::RULong64>(0ULL));
        records.emplace_back("mc_Energy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VertexX", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VertexY", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_VertexZ", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_IsVertexFiducial", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_PdgCode", static_cast<LArNtupleRecord::RInt>(0));
        records.emplace_back("mc_NuanceCode", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("mc_MomentumX", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_MomentumY", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_MomentumZ", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_DirectionCosineX", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_DirectionCosineY", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_DirectionCosineZ", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_EnergyWeightedContainedPfoFraction", static_cast<LArNtupleRecord::RFloat>(0.f));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> CommonMCNtupleTool::ProduceNonNeutrinoPfoMCRecords(
    const ParticleFlowObject *const pPfo, const PfoList &, const MCParticle *const pMCParticle, const MCParticleList *const) const
{
    std::vector<LArNtupleRecord> records;

    if (pMCParticle)
    {
        const bool isShower = LArAnalysisHelper::IsTrueShower(pMCParticle);

        if (pPfo)
        {
            const auto [matchPurity, matchCompleteness]             = this->GetRecoMcMatchQuality(pPfo, pMCParticle);
            const auto [wPlaneMatchPurity, wPlaneMatchCompleteness] = this->GetCollectionPlaneRecoMcMatchQuality(pPfo, pMCParticle);

            records.emplace_back("mc_RecoMcMatchCompleteness", static_cast<LArNtupleRecord::RFloat>(matchCompleteness));
            records.emplace_back("mc_RecoMcMatchPurity", static_cast<LArNtupleRecord::RFloat>(matchPurity));
            records.emplace_back("mc_RecoMcMatchCollectionPlaneCompleteness", static_cast<LArNtupleRecord::RFloat>(wPlaneMatchPurity));
            records.emplace_back("mc_RecoMcMatchCollectionPlanePurity", static_cast<LArNtupleRecord::RFloat>(wPlaneMatchCompleteness));
        }

        else
        {
            records.emplace_back("mc_RecoMcMatchCompleteness", static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("mc_RecoMcMatchPurity", static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("mc_RecoMcMatchCollectionPlaneCompleteness", static_cast<LArNtupleRecord::RFloat>(0.f));
            records.emplace_back("mc_RecoMcMatchCollectionPlanePurity", static_cast<LArNtupleRecord::RFloat>(0.f));
        }

        records.emplace_back("mc_IsShower", static_cast<LArNtupleRecord::RBool>(isShower));
        records.emplace_back("mc_IsTrack", static_cast<LArNtupleRecord::RBool>(!isShower));
        records.emplace_back("mc_KineticEnergy", static_cast<LArNtupleRecord::RFloat>(LArAnalysisHelper::GetTrueKineticEnergy(pMCParticle)));
        records.emplace_back("mc_Mass", static_cast<LArNtupleRecord::RFloat>(LArAnalysisHelper::GetTrueMass(pMCParticle)));
        records.emplace_back("mc_Momentum", static_cast<LArNtupleRecord::RFloat>(pMCParticle->GetMomentum().GetMagnitude()));
        records.emplace_back("mc_IsPrimary", static_cast<LArNtupleRecord::RBool>(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle)));
        records.emplace_back("mc_IsCosmicRay", static_cast<LArNtupleRecord::RBool>(LArMCParticleHelper::IsCosmicRay(pMCParticle)));
    }

    else // null values for size consistency
    {
        records.emplace_back("mc_RecoMcMatchCompleteness", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_RecoMcMatchPurity", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_RecoMcMatchCompleteness", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_RecoMcMatchPurity", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_IsShower", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_IsTrack", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_KineticEnergy", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_Mass", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_Momentum", static_cast<LArNtupleRecord::RFloat>(0.f));
        records.emplace_back("mc_IsPrimary", static_cast<LArNtupleRecord::RBool>(false));
        records.emplace_back("mc_IsCosmicRay", static_cast<LArNtupleRecord::RBool>(false));
    }

    return records;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CommonMCNtupleTool::RecursivelyGetContainedEnergy(const MCParticle *const pMCParticle, float &containedEnergy, float &totalEnergy) const
{
    if (LArMCParticleHelper::IsVisible(pMCParticle))
    {
        const float kineticEnergy = LArAnalysisHelper::GetTrueKineticEnergy(pMCParticle);
        totalEnergy += kineticEnergy;

        // Approximate strict containment as both the vertex and endpoint being fiducial
        const bool isVertexFiducial   = this->IsPointFiducial(pMCParticle->GetVertex());
        const bool isEndpointFiducial = this->IsPointFiducial(pMCParticle->GetEndpoint());

        if (isVertexFiducial && isEndpointFiducial)
            containedEnergy += kineticEnergy;
    }

    for (const MCParticle *const pDaughter : pMCParticle->GetDaughterList())
        this->RecursivelyGetContainedEnergy(pDaughter, containedEnergy, totalEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CommonMCNtupleTool::RecursivelySumVisibleMomentum(const MCParticle *const pMCParticle, CartesianVector &visibleMomentum) const
{
    if (LArMCParticleHelper::IsVisible(pMCParticle))
        visibleMomentum += pMCParticle->GetMomentum();

    for (const MCParticle *const pDaughter : pMCParticle->GetDaughterList())
        this->RecursivelySumVisibleMomentum(pDaughter, visibleMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float> CommonMCNtupleTool::GetRecoMcMatchQuality(const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) const
{
    return this->GetRecoMcMatchQualityImpl(pMCParticle, [&]() { return this->GetAllDownstreamTwoDHits(pPfo); },
        [](const CaloHit *const pCaloHit) {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                case TPC_VIEW_V:
                case TPC_VIEW_W:
                    return true;
                default:
                    break;
            }

            return false;
        });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float> CommonMCNtupleTool::GetCollectionPlaneRecoMcMatchQuality(
    const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle) const
{
    return this->GetRecoMcMatchQualityImpl(pMCParticle, [&]() { return this->GetAllDownstreamWHits(pPfo); },
        [](const CaloHit *const pCaloHit) { return pCaloHit->GetHitType() == TPC_VIEW_W; });
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float> CommonMCNtupleTool::GetRecoMcMatchQualityImpl(
    const MCParticle *const pMCParticle, const HitGetter &pfoHitGetter, const HitSelector &hitSelector) const
{
    const float hitsAssociatedWithMcParticle = this->GetHitWeightAssociatedWithMcParticle(pMCParticle, hitSelector);
    float       hitsAssociatedWithPfo(0.f), hitsAssociatedWithBoth(0.f);

    for (const CaloHit *const pCaloHit : pfoHitGetter())
    {
        const MCParticleWeightMap &weightMap = pCaloHit->GetMCParticleWeightMap();

        for (const auto &entry : weightMap)
            hitsAssociatedWithPfo += entry.second;

        const auto findIter = weightMap.find(pMCParticle);

        if (findIter != weightMap.end())
            hitsAssociatedWithBoth += findIter->second;
    }

    // Calculate the purity and completeness
    float purity(0.f), completeness(0.f);

    if (hitsAssociatedWithPfo > std::numeric_limits<float>::epsilon())
        purity = hitsAssociatedWithBoth / hitsAssociatedWithPfo;

    if (hitsAssociatedWithMcParticle > std::numeric_limits<float>::epsilon())
        completeness = hitsAssociatedWithBoth / hitsAssociatedWithMcParticle;

    return {purity, completeness};
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CommonMCNtupleTool::GetHitWeightAssociatedWithMcParticle(const MCParticle *const pMCParticle, const HitSelector &hitSelector) const
{
    if (!m_pTwoDCaloHitList)
    {
        std::cerr << "CommonMCNtupleTool: Could not get reco/MC match quality because there was no hit list" << std::endl;
        throw STATUS_CODE_NOT_FOUND;
    }

    float hitsAssociatedWithMcParticle(0.f);

    for (const CaloHit *const pCaloHit : *m_pTwoDCaloHitList)
    {
        if (!hitSelector(pCaloHit))
            continue;

        const MCParticleWeightMap &weightMap = pCaloHit->GetMCParticleWeightMap();
        const auto                 findIter  = weightMap.find(pMCParticle);

        if (findIter != weightMap.end())
            hitsAssociatedWithMcParticle += findIter->second;
    }

    return hitsAssociatedWithMcParticle;
}

} // namespace lar_physics_content
