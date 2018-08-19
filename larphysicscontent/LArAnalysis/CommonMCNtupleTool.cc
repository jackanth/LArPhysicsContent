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
CommonMCNtupleTool::CommonMCNtupleTool() :
    NtupleVariableBaseTool(),
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minFiducialCoordinates(0.f, 0.f, 0.f),
    m_maxFiducialCoordinates(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CommonMCNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
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
        records.emplace_back("mc_IsVertexFiducial", static_cast<LArNtupleRecord::RBool>(LArAnalysisHelper::IsPointFiducial(
                                                        vertexPosition, m_minFiducialCoordinates, m_maxFiducialCoordinates)));
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
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const pMCParticle, const MCParticleList *const) const
{
    std::vector<LArNtupleRecord> records;

    if (pMCParticle)
    {
        const bool isShower = LArAnalysisHelper::IsTrueShower(pMCParticle);

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
        const bool isVertexFiducial = LArAnalysisHelper::IsPointFiducial(pMCParticle->GetVertex(), m_minFiducialCoordinates, m_maxFiducialCoordinates);
        const bool isEndpointFiducial =
            LArAnalysisHelper::IsPointFiducial(pMCParticle->GetEndpoint(), m_minFiducialCoordinates, m_maxFiducialCoordinates);

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

} // namespace lar_physics_content
