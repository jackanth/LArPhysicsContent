/**
 *  @file   larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.cc
 *
 *  @brief  Implementation of the energy estimator ntuple tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.h"

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
EnergyEstimatorNtupleTool::EnergyEstimatorNtupleTool() :
    NtupleVariableBaseTool(),
    m_writeEnergiesToNtuple(true),
    m_useParticleId(true),
    m_trainingSetMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyEstimatorNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteEnergiesToNtuple", m_writeEnergiesToNtuple));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseParticleId", m_useParticleId));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingSetMode", m_trainingSetMode));

    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessEvent(const PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessNeutrino(
    const ParticleFlowObject *const, const PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const pMcParticle, const MCParticleList *const mcParticleList)
{
    if (m_trainingSetMode)
        return this->ProduceTrainingRecords(pPfo, pfoList, pMcParticle, mcParticleList);

    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const pPfo, const PfoList &pfoList, const MCParticle *const pMcParticle, const MCParticleList *const mcParticleList)
{
    if (m_trainingSetMode)
        return this->ProduceTrainingRecords(pPfo, pfoList, pMcParticle, mcParticleList);

    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::CaloHitToThreeDDistance(const CaloHit *const pCaloHit, const ThreeDSlidingFitResult &trackFit) const
{
    if (pCaloHit->GetHitType() != TPC_VIEW_W)
    {
        std::cerr << "EnergyEstimatorNtupleTool: Can only calculate CaloHit 3D distance for W-view hits" << std::endl;
        throw STATUS_CODE_NOT_ALLOWED;
    }

    const CartesianVector fitDirection = LArAnalysisHelper::GetFittedDirectionAtPosition(trackFit, pCaloHit->GetPositionVector());
    const float           wirePitch    = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W);
    return this->CellToThreeDDistance(pCaloHit->GetCellSize1(), wirePitch, fitDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyEstimatorNtupleTool::CellToThreeDDistance(const float hitWidth, const float wirePitch, const CartesianVector &fitDirection) const
{
    float polarAngle = 0.f, azimuthalAngle = 0.f;
    std::tie(polarAngle, azimuthalAngle) = this->GetPolarAnglesFromDirection(fitDirection);

    if (polarAngle <= std::numeric_limits<float>::epsilon()) // negligible polar angle, so the wire pitch is the hit separation
        return wirePitch;

    const float cosPhi_sinTheta = std::fabs(std::cos(azimuthalAngle) * std::sin(polarAngle));
    const float sinPhi_sinTheta = std::fabs(std::sin(azimuthalAngle) * std::sin(polarAngle));

    float dx_p = std::numeric_limits<float>::max();
    float dx_w = std::numeric_limits<float>::max();

    if (cosPhi_sinTheta > std::numeric_limits<float>::epsilon())
        dx_p = wirePitch / cosPhi_sinTheta;

    if (sinPhi_sinTheta > std::numeric_limits<float>::epsilon())
        dx_w = hitWidth / sinPhi_sinTheta;

    if ((dx_p < std::numeric_limits<float>::max()) || (dx_w < std::numeric_limits<float>::max()))
        return std::min(dx_p, dx_w);

    // For non-negligible polar angle, this shouldn't happen
    std::cerr << "EnergyEstimatorNtupleTool: Failed to calculate hit's 3D distance" << std::endl;
    throw STATUS_CODE_FAILURE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<float, float> EnergyEstimatorNtupleTool::GetPolarAnglesFromDirection(const CartesianVector &direction) const
{
    const float polarAngle    = std::acos(std::fabs(direction.GetY()));
    const float sinPolarAngle = std::sin(polarAngle);

    if (sinPolarAngle <= std::numeric_limits<float>::epsilon())
        return {0.f, 0.f}; // negligible polar angle means azimuthal angle is undefined

    const float azimuthalAngle = std::asin(std::fabs(direction.GetX() / sinPolarAngle));
    return {polarAngle, azimuthalAngle};
}


//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> EnergyEstimatorNtupleTool::ProduceTrainingRecords(
    const ParticleFlowObject *const pPfo, const PfoList &, const MCParticle *const pMcParticle, const MCParticleList *const)
{
    std::vector<LArNtupleRecord> records;

    if (pPfo && pMcParticle)
    {
        LArNtupleRecord::RFloatVector threeDDistances, integratedAdcCounts;
        LArNtupleRecord::RFloat       showerAdcCount(0.f);
        LArNtupleRecord::RUInt        numVectorEntries(0U);

        for (const ParticleFlowObject *const pDownstreamPfo : this->GetAllDownstreamPfos(pPfo))
        {
            CaloHitList collectionPlaneHits;
            LArPfoHelper::GetCaloHits(pDownstreamPfo, TPC_VIEW_W, collectionPlaneHits);

            if (LArPfoHelper::IsShower(pDownstreamPfo))
            {
                for (const CaloHit *const pCaloHit : collectionPlaneHits)
                    showerAdcCount += pCaloHit->GetInputEnergy();
            }

            else
            {
                const LArNtupleHelper::TrackFitSharedPtr &spTrackFit = this->GetTrackFit(pDownstreamPfo);

                if (!spTrackFit)
                    continue;

                for (const CaloHit *const pCaloHit : collectionPlaneHits)
                {
                    try // failure to calculate a hit's 3D distance shouldn't be the end of the world
                    {
                        const float threeDDistance = this->CaloHitToThreeDDistance(pCaloHit, *spTrackFit);
                        threeDDistances.push_back(threeDDistance);
                        integratedAdcCounts.push_back(pCaloHit->GetInputEnergy());
                        ++numVectorEntries;
                    }

                    catch (...)
                    {
                        continue;
                    }
                }
            }
        }

        records.emplace_back("ThreeDDistances", threeDDistances);
        records.emplace_back("IntegratedAdcCounts", integratedAdcCounts);
        records.emplace_back("NumVectorEntries", numVectorEntries);
        records.emplace_back("ShowerAdcCount", showerAdcCount);
    }

    else
    {
        records.emplace_back("ThreeDDistances", LArNtupleRecord::RFloatVector());
        records.emplace_back("IntegratedAdcCounts", LArNtupleRecord::RFloatVector());
        records.emplace_back("NumVectorEntries", static_cast<LArNtupleRecord::RUInt>(0U));
        records.emplace_back("ShowerAdcCount", static_cast<LArNtupleRecord::RFloat>(0.f));
    }

    return records;
}

} // namespace lar_physics_content
