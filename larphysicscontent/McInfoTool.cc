/**
 *  @file   larphysicscontent/McInfoTool.cc
 *
 *  @brief  Implementation of the MC info tool class.
 *
 *  $Log: $
 */

#include "larphysicscontent/McInfoTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{

McInfoTool::McInfoTool() :
    m_fiducialCutLowMargins(10.f, 20.f, 10.f),
    m_fiducialCutHighMargins(10.f, 20.f, 10.f),
    m_minCoordinates(0.f, 0.f, 0.f),
    m_maxCoordinates(0.f, 0.f, 0.f),
    m_mcContainmentFractionLowerBound(0.9f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool McInfoTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const pMCParticle, LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    pfoMcInfo = this->GetMcInformation(pMCParticle);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticleHelper::PfoMcInfo McInfoTool::GetMcInformation(const MCParticle *const pMCParticle) const
{
    if (!pMCParticle)
    {
        std::cout << "LArAnalysisParticleHelper: could not get MC information because there was no MC particle" << std::endl;
        throw STATUS_CODE_INVALID_PARAMETER;
    }

    LArAnalysisParticleHelper::PfoMcInfo pfoMcInfo;
    pfoMcInfo.m_pMCParticle = pMCParticle;
    this->CreateMcTypeTree(pMCParticle, pfoMcInfo.m_mcTypeTree);

    pfoMcInfo.m_mcEnergy         = pMCParticle->GetEnergy();
    pfoMcInfo.m_mcKineticEnergy  = pMCParticle->GetEnergy() - PdgTable::GetParticleMass(pMCParticle->GetParticleId());
    pfoMcInfo.m_mcMass           = PdgTable::GetParticleMass(pMCParticle->GetParticleId());
    pfoMcInfo.m_mcType           = pfoMcInfo.m_mcTypeTree.Type();
    pfoMcInfo.m_mcVertexPosition = pMCParticle->GetVertex();
    pfoMcInfo.m_mcMomentum       = pMCParticle->GetMomentum();

    if (pfoMcInfo.m_mcMomentum.GetMagnitude() < std::numeric_limits<float>::epsilon())
            std::cout << "LArAnalysisParticleHelper: could not get direction from MC momentum as it was too small" << std::endl;
    else
        pfoMcInfo.m_mcDirectionCosines = pfoMcInfo.m_mcMomentum.GetUnitVector();

    pfoMcInfo.m_mcPdgCode = pMCParticle->GetParticleId();

    float escapedEnergy(0.f), totalEnergy(0.f);
    this->RecursivelyAddEscapedEnergy(pMCParticle, escapedEnergy, totalEnergy);

    if (totalEnergy > std::numeric_limits<float>::epsilon())
        pfoMcInfo.m_mcContainmentFraction = (totalEnergy - escapedEnergy) / totalEnergy;

    else
        pfoMcInfo.m_mcContainmentFraction = 0.f;

    pfoMcInfo.m_mcIsContained  = (pfoMcInfo.m_mcContainmentFraction >= m_mcContainmentFractionLowerBound);
    pfoMcInfo.m_mcIsShower     = (pfoMcInfo.m_mcType == LArAnalysisParticle::TYPE::SHOWER);
    pfoMcInfo.m_mcIsProton     = (pfoMcInfo.m_mcType == LArAnalysisParticle::TYPE::PROTON);
    pfoMcInfo.m_mcIsPionOrMuon = (pfoMcInfo.m_mcType == LArAnalysisParticle::TYPE::PION_MUON);
    pfoMcInfo.m_mcIsCosmicRay  = (pfoMcInfo.m_mcType == LArAnalysisParticle::TYPE::COSMIC_RAY);

    return pfoMcInfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void McInfoTool::RecursivelyAddEscapedEnergy(const MCParticle *const pCurrentMCParticle, float &escapedEnergy,
    float &totalEnergy) const
{
    bool allEnergyContained(true);

    switch (pCurrentMCParticle->GetParticleId())
    {
        case E_MINUS:
        case E_PLUS:
        case MU_MINUS:
        case MU_PLUS:
        case PI_PLUS:
        case PI_MINUS:
        case PROTON:
        {
            const float mcParticleEnergy = pCurrentMCParticle->GetEnergy() - PdgTable::GetParticleMass(pCurrentMCParticle->GetParticleId());
            const CartesianVector displacementVector = pCurrentMCParticle->GetEndpoint() - pCurrentMCParticle->GetVertex();

            if (displacementVector.GetMagnitude() < std::numeric_limits<float>::epsilon())
                break;

            float muMin(0.f), muMax(1.f);
            bool forceZeroContainment(false);

            this->AdjustMusForContainmentFraction(CartesianVector(m_minCoordinates.GetX(), 0.f, 0.f),
                CartesianVector(-1.f, 0.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);

            this->AdjustMusForContainmentFraction(CartesianVector(m_maxCoordinates.GetX(), 0.f, 0.f),
                CartesianVector(1.f, 0.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);

            this->AdjustMusForContainmentFraction(CartesianVector(0.f, m_minCoordinates.GetY(), 0.f),
                CartesianVector(0.f, -1.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);

            this->AdjustMusForContainmentFraction(CartesianVector(0.f, m_maxCoordinates.GetY(), 0.f),
                CartesianVector(0.f, 1.f, 0.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);

            this->AdjustMusForContainmentFraction(CartesianVector(0.f, 0.f, m_minCoordinates.GetZ()),
                CartesianVector(0.f, 0.f, -1.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);

            this->AdjustMusForContainmentFraction(CartesianVector(0.f, 0.f, m_maxCoordinates.GetZ()),
                CartesianVector(0.f, 0.f, 1.f), pCurrentMCParticle->GetVertex(), displacementVector, muMin, muMax, forceZeroContainment);

            if (forceZeroContainment)
            {
                escapedEnergy += mcParticleEnergy;
                allEnergyContained = false;
            }

            else
            {
                const float containmentFraction = std::max(0.f, muMax - muMin);

                if (containmentFraction < 1.f)
                {
                    escapedEnergy += (1.f - containmentFraction) * mcParticleEnergy;
                    allEnergyContained = false;
                }
            }

            totalEnergy += mcParticleEnergy;
        }

        default: break;
    }

    if (allEnergyContained)
    {
        for (const MCParticle *const pDaughterParticle : pCurrentMCParticle->GetDaughterList())
            this->RecursivelyAddEscapedEnergy(pDaughterParticle, escapedEnergy, totalEnergy);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void McInfoTool::AdjustMusForContainmentFraction(const CartesianVector &planePoint, const CartesianVector &planeNormal,
    const CartesianVector &vertexPosition, const CartesianVector &originalDisplacementVector, float &muMin, float &muMax,
    bool &forceZeroContainment) const
{
    const float projectedDisplacement = originalDisplacementVector.GetDotProduct(planeNormal);

    if (std::fabs(projectedDisplacement) < std::numeric_limits<float>::epsilon())
        return;

    const float muIntercept = (planePoint - vertexPosition).GetDotProduct(planeNormal) / projectedDisplacement;
    const bool isAligned = (projectedDisplacement > 0.f);

    if (isAligned)
    {
        if (muIntercept > std::numeric_limits<float>::min() && muIntercept < std::numeric_limits<float>::max())
            muMax = std::min(muMax, muIntercept);

        else if (muIntercept < 0.f)
            forceZeroContainment = true;
    }

    else
    {
        if (muIntercept > std::numeric_limits<float>::min() && muIntercept < std::numeric_limits<float>::max())
            muMin = std::max(muMin, muIntercept);

        else if (muIntercept > 0.f)
            forceZeroContainment = true;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool McInfoTool::CreateMcTypeTree(const MCParticle *const pMCParticle, LArAnalysisParticle::TypeTree &typeTree) const
{
    LArAnalysisParticle::TypeTree::List daughterTypeTrees;

    const LArAnalysisParticle::TYPE type = this->GetMcParticleType(pMCParticle);

    if (type == LArAnalysisParticle::TYPE::UNKNOWN)
        return false;

    if (type != LArAnalysisParticle::TYPE::SHOWER)
    {
        for (const MCParticle *const pDaughterParticle : pMCParticle->GetDaughterList())
        {
            if (pDaughterParticle->GetEnergy() > 0.05f)
            {
                LArAnalysisParticle::TypeTree daughterTypeTree;

                if (this->CreateMcTypeTree(pDaughterParticle, daughterTypeTree))
                    daughterTypeTrees.push_back(daughterTypeTree);
            }
        }
    }

    typeTree = LArAnalysisParticle::TypeTree(type, daughterTypeTrees);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArAnalysisParticle::TYPE McInfoTool::GetMcParticleType(const MCParticle *const pMCParticle) const
{
    if (LArMCParticleHelper::IsNeutrino(pMCParticle))
        return LArAnalysisParticle::TYPE::NEUTRINO;

    if (LArMCParticleHelper::IsCosmicRay(pMCParticle))
        return LArAnalysisParticle::TYPE::COSMIC_RAY;

    switch (pMCParticle->GetParticleId())
    {
        case PROTON:   return LArAnalysisParticle::TYPE::PROTON;
        case MU_MINUS:
        case MU_PLUS:
        case PI_MINUS:
        case PI_PLUS:  return LArAnalysisParticle::TYPE::PION_MUON;
        case PHOTON:
        case E_MINUS:
        case E_PLUS:   return LArAnalysisParticle::TYPE::SHOWER;
        case NEUTRON:
        default: break;
    }

    return LArAnalysisParticle::TYPE::UNKNOWN;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode McInfoTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutLowMargins", m_fiducialCutLowMargins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutHighMargins", m_fiducialCutHighMargins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "McContainmentFractionLowerBound",
        m_mcContainmentFractionLowerBound));

    LArAnalysisParticleHelper::GetFiducialCutParameters(this->GetPandora(), m_fiducialCutLowMargins, m_fiducialCutHighMargins, m_minCoordinates, m_maxCoordinates);
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
