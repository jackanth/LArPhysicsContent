/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.cc
 *
 *  @brief  Implementation of the analysis ntuple algorithm class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h"
#include "larphysicscontent/LArHelpers/LArAnalysisNtupleHelper.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{

AnalysisNtupleAlgorithm::AnalysisNtupleAlgorithm() : m_pfoListName()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AnalysisNtupleAlgorithm::~AnalysisNtupleAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisNtupleAlgorithm::Run()
{
    // Get input Pfo List
    const PfoList *pPfoList(NULL);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_pfoListName, pPfoList) || !pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "AnalysisNtupleAlgorithm: cannot find pfo list " << m_pfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    std::cout << "The PFO list name is " << m_pfoListName << " with size " << pPfoList->size() << std::endl;

    const PfoList &pfoList = *pPfoList;

    for (const ParticleFlowObject *const pPfo : pfoList)
        this->ProcessPfo(pPfo, pfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisNtupleAlgorithm::ProcessPfo(const ParticleFlowObject *const pPfo, const PfoList &pfoList) const
{
    (void) pfoList;
    
    // Work out what kind of PFO we're dealing with.
    const LArAnalysisNtupleHelper::PARTICLE_CLASS particleClass = LArAnalysisNtupleHelper::GetParticleClass(pPfo);
    std::cout << "Particle type = " << LArAnalysisNtupleHelper::ToString(particleClass) << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AnalysisNtupleAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    return STATUS_CODE_SUCCESS;
}
} // namespace lar_physics_content
