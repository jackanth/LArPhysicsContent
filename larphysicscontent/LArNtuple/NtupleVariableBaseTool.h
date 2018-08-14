/**
 *  @file   larphysicscontent/LArNtuple/NtupleVariableBaseTool.h
 *
 *  @brief  Header file for the ntuple variable base tool class.
 *
 *  $Log: $
 */
#ifndef LAR_NTUPLE_VARIABLE_BASE_TOOL_H
#define LAR_NTUPLE_VARIABLE_BASE_TOOL_H 1

#include "Api/PandoraContentApi.h"
#include "Objects/ParticleFlowObject.h"
#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"
#include "larphysicscontent/LArNtuple/LArNtupleRecord.h"

#include <functional>

using namespace pandora;

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the AnalysisNtupleAlgorithm class
 */
class AnalysisNtupleAlgorithm;

/**
 *  @brief  NtupleVariableBaseTool class
 */
class NtupleVariableBaseTool : public AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    NtupleVariableBaseTool() noexcept;

    /**
     * @brief  Default copy constructor
     */
    NtupleVariableBaseTool(const NtupleVariableBaseTool &) = default;

    /**
     * @brief  Default move constructor
     */
    NtupleVariableBaseTool(NtupleVariableBaseTool &&) = default;

    /**
     * @brief  Default copy assignment operator
     */
    NtupleVariableBaseTool &operator=(const NtupleVariableBaseTool &) = default;

    /**
     * @brief  Default move assignment operator
     */
    NtupleVariableBaseTool &operator=(NtupleVariableBaseTool &&) = default;

    /**
     * @brief  Default destructor
     */
    ~NtupleVariableBaseTool() = default;

protected:
    StatusCode ReadSettings(const TiXmlHandle);

    /**
     *  @brief  Process an event
     *
     *  @param  pfoList The list of all PFOs
     */
    virtual std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process any particle (including non-primary daughters)
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    virtual std::vector<LArNtupleRecord> ProcessParticle(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a neutrino
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    virtual std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a primary neutrino daughter
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    virtual std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a cosmic ray
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    virtual std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const MCParticle *const pMCParticle, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process an event (wrapper method)
     *
     *  @param  pfoList The list of all PFOs
     */
    virtual std::vector<LArNtupleRecord> ProcessEventWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
        const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get an MC particle
     *
     *  @param  pfoList The list of all PFOs
     */
    virtual const MCParticle * GetMCParticle(const pandora::ParticleFlowObject *const pPfo,
    const MCParticleList *const pMCParticleList);

    friend class AnalysisNtupleAlgorithm;

private:
    using Processor = std::function<std::vector<LArNtupleRecord>(const MCParticle *const)>; ///< Alias for a function to process a PFO

    /**
     *  @brief  Process any particle (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    std::vector<LArNtupleRecord> ProcessParticleWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a neutrino (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    std::vector<LArNtupleRecord> ProcessNeutrinoWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a primary neutrino daughter (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    std::vector<LArNtupleRecord> ProcessPrimaryWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a cosmic ray (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    std::vector<LArNtupleRecord> ProcessCosmicRayWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const MCParticleList *const pMCParticleList);

    /**
     *  @brief  Implementation of PFO processing wrapper
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOs
     */
    std::vector<LArNtupleRecord> ProcessImpl(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const MCParticleList *const pMCParticleList, const std::string &prefix, const Processor &processor);


    MCParticleWeightMap GetMCParticleWeightMap(const CaloHitList &caloHitList) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline StatusCode NtupleVariableBaseTool::ReadSettings(const TiXmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessEvent(const pandora::PfoList &, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessParticle(
    const pandora::ParticleFlowObject *const, const pandora::PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessNeutrino(
    const pandora::ParticleFlowObject *const, const pandora::PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessPrimary(
    const pandora::ParticleFlowObject *const, const pandora::PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessCosmicRay(
    const pandora::ParticleFlowObject *const, const pandora::PfoList &, const MCParticle *const, const MCParticleList *const)
{
    return {};
}
} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_VARIABLE_BASE_TOOL_H
