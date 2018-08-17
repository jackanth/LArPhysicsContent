/**
 *  @file   larphysicscontent/LArNtuple/NtupleVariableBaseTool.h
 *
 *  @brief  Header file for the ntuple variable base tool class.
 *
 *  $Log: $
 */
#ifndef LAR_NTUPLE_VARIABLE_BASE_TOOL_H
#define LAR_NTUPLE_VARIABLE_BASE_TOOL_H 1

#include "larphysicscontent/LArNtuple/LArNtupleRecord.h"

#include "Api/PandoraContentApi.h"
#include "Objects/ParticleFlowObject.h"
#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include <functional>
#include <memory>

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the AnalysisNtupleAlgorithm class
 */
class AnalysisNtupleAlgorithm;

/**
 *  @brief  Forward declaration of the LArNtuple class
 */
class LArNtuple;

/**
 *  @brief  NtupleVariableBaseTool class
 */
class NtupleVariableBaseTool : public pandora::AlgorithmTool
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
     * @brief  Default virtual destructor
     */
    virtual ~NtupleVariableBaseTool() = default;

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle);

    /**
     *  @brief  Process an event - to be overriden
     *
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the event records
     */
    virtual std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process any particle (including non-primary daughters) - to be overriden
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional pointer to the corresponding MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the particle records
     */
    virtual std::vector<LArNtupleRecord> ProcessParticle(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a neutrino - to be overriden
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional pointer to the corresponding MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the neutrino records
     */
    virtual std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a primary neutrino daughter - to be overriden
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional pointer to the corresponding MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the primary records
     */
    virtual std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a cosmic ray - to be overriden
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional pointer to the corresponding MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the cosmic ray records
     */
    virtual std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get all the downstream 3D hits of a PFO (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamThreeDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream 2D hits of a PFO (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    pandora::CaloHitList GetAllDownstreamTwoDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream U hits of a PFO (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamUHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream V hits of a PFO (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamVHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream W hits of a PFO (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamWHits(const pandora::ParticleFlowObject *const pPfo) const;

    friend class AnalysisNtupleAlgorithm;

private:
    using Processor = std::function<std::vector<LArNtupleRecord>()>; ///< Alias for a function to process a PFO

    std::shared_ptr<LArNtuple> m_spNtuple; ///< Shared pointer to the ntuple

    /**
     *  @brief  Process an event (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the event records
     */
    std::vector<LArNtupleRecord> ProcessEventWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::PfoList &pfoList,
        const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process any particle (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional address of the MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the particle records
     */
    std::vector<LArNtupleRecord> ProcessParticleWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::PfoList &pfoList, const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a neutrino (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional address of the MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the neutrino records
     */
    std::vector<LArNtupleRecord> ProcessNeutrinoWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::PfoList &pfoList, const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a primary neutrino daughter (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional addresss of the MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the primary records
     */
    std::vector<LArNtupleRecord> ProcessPrimaryWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::PfoList &pfoList, const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Process a cosmic ray (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional addresss of the MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the cosmic ray records
     */
    std::vector<LArNtupleRecord> ProcessCosmicRayWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::PfoList &pfoList, const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Implementation of PFO processing wrapper
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  prefix the prefix to apply to branch names
     *  @param  processor the PFO processor method
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProcessImpl(const AnalysisNtupleAlgorithm *const pAlgorithm, const std::string &prefix, const Processor &processor);

    /**
     *  @brief  Set the ntuple shared pointer
     *
     *  @param  spNtuple shared pointer to the ntuple
     */
    void SetNtuple(std::shared_ptr<LArNtuple> spNtuple) noexcept;

    /**
     *  @brief  Get an MC particle
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    const pandora::MCParticle *GetMCParticle(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode NtupleVariableBaseTool::ReadSettings(const pandora::TiXmlHandle)
{
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessEvent(const pandora::PfoList &, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessParticle(const pandora::ParticleFlowObject *const,
    const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessNeutrino(const pandora::ParticleFlowObject *const,
    const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessPrimary(const pandora::ParticleFlowObject *const,
    const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessCosmicRay(const pandora::ParticleFlowObject *const,
    const pandora::PfoList &, const pandora::MCParticle *const, const pandora::MCParticleList *const)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NtupleVariableBaseTool::SetNtuple(std::shared_ptr<LArNtuple> spNtuple) noexcept
{
    m_spNtuple = std::move_if_noexcept(spNtuple);
}
} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_VARIABLE_BASE_TOOL_H
