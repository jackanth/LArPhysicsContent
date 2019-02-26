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

#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"
#include "larphysicscontent/LArHelpers/LArNtupleHelper.h"
#include "larphysicscontent/LArNtuple/LArBranchPlaceholder.h"
#include "larphysicscontent/LArObjects/LArInteractionValidationInfo.h"
#include "larphysicscontent/LArObjects/LArRootRegistry.h"

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
     *  @brief  Prepare an event - to be overriden
     *
     *  @param  pfoList the list of all PFOs
     *  @param  eventValidationInfo vector of shared pointers to the interaction validation info objects
     *
     *  @return the event records
     */
    virtual void PrepareEvent(const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo);

    /**
     *  @brief  Process an event - to be overriden
     *
     *  @param  pfoList the list of all PFOs
     *  @param  eventValidationInfo vector of shared pointers to the interaction validation info objects
     *
     *  @return the event records
     */
    virtual std::vector<LArNtupleRecord> ProcessEvent(
        const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo);

    /**
     *  @brief  Process a neutrino - to be overriden
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  spInteractionInfo shared plinter to the interaction validation info object
     *
     *  @return the neutrino records
     */
    virtual std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo);

    /**
     *  @brief  Process a primary neutrino daughter - to be overriden
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  spMcTarget shared pointer to the MC target
     *
     *  @return the primary records
     */
    virtual std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget);

    /**
     *  @brief  Process a cosmic ray - to be overriden
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  spMcTarget shared pointer to the MC target
     *
     *  @return the cosmic ray records
     */
    virtual std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget);

    /**
     *  @brief  Get all the downstream 3D hits of a PFO, including from the PFO itself (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamThreeDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream 2D hits of a PFO, including from the PFO itself (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    pandora::CaloHitList GetAllDownstreamTwoDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream U hits of a PFO, including from the PFO itself (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamUHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream V hits of a PFO, including from the PFO itself (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamVHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream W hits of a PFO, including from the PFO itself (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamWHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all the downstream PFOs of a PFO, including the PFO itself (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the downstream PFOs
     */
    const pandora::PfoList &GetAllDownstreamPfos(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get the track fit for a PFO (from the cache if possible)
     *
     *  @param  pPfo address of the PFO
     *
     *  @return shared pointer to the track fit, or nullptr if fit is not possible
     */
    const LArNtupleHelper::TrackFitSharedPtr &GetTrackFit(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Retrieve an event record
     *
     *  @param  branchName the unprefixed branch name
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetEventRecord(const std::string &branchName) const;

    /**
     *  @brief  Retrieve a particle record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pPfo address of the PFO
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetParticleRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Retrieve a particle record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetParticleRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Retrieve a primary record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pPfo address of the PFO
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetPrimaryRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Retrieve a primary record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetPrimaryRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Retrieve a cosmic record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pPfo address of the PFO
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetCosmicRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Retrieve a cosmic record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetCosmicRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Retrieve a neutrino record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pPfo address of the PFO
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetNeutrinoRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Retrieve a neutrino record
     *
     *  @param  branchName the unprefixed branch name
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the record value
     */
    template <typename T>
    const std::decay_t<T> &GetNeutrinoRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Check whether a point is fiducial
     *
     *  @param  point the point
     *
     *  @return whether it is fiducial
     */
    bool IsPointFiducial(const pandora::CartesianVector &point) const;

    /**
     *  @brief  Get the algorithm address
     *
     *  @return the algorithm address
     */
    const pandora::Algorithm *GetAlgorithm() const noexcept;

    /**
     *  @brief  Get the plots ROOT registry.
     *
     *  @return the plots ROOT registry
     */
    const std::shared_ptr<LArRootRegistry> &GetPlotsRegistry() const noexcept;

    /**
     *  @brief  Get the tmp ROOT registry.
     *
     *  @return the tmp ROOT registry
     */
    const std::shared_ptr<LArRootRegistry> &GetTmpRegistry() const noexcept;

    friend class AnalysisNtupleAlgorithm;

private:
    using Processor = std::function<std::vector<LArNtupleRecord>()>; ///< Alias for a function to process a PFO

    std::shared_ptr<LArNtuple>       m_spNtuple;                 ///< Shared pointer to the ntuple
    std::string                      m_eventPrefix;              ///< The event prefix
    std::string                      m_neutrinoPrefix;           ///< The neutrino prefix
    std::string                      m_primaryPrefix;            ///< The primary prefix
    std::string                      m_particlePrefix;           ///< The particle prefix
    std::string                      m_cosmicPrefix;             ///< The cosmic prefix
    pandora::CartesianVector         m_fiducialRegion1MinCoords; ///< The minimum fiducial coordinates of region 1
    pandora::CartesianVector         m_fiducialRegion1MaxCoords; ///< The maximum fiducial coordinates of region 1
    pandora::CartesianVector         m_fiducialRegion2MinCoords; ///< The minimum fiducial coordinates of region 2
    pandora::CartesianVector         m_fiducialRegion2MaxCoords; ///< The maximum fiducial coordinates of region 2
    const pandora::Algorithm *       m_pAlgorithm;               ///< The address of the calling algorithm
    std::shared_ptr<LArRootRegistry> m_spPlotsRegistry;          ///< The plots ROOT registry
    std::shared_ptr<LArRootRegistry> m_spTmpRegistry;            ///< The tmp ROOT registry
    bool                             m_isSetup;                  ///< Whether the tool has been set up.

    /**
     *  @brief  Prepare an event (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pfoList the list of all PFOs
     *  @param  eventValidationInfo vector of shared pointers to the interaction validation info objects
     */
    void PrepareEventWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::PfoList &pfoList,
        const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo);

    /**
     *  @brief  Process an event (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pfoList the list of all PFOs
     *  @param  eventValidationInfo vector of shared pointers to the interaction validation info objects
     *
     *  @return the event records
     */
    std::vector<LArNtupleRecord> ProcessEventWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::PfoList &pfoList,
        const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo);

    /**
     *  @brief  Process a neutrino (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  spInteractionInfo shared pointer to the interaction validation info object
     *
     *  @return the neutrino records
     */
    std::vector<LArNtupleRecord> ProcessNeutrinoWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::PfoList &pfoList, const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo);

    /**
     *  @brief  Process a primary neutrino daughter (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  spMcTarget shared pointer to the MC target
     *
     *  @return the primary records
     */
    std::vector<LArNtupleRecord> ProcessPrimaryWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget);

    /**
     *  @brief  Process a cosmic ray (wrapper method)
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  spMcTarget shared pointer to the MC target
     *
     *  @return the cosmic ray records
     */
    std::vector<LArNtupleRecord> ProcessCosmicRayWrapper(const AnalysisNtupleAlgorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget);

    /**
     *  @brief  Implementation of PFO processing wrapper
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  prefix the prefix to apply to branch names
     *  @param  pPfo optional address of the PFO
     *  @param  pMcParticle optional address of the MC particle
     *  @param  processor the PFO processor method
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProcessImpl(const AnalysisNtupleAlgorithm *const pAlgorithm, const std::string &prefix,
        const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMcParticle, const Processor &processor);

    /**
     *  @brief  Set the ntuple shared pointer
     *
     *  @param  spNtuple shared pointer to the ntuple
     *  @param  pAlgorithm the algorithm address
     *  @param  fiducialRegion1MinCoords the minimum fiducial coordinates of region 1
     *  @param  fiducialRegion1MaxCoords the maximum fiducial coordinates of region 1
     *  @param  fiducialRegion2MinCoords the minimum fiducial coordinates of region 2
     *  @param  fiducialRegion2MaxCoords the maximum fiducial coordinates of region 2
     *  @param  spPlotsRegistry shared pointer to the plots ROOT registry
     *  @param  spTmpRegistry shared pointer to the tmp ROOT registry
     */
    void Setup(std::shared_ptr<LArNtuple> spNtuple, const pandora::Algorithm *const pAlgorithm,
        pandora::CartesianVector fiducialRegion1MinCoords, pandora::CartesianVector fiducialRegion1MaxCoords,
        pandora::CartesianVector fiducialRegion2MinCoords, pandora::CartesianVector fiducialRegion2MaxCoords,
        std::shared_ptr<LArRootRegistry> spPlotsRegistry, std::shared_ptr<LArRootRegistry> spTmpRegistry);

    /**
     *  @brief  Get a scalar record
     *
     *  @param  branchName the branch name
     *
     *  @return shared pointer to the record
     */
    LArBranchPlaceholder::NtupleRecordSPtr GetScalarRecord(const std::string &branchName) const;

    /**
     *  @brief  Get a vector record element
     *
     *  @param  type the vector type
     *  @param  branchName the branch name
     *  @param  pPfo address of the PFO
     *
     *  @return shared pointer to the record
     */
    LArBranchPlaceholder::NtupleRecordSPtr GetVectorRecordElement(
        LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get a vector record element
     *
     *  @param  type the vector type
     *  @param  branchName the branch name
     *  @param  pMCParticle address of the MC particle
     *
     *  @return shared pointer to the record
     */
    LArBranchPlaceholder::NtupleRecordSPtr GetVectorRecordElement(
        LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const pandora::MCParticle *const pMCParticle) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode NtupleVariableBaseTool::ReadSettings(const pandora::TiXmlHandle)
{
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void NtupleVariableBaseTool::PrepareEvent(const pandora::PfoList &, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessEvent(
    const pandora::PfoList &, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessNeutrino(
    const pandora::ParticleFlowObject *const, const pandora::PfoList &, const std::shared_ptr<LArInteractionValidationInfo> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessPrimary(
    const pandora::ParticleFlowObject *const, const pandora::PfoList &, const std::shared_ptr<LArMCTargetValidationInfo> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<LArNtupleRecord> NtupleVariableBaseTool::ProcessCosmicRay(
    const pandora::ParticleFlowObject *const, const pandora::PfoList &, const std::shared_ptr<LArMCTargetValidationInfo> &)
{
    return {};
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetEventRecord(const std::string &branchName) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetScalarRecord(m_eventPrefix + branchName)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetParticleRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::PARTICLE, m_particlePrefix + branchName, pPfo)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetParticleRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::PARTICLE, m_particlePrefix + branchName, pMCParticle)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetPrimaryRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::PRIMARY, m_primaryPrefix + branchName, pPfo)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetPrimaryRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::PRIMARY, m_primaryPrefix + branchName, pMCParticle)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetCosmicRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::COSMIC_RAY, m_cosmicPrefix + branchName, pPfo)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetCosmicRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::COSMIC_RAY, m_cosmicPrefix + branchName, pMCParticle)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetNeutrinoRecord(const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::NEUTRINO, m_neutrinoPrefix + branchName, pPfo)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &NtupleVariableBaseTool::GetNeutrinoRecord(const std::string &branchName, const pandora::MCParticle *const pMCParticle) const
{
    if (!m_spNtuple)
    {
        std::cerr << "NtupleVariableBaseTool: Could not call ntuple method because no ntuple was set" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return this->GetVectorRecordElement(LArNtupleHelper::VECTOR_BRANCH_TYPE::NEUTRINO, m_neutrinoPrefix + branchName, pMCParticle)->Value<std::decay_t<T>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool NtupleVariableBaseTool::IsPointFiducial(const pandora::CartesianVector &point) const
{
    return LArAnalysisHelper::IsPointFiducial(point, m_fiducialRegion1MinCoords, m_fiducialRegion1MaxCoords) ||
           LArAnalysisHelper::IsPointFiducial(point, m_fiducialRegion2MinCoords, m_fiducialRegion2MaxCoords);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Algorithm *NtupleVariableBaseTool::GetAlgorithm() const noexcept
{
    return m_pAlgorithm;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::shared_ptr<LArRootRegistry> &NtupleVariableBaseTool::GetPlotsRegistry() const noexcept
{
    return m_spPlotsRegistry;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::shared_ptr<LArRootRegistry> &NtupleVariableBaseTool::GetTmpRegistry() const noexcept
{
    return m_spTmpRegistry;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_VARIABLE_BASE_TOOL_H
