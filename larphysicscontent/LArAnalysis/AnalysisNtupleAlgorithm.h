/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h
 *
 *  @brief  Header file for the analysis ntuple algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
#define LAR_ANALYSIS_NTUPLE_ALGORITHM_H 1

#include "larphysicscontent/LArAnalysis/EventValidationTool.h"
#include "larphysicscontent/LArNtuple/LArNtuple.h"
#include "larphysicscontent/LArObjects/LArRootRegistry.h"

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

namespace lar_physics_content
{
/**
 *  @brief  AnalysisNtupleAlgorithm class
 */
class AnalysisNtupleAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Constructor
     */
    AnalysisNtupleAlgorithm();

    /**
     *  @brief  Default copy constructor
     */
    AnalysisNtupleAlgorithm(const AnalysisNtupleAlgorithm &) = default;

    /**
     *  @brief  Default move constructor
     */
    AnalysisNtupleAlgorithm(AnalysisNtupleAlgorithm &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    AnalysisNtupleAlgorithm &operator=(const AnalysisNtupleAlgorithm &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    AnalysisNtupleAlgorithm &operator=(AnalysisNtupleAlgorithm &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~AnalysisNtupleAlgorithm() = default;

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode Run();

private:
    using PfoHypothesisMap = std::unordered_map<unsigned int, pandora::PfoVector>;    ///< Alias for a map from ID to PFO hypothesis
    using McTargetVector   = std::vector<std::shared_ptr<LArMCTargetValidationInfo>>; ///< Alias for a vector of MC target shared pointers
    using McInteractionVector = std::vector<std::shared_ptr<LArInteractionValidationInfo>>; ///< Alias for a vector of MC interaction shared pointers

    template <typename T>
    using VectorRecordProcessor = std::function<std::vector<LArNtupleRecord>(NtupleVariableBaseTool *const,
        const pandora::ParticleFlowObject *const, const std::shared_ptr<std::decay_t<T>> &)>; ///< Alias for a vector record processor

    unsigned int                          m_eventNumber;            ///< The current event number
    EventValidationTool *                 m_pEventValidationTool;   ///< Address of the event validation tool
    std::string                           m_caloHitListName;        ///< The CaloHit list name
    std::string                           m_mcParticleListName;     ///< The MCParticle list name
    bool                                  m_printValidation;        ///< Whether to print the validation
    bool                                  m_produceAllOutcomes;     ///< Whether to produce all outcomes
    std::string                           m_pfoListName;            ///< If not all outcomes, the PFO list to use
    std::string                           m_ntupleOutputFile;       ///< The ntuple ROOT tree output file
    std::string                           m_ntupleTreeName;         ///< The ntuple ROOT tree name
    std::string                           m_ntupleTreeTitle;        ///< The ntuple ROOT tree title
    std::string                           m_plotsOutputFile;        ///< The plots ROOT output file
    std::string                           m_tmpOutputFile;          ///< The tmp ROOT output file
    std::shared_ptr<LArNtuple>            m_spNtuple;               ///< Shared pointer to the ntuple
    int                                   m_fileIdentifier;         ///< The input file identifier
    bool                                  m_appendNtuple;           ///< Whether to append to an existing ntuple
    pandora::CartesianVector              m_fiducialCutLowMargins;  ///< The low-coordinate margins for the fiducial cut
    pandora::CartesianVector              m_fiducialCutHighMargins; ///< The high-coordinate margins for the fiducial cut
    pandora::CartesianVector              m_minFiducialCoordinates; ///< The minimum fiducial coordinates
    pandora::CartesianVector              m_maxFiducialCoordinates; ///< The maximum fiducial coordinates
    std::shared_ptr<LArRootRegistry>      m_spTmpRegistry;          ///< Shared pointer to the tmp ROOT registry
    std::shared_ptr<LArRootRegistry>      m_spPlotsRegistry;        ///< Shared pointer to the plots ROOT registry
    std::vector<NtupleVariableBaseTool *> m_ntupleVariableTools;    ///< The ntuple variable tools

    /**
     *  @brief  Collect all possible PFO outcomes
     *
     *  @param  clearCosmics the clear cosmics (to populate)
     *  @param  pfoHypotheses the PFO hypotheses (to populate)
     */
    void CollectAllPfoOutcomes(pandora::PfoVector &clearCosmics, PfoHypothesisMap &pfoHypotheses) const;

    /**
     *  @brief  Get the Pandora workers
     *
     *  @return addresses of the slicing worker, the slice nu worker, and the slice CR worker
     */
    std::tuple<const pandora::Pandora *, const pandora::Pandora *, const pandora::Pandora *> GetPandoraWorkers() const;

    /**
     *  @brief  Build the PFO hypotheses
     *
     *  @param  pfoHypotheses the PFO hypotheses (to populate)
     *  @param  pSlicingWorker address of the slicing worker
     *  @param  pSliceNuWorker address of the slice nu worker
     *  @param  pSliceCRWorker address of the slice CR worker
     */
    void BuildHypotheses(PfoHypothesisMap &pfoHypotheses, const pandora::Pandora *pSlicingWorker, const pandora::Pandora *pSliceNuWorker,
        const pandora::Pandora *pSliceCRWorker) const;

    /**
     *  @brief  Collect all connected PFOs into a vector
     *
     *  @param  parentPfoList the parent PFO list
     *  @param  pfoVector the PFO vector (to populate)
     */
    void CollectPfos(const pandora::PfoList &parentPfoList, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Process an event hypothesis
     *
     *  @param  hypothesisId the hypothesis ID
     *  @param  allPfos the list of all PFOs
     *  @param  caloHitList the CaloHit list
     *  @param  pMCParticleList address of the MC particle list
     */
    void ProcessEventHypothesis(const int hypothesisId, const pandora::PfoList &allPfos, const pandora::CaloHitList &caloHitList,
        const pandora::MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Print the validation info
     *
     *  @param  eventValidationInfo the vector of shared pointers to LArInteractionValidationInfo objects
     */
    void PrintValidation(const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const;

    /**
     *  @brief  Get the particle lists
     *
     *  @param  pfoList the list of all PFOs
     *
     *  @return the list of neutrinos, the list of cosmics, and the list of primaries
     */
    std::tuple<pandora::PfoList, pandora::PfoList, pandora::PfoList> GetParticleLists(const pandora::PfoList &pfoList) const;

    /**
     *  @brief  Get the MC particle lists
     *
     *  @param  eventValidationInfo the vector of LArInteractionValidationInfo shared pointers
     *
     *  @return the list of MC neutrino interactions, the list of MC cosmic targets, and the list of MC primary targets
     */
    std::tuple<McInteractionVector, McTargetVector, McTargetVector> GetMCParticleLists(
        const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const;

    /**
     *  @brief  Get the map from PFOs to MC targets
     *
     *  @param  eventValidationInfo the vector of LArInteractionValidationInfo shared pointers
     *
     *  @return the map from PFOs to best MC targets and the map from neutrino PFOs to interactions
     */
    std::tuple<LArAnalysisHelper::PfoToTargetMap, LArAnalysisHelper::PfoToInteractionMap> GetPfoToMcObjectMaps(
        const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const;

    /**
     *  @brief  Register the vector records for a given processor
     *
     *  @param  particles the list of all PFOs
     *  @param  pfoToMcObjectMap the map from PFOs to MC objects
     *  @param  allMcObjects the list of all MC objects
     *  @param  type the vector type
     *  @param  processor the PFO processor
     *
     *  @return the size of the vector records registered
     */
    template <typename T>
    std::size_t RegisterVectorRecords(const pandora::PfoList &                                           particles,
        const std::unordered_map<const pandora::ParticleFlowObject *, std::shared_ptr<std::decay_t<T>>> &pfoToMcObjectMap,
        const std::vector<std::shared_ptr<std::decay_t<T>>> &allMcObjects, const LArNtupleHelper::VECTOR_BRANCH_TYPE type,
        const VectorRecordProcessor<std::decay_t<T>> &processor) const;

    /**
     *  @brief  Register the ntuple records
     *
     *  @param  neutrinos the neutrinos
     *  @param  cosmicRays the cosmic rays
     *  @param  primaries the primaries
     *  @param  pfoList the list of all PFOs
     *  @param  mcNeutrinoInts the MC neutrino interactions
     *  @param  mcCosmicRayTargets the MC cosmic ray targets
     *  @param  mcPrimaryTargets the MC primary targets
     *  @param  pfoToTargetMap the map from PFOs to MC targets
     *  @param  pfoToInteractionMap the map from PFOs to MC interactions
     *  @param  eventValidationInfo the vector of LArInteractionValidationInfo shared pointers
     */
    void RegisterNtupleRecords(const int hypothesisId, const pandora::PfoList &neutrinos, const pandora::PfoList &cosmicRays,
        const pandora::PfoList &primaries, const pandora::PfoList &pfoList, const McInteractionVector &mcNeutrinoInts,
        const McTargetVector &mcCosmicRayTargets, const McTargetVector &mcPrimaryTargets,
        const LArAnalysisHelper::PfoToTargetMap &pfoToTargetMap, const LArAnalysisHelper::PfoToInteractionMap &pfoToInteractionMap,
        const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
std::size_t AnalysisNtupleAlgorithm::RegisterVectorRecords(const pandora::PfoList &                  particles,
    const std::unordered_map<const pandora::ParticleFlowObject *, std::shared_ptr<std::decay_t<T>>> &pfoToMcObjectMap,
    const std::vector<std::shared_ptr<std::decay_t<T>>> &allMcObjects, const LArNtupleHelper::VECTOR_BRANCH_TYPE type,
    const VectorRecordProcessor<std::decay_t<T>> &processor) const
{
    using T_D        = std::decay_t<T>;
    using TSharedPtr = std::shared_ptr<T_D>;

    // Run over the reco PFOs, matching to MC particles where possible
    std::unordered_set<TSharedPtr> encounteredMcObjects;

    for (const pandora::ParticleFlowObject *const pPfo : particles)
    {
        const auto        findIter   = pfoToMcObjectMap.find(pPfo);
        const TSharedPtr &spMcObject = (findIter == pfoToMcObjectMap.end()) ? nullptr : findIter->second;

        for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
        {
            for (const LArNtupleRecord &record : processor(pNtupleTool, pPfo, spMcObject))
                m_spNtuple->AddVectorRecordElement(record, type);

            if (spMcObject)
                encounteredMcObjects.insert(spMcObject);
        }

        m_spNtuple->FillVectors(type);
    }

    // Find all the MC particles that are in our main classes but not matched to a PFO
    std::size_t extraMCRecords(0UL);

    for (const TSharedPtr &spMcObject : allMcObjects)
    {
        if (encounteredMcObjects.find(spMcObject) != encounteredMcObjects.end())
            continue;

        for (NtupleVariableBaseTool *const pNtupleTool : m_ntupleVariableTools)
        {
            for (const LArNtupleRecord &record : processor(pNtupleTool, nullptr, spMcObject))
                m_spNtuple->AddVectorRecordElement(record, type);
        }

        ++extraMCRecords;
        m_spNtuple->FillVectors(type);
    }

    m_spNtuple->PushVectors(type);
    return particles.size() + extraMCRecords;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
