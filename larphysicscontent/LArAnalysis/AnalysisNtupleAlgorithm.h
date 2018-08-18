/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h
 *
 *  @brief  Header file for the analysis ntuple algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
#define LAR_ANALYSIS_NTUPLE_ALGORITHM_H 1

#include "larphysicscontent/LArNtuple/LArNtuple.h"
#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

#include "Objects/ParticleFlowObject.h"
#include "Pandora/Algorithm.h"

#include <functional>
#include <memory>

namespace lar_physics_content
{
/**
 *  @brief  AnalysisNtupleAlgorithm class
 */
class AnalysisNtupleAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
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
    using MCParticleRetriever =
        std::function<const pandora::MCParticle *(const pandora::ParticleFlowObject *const pPfo)>; ///< Alias for an MCParticle retriever
    using VectorRecordProcessor = std::function<std::vector<LArNtupleRecord>(NtupleVariableBaseTool *const pNtupleTool,
        const pandora::ParticleFlowObject *const, const pandora::MCParticle *const)>; ///< Alias for vector record processing function

    std::string                           m_pfoListName;                  ///< The PFO list name
    std::string                           m_mcParticleListName;           ///< The MC particle list name
    std::vector<NtupleVariableBaseTool *> m_ntupleVariableTools;          ///< The ntuple variable tools
    std::string                           m_ntupleOutputFile;             ///< The ntuple ROOT tree output file
    std::string                           m_ntupleTreeName;               ///< The ntuple ROOT tree name
    std::string                           m_ntupleTreeTitle;              ///< The ntuple ROOT tree title
    std::shared_ptr<LArNtuple>            m_spNtuple;                     ///< Shared pointer to the ntuple
    int                                   m_fileIdentifier;               ///< The input file identifier
    int                                   m_eventNumber;                  ///< The event number
    bool                                  m_appendNtuple;                 ///< Whether to append to an existing ntuple
    float                                 m_minUnmatchedMcParticleEnergy; ///< The minimum unmatched MCParticle energy that passes the cut

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
     *  @param  pMCParticleList optional address of the MC particle list
     *
     *  @return the list of MC neutrinos, the list of MC cosmics, and the list of MC primaries
     */
    std::tuple<pandora::MCParticleList, pandora::MCParticleList, pandora::MCParticleList> GetMCParticleLists(
        const pandora::MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Register the vector records for a given processor
     *
     *  @param  particles the list of all PFOs
     *  @param  mcParticleList the list of all relevant MC particles
     *  @param  pMCParticleList optional address of the list of all MC particles
     *  @param  type the vector type
     *  @param  processor the PFO processor
     *  @param  mcParticleRetriever the MCParticle retriever
     *
     *  @return the size of the vector records registered
     */
    std::size_t RegisterVectorRecords(const pandora::PfoList &particles, const pandora::MCParticleList &mcParticleList,
        const pandora::MCParticleList *const pMCParticleList, const LArNtupleHelper::VECTOR_BRANCH_TYPE type,
        const VectorRecordProcessor &processor, const MCParticleRetriever &mcParticleRetriever) const;

    /**
     *  @brief  Register the ntuple records
     *
     *  @param  neutrinos the neutrinos
     *  @param  cosmicRays the cosmic rays
     *  @param  primaries the primaries
     *  @param  pfoList the list of all PFOs
     *  @param  mcNeutrinos the MC neutrinos
     *  @param  mcCosmicRays the MC cosmic rays
     *  @param  mcPrimaries the MC primaries
     *  @param  pMCParticleList the optional address of the MC particle list
     */
    void RegisterNtupleRecords(const pandora::PfoList &neutrinos, const pandora::PfoList &cosmicRays, const pandora::PfoList &primaries,
        const pandora::PfoList &pfoList, const pandora::MCParticleList &mcNeutrinos, const pandora::MCParticleList &mcCosmicRays,
        const pandora::MCParticleList &mcPrimaries, const pandora::MCParticleList *const pMCParticleList) const;

    /**
     *  @brief  Test the quality of an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *
     *  @return whether it passes the quality checks
     */
    bool TestMCParticleQuality(const pandora::MCParticle *const pMCParticle) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool AnalysisNtupleAlgorithm::TestMCParticleQuality(const pandora::MCParticle *const pMCParticle) const
{
    return pMCParticle->GetEnergy() > m_minUnmatchedMcParticleEnergy;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
