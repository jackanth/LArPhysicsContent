/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h
 *
 *  @brief  Header file for the analysis ntuple algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
#define LAR_ANALYSIS_NTUPLE_ALGORITHM_H 1

#include "Objects/ParticleFlowObject.h"
#include "Pandora/Algorithm.h"
#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

#include <functional>
#include <memory>

namespace lar_physics_content
{
/**
 *  @brief  AnalysisNtupleAlgorithm class
 */
class AnalysisNtupleAlgorithm : public Algorithm
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
    StatusCode ReadSettings(const TiXmlHandle xmlHandle);
    StatusCode Run();

private:
    using VectorRecordProcessor =
        std::function<std::vector<LArNtupleRecord>(const pandora::ParticleFlowObject *const)>; ///< Alias for vector record processing function

    std::string m_pfoListName;                                   ///< The PFO list name
    std::string m_mcParticleListName;                            ///< The MC particle list name
    std::vector<NtupleVariableBaseTool *> m_ntupleVariableTools; ///< The ntuple variable tools
    std::string m_ntupleOutputFile;                              ///< The ntuple ROOT tree output file
    std::string m_ntupleTreeName;                                ///< The ntuple ROOT tree name
    std::string m_ntupleTreeTitle;                               ///< The ntuple ROOT tree title
    std::shared_ptr<LArNtuple> m_spNtuple;                       ///< Shared pointer to the ntuple
    int m_fileIdentifier;                                        ///< The input file identifier
    int m_eventNumber;                                           ///< The event number
    bool m_appendNtuple;                                         ///< Whether to append to an existing ntuple

    /**
     *  @brief  Get the particle lists.
     *
     *  @param  pfoList the list of all PFOs
     *
     *  @return the neutrinos, the list of cosmics and the list of primaries
     */
    std::tuple<PfoList, PfoList, PfoList> GetParticleLists(const PfoList &pfoList) const;

    /**
     *  @brief  Register the vector records for a given processor
     *
     *  @param  particles the list of all PFOs
     *  @param  processor the PFO processor
     */
    void RegisterVectorRecords(const PfoList &particles, const VectorRecordProcessor &processor) const;

    /**
     *  @brief  Register the ntuple records
     *
     *  @param  neutrinos the neutrinos
     *  @param  particles the list of all PFOs
     *  @param  cosmicRays the cosmic rays
     *  @param  primaries the primaries
     *  @param  pfoList the list of all PFOs
     *  @param  optional address of the MC particle list
     */
    void RegisterNtupleRecords(const PfoList &neutrinos, const PfoList &cosmicRays, const PfoList &primaries, const PfoList &pfoList,
        const MCParticleList *const pMCParticleList) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
