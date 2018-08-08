/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisNtupleAlgorithm.h
 *
 *  @brief  Header file for the write analysis ntuple algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
#define LAR_ANALYSIS_NTUPLE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Objects/ParticleFlowObject.h"

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief  AnalysisNtupleAlgorithm class
 *
 */
class AnalysisNtupleAlgorithm : public Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    AnalysisNtupleAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~AnalysisNtupleAlgorithm();

protected:
    StatusCode ReadSettings(const TiXmlHandle xmlHandle);
    StatusCode Run();

private:
    std::string m_pfoListName; ///< The PFO list name.

    /**
     *  @brief  Process a PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList The list of all PFOS
     */
    void ProcessPfo(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_NTUPLE_ALGORITHM_H
