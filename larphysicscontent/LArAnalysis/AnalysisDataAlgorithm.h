/**
 *  @file   larphysicscontent/LArAnalysis/AnalysisDataAlgorithm.h
 *
 *  @brief  Header file for the lar analysis data algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_DATA_ALGORITHM_H
#define LAR_ANALYSIS_DATA_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Objects/MCParticle.h"
#include "Objects/ParticleFlowObject.h"

#include "larphysicscontent/LArAnalysis/HitPurityTool.h"
#include "larphysicscontent/LArAnalysis/McInfoTool.h"
#include "larphysicscontent/LArAnalysis/TrackHitEnergyTool.h"

#include "larphysicscontent/LArHelpers/LArAnalysisParticleHelper.h"

#include "TNtuple.h"
#include "TTree.h"

#include <unordered_map>
#include <unordered_set>

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief  AnalysisDataAlgorithm class
 */
class AnalysisDataAlgorithm : public Algorithm
{
public:
    /**
     *  @brief Default onstructor
     */
    AnalysisDataAlgorithm();

    /**
     *  @brief Destructor
     */
    ~AnalysisDataAlgorithm();

protected:
    StatusCode ReadSettings(const TiXmlHandle xmlHandle);
    StatusCode Run();

private:
    using MCParticleMap = std::unordered_map<const ParticleFlowObject *, const MCParticle *>; ///< Alias for a map from PFOs to MC particles
    using PdgCodeSet = std::unordered_set<int>;                                               ///< Alias for a set of PDG codes

    CartesianVector m_fiducialCutLowMargins;  ///< The low fiducial margins
    CartesianVector m_fiducialCutHighMargins; ///< The high fiducial cut margins
    CartesianVector m_minCoordinates;         ///< The detector's minimum fiducial coordinates
    CartesianVector m_maxCoordinates;         ///< The detector's maximum fiducial coordinates

    bool m_produceBirksFitData;                    ///< Whether to produce the Birks fit data
    bool m_produceEnergyFromRangeData;             ///< Whether to produce energy-from-range data
    bool m_producePidData;                         ///< Whether to produce proton PID data
    std::string m_pfoListName;                     ///< The PFO list name
    std::string m_rootDataFileName;                ///< The file path to save the fit data
    std::string m_caloHitListName;                 ///< The CaloHit list name
    std::string m_pidDataProtonsTTreeName;         ///< The name of the TTree to use for the proton PID data
    std::string m_pidDataMuonsPionsTTreeName;      ///< The name of the TTree to use for the pion/muon PID data
    TFile *m_pRootDataFile;                        ///< Address of the fit data TFile object
    TTree *m_pBirksFitDataTree;                    ///< Address of the Birks fit data TTree object
    TNtuple *m_pEnergyFromRangeProtonDataNtuple;   ///< Address of the proton energy-from-range data TNtuple object
    TNtuple *m_pEnergyFromRangePionMuonDataNtuple; ///< Address of the pion/muon energy-from-range data TNtuple object
    TrackHitEnergyTool *m_pTrackHitEnergyTool;     ///< Address of the track hit energy tool
    McInfoTool *m_pMcInfoTool;                     ///< Address of the MC info tool
    HitPurityTool *m_pHitPurityTool;               ///< Address of the hit purity tool

    /**
     *  @brief  Recurse through the PFO hierarchy and append the MC particle map
     *
     *  @param  pPfo address of the current PFO
     *  @param  mcParticleMap the MC particle map to append
     */
    void RecursivelyAppendMCParticleMap(const ParticleFlowObject *const pPfo, MCParticleMap &mcParticleMap) const;

    /**
     *  @brief  Produce Birks fit data for a given PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFitMap the track fit map
     *  @param  mcParticleMap the MC particle map
     *  @param  trackHitEnergyMap the track hit energy map
     *
     *  @return success
     */
    bool ProduceBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                             const MCParticleMap &mcParticleMap) const;

    /**
     *  @brief  Get the Birks fit data for a given PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFitMap the track fit map
     *  @param  trackHitEnergyMap the track hit energy map
     *  @param  totalNoBirksAdcIntegral the sum of the charges not to be corrected
     *  @param  birksAdcIntegrals the vector of charges to be corrected
     *  @param  threeDDistances the vector of 3D distances over which the charges were deposited
     *
     *  @return success
     */
    bool GetBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                         float &totalNoBirksAdcIntegral, FloatVector &birksAdcIntegrals, FloatVector &threeDDistances) const;

    /**
     *  @brief  Add up the shower ADCs for a PFO
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the total shower ADC value
     */
    float AddUpShowerAdcs(const ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get the track ADCs and distances for the hits of a given track-like PFO
     *
     *  @param  trackHitEnergies the track hit energy vector for the PFO
     *  @param  totalNoBirksAdcIntegral the total ADC integral not to be corrected for recombination (to be populated)
     *  @param  birksAdcIntegrals the vector of ADC integrals to be corrected for recombination (to be populated)
     *  @param  threeDDistances the vector of 3D distances over which the charges were deposited (to be populated)
     */
    void GetTrackAdcsAndDistances(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergies, float &totalNoBirksAdcIntegral,
                                  FloatVector &birksAdcIntegrals, FloatVector &threeDDistances) const;

    /**
     *  @brief  Recurse through the PFO hierarchy and produce energy from range data for a given set of PDG codes
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFitMap the track fit map
     *  @param  mcParticleMap the MC particle map
     *  @param  pNtuple address of the TNtuple object
     *  @param  pdgCodeSet The set of PDG codes
     */
    void RecursivelyProduceEnergyFromRangeData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                                               const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap, TNtuple *const pNtuple,
                                               const PdgCodeSet &pdgCodeSet) const;

    /**
     *  @brief  Get the kinetic energy of an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the kinetic energy
     */
    float GetTrueEnergy(const MCParticle *pMCParticle) const;

    /**
     *  @brief  Recurse through the PFO hierarchy and produce proton PID data
     *
     *  @param  pPfo address of the current PFO
     *  @param  trackFitMap the track fit map
     *  @param  mcParticleMap the MC particle map
     *  @param  isCosmicRay whether the particle is a cosmic ray
     */
    void RecursivelyProducePidData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
                                   const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap, const bool isCosmicRay) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_DATA_ALGORITHM_H
