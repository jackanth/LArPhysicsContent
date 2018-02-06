/**
 *  @file LArPhysicsContent/include/AnalysisDataAlgorithm.h
 *
 *  @brief Header file for the LEE analysis data algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LEE_ANALYSIS_DATA_ALGORITHM_H
#define LAR_LEE_ANALYSIS_DATA_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Objects/ParticleFlowObject.h"
#include "Objects/MCParticle.h"

#include "BirksHitSelectionTool.h"

#include "LArAnalysisParticleHelper.h"
#include "TNtuple.h"
#include "TTree.h"

#include <unordered_map>
#include <unordered_set>

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief AnalysisDataAlgorithm class.
 * 
 */
class AnalysisDataAlgorithm : public Algorithm
{
public:
    /**
     *  @brief Constructor.
     * 
     */
    AnalysisDataAlgorithm();
    
    /**
     *  @brief Destructor.
     * 
     */
    ~AnalysisDataAlgorithm();

protected:
    StatusCode ReadSettings(const TiXmlHandle xmlHandle);
    StatusCode Run();
    
private:
    using MCParticleMap = std::unordered_map<const ParticleFlowObject *, const MCParticle *>; /// Alias for a map to MC particles.
    using PdgCodeSet    = std::unordered_set<int>;                                            ///< ...

    bool          m_produceBirksFitData;                 ///< 
    bool          m_produceEnergyFromRangeData;          ///< 
    bool          m_producePidData;                      ///< 
    std::string   m_pfoListName;                         ///<
    float         m_fiducialCutLowXMargin;                 ///< 
    float         m_fiducialCutHighXMargin;                 ///< 
    float         m_fiducialCutLowYMargin;                 ///< 
    float         m_fiducialCutHighYMargin;                 ///< 
    float         m_fiducialCutLowZMargin;                 ///< 
    float         m_fiducialCutHighZMargin;                 ///< 
    unsigned int  m_trackSlidingFitWindow;               ///< 
    std::string   m_rootDataFileName;                    ///<
    TFile        *m_pRootDataFile;                       ///<
    TTree        *m_pBirksFitDataTree;                   ///<
    TNtuple      *m_pEnergyFromRangeProtonDataNtuple;    ///<
    TNtuple      *m_pEnergyFromRangePionMuonDataNtuple;  ///<
    std::string   m_pPidDataProtonsTTreeName;            ///<
    std::string   m_pPidDataMuonsPionsTTreeName;         ///<
    mutable int   m_uniquePlotIdentifier;                ///<
    BirksHitSelectionTool *m_pBirksHitSelectionTool;
    std::string m_caloHitListName;
    
    /**
     *  @brief ...
     * 
     */
    void RecursivelyAppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo, LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap, 
                                            const LArAnalysisParticleHelper::TrackFitMap &trackFitMap) const;
    
    /**
     *  @brief ...
     * 
     */
    LArTrackHitEnergy::Vector AppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit) const;
    
    
    /**
     *  @brief ...
     * 
     */
    void RecursivelyAppendMCParticleMap(const ParticleFlowObject *const pPfo, MCParticleMap &mcParticleMap) const;
    
    /**
     *  @brief ...
     * 
     */
    bool ProduceBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                             const MCParticleMap &mcParticleMap, const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap) const;
    
    /**
     *  @brief ...
     * 
     */
    bool GetBirksFitData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                         const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap, float &totalNoBirksAdcIntegral,
                         FloatVector &birksAdcIntegrals, FloatVector &threeDDistances) const;
    
    /**
     *  @brief ...
     * 
     */
    float AddUpShowerAdcs(const ParticleFlowObject *const pPfo) const;
    
    /**
     *  @brief ...
     * 
     */
    void GetTrackAdcsAndDistances(const LArTrackHitEnergy::Vector &trackHitEnergies, float &totalNoBirksAdcIntegral,
                                  FloatVector &birksAdcIntegrals, FloatVector &threeDDistances) const;
    
    /**
     *  @brief ...
     * 
     */
    void RecursivelyProduceEnergyFromRangeData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                               const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap, TNtuple *const pNtuple,
                                               const PdgCodeSet &pdgCodeSet) const;
    
    /**
     *  @brief ...
     * 
     */
    float GetTrueEnergy(const MCParticle *pMCParticle) const;
       
    /**
     *  @brief ...
     * 
     */
    void RecursivelyProducePidData(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                   const AnalysisDataAlgorithm::MCParticleMap &mcParticleMap, const bool isCosmicRay) const;
};
} // namespace lar_physics_content

#endif // #ifndef LAR_LEE_ANALYSIS_DATA_ALGORITHM_H
