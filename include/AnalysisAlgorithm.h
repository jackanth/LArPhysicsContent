/**
 *  @file LArPhysicsContent/include/AnalysisAlgorithm.h
 *
 *  @brief Header file for the LEE analysis algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LEE_ANALYSIS_ALGORITHM_H
#define LAR_LEE_ANALYSIS_ALGORITHM_H 1

#include "Objects/ParticleFlowObject.h"

#include "larpandoracontent/LArCustomParticles/CustomParticleCreationAlgorithm.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "LArAnalysisParticleHelper.h"
#include "LArAnalysisParticle.h"
#include "DebugDefinitions.h"
#include "LArTrackHitEnergy.h"
#include "BirksHitSelectionTool.h"

#include "TNtuple.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

/**
 *  @brief AnalysisAlgorithm class.
 * 
 */
class AnalysisAlgorithm : public CustomParticleCreationAlgorithm
{
public:
    /**
     *  @brief Constructor.
     * 
     */
    AnalysisAlgorithm();
    
    /**
     *  @brief Destructor.
     * 
     */
    ~AnalysisAlgorithm();

protected:
    StatusCode ReadSettings(const TiXmlHandle xmlHandle) override;
    void CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject*& pOutputPfo) const override;

private:
    class EnergyFromRangeData
    {
    public:
        EnergyFromRangeData(const float rangeMin, const float rangeMax, const float energy) noexcept :
            m_rangeMin{rangeMin}, m_rangeMax{rangeMax}, m_energy{energy}
        {}
        
        float m_rangeMin;
        float m_rangeMax;
        float m_energy;
    };
    
    using EnergyFromRangeDataVector = std::vector<EnergyFromRangeData>;
    
    std::string   m_mcParticleListName;
    float         m_fiducialCutLowXMargin;                 ///< 
    float         m_fiducialCutHighXMargin;                 ///< 
    float         m_fiducialCutLowYMargin;                 ///< 
    float         m_fiducialCutHighYMargin;                 ///< 
    float         m_fiducialCutLowZMargin;                 ///< 
    float         m_fiducialCutHighZMargin;                 ///< 
    unsigned int  m_trackSlidingFitWindow;              ///< 
    TFile        *m_pLArAnalysisParticleFile;              ///<
    TNtuple      *m_pLArAnalysisParticleNtuple;            ///< 
    std::string   m_parametersFile;                     ///<
    std::string   m_birksFitNtupleName;                 ///<
    std::string   m_protonEnergyFromRangeNtupleName;    ///<
    std::string   m_pionMuonEnergyFromRangeNtupleName;  ///<
    float         m_birksFitAlpha;                      ///<
    float         m_birksFitBeta;                       ///< 
    float         m_birksFitPole;                       ///< 
    float         m_birksSelectionMaxdEdX;              ///<
    mutable int   m_uniquePlotIdentifier;               ///< 
    bool          m_addMcInformation;                   ///< 
    EnergyFromRangeDataVector m_protonEnergyFromRangeDataVector;
    EnergyFromRangeDataVector m_pionMuonEnergyFromRangeDataVector;
    TMVA::Reader *m_pTmvaReader;
    mutable float m_tmvaTrackLength;
    mutable float m_tmvaAvgEnergyDeposition;
    CartesianVector m_minCoordinates;
    CartesianVector m_maxCoordinates;
    BirksHitSelectionTool *m_pBirksHitSelectionTool;
    float m_mcContainmentFractionLowerBound;
    std::string m_caloHitListName; 
    std::string m_tmvaWeights;

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
    void RecursivelyAppendParticleTypeMap(const ParticleFlowObject *const pPfo, LArAnalysisParticle::PfoTypeMap &pfoTypeMap,
                                          const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap,
                                          const LArAnalysisParticleHelper::TrackFitMap &trackFitMap) const;
    
    /**
     *  @brief ...
     * 
     */
    LArAnalysisParticle::TYPE EstimateParticleType(const ParticleFlowObject *const pPfo,
                                                const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap,
                                                const LArAnalysisParticleHelper::TrackFitMap &trackFitMap) const;
    
    /**
     *  @brief ...
     * 
     */
    void EstimateParticleEnergy(const ParticleFlowObject *const pPfo, const LArAnalysisParticle::PfoTypeMap &typeMap,
                                 const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
                                 const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap, float &particleEnergy,
    float &energySourcedFromRange, float &energySourcedFromShowerCharge, 
    float &energySourcedFromTrackCharge, float &energySourcedFromCorrectedTrackCharge) const;
    
    /**
     *  @brief ...
     * 
     */
    float EstimateShowerEnergy(const ParticleFlowObject *const pPfo) const;
    
    /**
     *  @brief ...
     * 
     */
    float EstimateEnergyFromCharge(const ParticleFlowObject *const pPfo) const;
    
    /**
     *  @brief ...
     * 
     */
    float CaloHitToScaledEnergy(const CaloHit *const pCaloHit) const;
        
    /**
     *  @brief ...
     * 
     */                                          
    float RecombinationCorrectEnergy(const float scaledEnergy, const float threeDDistance) const;
    
    /**
     *  @brief ...
     * 
     */
    float EstimateTrackEnergyFromCharge(const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;
    
    /**
     *  @brief ...
     * 
     */
    void DecideWhichTrackHitsToCorrect(LArTrackHitEnergy::Vector &trackHitEnergyVector) const;
                                           
    /**
     *  @brief ...
     * 
     */
    float EstimateTrackEnergyFromRange(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit,
                                       const LArAnalysisParticle::TYPE particleType, const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;
    
    /**
     *  @brief ...
     * 
     */
    LArAnalysisParticle::TypeTree CreateTypeTree(const ParticleFlowObject *const pPfo, const LArAnalysisParticle::PfoTypeMap &typeMap) const;
    
    /**
     *  @brief ...
     * 
     */
    CartesianVector GetCosmicRayDirectionAtVertex(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap, 
                                         const Vertex *const pVertex) const;
                                         
    /**
     *  @brief ...
     * 
     */
    CartesianVector GetPrimaryDirectionAtVertex(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap, 
                                         const Vertex *const pVertex) const;
    
    /**
     *  @brief ...
     * 
     */
    bool GetMcInformation(const ParticleFlowObject *const pPfo, float &mcEnergy, float &mcKineticEnergy, float &mcMass, LArAnalysisParticle::TypeTree &typeTree, 
        LArAnalysisParticle::TYPE &mcType, CartesianVector &mcVertexPosition, CartesianVector &mcMomentum, int &mcPdgCode, const bool isNeutrino,
        float &mcContainmentFraction, const MCParticle * &pMcMainMCParticle, float &mcHitPurity, float &mcHitCompleteness,
        float &mcCollectionPlaneHitPurity, float &mcCollectionPlaneHitCompleteness) const;
    
    /**
     *  @brief ...
     * 
     */
    void CalculateHitPurityAndCompleteness(const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle, 
        const CaloHitList *const pCaloHitList, const bool isNeutrino, float &hitPurity, float &hitCompleteness, float &mcCollectionPlaneHitPurity,
        float &mcCollectionPlaneHitCompleteness) const;
    
    void CalculateHitPurityAndCompleteness(const CaloHitList &pfoAssociatedCaloHits, const MCParticle *const pMCParticle, 
        const CaloHitList *const pCaloHitList, const bool isNeutrino, float &hitPurity, float &hitCompleteness, 
        const bool useCollectionPlaneOnly) const;
    
    /**
     *  @brief ...
     * 
     */
    void CountNumberOfDownstreamParticles(const ParticleFlowObject *const pPfo, unsigned &numberOfParticles) const;
    
    /**
     *  @brief ...
     * 
     */
    void WriteLArAnalysisParticleToFile(const LArAnalysisParticle *const pLArAnalysisParticle) const;
};
} // namespace lar_physics_content

#endif // #ifndef LAR_LEE_ANALYSIS_ALGORITHM_H
