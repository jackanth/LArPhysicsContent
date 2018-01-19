/**
 *  @file   LArPhysicsContent/include/LArAnalysisParticleHelper.h
 *
 *  @brief  Header file for the LEE analysis helper class.
 *
 *  $Log: $
 */
#ifndef LEE_ANALYSIS_HELPER_H
#define LEE_ANALYSIS_HELPER_H 1

#include "Pandora/Pandora.h"
#include "Objects/ParticleFlowObject.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "TNtuple.h"
#include "LArTrackHitEnergy.h"

#include <tuple>
#include <unordered_map>

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{
/**
 *  @brief LArAnalysisParticleHelper class.
 * 
 */
class LArAnalysisParticleHelper
{
public:
    using LArTrackHitEnergyMap   = std::unordered_map<const ParticleFlowObject *, LArTrackHitEnergy::Vector>; ///< 
    using HitProjectionPair   = std::pair<const CaloHit *, float>;
    using HitProjectionVector = std::vector<HitProjectionPair>;
    using TrackFitMap = std::unordered_map<const ParticleFlowObject *, ThreeDSlidingFitResult>; ///< Alias for map from PFOs to track fits.

    /**
     *  @brief  ...
     * 
     */
    static void GetFiducialCutParameters(const Pandora &pandoraInstance, const float fiducialCutXMargin, const float fiducialCutYMargin, const float fiducialCutZMargin,
        CartesianVector &minCoordinates, CartesianVector &maxCoordinates);
    
    /**
     *  @brief  ...
     * 
     */
    static void RecursivelyAppendTrackFitMap(const Pandora &pandoraInstance, const ParticleFlowObject *const pPfo, TrackFitMap &trackFitMap,
                                             const unsigned int slidingFitWindow);
    
    /**
     *  @brief  ...
     * 
     */
    static ThreeDSlidingFitResult PerformSlidingTrackFit(const Pandora &pandoraInstance, const ParticleFlowObject *const pPfo,
                                                         const unsigned int slidingFitWindow);
                                                         
    /**
     *  @brief ...
     * 
     */
    static PfoList GetRecoNeutrinoList(const Algorithm &algorithm, const std::string &pfoListName);
    
    /**
     *  @brief ...
     * 
     */
    static float CaloHitToThreeDDistance(const Pandora &pandoraInstance, const CaloHit *const pCaloHit,
                                         const ThreeDSlidingFitResult &trackFit);
    
    /**
     *  @brief ...
     * 
     */
    static float CaloHitToThreeDDistance(const Pandora &pandoraInstance, const float hitWidth, const CartesianVector &fitDirection);
    
    /**
     *  @brief ...
     * 
     */
    static std::tuple<float, float> GetPolarAnglesFromDirection(const CartesianVector &direction);
        
    /**
     *  @brief ...
     * 
     */
    static CaloHitList GetHitsOfType(const ParticleFlowObject *const pPfo, const HitType hitType, const bool recurseOverDaughters);
    
    /**
     *  @brief ...
     * 
     */
    static CartesianVector GetFittedDirectionAtPosition(const ThreeDSlidingFitResult &trackFit, const CartesianVector &position);
    
    /**
     *  @brief ...
     * 
     */
    static bool RecursivelyCheckFiducialCut(const ParticleFlowObject *const pPfo, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates);
    
    /**
     *  @brief ...
     * 
     */
    static bool CheckFiducialCut(const ParticleFlowObject *const pPfo, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates);
    
    /**
     *  @brief ...
     * 
     */
    static bool IsPointFiducial(const CartesianVector &point, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates);
    
    /**
     *  @brief ...
     * 
     */
    static float GetParticleRange(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit);
    
    /**
     *  @brief ...
     * 
     */
    static HitProjectionVector OrderHitsByProjectionOnToTrackFit(const CaloHitList &caloHitList, const ThreeDSlidingFitResult &trackFit);
    
    /**
     *  @brief ...
     * 
     */
    static void WriteNTuple(TNtuple *const pNtuple, const std::string &fileName, const bool verboseMode);
    
    /**
     *  @brief ...
     * 
     */
    static TNtuple * LoadNTupleFromFile(const std::string &filePath, const std::string &nTupleName);
    
    /**
     *  @brief ...
     * 
     */
    static const MCParticle *GetMainMCParticle(const ParticleFlowObject *const pPfo);
    
private:
    using McParticleVotingMap    = std::unordered_map<const MCParticle *, unsigned int>; ///< ...
    using McParticleVotingPair   = std::pair<const MCParticle *, unsigned int>;          ///< ...
    using McParticleVotingVector = std::vector<McParticleVotingPair>;                    ///< ...

    /**
     *  @brief ...
     * 
     */
    static bool SortRecoNeutrinos(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs);
        
    /**
     *  @brief ...
     * 
     */
    static float CellToThreeDDistance(const float hitWidth, const float wirePitch, const CartesianVector &fitDirection);
};
} // namespace lar_physics_content

#endif // #ifndef LEE_ANALYSIS_HELPER_H
