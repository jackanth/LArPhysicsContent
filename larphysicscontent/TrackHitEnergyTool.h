/**
 *  @file   larphysicscontent/TrackHitEnergyTool.h
 *
 *  @brief  Header file for the track hit energy tool class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_HIT_ENERGY_TOOL_H
#define LAR_TRACK_HIT_ENERGY_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

#include "larphysicscontent/LArTrackHitValue.h"
#include "larphysicscontent/LArAnalysisParticleHelper.h"
#include "larphysicscontent/LArFittedTrackInfo.h"

#include <functional>

using namespace pandora;

namespace lar_physics_content
{

/**
 *  @brief  TrackHitEnergyTool class
 */
class TrackHitEnergyTool : public AlgorithmTool
{
public:
    using HitPurityToolCallback = std::function<bool(LArFittedTrackInfo::TrackHitValueVector &, float &)>; ///< Alias for hit purity tool callback
    
    /**
     *  @brief  Default constructor
     */
    TrackHitEnergyTool();

    bool Run(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap,
        float &excessCharge, const HitPurityToolCallback &hitPurityToolCallback);

private:
    using HitProjectionPair   = std::pair<const CaloHit *, float>; ///< Alias for a map from CaloHits to their projected track coordinate
    using HitProjectionVector = std::vector<HitProjectionPair>;    ///< Alias for a vector of hit projection pairs
    
    unsigned    m_trackSlidingFitWindow;    ///< The track sliding fit window to use.

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    /**
     *  @brief  Recurse through the PFO hierarchy and append the track fit map
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pPfo address of the current PFO
     *  @param  trackFitMap the track fit map to append
     *  @param  slidingFitWindow the sliding fit window size
     */
    void RecursivelyAppendMap(const ParticleFlowObject *const pPfo, float &excessCharge, HitPurityToolCallback hitPurityToolCallback,
        LArAnalysisParticleHelper::FittedTrackInfoMap &fittedTrackInfoMap) const;
    
    /**
     *  @brief  Perform a 3D sliding track fit for a given PFO
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pPfo address of the PFO
     *  @param  slidingFitWindow the sliding fit window size
     *
     *  @return the 3D track fit
     */
    ThreeDSlidingFitResult PerformSlidingTrackFit(const ParticleFlowObject *const pPfo) const;
    
    /**
     *  @brief  Append the track hit energy map for a given PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFit the PFO's 3D fit
     *
     *  @return the vector of track hit energies
     */
    LArFittedTrackInfo::TrackHitValueVector AppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo,
        const ThreeDSlidingFitResult &trackFit, float &excessCharge, HitPurityToolCallback hitPurityToolCallback) const;
        
    /**
     *  @brief  Get a 3D distance from a cell
     *
     *  @param  hitWidth the hit width
     *  @param  wirePitch the wire pitch
     *  @param  fitDirection the fit direction at the hit
     *
     *  @return the 3D distance
     */
    float CellToThreeDDistance(const float hitWidth, const float wirePitch, const CartesianVector &fitDirection) const;
    
    /**
     *  @brief  Get the 3D distance from a CaloHit
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pCaloHit address of the CaloHit
     *  @param  trackFit the 3D track fit to which the hit belongs
     *
     *  @return the 3D distance
     */
    float CaloHitToThreeDDistance(const CaloHit *const pCaloHit, const ThreeDSlidingFitResult &trackFit) const;

    /**
     *  @brief  Turn a direction in polar and azimuthal angles
     *
     *  @param  direction the Cartesian direction
     *
     *  @return the polar and azimuthal angles
     */
    std::pair<float, float> GetPolarAnglesFromDirection(const CartesianVector &direction) const;
    
        /**
     *  @brief  Get the range of a track-like PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFit the track fit for the PFO
     *
     *  @return the range
     */
    float GetParticleRange(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit) const;

    /**
     *  @brief  Order a set of CaloHits by their projection onto a track fit and return the ordered projections
     *
     *  @param  caloHitList the list of CaloHits
     *  @param  trackFit the track fit object
     *
     *  @return the ordered list of hit projections
     */
    HitProjectionVector OrderHitsByProjectionOnToTrackFit(const CaloHitList &caloHitList, const ThreeDSlidingFitResult &trackFit) const;

};

} // namespace lar_physics_content

#endif // #ifndef LAR_TRACK_HIT_ENERGY_TOOL_H
