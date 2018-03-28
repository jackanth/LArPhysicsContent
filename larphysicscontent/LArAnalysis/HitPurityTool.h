/**
 *  @file   larphysicscontent/LArAnalysis/HitPurityTool.h
 *
 *  @brief  Header file for the hit purity tool class.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_PURITY_TOOL_H
#define LAR_HIT_PURITY_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

#include "larphysicscontent/LArObjects/LArTrackHitValue.h"
#include "larphysicscontent/LArObjects/LArFittedTrackInfo.h"

#include "larphysicscontent/LArHelpers/LArAnalysisParticleHelper.h"

using namespace pandora;

namespace lar_physics_content
{

/**
 *  @brief  HitPurityTool class
 */
class HitPurityTool : public AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    HitPurityTool();

    bool Run(const Algorithm *const pAlgorithm, LArFittedTrackInfo::TrackHitValueVector &trackHitValueVector, float &excessCaloValue);

private:
    using FloatMatrix = std::vector<FloatVector>; ///< Alias for a matrix of inter-datapoint distances

    float          m_maxImpurityScore;            ///< The maximum impurity score before we decide to adjust
    std::size_t    m_valueAverageSearchRadius;    ///< The search radius for calculating value averages
    std::size_t    m_nearestNeighbourNumber;      ///< The number of nearest neighbours to consider
    bool           m_makePlots;                   ///< ATTN temporary

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the calo value averages about each hit (excluding the hit)
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  nHits the number of hits
     *
     *  @return the vector of averages
     */
    FloatVector GetValueAverages(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits) const;

    /**
     *  @brief  Get the calo value/coordinate scale factor, inter-datapoint distance mean and standard deviation
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  nHits the number of hits
     *  @param  scaleFactor the scale factor (to populate)
     *  @param  mean the mean (to populate)
     *  @param  sigma the standard deviation (to populate)
     *
     *  @return success
     */
    bool GetStatistics(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits, float &scaleFactor,
        float &mean, float &sigma) const;

    /**
     *  @brief  Get the calo value/coordinate ranges
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  nHits the number of hits
     *  @param  coordinateRange the coordinate range (to populate)
     *  @param  caloValueRange the calo value range (to populate)
     */
    void CalculateRanges(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits,
        float &coordinateRange, float &caloValueRange) const;

    /**
     *  @brief  Get the inter-datapoint distance mean
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  nHits the number of hits
     *  @param  scaleFactor the scale factor
     *
     *  @return the mean
     */
    float CalculateMean(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits, const float scaleFactor) const;

    /**
     *  @brief  Get the inter-datapoint distance standard deviation
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  nHits the number of hits
     *  @param  scaleFactor the scale factor
     *  @param  mean the mean
     *
     *  @return the standard deviation
     */
    float CalculateStandardDeviation(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector, const std::size_t nHits,
        const float scaleFactor, const float mean) const;

    /**
     *  @brief  Calculate the impurity scores
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  scaleFactor the scale factor
     *  @param  nHits the number of hits
     *  @param  mean the mean
     *  @param  sigma the standard deviation
     *
     *  @return the impurity scores
     */
    FloatVector CalculateImpurityScores(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector, const float scaleFactor,
        const std::size_t nHits, const float mean, const float sigma) const;

    /**
     *  @brief  Get the vector of inter-datapoint distances, sorted smallest-to-largest
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  currentTrackHitValue the current track hit value
     *  @param  scaleFactor the scale factor
     *  @param  nHits the number of hits
     *
     *  @return the distance vector
     */
    FloatVector GetSortedDistanceVector(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector,
        const LArTrackHitValue &currentTrackHitValue, const float scaleFactor, const std::size_t nHits) const;

    // ATTN temporary
    void MakePlots(const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergyVector, const LArFittedTrackInfo::TrackHitValueVector &trackHitEnergiesChanged) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_HIT_PURITY_TOOL_H
