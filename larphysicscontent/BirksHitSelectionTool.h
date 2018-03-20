/**
 *  @file   LArPhysicsContent/include/BirksHitSelectionTool.h
 *
 *  @brief  Header file for the Birks' hit selection tool class.
 *
 *  $Log: $
 */
#ifndef BIRKS_HIT_SELECTION_TOOL_H
#define BIRKS_HIT_SELECTION_TOOL_H 1

#include "Pandora/AlgorithmTool.h"
#include "LArTrackHitEnergy.h"

using namespace pandora;

namespace lar_physics_content
{

/**
 *  @brief  BirksHitSelectionTool class
 */
class BirksHitSelectionTool : public AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    BirksHitSelectionTool();

    bool Run(const Algorithm *const pAlgorithm, LArTrackHitEnergy::Vector &trackHitEnergyVector, int &uniquePlotIdentifier);

private:
    using PointDistancePair   = std::pair<float, float>;          ///< Alias for a pair of datapoint distances
    using PointDistanceVector = std::vector<PointDistancePair>;   ///< Alias for a vector of datapoint distance pairs
    using PointDistanceMatrix = std::vector<PointDistanceVector>; ///< Alias for a matrix of datapoint distance pairs

    float    m_birksSelectionBinWidth;                ///< The data bin width
    float    m_birksSelectionSearchXWidthFraction;    ///< The x-width search fraction
    float    m_birksSelectionSearchYWidth;            ///< The y-width search fraction
    float    m_birksSelectionMinDensityFraction;      ///< The min density fraction
    bool     m_makePlots;                             ///< Whether to produce plots (ATTN: temporary)

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Decide which hits in the track hit energy vector should be corrected for recombination
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  uniquePlotIdentifier the unique plot identifier (ATTN: temporary)
     */
    void DecideWhichTrackHitsToCorrect(LArTrackHitEnergy::Vector &trackHitEnergyVector, int &uniquePlotIdentifier) const;

    /**
     *  @brief  Get the extremal track coordinates and energies from the track hit energy vector
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *
     *  @return the minimum and maximum track coordinates, and the minimum and maximum hit energies
     */
    std::tuple<float, float, float, float> GetExtremalBirksSelectionValues(const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;

    /**
     *  @brief  Create the matrix of inter-datapoint distances from the track hit energy vector
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *
     *  @return the matrix of distances
     */
    PointDistanceMatrix CreateMatrixOfPointDistances(const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;

    /**
     *  @brief  Get the point densities and average energies for each bin
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  searchRadiusX the search radius for the x-coordinate
     *  @param  matrixOfPointDistances the matrix of inter-datapoint distances
     *  @param  numBins the number of bins
     *  @param  minCoordinate the minimum track coordinate
     *
     *  @return the point densities and average energies in each bin
     */
    std::tuple<FloatVector, FloatVector> GetBinAverages(const LArTrackHitEnergy::Vector &trackHitEnergyVector, const float searchRadiusX,
        const PointDistanceMatrix &matrixOfPointDistances, const int numBins, const float minCoordinate) const;

    /**
     *  @brief  Get the point density and average energy for a given bin
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  xLower the lower x-coordinate
     *  @param  xUpper the upper x-coordinate
     *  @param  searchRadiusX the search radius for the x-coordinate
     *  @param  matrixOfPointDistances the matrix of inter-datapoint distances
     *
     *  @return the point density and average energy
     */
    std::tuple<float, float> GetBinAverage(const LArTrackHitEnergy::Vector &trackHitEnergyVector, const float xLower, const float xUpper,
        const float searchRadiusX, const PointDistanceMatrix &matrixOfPointDistances) const;

    /**
     *  @brief  Iteratively make the selection changes to the track hit energy vector until no more changes are made
     *
     *  @param  trackHitEnergyVector the track hit energy vector
     *  @param  searchRadiusX the search radius for the x-coordinate
     *  @param  matrixOfPointDistances the matrix of inter-datapoint distances
     *  @param  numBins the number of bins
     *  @param  changeMade whether a change has been made this iteration
     *  @param  minCoordinate the minimum x-coordinate
     *  @param  averageBinPointDensity the set of point densities for each bin
     *  @param  averageBinEnergy the set of average energies for each bin
     *  @param  trackHitEnergiesChanged the track hit energies that have been changed this iteration
     *  @param  allTrackHitEnergiesChanged all the track hit energies that have been changed
     */
    void MakeBirksSelectionChanges(LArTrackHitEnergy::Vector &trackHitEnergyVector, const float searchRadiusX,
        const PointDistanceMatrix &matrixOfPointDistances, const int numBins, bool &changeMade, const float minCoordinate,
        const FloatVector &averageBinPointDensity, const FloatVector &averageBinEnergy, LArTrackHitEnergy::Vector &trackHitEnergiesChanged,
        LArTrackHitEnergy::Vector &allTrackHitEnergiesChanged) const;

    /**
     *  @brief  Produce plots (ATTN: temporary)
     *
     *  @param  trackHitEnergyVector the vector of track hit energies
     *  @param  trackHitEnergiesChanged the set of track hit energies that were changed
     *  @param  birksPointsOnly whether to plot the Birks-corrected points only
     *  @param  pause whether to pause
     *  @param  plotUncorrectedEnergiesOnly whether to plot the uncorrected energies only
     *  @param  uniquePlotIdentifier the unique plot identifier
     */
    void ProduceHitSelectionPlots(const LArTrackHitEnergy::Vector &trackHitEnergyVector, const LArTrackHitEnergy::Vector &trackHitEnergiesChanged,
        const bool birksPointsOnly, const bool pause, const bool plotUncorrectedEnergiesOnly, int &uniquePlotIdentifier) const;
};

} // namespace lar_physics_content

#endif // #ifndef BIRKS_HIT_SELECTION_TOOL_H
