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
    using PointDistancePair   = std::pair<float, float>;          ///< 
    using PointDistanceVector = std::vector<PointDistancePair>;   ///< 
    using PointDistanceMatrix = std::vector<PointDistanceVector>; ///< 
    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    /**
     *  @brief ...
     * 
     */
    void DecideWhichTrackHitsToCorrect(LArTrackHitEnergy::Vector &trackHitEnergyVector, int &uniquePlotIdentifier) const;
    
    /**
     *  @brief ...
     * 
     */
    std::tuple<float, float, float, float> GetExtremalBirksSelectionValues(const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;
    
    /**
     *  @brief ...
     * 
     */
    PointDistanceMatrix CreateMatrixOfPointDistances(const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;
    
    /**
     *  @brief ...
     * 
     */
    std::tuple<FloatVector, FloatVector> GetBinAverages(const LArTrackHitEnergy::Vector &trackHitEnergyVector, const float searchRadiusX,
                                                        const PointDistanceMatrix &matrixOfPointDistances, const int numBins,
                                                        const float minCoordinate) const;
    
    /**
     *  @brief ...
     * 
     */
    std::tuple<float, float> GetBinAverage(const LArTrackHitEnergy::Vector &trackHitEnergyVector, const float xLower,
                                           const float xUpper, const float searchRadiusX,
                                           const PointDistanceMatrix &matrixOfPointDistances) const;
    
    /**
     *  @brief  ...
     * 
     */
    void MakeBirksSelectionChanges(LArTrackHitEnergy::Vector &trackHitEnergyVector, const float searchRadiusX,
                                          const PointDistanceMatrix &matrixOfPointDistances, const int numBins, bool &changeMade,
                                          const float minCoordinate, const FloatVector &averageBinPointDensity,
                                          const FloatVector &averageBinEnergy, LArTrackHitEnergy::Vector &trackHitEnergiesChanged,
                                          LArTrackHitEnergy::Vector &allTrackHitEnergiesChanged) const;
    
    /**
     *  @brief  ...
     * 
     */
    void ProduceHitSelectionPlots(const LArTrackHitEnergy::Vector &trackHitEnergyVector,
                                         const LArTrackHitEnergy::Vector &trackHitEnergiesChanged, const bool birksPointsOnly,
                                         const bool pause, const bool plotUncorrectedEnergiesOnly, int &uniquePlotIdentifier) const;
    
    float m_birksSelectionBinWidth;             ///< 
    float m_birksSelectionSearchXWidthFraction; ///< 
    float m_birksSelectionSearchYWidth;         ///< 
    float m_birksSelectionMinDensityFraction;   ///<
    bool  m_makePlots;
};

} // namespace lar_physics_content

#endif // #ifndef BIRKS_HIT_SELECTION_TOOL_H
