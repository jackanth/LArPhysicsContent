/**
 *  @file   larphysicscontent/LArAnalysis/EnergyEstimatorNtupleTool.h
 *
 *  @brief  Header file for the energy estimator ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
#define LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArHelpers/LArRootHelper.h"
#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

namespace lar_physics_content
{
/**
 *  @brief  EnergyEstimatorNtupleTool class
 */
class EnergyEstimatorNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    EnergyEstimatorNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    EnergyEstimatorNtupleTool(const EnergyEstimatorNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    EnergyEstimatorNtupleTool(EnergyEstimatorNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    EnergyEstimatorNtupleTool &operator=(const EnergyEstimatorNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    EnergyEstimatorNtupleTool &operator=(EnergyEstimatorNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~EnergyEstimatorNtupleTool() = default;

protected:
    std::vector<LArNtupleRecord> ProcessEvent(
        const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

private:
    /**
     *  @brief  Struct containing hit calorimetry info
     */
    struct HitCalorimetryInfo
    {
        /**
         *  @brief  Constructor
         */
        HitCalorimetryInfo();

        bool                     m_projectionSuccessful; ///< Whether the projection was successful
        pandora::CartesianVector m_threeDPosition;       ///< If projection successful, the 3D position
        float                    m_projectionError;      ///< If projection successful, the projection error
        float                    m_coordinate;           ///< If projection successful, the 3D track coordinate
        float                    m_dQ;                   ///< The hit charge
        float                    m_dX;                   ///< The 3D dx
        bool                     m_isNoise;              ///< Whether this hit is tagged as noisy
        bool                     m_isFake;               ///< Whether this is a fake hit
    };

    using HitCalorimetryInfoPtr = std::shared_ptr<HitCalorimetryInfo>; ///< Alias for a shared to a HitCalorimetryInfo object
    using HitCalorimetryInfoMap =
        std::unordered_map<const pandora::CaloHit *, HitCalorimetryInfoPtr>; ///< Alias for a map from CaloHits to HitCalorimetryInfo shared pointers
    using CaloHitMap = std::unordered_map<const pandora::CaloHit *, const pandora::CaloHit *>; ///< Alias for a map between CaloHits

    template <typename T>
    using HitDataGetter = std::function<std::optional<std::decay_t<T>>(const pandora::CaloHit *const, const HitCalorimetryInfo &hitInfo)>; ///< Alias for a CaloHit info getter function

    bool        m_writeEnergiesToNtuple;           ///< Whether to write the energies to the ntuple
    bool        m_useParticleId;                   ///< Whether to use particle ID
    bool        m_trainingSetMode;                 ///< Whether to run in training set mode
    bool        m_makePlots;                       ///< Make plots
    std::string m_recombinationCorrectionDataFile; ///< The recombination correction data file path
    float       m_birksAlpha;                      ///< The Birks' alpha fit parameter
    float       m_birksBeta;                       ///< The Birks' beta fit parameter
    float       m_birksdQdXPole;                   ///< The Birks' dQdX pole
    bool        m_birksFitParametersSet;           ///< Whether the Birks fit parameters have been set

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get energy estimator records for a PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticle optional address of the MC particle
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> GetEnergyEstimatorRecords(
        const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Get CaloHit calorimetry info
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticle optional address of the MC particle
     *
     *  @return the dQ/dx vector, the dx vector, the shower charge, and the number of hits lost to errors
     */
    std::tuple<LArNtupleRecord::RFloatVector, LArNtupleRecord::RFloatVector, LArNtupleRecord::RFloat, LArNtupleRecord::RUInt> GetHitCalorimetryInfo(
        const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Get a 3D distance from a cell
     *
     *  @param  hitWidth the hit width
     *  @param  wirePitch the wire pitch
     *  @param  fitDirection the fit direction at the hit
     *
     *  @return the 3D distance
     */
    float CellToThreeDDistance(const float hitWidth, const float wirePitch, const pandora::CartesianVector &fitDirection) const;

    /**
     *  @brief  Get the 3D distance from CaloHit info
     *
     *  @param  hitWidth the instance of Pandora
     *  @param  trackFit the 3D track fit to which the hit belongs
     *  @param  threeDPosition the 3D position
     *  @param  threeDDistance the 3D distance (to populate)
     *
     *  @return the status code
     */
    pandora::StatusCode CaloHitToThreeDDistance(const float hitWidth, const lar_content::ThreeDSlidingFitResult &trackFit,
        const pandora::CartesianVector &threeDPosition, float &threeDDistance) const;

    /**
     *  @brief  Turn a direction in polar and azimuthal angles
     *
     *  @param  direction the Cartesian direction
     *
     *  @return the polar and azimuthal angles
     */
    std::tuple<float, float> GetPolarAnglesFromDirection(const pandora::CartesianVector &direction) const;

    /**
     *  @brief  Produce training records for a cosmic or primary PFO
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional pointer to the corresponding MCParticle
     *
     *  @return the training records
     */
    std::vector<LArNtupleRecord> ProduceTrainingRecords(
        const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const pandora::MCParticle *const pMcParticle);

    /**
     *  @brief  Select noisy hits
     *
     *  @param  trackFit the 3D track fit
     *  @param  caloHitList the CaloHit list
     *  @param  pMCParticle optional address of the MCParticle
     *  @param  caloHitMap the map from 2D to 3D CaloHits
     *
     *  @return the vector of dQ/dx values, the vector of dx values, and the total shower charge
     */
    std::tuple<pandora::FloatVector, pandora::FloatVector, float> SelectNoisyHits(const lar_content::ThreeDSlidingFitResult &trackFit,
        pandora::CaloHitList caloHitList, const pandora::MCParticle *const pMCParticle, const CaloHitMap &caloHitMap) const;

    /**
     *  @brief  Identify noisy hits that weren't projected
     *
     *  @param  hitInfoMap the hit info map
     */
    void IdentifyNoisyUnprojectedHits(const HitCalorimetryInfoMap &hitInfoMap) const;

    /**
     *  @brief  Sum up the noisy charge
     *
     *  @param  hitInfoMap the hit info map
     *  @param  useAllHits whether to use all hits
     *
     *  @return the summed charge
     */
    float SumNoisyCharge(const HitCalorimetryInfoMap &hitInfoMap, const bool useAllHits) const;

    /**
     *  @brief  Generate the hit calorimetry info map
     *
     *  @param  trackFit the track fit object
     *  @param  caloHitList the CaloHit list
     *  @param  caloHitMap the map from 2D to 3D CaloHits
     *
     *  @return the hit calorimetry info map
     */
    HitCalorimetryInfoMap CalculateHitCalorimetryInfo(
        const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CaloHitList &caloHitList, const CaloHitMap &caloHitMap) const;

    /**
     *  @brief  Calculate the hit calorimetry info
     *
     *  @param  trackFit the track fit object
     *  @param  twoDPositionVector the 2D position
     *  @param  dQ the hit charge
     *  @param  hitWidth the hit width
     *  @param  caloHitMap the map from 2D to 3D CaloHits
     *  @param  pCaloHit address of the CaloHit
     *
     *  @return shared pointer to the hit calorimetry info object
     */
    HitCalorimetryInfoPtr CalculateHitCalorimetryInfo(const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CartesianVector &twoDPositionVector,
        const float dQ, const float hitWidth, const CaloHitMap &caloHitMap, const pandora::CaloHit *const pCaloHit) const;

    /**
     *  @brief  Identify noisy hits
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  pdQdXFunction address of the dQ/dx fit function
     *  @param  resetMode whether to run in reset mode
     *
     *  @return whether a change was made
     */
    bool IdentifyNoisyHits(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdQdXFunction, const bool resetMode) const;

    /**
     *  @brief  Fit the dQ/dX function
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  iteration the iteration number
     *
     *  @return address of the fitted TF1 object
     */
    TF1 *FitdQdXFunction(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, const std::size_t iteration) const;

    /**
     *  @brief  Get the gap parameters
     *
     *  @param  trackFit the track fit objects
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *
     *  @return the live length, the dead length, the min coordinate, and the max coordinate
     */
    std::tuple<float, float, float, float> GetGapParameters(
        const lar_content::ThreeDSlidingFitResult &trackFit, pandora::CaloHitList &caloHitList, HitCalorimetryInfoMap &hitInfoMap) const;

    /**
     *  @brief  Calculate the coordinate range
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  minCoordinate the min coordinate (to populate)
     *  @param  maxCoordinate the max coordinate (to populate)
     *
     *  @return the status code
     */
    pandora::StatusCode CalculateCoordinateRange(
        const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, float &minCoordinate, float &maxCoordinate) const;

    /**
     *  @brief  Make a 1D hit histogram
     *
     *  @param  options the plot options
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  hitDataGetter the function for getting the data to plot
     *
     *  @return address of the TH1F object
     */
    TH1F *MakeOneDHitHistogram(const LArRootHelper::PlotOptions &options, const pandora::CaloHitList &caloHitList,
        const HitCalorimetryInfoMap &hitInfoMap, const HitDataGetter<float> &hitDataGetter) const;

    /**
     *  @brief  Make a 2D hit histogram
     *
     *  @param  options the plot options
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  hitDataGetter the function for getting the data to plot
     *
     *  @return address of the TH2F object
     */
    TH2F *MakeTwoDHitHistogram(const LArRootHelper::PlotOptions &options, const pandora::CaloHitList &caloHitList,
        const HitCalorimetryInfoMap &hitInfoMap, const HitDataGetter<std::pair<float, float>> &hitDataGetter) const;

    /**
     *  @brief  Get the standard deviation to the fit
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  pdQdXFunction address of the fitted dQ/dx function
     *
     *  @return the standard deviation
     */
    float GetFitStandardDeviation(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdQdXFunction) const;

    /**
     *  @brief  Make a colour-coded hit plot
     *
     *  @param  title the plot title
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  pdQdXFunction address of the dQ/dx function
     */
    void MakeColourCodedHitPlot(const std::string &title, const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap,
        TF1 *pdQdXFunction) const;

    /**
     *  @brief  Get the number of good hits
     *
     *  @param  hitInfoMap the hit info map
     *
     *  @return the number of good hits
     */
    std::size_t GetNumGoodHits(const HitCalorimetryInfoMap &hitInfoMap) const;

    /**
     *  @brief  Plot the hit dQ values
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     */
    void PlotHitdQ(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const;

    /**
     *  @brief  Plot the hit dx values
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     */
    void PlotHitdX(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const;

    /**
     *  @brief  Plot the hit dQ/dx values
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     */
    void PlotHitdQdX(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const;

    /**
     *  @brief  Plot the true matched hit dQ/dx
     *
     *  @param  caloHitList the CaloHit list
     *  @param  hitInfoMap the hit info map
     *  @param  pMCParticle address of the MCParticle
     */
    void PlotTrueMatchedHitdQdX(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap,
        const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Plot a histogram and a function on the same canvas
     *
     *  @param  pHistogram address of the histogram
     *  @param  pFunction address of the function
     */
    void PlotHistogramAndFunction(TH1 *pHistogram, TF1 *pFunction) const;

    /**
     *  @brief  Esimate the energy of a track hit
     *
     *  @param  dQdX the dQ/dx value
     *  @param  dX the dx value
     *
     *  @return the estimated energy
     */
    float EstimateTrackHitEnergy(const float dQdX, const float dX) const;

    /**
     *  @brief  Estimate the energy of some shower charge
     *
     *  @param  showerCharge the shower charge
     *
     *  @return the estimated energy
     */
    float EstimateShowerEnergy(const float showerCharge) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
