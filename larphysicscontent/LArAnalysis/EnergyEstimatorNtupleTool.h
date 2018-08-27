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
    std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

private:
    struct HitCalorimetryInfo
    {
        HitCalorimetryInfo() :
            m_projectionSuccessful(false),
            m_threeDPosition(0.f, 0.f, 0.f),
            m_projectionError(0.f),
            m_coordinate(0.f),
            m_dQ(0.f),
            m_dX(0.f),
            m_isNoise(false),
            m_isFake(false)
        {
        }

        bool                     m_projectionSuccessful;
        pandora::CartesianVector m_threeDPosition;
        float                    m_projectionError;
        float                    m_coordinate;
        float                    m_dQ;
        float                    m_dX;
        bool                     m_isNoise;
        bool                     m_isFake;
    };

    using HitCalorimetryInfoPtr = std::shared_ptr<HitCalorimetryInfo>;
    using HitCalorimetryInfoMap = std::unordered_map<const pandora::CaloHit *, HitCalorimetryInfoPtr>;

    template <typename T>
    using HitDataGetter = std::function<std::optional<std::decay_t<T>>(const pandora::CaloHit *const, const HitCalorimetryInfo &hitInfo)>;

    bool                             m_writeEnergiesToNtuple; ///< Whether to write the energies to the ntuple
    bool                             m_useParticleId;         ///< Whether to use particle ID
    bool                             m_trainingSetMode;       ///< Whether to run in training set mode
    bool                             m_makePlots;             ///< Make plots
    std::shared_ptr<LArRootRegistry> m_spTmpRegistry;         ///< Shared pointer to the tmp registry
    std::shared_ptr<LArRootRegistry> m_spPlotRegistry;        ///< Shared pointer to the plot registry

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

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
     *  @brief  Get the 3D distance from a CaloHit
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pCaloHit address of the CaloHit
     *  @param  trackFit the 3D track fit to which the hit belongs
     *
     *  @return the 3D distance
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
     *  @param  pMCParticle optional pointer to the corresponding MC particle
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return the training records
     */
    std::vector<LArNtupleRecord> ProduceTrainingRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &,
        const pandora::MCParticle *const pMcParticle, const pandora::MCParticleList *const);

    std::tuple<pandora::CaloHitList, pandora::CaloHitList> SelectNoisyHits(const lar_content::ThreeDSlidingFitResult &trackFit,
        pandora::CaloHitList caloHitList, const pandora::MCParticle *const pMCParticle) const;

    TF1 *FitHitGenerationFunction(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const;

    HitCalorimetryInfoMap CalculateHitCalorimetryInfo(const lar_content::ThreeDSlidingFitResult &trackFit, const pandora::CaloHitList &caloHitList) const;

    HitCalorimetryInfoPtr CalculateHitCalorimetryInfo(const lar_content::ThreeDSlidingFitResult &trackFit,
        const pandora::CartesianVector &twoDPositionVector, const float dQ, const float hitWidth) const;

    bool IdentifyNoisyHits(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdEdXFunction) const;

    TF1 *FitdQdXFunction(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap,
        const std::size_t iteration) const;

    void GenerateNewHits(const lar_content::ThreeDSlidingFitResult &trackFit,
        pandora::CaloHitList &caloHitList, HitCalorimetryInfoMap &hitInfoMap, TF1 *pHitGenerationFunction) const;

    pandora::StatusCode CalculateCoordinateRange(
        const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, float &minCoordinate, float &maxCoordinate) const;

    TH1F *MakeOneDHitHistogram(const LArRootHelper::PlotOptions &options,
        const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, const HitDataGetter<float> &hitDataGetter) const;

    TH2F *MakeTwoDHitHistogram(const LArRootHelper::PlotOptions &options, const pandora::CaloHitList &caloHitList,
        const HitCalorimetryInfoMap &hitInfoMap, const HitDataGetter<std::pair<float, float>> &hitDataGetter) const;

    float GetFitStandardDeviation(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdQdXFunction) const;

    void MakeColourCodedHitPlot(const std::string &title, const pandora::CaloHitList &caloHitList,
        const HitCalorimetryInfoMap &hitInfoMap, TF1 *pdQdXFunction) const;

    std::size_t GetNumGoodHits(const HitCalorimetryInfoMap &hitInfoMap) const;

    void PlotHitdQ(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const;

    void PlotHitdX(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const;

    void PlotHitdQdX(const pandora::CaloHitList &caloHitList, const HitCalorimetryInfoMap &hitInfoMap) const;

    void PlotTrueMatchedHitdQdX(const pandora::CaloHitList &caloHitList,
        const HitCalorimetryInfoMap &hitInfoMap, const pandora::MCParticle *const pMCParticle) const;

    void PlotHistogramAndFunction(TH1 *pHistogram, TF1 *pFunction) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ENERGY_ESTIMATOR_NTUPLE_TOOL_H
