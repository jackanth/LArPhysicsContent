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
 *  @brief  AnalysisAlgorithm class
 */
class AnalysisAlgorithm : public CustomParticleCreationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    AnalysisAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~AnalysisAlgorithm();

protected:
    StatusCode ReadSettings(const TiXmlHandle xmlHandle) override;
    void CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject*& pOutputPfo) const override;

private:
    /**
     *  @brief  Energy-from-range data entry class
     */
    class EnergyFromRangeData
    {
        public:
        /**
         *  @brief  Constructor
         *
         *  @param  rangeMin the lower range bound
         *  @param  rangeMax the upper range bound
         *  @param  energy the energy for this range bin
         */
        EnergyFromRangeData(const float rangeMin, const float rangeMax, const float energy) noexcept;

        float    m_rangeMin;    ///< The lower range bound
        float    m_rangeMax;    ///< The upper range bound
        float    m_energy;      ///< The energy for the range bin
    };

    using EnergyFromRangeDataVector = std::vector<EnergyFromRangeData>; ///< Alias for a vector of energy-from-range data entries

    float                     m_fiducialCutLowXMargin;             ///< The low-x fiducial margin
    float                     m_fiducialCutHighXMargin;            ///< The high-x fiducial cut margin
    float                     m_fiducialCutLowYMargin;             ///< The low-y fiducial cut margin
    float                     m_fiducialCutHighYMargin;            ///< The high-y fiducial cut margin
    float                     m_fiducialCutLowZMargin;             ///< The low-z fiducial cut margin
    float                     m_fiducialCutHighZMargin;            ///< The high-z fiducial cut margin
    float                     m_birksSelectionMaxdEdX;             ///< The maximum corrected dEdX value allowed.
    float                     m_mcContainmentFractionLowerBound;   ///< The lower containment fraction bound for MC containment
    unsigned int              m_trackSlidingFitWindow;             ///< The sliding fit window for 3D track fits
    std::string               m_mcParticleListName;                ///< The name of the MC particle list
    std::string               m_parametersFile;                    ///< The path to the file containing the fit data
    std::string               m_birksFitNtupleName;                ///< The name of the Birks fit ntuple
    std::string               m_protonEnergyFromRangeNtupleName;   ///< The name of the proton energy-from-range ntuple
    std::string               m_pionMuonEnergyFromRangeNtupleName; ///< The name of the pion/muon energy-from-range ntuple
    std::string               m_caloHitListName;                   ///< The name of the CaloHit list
    std::string               m_tmvaWeights;                       ///< The path to the file containing the TMVA weights for proton ID
    bool                      m_addMcInformation;                  ///< Whether to add MC information to the analysis particles

    float                     m_birksFitAlpha;                     ///< The Birks fit alpha parameter
    float                     m_birksFitBeta;                      ///< The Birks fit beta parameter
    float                     m_birksFitPole;                      ///< The Birks fit dEdX pole value
    EnergyFromRangeDataVector m_protonEnergyFromRangeDataVector;   ///< The vector of proton energy-from-range data entries
    EnergyFromRangeDataVector m_pionMuonEnergyFromRangeDataVector; ///< The vector of pion/muon energy-from-range data entries
    CartesianVector           m_minCoordinates;                    ///< The detector's minimum fiducial coordinates
    CartesianVector           m_maxCoordinates;                    ///< The detector's maximum fiducial coordinates
    BirksHitSelectionTool    *m_pBirksHitSelectionTool;            ///< Address of the Birks hit-selection tool
    TMVA::Reader             *m_pTmvaReader;                       ///< Address of the TMVA Reader object

    mutable float             m_tmvaTrackLength;                   ///< Mutable track length member variable for TMVA to use
    mutable float             m_tmvaAvgEnergyDeposition;           ///< Mutable average energy deposition member variable for TMVA to use
    mutable int               m_uniquePlotIdentifier;              ///< Unique plot identifier (ATTN temporary)

    /**
     *  @brief  Recurse through the PFO hierarchy and append the track hit energy map
     *
     *  @param  pPfo address of the current PFO
     *  @param  trackHitEnergyMap the track hit energy map to append
     *  @param  trackFitMap the track fit map
     */
    void RecursivelyAppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo,
        LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap) const;

    /**
     *  @brief  Append the track hit energy map for a given PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFit the PFO's 3D fit
     *
     *  @return the vector of track hit energies
     */
    LArTrackHitEnergy::Vector AppendLArTrackHitEnergyMap(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit) const;

    /**
     *  @brief  Recurse through the PFO hierarchy and append the particle type map
     *
     *  @param  pPfo address of the current PFO
     *  @param  pfoTypeMap the type map to append
     *  @param  trackHitEnergyMap the track hit energy map
     *  @param  trackFitMap the track fit map
     */
    void RecursivelyAppendParticleTypeMap(const ParticleFlowObject *const pPfo, LArAnalysisParticle::PfoTypeMap &pfoTypeMap,
        const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap) const;

    /**
     *  @brief  Estimate the particle type for a given PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  trackHitEnergyMap the track hit energy map
     *  @param  trackFitMap the track fit map
     *
     *  @return the estimated particle type
     */
    LArAnalysisParticle::TYPE EstimateParticleType(const ParticleFlowObject *const pPfo,
        const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap) const;

    /**
     *  @brief  Estimate the energy of a PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  typeMap the particle type map
     *  @param  trackFitMap the particle track fit map
     *  @param  trackHitEnergyMap the track hit energy map
     *  @param  particleEnergy the particle energy value to populate
     *  @param  energySourcedFromRange the amount of energy sourced from range value to populate
     *  @param  energySourcedFromShowerCharge the amount of energy sourced from shower charge value to populate
     *  @param  energySourcedFromTrackCharge the amount of energy sourced from uncorrected track charge value to populate
     *  @param  energySourcedFromCorrectedTrackCharge the amount of energy sourced from Birks-corrected track charge value to populate
     */
    void EstimateParticleEnergy(const ParticleFlowObject *const pPfo, const LArAnalysisParticle::PfoTypeMap &typeMap,
        const LArAnalysisParticleHelper::TrackFitMap &trackFitMap, const LArAnalysisParticleHelper::LArTrackHitEnergyMap &trackHitEnergyMap,
        float &particleEnergy, float &energySourcedFromRange, float &energySourcedFromShowerCharge, float &energySourcedFromTrackCharge,
        float &energySourcedFromCorrectedTrackCharge) const;

    /**
     *  @brief  Estimate the energy of a shower-like PFO
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the estimated energy
     */
    float EstimateShowerEnergy(const ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Estimate the energy of a PFO from its charge
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the estimated energy
     */
    float EstimateEnergyFromCharge(const ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get an energy value from a CaloHit (not Birks-corrected)
     *
     *  @param  pCaloHit address of the CaloHit
     *
     *  @return the hit energy
     */
    float CaloHitToScaledEnergy(const CaloHit *const pCaloHit) const;

    /**
     *  @brief Correct an energy value for recombination using Birks' correction
     *
     *  @param scaledEnergy the uncorrected energy value
     *  @param threeDDistance the 3D distance over which the energy was deposited
     *
     *  @return the recombination-corrected energy value
     */
    float RecombinationCorrectEnergy(const float scaledEnergy, const float threeDDistance) const;

    /**
     *  @brief  Estimate the energy of a track from its charge
     *
     *  @param  trackHitEnergyVector the track hit energy vector for the track
     *
     *  @return the estimated track energy
     */
    float EstimateTrackEnergyFromCharge(const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;

    /**
     *  @brief  Decide which track hits to correct for recombination for a track-like PFO
     *
     *  @param  trackHitEnergyVector the track hit energy vector for the track-like PFO
     */
    void DecideWhichTrackHitsToCorrect(LArTrackHitEnergy::Vector &trackHitEnergyVector) const;

    /**
     *  @brief  Estimate the energy of a track-like PFO using its range
     *
     *  @param  pPfo address of the track-like PFO
     *  @param  trackFit the PFO's 3D fit
     *  @param  particleType the particle's estimated type
     *  @param  trackHitEnergyVector the particle's track hit energy vector
     *
     *  @return the estimated energy
     */
    float EstimateTrackEnergyFromRange(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit,
        const LArAnalysisParticle::TYPE particleType, const LArTrackHitEnergy::Vector &trackHitEnergyVector) const;

    /**
     *  @brief  Create a type tree for a given PFO and its descendents using the type map
     *
     *  @param  pPfo address of the PFO
     *  @param  typeMap the type map
     *
     *  @return the type tree
     */
    LArAnalysisParticle::TypeTree CreateTypeTree(const ParticleFlowObject *const pPfo, const LArAnalysisParticle::PfoTypeMap &typeMap) const;

    /**
     *  @brief  Get the direction of a PFO at a vertex
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFitMap the track fit map
     *  @param  pVertex address of the vertex
     *  @param  isCosmicRay Whether the PFO is a cosmic ray
     *
     *  @return the direction of the PFO at the vertex
     */
    CartesianVector GetDirectionAtVertex(const ParticleFlowObject *const pPfo, const LArAnalysisParticleHelper::TrackFitMap &trackFitMap,
        const Vertex *const pVertex, const bool isCosmicRay) const;

    /**
     *  @brief  Get MC information for a PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  mcEnergy the MC energy value to populate
     *  @param  mcKineticEnergy the MC kinetic energy value to populate
     *  @param  mcMass the MC mass value to populate
     *  @param  typeTree the MC type tree value to populate
     *  @param  mcType the MC type value to populate
     *  @param  mcVertexPosition the MC vertex position value to populate
     *  @param  mcMomentum the MC momentum value to populate
     *  @param  mcPdgCode the MC PDG code value to populate
     *  @param  isNeutrino whether the PFO is a neutrino
     *  @param  mcContainmentFraction the MC containment fraction value to populate
     *  @param  pMcMainMCParticle address of the main MC particle for the PFO
     *  @param  mcHitPurity the MC hit purity value to populate
     *  @param  mcHitCompleteness the MC hit completeness value to populate
     *  @param  mcCollectionPlaneHitPurity the MC collection plane hit purity value to populate
     *  @param  mcCollectionPlaneHitCompleteness the MC collection plane hit completeness value to populate
     *
     *  @return success
     */
    bool GetMcInformation(const ParticleFlowObject *const pPfo, float &mcEnergy, float &mcKineticEnergy, float &mcMass,
        LArAnalysisParticle::TypeTree &typeTree, LArAnalysisParticle::TYPE &mcType, CartesianVector &mcVertexPosition, CartesianVector &mcMomentum,
        int &mcPdgCode, const bool isNeutrino, float &mcContainmentFraction, const MCParticle * &pMcMainMCParticle, float &mcHitPurity,
        float &mcHitCompleteness, float &mcCollectionPlaneHitPurity, float &mcCollectionPlaneHitCompleteness) const;

    /**
     *  @brief  Recursively count the number of reconstructed downstream particles of a PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  numberOfParticles The number of downstream particles value to populate
     */
    void CountNumberOfDownstreamParticles(const ParticleFlowObject *const pPfo, unsigned &numberOfParticles) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline AnalysisAlgorithm::EnergyFromRangeData::EnergyFromRangeData(const float rangeMin, const float rangeMax,
    const float energy) noexcept :
    m_rangeMin(rangeMin),
    m_rangeMax(rangeMax),
    m_energy(energy)
{
}
} // namespace lar_physics_content

#endif // #ifndef LAR_LEE_ANALYSIS_ALGORITHM_H