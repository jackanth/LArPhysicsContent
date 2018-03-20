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

#include "LArAnalysisParticle.h"

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
 */
class LArAnalysisParticleHelper
{
public:
    using LArTrackHitEnergyMap = std::unordered_map<const ParticleFlowObject *, LArTrackHitEnergy::Vector>; ///< Alias for a map from PFOs to their track hit energies
    using TrackFitMap          = std::unordered_map<const ParticleFlowObject *, ThreeDSlidingFitResult>;    ///< Alias for map from PFOs to track fits.
    using HitProjectionPair    = std::pair<const CaloHit *, float>; ///< Alias for a map from CaloHits to their projected track coordinate
    using HitProjectionVector  = std::vector<HitProjectionPair>;    ///< Alias for a vector of hit projection pairs

    /**
     *  @brief  Get the particle type as a string
     *
     *  @param  type the type
     *
     *  @return the string type
     */
    static std::string TypeAsString(const LArAnalysisParticle::TYPE type);

    /**
     *  @brief  Get the particle type tree as a string
     *
     *  @param  typeTree the typeTree
     *
     *  @return the string type
     */
    static std::string TypeTreeAsString(const LArAnalysisParticle::TypeTree &typeTree);

    /**
     *  @brief  Get the fiducial cut parameters
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  fiducialCutLowXMargin the fiducial cut low-x margin
     *  @param  fiducialCutHighXMargin the fiducial cut high-x margin
     *  @param  fiducialCutLowYMargin the fiducial cut low-y margin
     *  @param  fiducialCutHighYMargin the fiducial cut high-y margin
     *  @param  fiducialCutLowZMargin the fiducial cut low-z margin
     *  @param  fiducialCutHighZMargin the fiducial cut high-z margin
     *  @param  minCoordinates the minimum fiducial coordinates (to populate)
     *  @param  maxCoordinates the maximum fiducial coordinates (to populate)
     */
    static void GetFiducialCutParameters(const Pandora &pandoraInstance, const float fiducialCutLowXMargin, const float fiducialCutHighXMargin,
        const float fiducialCutLowYMargin, const float fiducialCutHighYMargin, const float fiducialCutLowZMargin, const float fiducialCutHighZMargin,
        CartesianVector &minCoordinates, CartesianVector &maxCoordinates);

    /**
     *  @brief  Recurse through the PFO hierarchy and append the track fit map
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pPfo address of the current PFO
     *  @param  trackFitMap the track fit map to append
     *  @param  slidingFitWindow the sliding fit window size
     */
    static void RecursivelyAppendTrackFitMap(const Pandora &pandoraInstance, const ParticleFlowObject *const pPfo, TrackFitMap &trackFitMap,
        const unsigned int slidingFitWindow);

    /**
     *  @brief  Perform a 3D sliding track fit for a given PFO
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pPfo address of the PFO
     *  @param  slidingFitWindow the sliding fit window size
     *
     *  @return the 3D track fit
     */
    static ThreeDSlidingFitResult PerformSlidingTrackFit(const Pandora &pandoraInstance, const ParticleFlowObject *const pPfo,
        const unsigned int slidingFitWindow);

    /**
     *  @brief  Get the list of reco neutrinos
     *
     *  @param  algorithm the calling algorithm
     *  @param  pfoListName the PFO list name
     *
     *  @return the list of reco neutrinos
     */
    static PfoList GetRecoNeutrinoList(const Algorithm &algorithm, const std::string &pfoListName);

    /**
     *  @brief  Get the 3D distance from a CaloHit
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pCaloHit address of the CaloHit
     *  @param  trackFit the 3D track fit to which the hit belongs
     *
     *  @return the 3D distance
     */
    static float CaloHitToThreeDDistance(const Pandora &pandoraInstance, const CaloHit *const pCaloHit, const ThreeDSlidingFitResult &trackFit);

    /**
     *  @brief  Get the 3D distance from a CaloHit
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  hitWidth the CaloHit's width
     *  @param  fitDirection the direction of the fit at the CaloHit's position
     *
     *  @return the 3D distance
     */
    static float CaloHitToThreeDDistance(const Pandora &pandoraInstance, const float hitWidth, const CartesianVector &fitDirection);

    /**
     *  @brief  Turn a direction in polar and azimuthal angles
     *
     *  @param  direction the Cartesian direction
     *
     *  @return the polar and azimuthal angles
     */
    static std::tuple<float, float> GetPolarAnglesFromDirection(const CartesianVector &direction);

    /**
     *  @brief  Get all the hits of a given type from a PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  hitType the hit type
     *  @param  recurseOverDaughters whether to recurse through the hierarchy
     *
     *  @return the list of hits
     */
    static CaloHitList GetHitsOfType(const ParticleFlowObject *const pPfo, const HitType hitType, const bool recurseOverDaughters);

    /**
     *  @brief  Get the fitted direction at a given position
     *
     *  @param  trackFit the fit object
     *  @param  position the position
     *  @param  pointTowardsMiddle whether to ensure that the direction points towards the middle of the fit, otherwise it will point in the
     *          negative y direction
     *
     *  @return the direction
     */
    static CartesianVector GetFittedDirectionAtPosition(const ThreeDSlidingFitResult &trackFit, const CartesianVector &position,
        const bool pointTowardsMiddle);

    /**
     *  @brief  Get the fraction of a PFO's hits that lie in the fiducial region of the detector
     *
     *  @param  pPfo address of the PFO
     *  @param  minCoordinates the minimum fiducial volume coordinates
     *  @param  maxCoordinates the maximum fiducial volume coordinates
     *
     *  @return the fraction of fiducial hits
     */
    static float GetFractionOfFiducialHits(const ParticleFlowObject *const pPfo, const CartesianVector &minCoordinates,
        const CartesianVector &maxCoordinates);

    /**
     *  @brief  Find out whether a given point lies in the fiducial region of the detector
     *
     *  @param  point the point
     *  @param  minCoordinates the minimum fiducial volume coordinates
     *  @param  maxCoordinates the maximum fiducial volume coordinates
     *
     *  @return whether the point is fiducial
     */
    static bool IsPointFiducial(const CartesianVector &point, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates);

    /**
     *  @brief  Get the range of a track-like PFO
     *
     *  @param  pPfo address of the PFO
     *  @param  trackFit the track fit for the PFO
     *
     *  @return the range
     */
    static float GetParticleRange(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &trackFit);

    /**
     *  @brief  Order a set of CaloHits by their projection onto a track fit and return the ordered projections
     *
     *  @param  caloHitList the list of CaloHits
     *  @param  trackFit the track fit object
     *
     *  @return the ordered list of hit projections
     */
    static HitProjectionVector OrderHitsByProjectionOnToTrackFit(const CaloHitList &caloHitList, const ThreeDSlidingFitResult &trackFit);

    /**
     *  @brief  Write a TNtuple object to a file
     *
     *  @param  pNtuple address of the TNtuple object
     *  @param  fileName the path to the new ROOT file
     *  @param  verboseMode whether to save in verbose mode
     */
    static void WriteNTuple(TNtuple *const pNtuple, const std::string &fileName, const bool verboseMode);

    /**
     *  @brief  Load a TNtuple object from a file
     *
     *  @param  filePath the path to the ROOT file
     *  @param  nTupleName the name of the TNtuple object to load
     *
     *  @return address of the loaded TNtuple object
     */
    static TNtuple * LoadNTupleFromFile(const std::string &filePath, const std::string &nTupleName);

    /**
     *  @brief  Get a majority-wins main MC particle for a PFO
     *
     *  @param  pPfo address of the PFO
     *
     *  @return address of the MC particle
     */
    static const MCParticle *GetMainMCParticle(const ParticleFlowObject *const pPfo);

    /**
     *  @brief  Find out whether a PFO is a neutrino
     *
     *  @param  pPfo address of the PFO
     *
     *  @return whether it is a neutrino
     */
    static bool IsNeutrino(const ParticleFlowObject *const pPfo);

    /**
     *  @brief  Find out whether a PFO is a cosmic ray
     *
     *  @param  pPfo address of the PFO
     *
     *  @return whether it is a cosmic ray
     */
    static bool IsCosmicRay(const ParticleFlowObject *const pPfo);

    /**
     *  @brief  Find out whether a PFO is a primary daughter of the neutrino
     *
     *  @param  pPfo address of the PFO
     *
     *  @return whether it is a primary daughter
     */
    static bool IsPrimaryNeutrinoDaughter(const ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get MC information for a given MC particle
     *
     *  @param  pMCParticle
     *  @param  mcEnergy the energy of the MC particle (to populate)
     *  @param  mcKineticEnergy the kinetic energy of the MC particle (to populate)
     *  @param  mcMass the mass of the MC particle (to populate)
     *  @param  typeTree the type tree of the MC particle (to populate)
     *  @param  mcType the type of the MC particle (to populate)
     *  @param  mcVertexPosition the vertex position of the MC particle (to populate)
     *  @param  mcMomentum the momentum of the MC particle (to populate)
     *  @param  mcPdgCode the PDG code of the MC particle (to populate)
     *  @param  mcContainmentFraction the containment fraction of the MC particle (to populate)
     *  @param  minCoordinates the minimum fiducial coordinates of the detector
     *  @param  maxCoordinates the maximum fiducial coordinates of the detector
     *
     *  @return success
     */
    static bool GetMcInformation(const MCParticle *const pMCParticle, float &mcEnergy, float &mcKineticEnergy, float &mcMass,
        LArAnalysisParticle::TypeTree &typeTree, LArAnalysisParticle::TYPE &mcType, CartesianVector &mcVertexPosition, CartesianVector &mcMomentum,
        int &mcPdgCode, float &mcContainmentFraction, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates);

    /**
     *  @brief  Calculate hit-counting-based purity and completeness values
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticle address of the MC particle
     *  @param  pCaloHitList address of the CaloHit list
     *  @param  isNeutrino whether the PFO is a beam neutrino
     *  @param  hitPurity the purity (to populate)
     *  @param  hitCompleteness the completeness (to populate)
     *  @param  mcCollectionPlaneHitPurity the purity in the collection plane (to populate)
     *  @param  mcCollectionPlaneHitCompleteness the completeness in the collection plane (to populate)
     */
    static void CalculateHitPurityAndCompleteness(const ParticleFlowObject *const pPfo, const MCParticle *const pMCParticle,
        const CaloHitList *const pCaloHitList, const bool isNeutrino, float &hitPurity, float &hitCompleteness, float &mcCollectionPlaneHitPurity,
        float &mcCollectionPlaneHitCompleteness);

private:
    /**
     *  @brief  Sort reco neutrinos
     *
     *  @param  pLhs address of the LHS neutrino
     *  @param  pRhs address of the RHS neutrino
     *
     *  @return whether LHS < RHS
     */
    static bool SortRecoNeutrinos(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs);

    /**
     *  @brief  Get a 3D distance from a cell
     *
     *  @param  hitWidth the hit width
     *  @param  wirePitch the wire pitch
     *  @param  fitDirection the fit direction at the hit
     *
     *  @return the 3D distance
     */
    static float CellToThreeDDistance(const float hitWidth, const float wirePitch, const CartesianVector &fitDirection);

    /**
     *  @brief  Recurse through the MC particle hierarchy and add up the escaped energy for the containment fraction calculation
     *
     *  @param  pCurrentMCParticle address of the current MC particle
     *  @param  escapedEnergy the escaped energy (to populate)
     *  @param  totalEnergy the total energy (to populate)
     *  @param  minCoordinates the minimum fiducial volume coordinates
     *  @param  maxCoordinates the maximum fiducial volume coordinates
     */
    static void RecursivelyAddEscapedEnergy(const MCParticle *const pCurrentMCParticle, float &escapedEnergy, float &totalEnergy,
        const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates);

    /**
     *  @brief  Adjust the line equation mu values for a given detector face constraint
     *
     *  @param  planePoint a point on the detector face
     *  @param  planeNormal the normal to the detector face
     *  @param  vertexPosition the position of the MC vertex
     *  @param  originalDisplacementVector the original MC vertex-to-endpoint displacement vector
     *  @param  muMin the minimum mu value to be adjusted
     *  @param  muMax the maximum mu value to be adjusted
     *  @param  forceZeroContainment whether this face constraint tells us that the particle definitely has zero containment (to populate)
     */
    static void AdjustMusForContainmentFraction(const CartesianVector &planePoint, const CartesianVector &planeNormal,
        const CartesianVector &vertexPosition, const CartesianVector &originalDisplacementVector, float &muMin, float &muMax,
        bool &forceZeroContainment);

    /**
     *  @brief  Create a type tree for an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *  @param  typeTree the type tree (to populate)
     *
     *  @return success
     */
    static bool CreateMcTypeTree(const MCParticle *const pMCParticle, LArAnalysisParticle::TypeTree &typeTree);

    /**
     *  @brief  Get the type of an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the particle type
     */
    static LArAnalysisParticle::TYPE GetMcParticleType(const MCParticle *const pMCParticle);

    /**
     *  @brief  Calculate the hit-counting-based purity and completeness for a PFO
     *
     *  @param  pfoAssociatedCaloHits the hits associated with the PFO
     *  @param  pMCParticle address of the PFO's main MC particle
     *  @param  pCaloHitList address of the CaloHit list
     *  @param  isNeutrino whether the PFO is a beam neutrino
     *  @param  hitPurity the hit purity (to populate)
     *  @param  hitCompleteness the hit completeness (to populate)
     *  @param  useCollectionPlaneOnly whether to use only the collection plane
     */
    static void CalculateHitPurityAndCompleteness(const CaloHitList &pfoAssociatedCaloHits, const MCParticle *const pMCParticle,
        const CaloHitList *const pCaloHitList, const bool isNeutrino, float &hitPurity, float &hitCompleteness,
        const bool useCollectionPlaneOnly);

    /**
     *  @brief  Get the particle type tree as a string (implementation)
     *
     *  @param  typeTree the typeTree
     *  @param  printTrailingDelimiter whether to print the trailing delimiter
     *
     *  @return the string type
     */
    static std::string TypeTreeAsStringImpl(const LArAnalysisParticle::TypeTree &typeTree, const bool printTrailingDelimiter);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string LArAnalysisParticleHelper::TypeTreeAsString(const LArAnalysisParticle::TypeTree &typeTree)
{
    return LArAnalysisParticleHelper::TypeTreeAsStringImpl(typeTree, false);
}

} // namespace lar_physics_content

#endif // #ifndef LEE_ANALYSIS_HELPER_H
