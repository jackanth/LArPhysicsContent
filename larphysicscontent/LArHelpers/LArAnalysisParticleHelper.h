/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisParticleHelper.h
 *
 *  @brief  Header file for the lar analysis particle helper class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_PARTICLE_HELPER_H
#define LAR_ANALYSIS_PARTICLE_HELPER_H 1

#include "Objects/ParticleFlowObject.h"
#include "Pandora/Pandora.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larphysicscontent/LArObjects/LArAnalysisParticle.h"
#include "larphysicscontent/LArObjects/LArFittedTrackInfo.h"
#include "larphysicscontent/LArObjects/LArTrackHitValue.h"

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
    using FittedTrackInfoMap = std::unordered_map<const ParticleFlowObject *, LArFittedTrackInfo>; ///< Alias for a map to LArFittedTrackInfo objects

    /**
     *  @brief  PfoMcInfo class
     */
    class PfoMcInfo
    {
    public:
        /**
         *  @brief  Default constructor
         */
        PfoMcInfo() noexcept;

        const MCParticle *m_pMCParticle;
        float m_mcEnergy;
        float m_mcKineticEnergy;
        float m_mcMass;
        LArAnalysisParticle::TypeTree m_mcTypeTree;
        LArAnalysisParticle::TYPE m_mcType;
        CartesianVector m_mcVertexPosition;
        CartesianVector m_mcMomentum;
        CartesianVector m_mcDirectionCosines;
        int m_mcPdgCode;
        float m_mcContainmentFraction;
        bool m_mcIsShower;
        bool m_mcIsVertexFiducial;
        bool m_mcIsContained;
        bool m_mcIsProton;
        bool m_mcIsCosmicRay;
        bool m_mcIsPionOrMuon;
    };

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
    static void GetFiducialCutParameters(const Pandora &pandoraInstance, const CartesianVector &fiducialCutLowMargins,
        const CartesianVector &fiducialCutHighMargins, CartesianVector &minCoordinates, CartesianVector &maxCoordinates);

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
    static CartesianVector GetFittedDirectionAtPosition(
        const ThreeDSlidingFitResult &trackFit, const CartesianVector &position, const bool pointTowardsMiddle);

    /**
     *  @brief  Get the fraction of a PFO's hits that lie in the fiducial region of the detector
     *
     *  @param  pPfo address of the PFO
     *  @param  minCoordinates the minimum fiducial volume coordinates
     *  @param  maxCoordinates the maximum fiducial volume coordinates
     *
     *  @return the fraction of fiducial hits
     */
    static float GetFractionOfFiducialHits(
        const ParticleFlowObject *const pPfo, const CartesianVector &minCoordinates, const CartesianVector &maxCoordinates);

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

private:
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

inline LArAnalysisParticleHelper::PfoMcInfo::PfoMcInfo() noexcept :
    m_pMCParticle(nullptr),
    m_mcEnergy(0.f),
    m_mcKineticEnergy(0.f),
    m_mcMass(0.f),
    m_mcTypeTree(),
    m_mcType(LArAnalysisParticle::TYPE::UNKNOWN),
    m_mcVertexPosition(0.f, 0.f, 0.f),
    m_mcMomentum(0.f, 0.f, 0.f),
    m_mcDirectionCosines(0.f, 0.f, 0.f),
    m_mcPdgCode(0),
    m_mcContainmentFraction(0.f),
    m_mcIsShower(false),
    m_mcIsVertexFiducial(false),
    m_mcIsContained(false),
    m_mcIsProton(false),
    m_mcIsCosmicRay(false),
    m_mcIsPionOrMuon(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string LArAnalysisParticleHelper::TypeTreeAsString(const LArAnalysisParticle::TypeTree &typeTree)
{
    return LArAnalysisParticleHelper::TypeTreeAsStringImpl(typeTree, false);
}

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_PARTICLE_HELPER_H
