/**
 *  @file   larphysicscontent/McInfoTool.h
 *
 *  @brief  Header file for the track hit energy tool class.
 *
 *  $Log: $
 */
#ifndef LAR_MC_INFO_TOOL
#define LAR_MC_INFO_TOOL 1

#include "Pandora/AlgorithmTool.h"

#include "larphysicscontent/LArAnalysisParticleHelper.h"

using namespace pandora;

namespace lar_physics_content
{

/**
 *  @brief  McInfoTool class
 */
class McInfoTool : public AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    McInfoTool();

    bool Run(const Algorithm *const pAlgorithm, const MCParticle *const pMCParticle, LArAnalysisParticleHelper::PfoMcInfo &pfoMcInfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    CartesianVector    m_fiducialCutLowMargins;              ///< The low fiducial margins
    CartesianVector    m_fiducialCutHighMargins;             ///< The high fiducial cut margins
    CartesianVector              m_minCoordinates;                       ///< The detector's minimum fiducial coordinates
    CartesianVector              m_maxCoordinates;                       ///< The detector's maximum fiducial coordinates
    float              m_mcContainmentFractionLowerBound;    ///< The lower containment fraction bound for MC containment

    /**
     *  @brief  Get MC information for a given MC particle
     *
     *  @param  pMCParticle
     *  @param  minCoordinates the minimum fiducial coordinates of the detector
     *  @param  maxCoordinates the maximum fiducial coordinates of the detector
     *
     *  @return the mc information
     */
    LArAnalysisParticleHelper::PfoMcInfo GetMcInformation(const MCParticle *const pMCParticle) const;

      /**
     *  @brief  Recurse through the MC particle hierarchy and add up the escaped energy for the containment fraction calculation
     *
     *  @param  pCurrentMCParticle address of the current MC particle
     *  @param  escapedEnergy the escaped energy (to populate)
     *  @param  totalEnergy the total energy (to populate)
     *  @param  minCoordinates the minimum fiducial volume coordinates
     *  @param  maxCoordinates the maximum fiducial volume coordinates
     */
    void RecursivelyAddEscapedEnergy(const MCParticle *const pCurrentMCParticle, float &escapedEnergy, float &totalEnergy) const;

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
    void AdjustMusForContainmentFraction(const CartesianVector &planePoint, const CartesianVector &planeNormal,
        const CartesianVector &vertexPosition, const CartesianVector &originalDisplacementVector, float &muMin, float &muMax,
        bool &forceZeroContainment) const;

    /**
     *  @brief  Create a type tree for an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *  @param  typeTree the type tree (to populate)
     *
     *  @return success
     */
    bool CreateMcTypeTree(const MCParticle *const pMCParticle, LArAnalysisParticle::TypeTree &typeTree) const;

    /**
     *  @brief  Get the type of an MC particle
     *
     *  @param  pMCParticle address of the MC particle
     *
     *  @return the particle type
     */
    LArAnalysisParticle::TYPE GetMcParticleType(const MCParticle *const pMCParticle) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_MC_INFO_TOOL
