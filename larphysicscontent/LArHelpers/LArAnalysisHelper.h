/**
 *  @file   larphysicscontent/LArHelpers/LArAnalysisHelper.h
 *
 *  @brief  Header file for the lar analysis helper class.
 *
 *  $Log: $
 */
#ifndef LAR_ANALYSIS_HELPER_H
#define LAR_ANALYSIS_HELPER_H 1

#include "Objects/CartesianVector.h"
#include "Pandora/Pandora.h"

#include <tuple>

namespace lar_physics_content
{
/**
 *  @brief  LArAnalysisHelper class
 */
class LArAnalysisHelper
{
public:
    /**
     *  @brief  Deleted copy constructor
     */
    LArAnalysisHelper(const LArAnalysisHelper &) = delete;

    /**
     *  @brief  Deleted move constructor
     */
    LArAnalysisHelper(LArAnalysisHelper &&) = delete;

    /**
     *  @brief  Deleted copy assignment operator
     */
    LArAnalysisHelper &operator=(const LArAnalysisHelper &) = delete;

    /**
     *  @brief  Deleted move assignment operator
     */
    LArAnalysisHelper &operator=(LArAnalysisHelper &&) = delete;

    /**
     *  @brief  Deleted destructor
     */
    ~LArAnalysisHelper() = delete;

    /**
     *  @brief  Get the minimum and maximum fiducial cut coordinates
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  fiducialCutLowMargins the fiducial cut low margins
     *  @param  fiducialCutHighMargins the fiducial cut high margins
     *
     *  @return the minimum and maximum fiducial cut coordinates
     */
    static std::tuple<pandora::CartesianVector, pandora::CartesianVector> GetFiducialCutCoordinates(const pandora::Pandora &pandoraInstance,
        const pandora::CartesianVector &fiducialCutLowMargins, const pandora::CartesianVector &fiducialCutHighMargins);

    /**
     *  @brief  Find out whether a given point lies in the fiducial region of the detector
     *
     *  @param  point the point
     *  @param  minCoordinates the minimum fiducial volume coordinates
     *  @param  maxCoordinates the maximum fiducial volume coordinates
     *
     *  @return whether the point is fiducial
     */
    static bool IsPointFiducial(const pandora::CartesianVector &point, const pandora::CartesianVector &minCoordinates,
        const pandora::CartesianVector &maxCoordinates);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ANALYSIS_HELPER_H
