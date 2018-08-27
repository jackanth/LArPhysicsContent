/**
 *  @file   larphysicscontent/LArHelpers/LArRootHelper.h
 *
 *  @brief  Header file for the lar ROOT helper class.
 *
 *  $Log: $
 */
#ifndef LAR_ROOT_HELPER_H
#define LAR_ROOT_HELPER_H 1

#include "Pandora/Pandora.h"
#include "PandoraMonitoringApi.h"
#include "larphysicscontent/LArObjects/LArRootRegistry.h"

#include "TH1F.h"
#include "TH2F.h"

#include <vector>

namespace lar_physics_content
{
/**
 *  @brief  LArRootHelper class
 */
class LArRootHelper
{
public:
    using FloatVector = std::vector<float>; ///< Alias for a vector of floats

    struct PlotOptions
    {
        PlotOptions() noexcept;

        std::string m_name;
        std::string m_title;
        std::string m_fileName;
        unsigned    m_numXBins;
        unsigned    m_numYBins;
        float       m_xMin;
        float       m_xMax;
        float       m_yMin;
        float       m_yMax;
        Color_t     m_markerColour;
        std::string m_xTitle;
        std::string m_yTitle;
    };

    /**
     *  @brief  Deleted copy constructor
     */
    LArRootHelper(const LArRootHelper &) = delete;

    /**
     *  @brief  Deleted move constructor
     */
    LArRootHelper(LArRootHelper &&) = delete;

    /**
     *  @brief  Deleted copy assignment operator
     */
    LArRootHelper &operator=(const LArRootHelper &) = delete;

    /**
     *  @brief  Deleted move assignment operator
     */
    LArRootHelper &operator=(LArRootHelper &&) = delete;

    /**
     *  @brief  Deleted destructor
     */
    ~LArRootHelper() = delete;

    static TH1F *CreateOneDHistogram(const LArRootRegistry &registry, const FloatVector &vector, const PlotOptions &options);

    static TH2F *CreateTwoDHistogram(const LArRootRegistry &registry, const FloatVector &xVector, const FloatVector &yVector, const PlotOptions &options);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ROOT_HELPER_H
