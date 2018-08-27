/**
 *  @file   larphysicscontent/LArHelpers/LArRootHelper.cxx
 *
 *  @brief  Implementation of the lar ROOT helper class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArHelpers/LArRootHelper.h"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace pandora;

namespace lar_physics_content
{

LArRootHelper::PlotOptions::PlotOptions() noexcept :
    m_name("PandoraPlot"),
    m_title("Pandora plot"),
    m_fileName("plots.root"),
    m_numXBins(80U),
    m_numYBins(80U),
    m_xMin(0.f),
    m_xMax(0.f),
    m_yMin(0.f),
    m_yMax(0.f),
    m_markerColour(kBlack),
    m_xTitle(),
    m_yTitle()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TH1F *LArRootHelper::CreateOneDHistogram(const LArRootRegistry &registry, const FloatVector &vector, const PlotOptions &options)
{
    if (vector.empty())
    {
        std::cerr << "LArRootHelper: Could not create histogram because there were no datapoints" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    TH1F *pHistogram =
        registry.CreateWithUniqueName<TH1F>(options.m_name, options.m_title.c_str(), options.m_numXBins, options.m_xMin, options.m_xMax);

    pHistogram->SetXTitle(options.m_xTitle.c_str());
    pHistogram->SetYTitle(options.m_yTitle.c_str());
    pHistogram->SetMarkerStyle(6);
    pHistogram->SetMarkerColor(options.m_markerColour);

    for (const float value : vector)
        pHistogram->Fill(value);

    return pHistogram;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TH2F *LArRootHelper::CreateTwoDHistogram(const LArRootRegistry &registry, const FloatVector &xVector, const FloatVector &yVector, const PlotOptions &options)
{
    if (xVector.empty())
    {
        std::cerr << "LArRootHelper: Could not create histogram because there were no datapoints" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    if (xVector.size() != yVector.size())
    {
        std::cerr << "LArRootHelper: x-values vector and y-values vector were not the same size" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    TH2F *pHistogram = registry.CreateWithUniqueName<TH2F>(options.m_name, options.m_title.c_str(), options.m_numXBins, options.m_xMin,
        options.m_xMax, options.m_numYBins, options.m_yMin, options.m_yMax);

    pHistogram->SetXTitle(options.m_xTitle.c_str());
    pHistogram->SetYTitle(options.m_yTitle.c_str());
    pHistogram->SetMarkerStyle(6);
    pHistogram->SetMarkerColor(options.m_markerColour);

    for (std::size_t i = 0UL; i < xVector.size(); ++i)
        pHistogram->Fill(xVector.at(i), yVector.at(i));

    return pHistogram;
}

} // namespace lar_physics_content
