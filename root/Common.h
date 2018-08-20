#ifndef LAR_ANALYSIS_ROOT_COMMON
#define LAR_ANALYSIS_ROOT_COMMON 1

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TNtuple.h"

#include <iostream>

#define TEXT_NORMAL "\033[0m"
#define TEXT_BOLD "\033[1m"

#define TEXT_RED "\033[0;31m"
#define TEXT_GREEN "\033[0;32m"
#define TEXT_MAGENTA "\033[0;35m"
#define TEXT_WHITE "\033[0;37m"

#define TEXT_RED_BOLD "\033[1;31m"
#define TEXT_GREEN_BOLD "\033[1;32m"
#define TEXT_MAGENTA_BOLD "\033[1;35m"
#define TEXT_WHITE_BOLD "\033[1;37m"

#define COUT(a) std::cout << a << TEXT_NORMAL << std::endl
#define CERR(a) std::cerr << TEXT_RED << a << TEXT_NORMAL << std::endl

//------------------------------------------------------------------------------------------------------------------------------------------

inline TNtuple *LoadNTupleFromFile(const char *const dataFilePath, const char *const nTupleName)
{
    std::cout << "Using data file: " << dataFilePath << std::endl;

    TFile *pFile = new TFile(dataFilePath, "READ");

    if (!pFile->IsOpen())
    {
        CERR("Failed to open file at " << dataFilePath);
        return NULL;
    }

    if (!pFile->GetListOfKeys()->Contains(nTupleName))
    {
        CERR("Data file at " << dataFilePath << " did not contain key '" << nTupleName << "'");
        return NULL;
    }

    return (TNtuple *)pFile->Get(nTupleName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TTree *LoadTreeFromFile(const char *const dataFilePath, const char *const treeName)
{
    std::cout << "Using data file: " << dataFilePath << std::endl;

    TFile *pFile = new TFile(dataFilePath, "READ");

    if (!pFile->IsOpen())
    {
        CERR("Failed to open file at " << dataFilePath);
        return NULL;
    }

    if (!pFile->GetListOfKeys()->Contains(treeName))
    {
        CERR("Data file at " << dataFilePath << " did not contain key '" << treeName << "'");
        return NULL;
    }

    return (TTree *)pFile->Get(treeName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

enum PLOT_TYPE
{
    HISTOGRAM,
    SCATTER,
    LINE,
    SAME
};

//------------------------------------------------------------------------------------------------------------------------------------------

struct PlotSettings2D
{
    char      title[100];
    char      xTitle[100];
    char      yTitle[100];
    float     xMin;
    float     xMax;
    int       xNumBins;
    float     yMin;
    float     yMax;
    int       yNumBins;
    PLOT_TYPE plotType;
    bool      useMaximumRange;
    bool      newCanvas;
    EColor    lineColor;
    EColor    pointColor;
    bool      trimSmallHistogramValues;
};

struct PlotSettings1D
{
    char      title[100];
    char      xTitle[100];
    float     xMin;
    float     xMax;
    int       xNumBins;
    PLOT_TYPE plotType;
    bool      useMaximumRange;
    bool      newCanvas;
    EColor    lineColor;
    EColor    pointColor;
    bool      trimSmallHistogramValues;
};

const struct PlotSettings2D g_defaultPlotSettings2D = {"", "", "", 0.f, 0.f, 80, 0.f, 0.f, 80, HISTOGRAM, true, true, kBlack, kBlack, true};
const struct PlotSettings1D g_defaultPlotSettings1D = {"", "", 0.f, 0.f, 80, HISTOGRAM, false, true, kBlack, kBlack, true};

//------------------------------------------------------------------------------------------------------------------------------------------

inline TCanvas *PlotNtuple1D(TNtuple *const pNtuple, const char *const xName, const char *const identifier, const PlotSettings1D &plotSettings)
{
    float     x          = 0.f;
    const int numEntries = pNtuple->GetEntries();

    if (numEntries == 0)
    {
        CERR("Number of ntuple entries was 0 so not drawing plot");
        return NULL;
    }

    pNtuple->SetBranchAddress(xName, &x);

    float xMax = plotSettings.xMax;
    float xMin = plotSettings.xMin;

    if (xMax == 0.f && xMin == 0.f && plotSettings.useMaximumRange)
    {
        bool firstEntry = true;

        for (int i = 0; i < numEntries; ++i)
        {
            pNtuple->GetEntry(i);

            if (firstEntry)
            {
                xMax = x;
                xMin = x;

                firstEntry = false;
                continue;
            }

            if (x > xMax)
                xMax = x;
            if (x < xMin)
                xMin = x;
        }
    }

    TH1F *pHistogram = new TH1F(identifier, identifier, plotSettings.xNumBins, xMin, xMax);

    for (int i = 0; i < numEntries; ++i)
    {
        pNtuple->GetEntry(i);
        pHistogram->Fill(x);
    }

    TCanvas *pCanvas = NULL;

    if (plotSettings.newCanvas)
        pCanvas = new TCanvas(identifier, identifier, 900, 600);

    pHistogram->SetXTitle((strcmp(plotSettings.xTitle, "") == 0) ? xName : plotSettings.xTitle);
    pHistogram->SetTitle((strcmp(plotSettings.title, "") == 0) ? identifier : plotSettings.title);
    pHistogram->SetMarkerStyle(6);
    pHistogram->SetStats(kTRUE);
    pHistogram->SetMarkerColor(plotSettings.pointColor);

    switch (plotSettings.plotType)
    {
        case SAME:
            pHistogram->Draw("SAME");
            break;
        case SCATTER:
            pHistogram->Draw("P");
            break;
        case LINE:
            pHistogram->Draw("L");
            break;
        case HISTOGRAM:
        default:
        {
            if (plotSettings.trimSmallHistogramValues)
                pHistogram->SetMinimum(static_cast<int>(std::round(pHistogram->GetEntries() / 10000.f)));

            pHistogram->Draw("COLZ");
            break;
        }
    }

    pNtuple->ResetBranchAddresses();
    return pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TCanvas *PlotNtuple2D(TNtuple *const pNtuple, const char *const xName, const char *const yName, const char *const identifier,
    const PlotSettings2D &plotSettings)
{
    float     x = 0.f, y = 0.f;
    const int numEntries = pNtuple->GetEntries();

    if (numEntries == 0)
    {
        CERR("Number of ntuple entries was 0 so not drawing plot");
        return NULL;
    }

    pNtuple->SetBranchAddress(xName, &x);
    pNtuple->SetBranchAddress(yName, &y);

    float xMax = plotSettings.xMax, yMax = plotSettings.yMax;
    float xMin = plotSettings.xMin, yMin = plotSettings.yMin;

    if (xMax == 0.f && yMax == 0.f && xMin == 0.f && yMin == 0.f && plotSettings.useMaximumRange)
    {
        bool firstEntry = true;

        for (int i = 0; i < numEntries; ++i)
        {
            pNtuple->GetEntry(i);

            if (firstEntry)
            {
                xMax = x;
                yMax = y;

                xMin = x;
                yMax = y;

                firstEntry = false;
                continue;
            }

            if (x > xMax)
                xMax = x;
            if (y > yMax)
                yMax = y;
            if (x < xMin)
                xMin = x;
            if (y < yMin)
                yMin = y;
        }
    }

    TH2F *pHistogram = new TH2F(identifier, identifier, plotSettings.xNumBins, xMin, xMax, plotSettings.yNumBins, yMin, yMax);

    for (int i = 0; i < numEntries; ++i)
    {
        pNtuple->GetEntry(i);
        pHistogram->Fill(x, y);
    }

    TCanvas *pCanvas = NULL;

    if (plotSettings.newCanvas)
        pCanvas = new TCanvas(identifier, identifier, 900, 600);

    pHistogram->SetXTitle((strcmp(plotSettings.xTitle, "") == 0) ? xName : plotSettings.xTitle);
    pHistogram->SetYTitle((strcmp(plotSettings.yTitle, "") == 0) ? yName : plotSettings.yTitle);
    pHistogram->SetTitle((strcmp(plotSettings.title, "") == 0) ? identifier : plotSettings.title);
    pHistogram->SetMarkerStyle(6);
    pHistogram->SetStats(kTRUE);
    pHistogram->SetMarkerColor(plotSettings.pointColor);

    switch (plotSettings.plotType)
    {
        case SAME:
            pHistogram->Draw("SAME");
            break;
        case SCATTER:
            pHistogram->Draw("P");
            break;
        case LINE:
            pHistogram->Draw("L");
            break;
        case HISTOGRAM:
        default:
        {
            if (plotSettings.trimSmallHistogramValues)
                pHistogram->SetMinimum(static_cast<int>(std::round(pHistogram->GetEntries() / 10000.f)));

            pHistogram->Draw("COLZ");
            break;
        }
    }

    pNtuple->ResetBranchAddresses();
    return pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TCanvas *PlotArrays2D(double xArray[], double yArray[], const int numEntries, const char *const identifier, const PlotSettings2D &plotSettings)
{
    float xMax = plotSettings.xMax, yMax = plotSettings.yMax;
    float xMin = plotSettings.xMin, yMin = plotSettings.yMin;

    if (xMax == 0.f && yMax == 0.f && xMin == 0.f && yMin == 0.f && plotSettings.useMaximumRange)
    {
        bool firstEntry = true;

        for (int i = 0; i < numEntries; ++i)
        {
            if (firstEntry)
            {
                xMax = xArray[i];
                yMax = yArray[i];

                xMin = xArray[i];
                xMax = yArray[i];

                firstEntry = false;
                continue;
            }

            if (xArray[i] > xMax)
                xMax = xArray[i];
            if (yArray[i] > yMax)
                yMax = yArray[i];
            if (xArray[i] < xMin)
                xMin = xArray[i];
            if (yArray[i] < yMin)
                yMin = yArray[i];
        }
    }

    TGraph *pGraph = new TGraph(numEntries, xArray, yArray);

    TCanvas *pCanvas = NULL;

    if (plotSettings.newCanvas)
        pCanvas = new TCanvas(identifier, identifier, 900, 600);

    pGraph->GetXaxis()->SetTitle((strcmp(plotSettings.xTitle, "") == 0) ? "x" : plotSettings.xTitle);
    pGraph->GetXaxis()->SetRangeUser(plotSettings.xMin, plotSettings.xMax);
    pGraph->GetYaxis()->SetTitle((strcmp(plotSettings.yTitle, "") == 0) ? "y" : plotSettings.yTitle);
    pGraph->GetYaxis()->SetRangeUser(plotSettings.yMin, plotSettings.yMax);
    pGraph->SetTitle((strcmp(plotSettings.title, "") == 0) ? identifier : plotSettings.title);
    pGraph->SetMarkerStyle(6);
    pGraph->SetLineColor(plotSettings.lineColor);

    if (plotSettings.newCanvas)
    {
        switch (plotSettings.plotType)
        {
            case SAME:
                pGraph->Draw("SAME");
                break;
            case SCATTER:
                pGraph->Draw("AP");
                break;
            case LINE:
            default:
                pGraph->Draw("AL");
                break;
        }
    }

    else
    {
        switch (plotSettings.plotType)
        {
            case SAME:
                pGraph->Draw("SAME");
                break;
            case SCATTER:
                pGraph->Draw("P");
                break;
            case LINE:
            default:
                pGraph->Draw("L");
                break;
        }
    }

    return pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TCanvas *PlotFunction2D(float (*pFunction)(float, va_list), const float xArgMin, const float xArgMax, const int nSteps,
    const char *const identifier, const PlotSettings2D plotSettings, ...)
{
    if (nSteps < 2)
    {
        CERR("Number of steps must be 2 or more");
        return NULL;
    }

    const float stepSize = (xArgMax - xArgMin) / (float)(nSteps - 1);

    float xMax = plotSettings.xMax, yMax = plotSettings.yMax;
    float xMin = plotSettings.xMin, yMin = plotSettings.yMin;

    // Populate x and y with the values.
    float x[10000], y[10000];

    for (int i = 0; i < nSteps; ++i)
    {
        x[i] = xArgMin + (float)i * stepSize;

        va_list argptr;
        va_start(argptr, plotSettings);
        y[i] = (*pFunction)(x[i], argptr);
        va_end(argptr);
    }

    if (xMax == 0.f && yMax == 0.f && xMin == 0.f && yMin == 0.f && plotSettings.useMaximumRange)
    {
        bool firstEntry = true;

        for (int i = 0; i < nSteps; ++i)
        {
            if (firstEntry)
            {
                xMax = x[i];
                yMax = y[i];

                xMin = x[i];
                xMax = y[i];

                firstEntry = false;
                continue;
            }

            if (x[i] > xMax)
                xMax = x[i];
            if (y[i] > yMax)
                yMax = y[i];
            if (x[i] < xMin)
                xMin = x[i];
            if (y[i] < yMin)
                yMin = y[i];
        }
    }

    TGraph *pGraph = new TGraph(nSteps, x, y);

    TCanvas *pCanvas = NULL;

    if (plotSettings.newCanvas)
        pCanvas = new TCanvas(identifier, identifier, 900, 600);

    pGraph->GetXaxis()->SetTitle((strcmp(plotSettings.xTitle, "") == 0) ? "x" : plotSettings.xTitle);
    pGraph->GetYaxis()->SetTitle((strcmp(plotSettings.yTitle, "") == 0) ? "y" : plotSettings.yTitle);
    pGraph->SetTitle((strcmp(plotSettings.title, "") == 0) ? identifier : plotSettings.title);
    pGraph->SetMarkerStyle(6);
    pGraph->SetLineColor(plotSettings.lineColor);

    if (plotSettings.newCanvas)
    {
        switch (plotSettings.plotType)
        {
            case SCATTER:
                pGraph->Draw("AP");
                break;
            case LINE:
            default:
                pGraph->Draw("AL");
                break;
        }
    }

    else
    {
        switch (plotSettings.plotType)
        {
            case SCATTER:
                pGraph->Draw("P");
                break;
            case LINE:
            default:
                pGraph->Draw("L");
                break;
        }
    }

    return pCanvas;
}

#endif // #ifndef LAR_ANALYSIS_ROOT_COMMON
