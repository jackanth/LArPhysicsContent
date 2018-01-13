#include "Common.h"

//------------------------------------------------------------------------------------------------------------------------------------------

enum TYPE : int
{
    UNKNOWN   = 0,
    TRACK     = 1,
    SHOWER    = 2,
    PROTON    = 3,
    PION_MUON = 4
};

//------------------------------------------------------------------------------------------------------------------------------------------

void ReadPandoraNtuple(const char *const inputFilePath)
{
    TTree *const pTree = LoadTreeFromFile(inputFilePath, "PandoraTree");
    
    float nu_mc_VertexX      = 0.f;
    float nu_mc_VertexY      = 0.f;
    float tmp_nu_mc_SinTheta = 0.f;

    pTree->SetBranchAddress("nu_mc_VertexX",      &nu_mc_VertexX);
    pTree->SetBranchAddress("nu_mc_VertexY",      &nu_mc_VertexY);
    pTree->SetBranchAddress("tmp_nu_mc_SinTheta", &tmp_nu_mc_SinTheta);
    
    // Define some new ntuples to plot later.
    TFile *const pFile = new TFile("tmp.root", "RECREATE");
    TNtuple *const pSinThetaPlotNtuple = new TNtuple("SinThetaPlotNtuple", "SinThetaPlotNtuple", "X:Y:SinTheta:Theta");
    
    // Loop over the input primaries.
    const int nEntries = pTree->GetEntries();
    
    for (int i = 0; i < nEntries; ++i)
    {
        pTree->GetEntry(i);
        
        if (tmp_nu_mc_SinTheta != 0.f)
        {
            const float degrees = std::asin(tmp_nu_mc_SinTheta) * 180.f / M_PI;
            pSinThetaPlotNtuple->Fill(nu_mc_VertexX, nu_mc_VertexY, tmp_nu_mc_SinTheta, degrees);
        }
    }
    
    TCanvas *pCanvas1 = new TCanvas("1", "1", 900, 600);
    pSinThetaPlotNtuple->Draw("X:Y>>SinThetaDistribution", "SinTheta", "colz");
    TGraph *graph1 = (TGraph*)gDirectory->Get("SinThetaDistribution");
    graph1->SetTitle("Neutrino vertex sin theta distribution;X coordinate (cm);Y coordinate (cm)");
    
    TCanvas *pCanvas2 = new TCanvas("2", "2", 900, 600);
    pSinThetaPlotNtuple->Draw("X:Y>>ThetaDistribution", "Theta", "colz");
    TGraph *graph2 = (TGraph*)gDirectory->Get("ThetaDistribution");
    graph2->SetTitle("Neutrino vertex theta distribution (degrees);X coordinate (cm);Y coordinate (cm)");
 
 /*
    
    {
        TString identifier = "SinThetaHistogram";
        struct PlotSettings1D plotSettings = g_defaultPlotSettings1D;
        strcpy(plotSettings.title, "Sin theta");
        strcpy(plotSettings.xTitle, "Sin theta");
        //plotSettings.xMax = 1.2f;
        PlotNtuple1D(pSinThetaPlotNtuple, "SinTheta", identifier, plotSettings)->SaveAs(identifier.Append(".png"));
    }
    
    {
        TString identifier = "XHistogram";
        struct PlotSettings1D plotSettings = g_defaultPlotSettings1D;
        strcpy(plotSettings.title, "X coordinate");
        strcpy(plotSettings.xTitle, "X coordinate");
        //plotSettings.xMax = 1.2f;
        PlotNtuple1D(pSinThetaPlotNtuple, "X", identifier, plotSettings)->SaveAs(identifier.Append(".png"));
    }
    
    {
        TString identifier = "YHistogram";
        struct PlotSettings1D plotSettings = g_defaultPlotSettings1D;
        strcpy(plotSettings.title, "Y coordinate");
        strcpy(plotSettings.xTitle, "Y coordinate");
        //plotSettings.xMax = 1.2f;
        PlotNtuple1D(pSinThetaPlotNtuple, "Y", identifier, plotSettings)->SaveAs(identifier.Append(".png"));
    }
    
    {
        TString identifier = "RHistogram";
        struct PlotSettings1D plotSettings = g_defaultPlotSettings1D;
        strcpy(plotSettings.title, "R coordinate");
        strcpy(plotSettings.xTitle, "R coordinate");
        //plotSettings.xMax = 1.2f;
        PlotNtuple1D(pSinThetaPlotNtuple, "R", identifier, plotSettings)->SaveAs(identifier.Append(".png"));
    }
     
    {
        TString identifier = "SinTheta2";
        struct PlotSettings2D plotSettings = g_defaultPlotSettings2D;
        strcpy(plotSettings.title, "Energy from charge discrepancy: fiducial particles with >0 collection plane hits");
        strcpy(plotSettings.xTitle, "True energy (GeV)");
        strcpy(plotSettings.yTitle, "Fractional energy discrepancy");
//        plotSettings.xMax = 1.2f;
//        plotSettings.yMin = -1.f;
//        plotSettings.yMax = 5.f;
//        plotSettings.xNumBins = 100;
//        plotSettings.yNumBins = 100;
        plotSettings.trimSmallHistogramValues = false;
        PlotNtuple2D(pSinThetaPlotNtuple, "SinTheta", "R", identifier, plotSettings)->SaveAs(identifier.Append(".png"));
    }
    
    {
        TString identifier = "SinTheta3";
        struct PlotSettings2D plotSettings = g_defaultPlotSettings2D;
        strcpy(plotSettings.title, "Energy from charge discrepancy: fiducial particles with >0 collection plane hits");
        strcpy(plotSettings.xTitle, "True energy (GeV)");
        strcpy(plotSettings.yTitle, "Fractional energy discrepancy");
//        plotSettings.xMax = 1.2f;
//        plotSettings.yMin = -1.f;
//        plotSettings.yMax = 5.f;
//        plotSettings.xNumBins = 100;
//        plotSettings.yNumBins = 100;
        plotSettings.trimSmallHistogramValues = false;
        PlotNtuple2D(pSinThetaPlotNtuple, "X", "Y", identifier, plotSettings)->SaveAs(identifier.Append(".png"));
    }
     */
}
