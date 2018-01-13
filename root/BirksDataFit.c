#include "Common.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TFitter.h"
#include "TMinuit.h"

#include <vector>

std::vector<float> g_trueEnergy;
std::vector<float> g_nonBirksHitSummeddQ;
std::vector<std::vector<float>> g_birksHitdQs;
std::vector<std::vector<float>> g_birksHitdXs;

//------------------------------------------------------------------------------------------------------------------------------------------

float GetEstimatedEnergy(const float alpha, const float beta, const int index)
{
    float estimatedEnergy = alpha * g_nonBirksHitSummeddQ.at(index);
    const int numBirksHits = g_birksHitdQs.at(index).size();
    
    for (int i = 0; i < numBirksHits; ++i)
    {
        const float unscaledValue = alpha * g_birksHitdQs.at(index).at(i);
        const float scaledValue   = unscaledValue / (1.f - beta * g_birksHitdQs.at(index).at(i) / g_birksHitdXs.at(index).at(i));
        
        if (scaledValue > 0.f && scaledValue < 1000.f)
            estimatedEnergy += scaledValue;
            
        else
            estimatedEnergy += unscaledValue;
    }
    
    return estimatedEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double CalculateMeanSquaredError(const double alpha, const double beta)
{
    float meanSquareError = 0.f;
    const int numDataPoints = g_trueEnergy.size();

    for (int i = 0; i < numDataPoints; i++)
    {
        const float trueEnergy = g_trueEnergy.at(i);
        const float estimatedEnergy = GetEstimatedEnergy(alpha, beta, i);
        meanSquareError += (trueEnergy - estimatedEnergy) * (trueEnergy - estimatedEnergy);
    }

    return meanSquareError / numDataPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MinuitFunction(Int_t &npar, Double_t *gin, Double_t &result, Double_t par[], Int_t iflag)
{
    result = CalculateMeanSquaredError(par[0], par[1]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const char * MinuitStatusToString(const int status)
{
    switch (status)
    {
        case 0:  return "no error";
        case 1:  return "covariance matrix was positive definite";
        case 2:  return "invalid Hesse matrix";
        case 3:  return "expected distance from the minimum above the maximum value";
        case 4:  return "reached call limit";
        case 5:  return "other failure";
        default: break;
    }
    
    return "unknown error";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RunMinuitFit(double &alpha, double &beta, const double alpha_initial, const double beta_initial, const double stepSizeFraction)
{
    TFitter *pTFitter = new TFitter(2);
    
    {
        double arg = -1.;
        pTFitter->ExecuteCommand("SET PRINTOUT", &arg, 1);
    }
    
    pTFitter->SetFCN(MinuitFunction);
    
    pTFitter->SetParameter(0, "alpha", alpha_initial, alpha_initial * stepSizeFraction, 0, 0);
    pTFitter->SetParameter(1, "beta",  beta_initial,  beta_initial * stepSizeFraction,  0, 0);
    
    double maxIterations[] = {500000000.};
    
    int status = pTFitter->ExecuteCommand("SIMPLEX", maxIterations, 1);
    
    if (status != 0)
        CERR("Simplex fit returned an error: " << MinuitStatusToString(status));
        
    TMinuit *const pMinuit = pTFitter->GetMinuit();
    int errFlag = 0;
    pMinuit->mnexcm("MIGRAD", maxIterations, 1, errFlag);
        
    if (status != 0)
        CERR("Migrad fit returned an error: " << MinuitStatusToString(status));
    
    alpha = pTFitter->GetParameter(0);
    beta  = pTFitter->GetParameter(1);
    
    double minFunctionError = CalculateMeanSquaredError(alpha, beta);
    
    COUT("Chi squared value at best parameters is " << minFunctionError);
    
    delete pTFitter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MakeDebugPlots(const double alpha, const double beta)
{
    TFile *pFile = new TFile("tmp.root", "RECREATE");
    TNtuple *const pNtuple1 = new TNtuple("BirksFit", "BirksFit", "EstimatedEnergy:TrueEnergy");
    TNtuple *const pNtuple2 = new TNtuple("BirksFit", "BirksFit", "Discrepancy:TrueEnergy");
    
    const int numDataPoints = g_trueEnergy.size();

    for (int i = 0; i < numDataPoints; i++)
    {
        const float trueEnergy = g_trueEnergy.at(i);
        const float estimatedEnergy = GetEstimatedEnergy(alpha, beta, i);
        
        if (estimatedEnergy > 0.f)
        {
            const float fractionalDiscrepancy = (estimatedEnergy - trueEnergy) / trueEnergy;
            pNtuple1->Fill(estimatedEnergy, trueEnergy);
            pNtuple2->Fill(fractionalDiscrepancy, trueEnergy);
        }
    }
    
    {
        struct PlotSettings2D plotSettings = g_defaultPlotSettings2D;
        plotSettings.xMax = 1000;
        plotSettings.yMax = 1000;
        plotSettings.xNumBins = 100;
        plotSettings.yNumBins = 100;
        plotSettings.trimSmallHistogramValues = false;
        strcpy(plotSettings.title, "Energy estimated from charge: with Birks' correction");
        strcpy(plotSettings.xTitle, "True energy (MeV)");
        strcpy(plotSettings.yTitle, "Estimated energy (MeV)");
        PlotNtuple2D(pNtuple1, "TrueEnergy", "EstimatedEnergy", "Estimator", plotSettings);
    }
    
    {
        struct PlotSettings2D plotSettings = g_defaultPlotSettings2D;
        plotSettings.xMax = 1000.f;
        plotSettings.yMin = -1.f;
        plotSettings.yMax = 2.f;
        plotSettings.xNumBins = 100;
        plotSettings.yNumBins = 100;
        plotSettings.trimSmallHistogramValues = false;
        strcpy(plotSettings.title, "Energy estimated from charge: with Birks' correction");
        strcpy(plotSettings.xTitle, "True energy (MeV)");
        strcpy(plotSettings.yTitle, "Fractional energy discrepancy");
        PlotNtuple2D(pNtuple2, "TrueEnergy", "Discrepancy", "Estimator_discrepancy", plotSettings);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OutputToNTuple(const char *const outputFilePath, const double alphaOut, const double betaOut, const double dQdxPole)
{
    TFile *pFile = new TFile(outputFilePath, "UPDATE");
    
    TNtuple *const pNtuple = new TNtuple("BirksFit", "BirksFit", "Alpha:Beta:dQdXPole");
    pNtuple->Fill(alphaOut, betaOut, dQdxPole);
    pNtuple->Write();
    
    pFile->Close();
    delete pFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BirksDataFit(const char *const inputFilePath, const char *const outputFilePath)
{
    TTree *const pTree = LoadTreeFromFile(inputFilePath, "BirksFitData");
    const int numEntries = pTree->GetEntries();
    
    float trueEnergy = 0.f;
    float nonBirksHitSummeddQ = 0.f;
    std::vector<float> *pBirksHitdQs = 0;
    std::vector<float> *pBirksHitdXs = 0;
    
    TBranch *pBranchBirksHitdQs = 0;
    TBranch *pBranchpBirksHitdXs = 0;
    TBranch *pBranchTrueEnergy = 0;
    TBranch *pBranchNonBirksHitSummeddQ = 0;
    
    pTree->SetBranchAddress("TruePrimaryEnergy", &trueEnergy, &pBranchTrueEnergy);
    pTree->SetBranchAddress("TotalNoBirksAdcIntegral", &nonBirksHitSummeddQ, &pBranchNonBirksHitSummeddQ);
    pTree->SetBranchAddress("BirksAdcIntegrals", &pBirksHitdQs, &pBranchBirksHitdQs);
    pTree->SetBranchAddress("ThreeDDistances", &pBirksHitdXs, &pBranchpBirksHitdXs);

    for (int i = 0; i < numEntries; ++i)
    {
        const auto tEntry = pTree->LoadTree(i);
        pBranchBirksHitdQs->GetEntry(tEntry);
        pBranchpBirksHitdXs->GetEntry(tEntry);
        pBranchTrueEnergy->GetEntry(tEntry);
        pBranchNonBirksHitSummeddQ->GetEntry(tEntry);
        
        g_trueEnergy.push_back(trueEnergy);
        g_nonBirksHitSummeddQ.push_back(nonBirksHitSummeddQ);

        std::vector<float> birksHitdQs;
        std::vector<float> birksHitdXs;

        for (int j = 0; j < (int)pBirksHitdQs->size(); ++j)
        {
            birksHitdQs.push_back(pBirksHitdQs->at(j));
            birksHitdXs.push_back(pBirksHitdXs->at(j));
        }
        
        g_birksHitdQs.push_back(birksHitdQs);
        g_birksHitdXs.push_back(birksHitdXs);
    }
    
    const double alphaInitial     = 0.003;
    const double betaInitial      = 0.000714286;
    const double stepSizeFraction = 0.1;
    
    double alpha = 0.005, beta = 0.00078;
    
    COUT("Running Minuit fit...");
    RunMinuitFit(alpha, beta, alphaInitial, betaInitial, stepSizeFraction);
    
    const double alphaOut = 1. / alpha;
    const double betaOut  = alpha / beta;
    const double dQdxPole = 1. / beta;
    
    COUT("Fitted values were:\n" << 
         "  - alpha'        = " << alpha    << " MeV/ADC\n" <<
         "  - beta'         = " << beta     << " cm/ADC\n" <<
         "  - alpha         = " << alphaOut << " ADC/MeV\n" <<
         "  - beta          = " << betaOut  << " MeV/cm\n" << 
         "  - Pole at dQ/dx = " << dQdxPole << " ADC/cm");
    
    MakeDebugPlots(alpha, beta);
    OutputToNTuple(outputFilePath, alphaOut, betaOut, dQdxPole);
}
