#include "Common.h"

#include "TFile.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TNtuple.h"

#include <vector>

std::vector<float>              g_trueEnergy;
std::vector<float>              g_showerChargeVector;
std::vector<std::vector<float>> g_dQdXVector;
std::vector<std::vector<float>> g_dXVector;

//------------------------------------------------------------------------------------------------------------------------------------------

float GetEstimatedEnergy(const float alpha, const float beta, const int index)
{
    // Scale up the charge induced by showers by alpha
    float estimatedEnergy = alpha * g_showerChargeVector.at(index);

    // Apply Birks' correction to get to dE = dx * dE/dx
    for (int i = 0, numBirksHits = g_dQdXVector.at(index).size(); i < numBirksHits; ++i)
    {
        const float dEdX_noBirks = alpha * g_dQdXVector.at(index).at(i);
        const float dEdX_Birks   = dEdX_noBirks / (1.f - beta * g_dQdXVector.at(index).at(i));

        if (dEdX_Birks > 0.f && dEdX_Birks < 1000.f)
            estimatedEnergy += dEdX_Birks * g_dXVector.at(index).at(i);

        else
            estimatedEnergy += dEdX_noBirks * g_dXVector.at(index).at(i);
    }

    return estimatedEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double CalculateMeanSquaredError(const double alpha, const double beta)
{
    float     meanSquareError = 0.f;
    const int numDataPoints   = g_trueEnergy.size();

    for (int i = 0; i < numDataPoints; i++)
    {
        const float trueEnergy      = g_trueEnergy.at(i);
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

const char *MinuitStatusToString(const int status)
{
    switch (status)
    {
        case 0:
            return "no error";
        case 1:
            return "covariance matrix was positive definite";
        case 2:
            return "invalid Hesse matrix";
        case 3:
            return "expected distance from the minimum above the maximum value";
        case 4:
            return "reached call limit";
        case 5:
            return "other failure";
        default:
            break;
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
    pTFitter->SetParameter(1, "beta", beta_initial, beta_initial * stepSizeFraction, 0, 0);

    double maxIterations[] = {500000000.};

    int status = pTFitter->ExecuteCommand("SIMPLEX", maxIterations, 1);

    if (status != 0)
        CERR("Simplex fit returned an error: " << MinuitStatusToString(status));

    TMinuit *const pMinuit = pTFitter->GetMinuit();
    int            errFlag = 0;
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

void MakeDebugPlots(const double alpha, const double beta, const char *const outputDir)
{
    TFile *        pFile    = new TFile(TString(outputDir) + "/Tmp.root", "RECREATE");
    TNtuple *const pNtuple1 = new TNtuple("BirksFit", "BirksFit", "EstimatedEnergy:TrueEnergy");
    TNtuple *const pNtuple2 = new TNtuple("BirksFit", "BirksFit", "Discrepancy:TrueEnergy");

    const int numDataPoints = g_trueEnergy.size();

    for (int i = 0; i < numDataPoints; i++)
    {
        const float trueEnergy      = g_trueEnergy.at(i);
        const float estimatedEnergy = GetEstimatedEnergy(alpha, beta, i);

        if (estimatedEnergy > 0.f)
        {
            const float fractionalDiscrepancy = (estimatedEnergy - trueEnergy) / trueEnergy;
            pNtuple1->Fill(estimatedEnergy, trueEnergy);
            pNtuple2->Fill(fractionalDiscrepancy, trueEnergy);
        }
    }

    {
        struct PlotSettings2D plotSettings    = g_defaultPlotSettings2D;
        plotSettings.xMax                     = 1;
        plotSettings.yMax                     = 1;
        plotSettings.xNumBins                 = 100;
        plotSettings.yNumBins                 = 100;
        plotSettings.trimSmallHistogramValues = false;
        strcpy(plotSettings.title, "Energy estimated from charge: with Birks' correction");
        strcpy(plotSettings.xTitle, "True energy (GeV)");
        strcpy(plotSettings.yTitle, "Estimated energy (GeV)");
        TCanvas *pCanvas = PlotNtuple2D(pNtuple1, "TrueEnergy", "EstimatedEnergy", "Estimator", plotSettings);
        pCanvas->SaveAs(TString(outputDir) + "/BirksFitPerformancePlot.root");
    }

    {
        struct PlotSettings2D plotSettings    = g_defaultPlotSettings2D;
        plotSettings.xMax                     = 1.f;
        plotSettings.yMin                     = -1.f;
        plotSettings.yMax                     = 2.f;
        plotSettings.xNumBins                 = 100;
        plotSettings.yNumBins                 = 100;
        plotSettings.trimSmallHistogramValues = false;
        strcpy(plotSettings.title, "Energy estimated from charge: with Birks' correction");
        strcpy(plotSettings.xTitle, "True energy (GeV)");
        strcpy(plotSettings.yTitle, "Fractional energy discrepancy");
        TCanvas *pCanvas = PlotNtuple2D(pNtuple2, "TrueEnergy", "Discrepancy", "Estimator_discrepancy", plotSettings);
        pCanvas->SaveAs(TString(outputDir) + "/BirksFitPerformanceDiscrepancyPlot.root");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OutputToNTuple(const TString &outputFilePath, const double alphaOut, const double betaOut, const double dQdxPole)
{
    TFile *pFile = new TFile(outputFilePath, "RECREATE");

    TNtuple *const pNtuple = new TNtuple("BirksFit", "BirksFit", "Alpha:Beta:dQdXPole");
    pNtuple->Fill(alphaOut, betaOut, dQdxPole);
    pNtuple->Write();

    pFile->Close();
    delete pFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FitBirksData(const char *const inputFilePath, const char *const ntupleName, const char *const outputDir, const char *const outputName)
{
    const float       minEnergyWeightedContainedPfoFraction = 0.9f;
    const float       minMcMatchCompleteness                = 0.9f;
    const float       minMcMatchPurity                      = 0.9f;
    const float       maxHitFracLostByFit                   = 0.1f;
    const std::size_t minNumCollectionPlaneHits             = 3UL;

    // Load the ntuple
    TFile       ntupleFile(inputFilePath);
    TTreeReader treeReader(ntupleName, &ntupleFile);

    // Primary fit data
    TTreeReaderValue<UInt_t>                            numPrimaryEntries(treeReader, "numPrimaryEntries");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> primary_dQdX(treeReader, "primary_dQdX");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> primary_dX(treeReader, "primary_dX");
    TTreeReaderValue<std::vector<Float_t>>              primary_ShowerCharge(treeReader, "primary_ShowerCharge");
    TTreeReaderValue<std::vector<Float_t>>              primary_mc_KineticEnergy(treeReader, "primary_mc_KineticEnergy");

    // Primary quality cut data
    TTreeReaderValue<std::vector<UInt_t>>  primary_NumVectorEntries(treeReader, "primary_NumVectorEntries");
    TTreeReaderValue<std::vector<UInt_t>>  primary_NumHitsLostToFittingErrors(treeReader, "primary_NumHitsLostToFittingErrors");
    TTreeReaderValue<std::vector<Bool_t>>  primary_WasReconstructedWithVertex(treeReader, "primary_WasReconstructedWithVertex");
    TTreeReaderValue<std::vector<Bool_t>>  primary_HasMCInfo(treeReader, "primary_HasMCInfo");
    TTreeReaderValue<std::vector<Bool_t>>  primary_IsVertexFiducial(treeReader, "primary_IsVertexFiducial");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_EnergyWeightedContainedPfoFraction(
        treeReader, "primary_mc_EnergyWeightedContainedPfoFraction");
    TTreeReaderValue<std::vector<Bool_t>>  primary_mc_IsGoodMatch(treeReader, "primary_mc_IsGoodMatch");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_MatchPurity(treeReader, "primary_mc_MatchPurity");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_MatchCompleteness(treeReader, "primary_mc_MatchCompleteness");

    // Cosmic fit data
    TTreeReaderValue<UInt_t>                            numCosmicRayEntries(treeReader, "numCosmicRayEntries");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> cr_dQdX(treeReader, "cr_dQdX");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> cr_dX(treeReader, "cr_dX");
    TTreeReaderValue<std::vector<Float_t>>              cr_ShowerCharge(treeReader, "cr_ShowerCharge");
    TTreeReaderValue<std::vector<Float_t>>              cr_mc_KineticEnergy(treeReader, "cr_mc_KineticEnergy");

    // Cosmic quality cut data
    TTreeReaderValue<std::vector<UInt_t>>  cr_NumVectorEntries(treeReader, "cr_NumVectorEntries");
    TTreeReaderValue<std::vector<UInt_t>>  cr_NumHitsLostToFittingErrors(treeReader, "cr_NumHitsLostToFittingErrors");
    TTreeReaderValue<std::vector<Bool_t>>  cr_WasReconstructedWithVertex(treeReader, "cr_WasReconstructedWithVertex");
    TTreeReaderValue<std::vector<Bool_t>>  cr_HasMCInfo(treeReader, "cr_HasMCInfo");
    TTreeReaderValue<std::vector<Bool_t>>  cr_IsVertexFiducial(treeReader, "cr_IsVertexFiducial");
    TTreeReaderValue<std::vector<Float_t>> cr_mc_EnergyWeightedContainedPfoFraction(treeReader, "cr_mc_EnergyWeightedContainedPfoFraction");
    TTreeReaderValue<std::vector<Bool_t>>  cr_mc_IsGoodMatch(treeReader, "cr_mc_IsGoodMatch");
    TTreeReaderValue<std::vector<Float_t>> cr_mc_MatchPurity(treeReader, "cr_mc_MatchPurity");
    TTreeReaderValue<std::vector<Float_t>> cr_mc_MatchCompleteness(treeReader, "cr_mc_MatchCompleteness");

    std::size_t numCosmicRayDatapoints(0UL), numPrimaryDatapoints(0UL);

    while (treeReader.Next())
    {
        for (std::size_t i = 0U; i < *numPrimaryEntries; ++i)
        {
            // Require reconstruction and good MC info
            if (!(*primary_WasReconstructedWithVertex)[i] || !(*primary_HasMCInfo)[i] || !(*primary_mc_IsGoodMatch)[i])
                continue;

            // Require a minimum purity and completeness of the MC match
            if ((*primary_mc_MatchPurity)[i] < minMcMatchPurity || (*primary_mc_MatchCompleteness)[i] < minMcMatchCompleteness)
                continue;

            // Require fiducial vertex and enough energy contained
            if (!(*primary_IsVertexFiducial)[i] || (*primary_mc_EnergyWeightedContainedPfoFraction)[i] < minEnergyWeightedContainedPfoFraction)
                continue;

            // Require at least n (and at least 1) collection plane hits
            if (minNumCollectionPlaneHits < 1UL || (*primary_NumVectorEntries)[i] < minNumCollectionPlaneHits)
                continue;

            // No more than a given fraction of hits may be lost due to fitting issues
            if ((*primary_NumHitsLostToFittingErrors)[i] / (*primary_NumVectorEntries)[i] > maxHitFracLostByFit)
                continue;

            g_trueEnergy.push_back((*primary_mc_KineticEnergy)[i]);
            g_showerChargeVector.push_back((*primary_ShowerCharge)[i]);
            g_dQdXVector.push_back((*primary_dQdX)[i]);
            g_dXVector.push_back((*primary_dX)[i]);

            ++numPrimaryDatapoints;
        }

        for (std::size_t i = 0U; i < *numCosmicRayEntries; ++i)
        {
            // Require reconstruction and good MC info
            if (!(*cr_WasReconstructedWithVertex)[i] || !(*cr_HasMCInfo)[i] || !(*cr_mc_IsGoodMatch)[i])
                continue;

            // Require a minimum purity and completeness of the MC match
            if ((*cr_mc_MatchPurity)[i] < minMcMatchPurity || (*cr_mc_MatchCompleteness)[i] < minMcMatchCompleteness)
                continue;

            // Require fiducial vertex and enough energy contained
            if (!(*cr_IsVertexFiducial)[i] || (*cr_mc_EnergyWeightedContainedPfoFraction)[i] < minEnergyWeightedContainedPfoFraction)
                continue;

            // Require at least n (and at least 1) collection plane hits
            if (minNumCollectionPlaneHits < 1UL || (*cr_NumVectorEntries)[i] < minNumCollectionPlaneHits)
                continue;

            // No more than a given fraction of hits may be lost due to fitting issues
            if ((*cr_NumHitsLostToFittingErrors)[i] / ((*cr_NumVectorEntries)[i] + (*cr_NumHitsLostToFittingErrors)[i]) > maxHitFracLostByFit)
                continue;

            g_trueEnergy.push_back((*cr_mc_KineticEnergy)[i]);
            g_showerChargeVector.push_back((*cr_ShowerCharge)[i]);
            g_dQdXVector.push_back((*cr_dQdX)[i]);
            g_dXVector.push_back((*cr_dX)[i]);

            ++numCosmicRayDatapoints;
        }
    }

    const std::size_t numDatapoints = g_trueEnergy.size();

    if (numDatapoints == 0UL)
    {
        COUT("There were no datapoints");
        return;
    }

    const double alphaInitial     = 0.000003;
    const double betaInitial      = 0.000714286;
    const double stepSizeFraction = 0.1;

    double alpha = alphaInitial, beta = betaInitial;

    COUT("Running Minuit fit...");
    RunMinuitFit(alpha, beta, alphaInitial, betaInitial, stepSizeFraction);

    const double alphaOut = 1. / alpha;
    const double betaOut  = alpha / beta;
    const double dQdxPole = 1. / beta;

    COUT("\nSuccessfully trained using " << numDatapoints << " datapoints (" << numPrimaryDatapoints << " primary, " << numCosmicRayDatapoints << " cosmic)");
    COUT("Fitted values were:\n"
         << "  - alpha'        = " << alpha << " GeV/ADC\n"
         << "  - beta'         = " << beta << " cm/ADC\n"
         << "  - alpha         = " << alphaOut << " ADC/GeV\n"
         << "  - beta          = " << betaOut << " GeV/cm\n"
         << "  - Pole at dQ/dx = " << dQdxPole << " ADC/cm");

    MakeDebugPlots(alpha, beta, outputDir);
    OutputToNTuple(TString(outputDir) + "/" + TString(outputName), alphaOut, betaOut, dQdxPole);
}
