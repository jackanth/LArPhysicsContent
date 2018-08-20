#include "Common.h"

#include "TFile.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TNtuple.h"

#include <vector>

std::vector<float>              g_trueEnergy;
std::vector<float>              g_nonBirksHitSummeddQ;
std::vector<std::vector<float>> g_birksHitdQs;
std::vector<std::vector<float>> g_birksHitdXs;

//------------------------------------------------------------------------------------------------------------------------------------------

float GetEstimatedEnergy(const float alpha, const float beta, const int index)
{
    float     estimatedEnergy = alpha * g_nonBirksHitSummeddQ.at(index);
    const int numBirksHits    = g_birksHitdQs.at(index).size();

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

void MakeDebugPlots(const double alpha, const double beta)
{
    TFile *        pFile    = new TFile("tmp.root", "RECREATE");
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
        plotSettings.xMax                     = 1000;
        plotSettings.yMax                     = 1000;
        plotSettings.xNumBins                 = 100;
        plotSettings.yNumBins                 = 100;
        plotSettings.trimSmallHistogramValues = false;
        strcpy(plotSettings.title, "Energy estimated from charge: with Birks' correction");
        strcpy(plotSettings.xTitle, "True energy (MeV)");
        strcpy(plotSettings.yTitle, "Estimated energy (MeV)");
        PlotNtuple2D(pNtuple1, "TrueEnergy", "EstimatedEnergy", "Estimator", plotSettings);
    }

    {
        struct PlotSettings2D plotSettings    = g_defaultPlotSettings2D;
        plotSettings.xMax                     = 1000.f;
        plotSettings.yMin                     = -1.f;
        plotSettings.yMax                     = 2.f;
        plotSettings.xNumBins                 = 100;
        plotSettings.yNumBins                 = 100;
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

void FitBirksData(const char *const inputFilePath, const char *const outputFilePath)
{
    const float minEnergyWeightedContainedPfoFraction     = 0.9f;
    const float minRecoMcMatchCompleteness                = 0.9f;
    const float minRecoMcMatchPurity                      = 0.9f;
    const float minRecoMcMatchCollectionPlaneCompleteness = 0.9f;
    const float minRecoMcMatchCollectionPlanePurity       = 0.9f;

    // Load the ntuple
    TFile       ntupleFile(inputFilePath);
    TTreeReader treeReader("PandoraNtuple_EnergyEstimator_TrainingSet_NoParticleId", &ntupleFile);

    // Prepare the values
    TTreeReaderValue<UInt_t> numPrimaryEntries(treeReader, "numPrimaryEntries");

    TTreeReaderValue<std::vector<std::vector<Float_t>>> primary_ThreeDDistances(treeReader, "primary_ThreeDDistances");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> primary_IntegratedAdcCounts(treeReader, "primary_IntegratedAdcCounts");
    TTreeReaderValue<std::vector<UInt_t>>               primary_NumVectorEntries(treeReader, "primary_NumVectorEntries");
    TTreeReaderValue<std::vector<Float_t>>              primary_ShowerAdcCount(treeReader, "primary_ShowerAdcCount");
    TTreeReaderValue<std::vector<Bool_t>>  primary_WasReconstructedWithVertex(treeReader, "primary_WasReconstructedWithVertex");
    TTreeReaderValue<std::vector<Bool_t>>  primary_HasMCInfo(treeReader, "primary_HasMCInfo");
    TTreeReaderValue<std::vector<Bool_t>>  primary_IsVertexFiducial(treeReader, "primary_IsVertexFiducial");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_EnergyWeightedContainedPfoFraction(
        treeReader, "primary_mc_EnergyWeightedContainedPfoFraction");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_KineticEnergy(treeReader, "primary_mc_KineticEnergy");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_RecoMcMatchCompleteness(treeReader, "primary_mc_RecoMcMatchCompleteness");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_RecoMcMatchPurity(treeReader, "primary_mc_RecoMcMatchPurity");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_RecoMcMatchCollectionPlaneCompleteness(
        treeReader, "primary_mc_RecoMcMatchCollectionPlaneCompleteness");
    TTreeReaderValue<std::vector<Float_t>> primary_mc_RecoMcMatchCollectionPlanePurity(
        treeReader, "primary_mc_RecoMcMatchCollectionPlanePurity");

    while (treeReader.Next())
    {
        for (std::size_t i = 0U; i < *numPrimaryEntries; ++i)
        {
            if (!(*primary_WasReconstructedWithVertex)[i] || !(*primary_HasMCInfo)[i])
                continue;

            if ((*primary_NumVectorEntries)[i] < 2)
                continue;

            if (!(*primary_IsVertexFiducial)[i] || ((*primary_mc_EnergyWeightedContainedPfoFraction)[i] < minRecoMcMatchCompleteness))
                continue;

            if (((*primary_mc_RecoMcMatchCompleteness)[i] < minRecoMcMatchCompleteness) || ((*primary_mc_RecoMcMatchPurity)[i] < minRecoMcMatchPurity) ||
                ((*primary_mc_RecoMcMatchCollectionPlaneCompleteness)[i] < minRecoMcMatchCollectionPlaneCompleteness) ||
                ((*primary_mc_RecoMcMatchCollectionPlanePurity)[i] < minRecoMcMatchCollectionPlanePurity))
            {
                continue;
            }

            g_trueEnergy.push_back((*primary_mc_KineticEnergy)[i]);
            g_nonBirksHitSummeddQ.push_back((*primary_ShowerAdcCount)[i]);
            g_birksHitdQs.push_back((*primary_IntegratedAdcCounts)[i]);
            g_birksHitdXs.push_back((*primary_ThreeDDistances)[i]);
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

    COUT("\nSuccessfully trained using " << numDatapoints << " datapoints");
    COUT("Fitted values were:\n"
         << "  - alpha'        = " << alpha << " GeV/ADC\n"
         << "  - beta'         = " << beta << " cm/ADC\n"
         << "  - alpha         = " << alphaOut << " ADC/GeV\n"
         << "  - beta          = " << betaOut << " GeV/cm\n"
         << "  - Pole at dQ/dx = " << dQdxPole << " ADC/cm");

    // MakeDebugPlots(alpha, beta);
    // OutputToNTuple(outputFilePath, alphaOut, betaOut, dQdxPole);
}
