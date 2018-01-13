#include "Common.h"


//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessEnergyFromRangeData(const char *const inputFilePath, const char *const outputFilePath, const bool isProton)
{
    TNtuple *const pNtuple = LoadNTupleFromFile(inputFilePath, isProton ? "EnergyFromRangeDataProtons" : "EnergyFromRangeDataPionsMuons");
    
    struct PlotSettings2D plotSettings = g_defaultPlotSettings2D;
    
    if (isProton)
    {
        plotSettings.xMax = 150.f;
        plotSettings.yMax = 1.f;
    }
    
    else
    {
        plotSettings.xMax = 600.f;
        plotSettings.yMax = 2.f;
    }
    
    plotSettings.xNumBins = 150;
    plotSettings.yNumBins = 150;
    
    strcpy(plotSettings.xTitle, "3D track length (cm)");
    strcpy(plotSettings.yTitle, "True energy (GeV)");
    strcpy(plotSettings.title, isProton ? "Energy from range: protons" : "Energy from range: pions/muons");
    
    const float binWidthX = 1.f;
    const float binWidthY = 0.01f;
    
    float range = 0.f, trueEnergy = 0.f;
    
    pNtuple->SetBranchAddress("Range", &range);
    pNtuple->SetBranchAddress("TrueEnergy", &trueEnergy);
    
    const int nEntries = pNtuple->GetEntries();
    
    float rangeMin = 0.f, rangeMax = 0.f, energyMin = 0.f, energyMax = 0.f;
    bool firstTime = true;
    
    for (int i = 0; i < nEntries; ++i)
    {
        pNtuple->GetEntry(i);
        
        if (firstTime)
        {
            rangeMin = range;
            rangeMax = range;
            energyMin = trueEnergy;
            energyMax = trueEnergy;
            firstTime = false;
            continue;
        }
        
        if (range < rangeMin) rangeMin = range;
        if (range > rangeMax) rangeMax = range;
        if (trueEnergy < energyMin) energyMin = trueEnergy;
        if (trueEnergy > energyMax) energyMax = trueEnergy;
    }
    
    int binPopulations[1000][1000] = {0};
    
    const int xNumBins = (int)round((rangeMax + binWidthX - rangeMin) / binWidthX);
    const int yNumBins = (int)round((energyMax + binWidthY - energyMin) / binWidthY);
    
    for (int i = 0; i < nEntries; ++i)
    {
        pNtuple->GetEntry(i);
        
        const int xBin = (int)round((float)xNumBins * (range - rangeMin) / (rangeMax + binWidthX - rangeMin));
        const int yBin = (int)round((float)yNumBins * (trueEnergy - energyMin) / (energyMax + binWidthY - energyMin));
        
        if (xBin < 1000 && yBin < 1000)
            ++binPopulations[xBin][yBin];
    }
    
    int mostPopulousBin[1000] = {0};
    
    for (int xBin = 0; xBin < 1000; ++xBin)
    {
        int biggestPopulation = 0;
        
        for (int yBin = 0; yBin < 1000; ++yBin)
        {
            if (binPopulations[xBin][yBin] > biggestPopulation)
            {
                biggestPopulation = binPopulations[xBin][yBin];
                mostPopulousBin[xBin] = yBin + 1;
            }
        }
    }
    
    TNtuple *const pNtuplePlot = new TNtuple("EnergyFromRange","EnergyFromRange", "AvgRange:Energy");
    
    TFile *const pFile = new TFile(outputFilePath, "UPDATE");
    TNtuple *const pNtupleBinned = new TNtuple(isProton ? "EnergyFromRangeProtons" : "EnergyFromRangePionsMuons","EnergyFromRange",
                                               "RangeMin:RangeMax:Energy");
    
    for (int xBin = 0; xBin < 1000; ++xBin)
    {
        if (mostPopulousBin[xBin] - 1 >= 0)
        {
            const float rangeBinMin = (float)xBin * binWidthX + rangeMin;
            const float rangeBinMax = rangeBinMin + binWidthX;
            const float averageRange = 0.5f * (rangeBinMin + rangeBinMax);
            
            if (averageRange < (isProton ? 100.f : 400.f))
            {
                const float energyBinMin = (float)(mostPopulousBin[xBin] - 1) * binWidthY + energyMin;
                const float energyBinMax = energyBinMin + binWidthY;
                const float averageEnergy = 0.5f * (energyBinMin + energyBinMax);
                
                pNtupleBinned->Fill(rangeBinMin, rangeBinMax, averageEnergy);
                pNtuplePlot->Fill(averageRange, averageEnergy);
            }
        }
    }
    
    
//    PlotNtuple2D(pNtuple, "Range", "TrueEnergy", "Histogram1", plotSettings);
//    
//    plotSettings.newCanvas  = false;
//    plotSettings.plotType   = SAME;
//    plotSettings.pointColor = kRed;
//    PlotNtuple2D(pNtuplePlot, "AvgRange", "Energy", "Histogram2", plotSettings);
    
    pNtupleBinned->Write();
    pFile->Close();
    delete pFile;
}
