/**
 *  @file   larphysicscontent/LArHelpers/LArRootHelper.cxx
 *
 *  @brief  Implementation of the lar root helper class.
 *
 *  $Log: $
 */

#include "larphysicscontent/LArHelpers/LArRootHelper.h"

#include "TFile.h"

#include <iostream>

namespace lar_physics_content
{
void LArRootHelper::WriteNTuple(TNtuple *const pNtuple, const std::string &fileName, const bool verboseMode)
{
    if (verboseMode)
        pNtuple->Print();

    TFile *pFile = new TFile(fileName.c_str(), "NEW");
    pNtuple->Write();

    pFile->Close();
    delete pFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TNtuple *LArRootHelper::LoadNTupleFromFile(const std::string &filePath, const std::string &nTupleName)
{
    TFile *pFile = new TFile(filePath.c_str(), "READ");

    if (!pFile->IsOpen())
    {
        std::cout << "LArRootHelper: failed to open file at " << filePath << std::endl;
        return NULL;
    }

    if (!pFile->GetListOfKeys()->Contains(nTupleName.c_str()))
    {
        std::cout << "LArRootHelper: data file at " << filePath << " did not contain key '" << nTupleName << "'" << std::endl;
        return NULL;
    }

    return (TNtuple *)pFile->Get(nTupleName.c_str());
}

} // namespace lar_physics_content
