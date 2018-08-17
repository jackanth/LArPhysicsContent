/**
 *  @file   test/ValidateNtuple.c
 *
 *  @brief  ROOT macro for validating ntuple output
 *
 *  $Log: $
 */

#include "Rtypes.h"
#include "TString.h"

#include <vector>

//------------------------------------------------------------------------------------------------------------------------------------------

#define TEXT_RED "\033[0;31m"
#define TEXT_GREEN "\033[0;32m"
#define TEXT_RED_BOLD "\033[1;31m"
#define TEXT_GREEN_BOLD "\033[1;32m"
#define TEXT_YELLOW_BOLD "\033[1;33m"
#define TEXT_WHITE_BOLD "\033[1;37m"
#define TEXT_NORMAL "\033[0m"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Macro for testing the value of an ntuple parameter
 * 
 *  @param  fn the function to retrieve the true value
 *  @param  valueExpr the expression giving the found value
 *  @param  counter the counter argument to the function
 */
#define TEST(fn, valueExpr, counter)                                                                                            \
{                                                                                                                               \
    if (fn(counter) == valueExpr)                                                                                               \
    {                                                                                                                           \
        ++successfulTests;                                                                                                      \
        std::cout << "[ " << TEXT_GREEN_BOLD << "SUCCESS" << TEXT_NORMAL << " ] "#fn"("#valueExpr", "#counter");" << std::endl; \
    }                                                                                                                           \
                                                                                                                                \
    else                                                                                                                        \
    {                                                                                                                           \
        ++failedTests;                                                                                                          \
        std::cerr << "[ " << TEXT_RED_BOLD << "FAILURE" << TEXT_NORMAL << " ] "#fn"("#valueExpr", "#counter");" << std::endl;   \
        std::cerr << "    - Got value     " << valueExpr << std::endl;                                                          \
        std::cerr << "    - Correct value " << fn(counter) << std::endl;                                                        \
    }                                                                                                                           \
}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Macro for testing the size of an ntuple vector parameter
 * 
 *  @param  theVector the vector
 *  @param  correctSize the correct size
 */
#define TEST_SIZE(theVector, correctSize)                                                                             \
{                                                                                                                     \
    if ((theVector).size() == correctSize)                                                                            \
    {                                                                                                                 \
        ++successfulTests;                                                                                            \
        std::cout << "[ " << TEXT_GREEN_BOLD << "SUCCESS" << TEXT_NORMAL << " ] ("#theVector").size();" << std::endl; \
    }                                                                                                                 \
                                                                                                                      \
    else                                                                                                              \
    {                                                                                                                 \
        ++failedTests;                                                                                                \
        std::cerr << "[ " << TEXT_RED_BOLD << "FAILURE" << TEXT_NORMAL << " ] ("#theVector").size();" << std::endl;   \
        std::cerr << "    - Got size     " << (theVector).size() << std::endl;                                        \
        std::cerr << "    - Correct size " << correctSize << std::endl;                                               \
    }                                                                                                                 \
}

//------------------------------------------------------------------------------------------------------------------------------------------

Float_t GetRFloatValue(const int counter);
Int_t GetRIntValue(const int counter);
Bool_t GetRBoolValue(const int counter);
UInt_t GetRUIntValue(const int counter);
ULong64_t GetRULong64Value(const int counter);
TString GetRTStringValue(const int counter);
std::vector<Float_t> GetRFloatVectorValue(const int counter);
std::vector<Int_t> GetRIntVectorValue(const int counter);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Function for printing vectors
 * 
 *  @param  os the ostream object
 *  @param  vec the vector to print
 * 
 *  @return the ostream object
 */
template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
{
    if (vec.empty())
        return os;

    os << '(';

    for (std::size_t i = 0UL; i < vec.size() - 1; ++i)
        os << vec[i] << ", ";

    os << vec.back() << ')';
    return os;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Validate the ntuple produced by the Pandora test ntuple tools
 * 
 *  @param  filePath path to the ntuple file
 *  @param  ntupleName the ntuple name
 */
void ValidateNtuple(const TString &filePath, const TString &ntupleName = "PandoraNtuple")
{
    // Load the ntuple
    TFile ntupleFile(filePath);
    TTreeReader treeReader(ntupleName, &ntupleFile);

    // Prepare the reserved values
    TTreeReaderValue<Int_t> fileId(treeReader, "fileId");
    TTreeReaderValue<Int_t> eventNum(treeReader, "eventNum");
    TTreeReaderValue<UInt_t> numNeutrinos(treeReader, "numNeutrinoEntries");
    TTreeReaderValue<UInt_t> numCosmicRays(treeReader, "numCosmicRayEntries");
    TTreeReaderValue<UInt_t> numPrimaries(treeReader, "numPrimaryEntries");
    TTreeReaderValue<UInt_t> numPfos(treeReader, "numPfoEntries");
    TTreeReaderValue<Bool_t> hasMcInfo(treeReader, "hasMcInfo");

    // Prepare the per-event values
    TTreeReaderValue<Float_t> evt_RFloat(treeReader, "evt_RFloat");
    TTreeReaderValue<Int_t> evt_RInt(treeReader, "evt_RInt");
    TTreeReaderValue<Bool_t> evt_RBool(treeReader, "evt_RBool");
    TTreeReaderValue<UInt_t> evt_RUInt(treeReader, "evt_RUInt");
    TTreeReaderValue<ULong64_t> evt_RULong64(treeReader, "evt_RULong64");
    TTreeReaderValue<TString> evt_RTString(treeReader, "evt_RTString");
    TTreeReaderValue<std::vector<Float_t>> evt_RFloatVector(treeReader, "evt_RFloatVector");
    TTreeReaderValue<std::vector<Int_t>> evt_RIntVector(treeReader, "evt_RIntVector");

    // Prepare the per-neutrino values
    TTreeReaderValue<std::vector<Float_t>> nu_RFloat(treeReader, "nu_RFloat");
    TTreeReaderValue<std::vector<Int_t>> nu_RInt(treeReader, "nu_RInt");
    TTreeReaderValue<std::vector<Bool_t>> nu_RBool(treeReader, "nu_RBool");
    TTreeReaderValue<std::vector<UInt_t>> nu_RUInt(treeReader, "nu_RUInt");
    TTreeReaderValue<std::vector<ULong64_t>> nu_RULong64(treeReader, "nu_RULong64");
    TTreeReaderValue<std::vector<TString>> nu_RTString(treeReader, "nu_RTString");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> nu_RFloatVector(treeReader, "nu_RFloatVector");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> nu_RIntVector(treeReader, "nu_RIntVector");

    // Prepare the per-PFO values
    TTreeReaderValue<std::vector<Float_t>> pfo_RFloat(treeReader, "pfo_RFloat");
    TTreeReaderValue<std::vector<Int_t>> pfo_RInt(treeReader, "pfo_RInt");
    TTreeReaderValue<std::vector<Bool_t>> pfo_RBool(treeReader, "pfo_RBool");
    TTreeReaderValue<std::vector<UInt_t>> pfo_RUInt(treeReader, "pfo_RUInt");
    TTreeReaderValue<std::vector<ULong64_t>> pfo_RULong64(treeReader, "pfo_RULong64");
    TTreeReaderValue<std::vector<TString>> pfo_RTString(treeReader, "pfo_RTString");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> pfo_RFloatVector(treeReader, "pfo_RFloatVector");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> pfo_RIntVector(treeReader, "pfo_RIntVector");

    // Prepare the per-primary values
    TTreeReaderValue<std::vector<Float_t>> primary_RFloat(treeReader, "primary_RFloat");
    TTreeReaderValue<std::vector<Int_t>> primary_RInt(treeReader, "primary_RInt");
    TTreeReaderValue<std::vector<Bool_t>> primary_RBool(treeReader, "primary_RBool");
    TTreeReaderValue<std::vector<UInt_t>> primary_RUInt(treeReader, "primary_RUInt");
    TTreeReaderValue<std::vector<ULong64_t>> primary_RULong64(treeReader, "primary_RULong64");
    TTreeReaderValue<std::vector<TString>> primary_RTString(treeReader, "primary_RTString");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> primary_RFloatVector(treeReader, "primary_RFloatVector");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> primary_RIntVector(treeReader, "primary_RIntVector");

    // Prepare the per-cosmic values
    TTreeReaderValue<std::vector<Float_t>> cr_RFloat(treeReader, "cr_RFloat");
    TTreeReaderValue<std::vector<Int_t>> cr_RInt(treeReader, "cr_RInt");
    TTreeReaderValue<std::vector<Bool_t>> cr_RBool(treeReader, "cr_RBool");
    TTreeReaderValue<std::vector<UInt_t>> cr_RUInt(treeReader, "cr_RUInt");
    TTreeReaderValue<std::vector<ULong64_t>> cr_RULong64(treeReader, "cr_RULong64");
    TTreeReaderValue<std::vector<TString>> cr_RTString(treeReader, "cr_RTString");
    TTreeReaderValue<std::vector<std::vector<Float_t>>> cr_RFloatVector(treeReader, "cr_RFloatVector");
    TTreeReaderValue<std::vector<std::vector<Int_t>>> cr_RIntVector(treeReader, "cr_RIntVector");

    std::cout << "Beginning ntuple validation" << std::endl;

    int successfulTests(0), failedTests(0), warnings(0);
    int evtCounter(0);

    while (treeReader.Next())
    {
        int nuCounter(0), pfoCounter(0), primaryCounter(0), cosmicCounter(0);

        // Output event info
        std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
        std::cout << TEXT_WHITE_BOLD << "Processing event #" << evtCounter << TEXT_NORMAL << std::endl;
        std::cout << "fileId        = " << *fileId << std::endl;
        std::cout << "eventNum      = " << *eventNum << std::endl;
        std::cout << "numNeutrinos  = " << *numNeutrinos << std::endl;
        std::cout << "numCosmicRays = " << *numCosmicRays << std::endl;
        std::cout << "numPrimaries  = " << *numPrimaries << std::endl;
        std::cout << "numPfos       = " << *numPfos << std::endl;
        std::cout << "hasMcInfo     = " << std::boolalpha << *hasMcInfo << std::endl;

        // Event numbers should get out of sync only if there are completely null events (i.e. all hits removed by CR tagger) or the ntuple
        // was composed through multiple pandora calls
        if (evtCounter != *eventNum)
        {
            std::cerr << "[ " << TEXT_YELLOW_BOLD << "WARNING" << TEXT_NORMAL << " ] Event numbers are misaligned" << std::endl;
            ++warnings;
        }

        // Per-event tests
        std::cout << std::endl << "Testing per-event parameters" << std::endl;

        TEST(GetRFloatValue, *evt_RFloat, *eventNum);
        TEST(GetRIntValue, *evt_RInt, *eventNum);
        TEST(GetRBoolValue, *evt_RBool, *eventNum);
        TEST(GetRUIntValue, *evt_RUInt, *eventNum);
        TEST(GetRULong64Value, *evt_RULong64, *eventNum);
        TEST(GetRTStringValue, *evt_RTString, *eventNum);
        TEST(GetRFloatVectorValue, *evt_RFloatVector, *eventNum);
        TEST(GetRIntVectorValue, *evt_RIntVector, *eventNum);

        // Per-neutrino tests
        std::cout << std::endl << "Testing per-neutrino parameters" << std::endl;

        TEST_SIZE(*nu_RFloat, *numNeutrinos);
        TEST_SIZE(*nu_RInt, *numNeutrinos);
        TEST_SIZE(*nu_RBool, *numNeutrinos);
        TEST_SIZE(*nu_RUInt, *numNeutrinos);
        TEST_SIZE(*nu_RULong64, *numNeutrinos);
        TEST_SIZE(*nu_RTString, *numNeutrinos);
        TEST_SIZE(*nu_RFloatVector, *numNeutrinos);
        TEST_SIZE(*nu_RIntVector, *numNeutrinos);

        for (std::size_t i = 0; i < *numNeutrinos; ++i)
        {
            if (i < (*nu_RFloat).size())
                TEST(GetRFloatValue, (*nu_RFloat).at(i), nuCounter);

            if (i < (*nu_RInt).size())
                TEST(GetRIntValue, (*nu_RInt).at(i), nuCounter);

            if (i < (*nu_RBool).size())
                TEST(GetRBoolValue, (*nu_RBool).at(i), nuCounter);

            if (i < (*nu_RUInt).size())
                TEST(GetRUIntValue, (*nu_RUInt).at(i), nuCounter);

            if (i < (*nu_RULong64).size())
                TEST(GetRULong64Value, (*nu_RULong64).at(i), nuCounter);

            if (i < (*nu_RTString).size())
                TEST(GetRTStringValue, (*nu_RTString).at(i), nuCounter);

            if (i < (*nu_RFloatVector).size())
                TEST(GetRFloatVectorValue, (*nu_RFloatVector).at(i), nuCounter);

            if (i < (*nu_RIntVector).size())
                TEST(GetRIntVectorValue, (*nu_RIntVector).at(i), nuCounter);

            ++nuCounter;
        }

        // Per-PFO tests
        std::cout << std::endl << "Testing per-PFO parameters" << std::endl;

        TEST_SIZE(*pfo_RFloat, *numPfos);
        TEST_SIZE(*pfo_RInt, *numPfos);
        TEST_SIZE(*pfo_RBool, *numPfos);
        TEST_SIZE(*pfo_RUInt, *numPfos);
        TEST_SIZE(*pfo_RULong64, *numPfos);
        TEST_SIZE(*pfo_RTString, *numPfos);
        TEST_SIZE(*pfo_RFloatVector, *numPfos);
        TEST_SIZE(*pfo_RIntVector, *numPfos);

        for (std::size_t i = 0; i < *numPfos; ++i)
        {
            if (i < (*pfo_RFloat).size())
                TEST(GetRFloatValue, (*pfo_RFloat).at(i), pfoCounter);

            if (i < (*pfo_RInt).size())
                TEST(GetRIntValue, (*pfo_RInt).at(i), pfoCounter);
            
            if (i < (*pfo_RBool).size())
                TEST(GetRBoolValue, (*pfo_RBool).at(i), pfoCounter);

            if (i < (*pfo_RUInt).size())
                TEST(GetRUIntValue, (*pfo_RUInt).at(i), pfoCounter);

            if (i < (*pfo_RULong64).size())
                TEST(GetRULong64Value, (*pfo_RULong64).at(i), pfoCounter);

            if (i < (*pfo_RTString).size())
                TEST(GetRTStringValue, (*pfo_RTString).at(i), pfoCounter);

            if (i < (*pfo_RFloatVector).size())
                TEST(GetRFloatVectorValue, (*pfo_RFloatVector).at(i), pfoCounter);

            if (i < (*pfo_RIntVector).size())
                TEST(GetRIntVectorValue, (*pfo_RIntVector).at(i), pfoCounter);

            ++pfoCounter;
        }

        // Per-primary tests
        std::cout << std::endl << "Testing per-primary parameters" << std::endl;

        TEST_SIZE(*primary_RFloat, *numPrimaries);
        TEST_SIZE(*primary_RInt, *numPrimaries);
        TEST_SIZE(*primary_RBool, *numPrimaries);
        TEST_SIZE(*primary_RUInt, *numPrimaries);
        TEST_SIZE(*primary_RULong64, *numPrimaries);
        TEST_SIZE(*primary_RTString, *numPrimaries);
        TEST_SIZE(*primary_RFloatVector, *numPrimaries);
        TEST_SIZE(*primary_RIntVector, *numPrimaries);

        for (std::size_t i = 0; i < *numPrimaries; ++i)
        {
            if (i < (*primary_RFloat).size())
                TEST(GetRFloatValue, (*primary_RFloat).at(i), primaryCounter);

            if (i < (*primary_RInt).size())
                TEST(GetRIntValue, (*primary_RInt).at(i), primaryCounter);

            if (i < (*primary_RBool).size())
                TEST(GetRBoolValue, (*primary_RBool).at(i), primaryCounter);

            if (i < (*primary_RUInt).size())
                TEST(GetRUIntValue, (*primary_RUInt).at(i), primaryCounter);

            if (i < (*primary_RULong64).size())
                TEST(GetRULong64Value, (*primary_RULong64).at(i), primaryCounter);

            if (i < (*primary_RTString).size())
                TEST(GetRTStringValue, (*primary_RTString).at(i), primaryCounter);

            if (i < (*primary_RFloatVector).size())
                TEST(GetRFloatVectorValue, (*primary_RFloatVector).at(i), primaryCounter);

            if (i < (*primary_RIntVector).size())
                TEST(GetRIntVectorValue, (*primary_RIntVector).at(i), primaryCounter);

            ++primaryCounter;
        }

        // Per-cosmic tests
        std::cout << std::endl << "Testing per-cosmic parameters" << std::endl;

        TEST_SIZE(*cr_RFloat, *numCosmicRays);
        TEST_SIZE(*cr_RInt, *numCosmicRays);
        TEST_SIZE(*cr_RBool, *numCosmicRays);
        TEST_SIZE(*cr_RUInt, *numCosmicRays);
        TEST_SIZE(*cr_RULong64, *numCosmicRays);
        TEST_SIZE(*cr_RTString, *numCosmicRays);
        TEST_SIZE(*cr_RFloatVector, *numCosmicRays);
        TEST_SIZE(*cr_RIntVector, *numCosmicRays);

        for (std::size_t i = 0; i < *numCosmicRays; ++i)
        {
            if (i < (*cr_RFloat).size())
                TEST(GetRFloatValue, (*cr_RFloat).at(i), cosmicCounter);

            if (i < (*cr_RInt).size())
                TEST(GetRIntValue, (*cr_RInt).at(i), cosmicCounter);

            if (i < (*cr_RBool).size())
                TEST(GetRBoolValue, (*cr_RBool).at(i), cosmicCounter);

            if (i < (*cr_RUInt).size())
                TEST(GetRUIntValue, (*cr_RUInt).at(i), cosmicCounter);

            if (i < (*cr_RULong64).size())
                TEST(GetRULong64Value, (*cr_RULong64).at(i), cosmicCounter);

            if (i < (*cr_RTString).size())
                TEST(GetRTStringValue, (*cr_RTString).at(i), cosmicCounter);

            if (i < (*cr_RFloatVector).size())
                TEST(GetRFloatVectorValue, (*cr_RFloatVector).at(i), cosmicCounter);

            if (i < (*cr_RIntVector).size())
                TEST(GetRIntVectorValue, (*cr_RIntVector).at(i), cosmicCounter);

            ++cosmicCounter;
        }

        ++evtCounter;
    }

    // Print summary
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Processed " << TEXT_WHITE_BOLD << evtCounter << " event(s) " << TEXT_NORMAL << "with " << TEXT_GREEN_BOLD << successfulTests
              << " passed test(s)" << TEXT_NORMAL << ", " << TEXT_RED_BOLD << failedTests << " failed test(s) " << TEXT_NORMAL << "and "
              << TEXT_YELLOW_BOLD << warnings <<  " warning(s)" << TEXT_NORMAL << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Float_t GetRFloatValue(const int counter)
{
    return static_cast<Float_t>(-1.234f) + static_cast<Float_t>(counter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

Int_t GetRIntValue(const int counter)
{
    return static_cast<Int_t>(-1234) + static_cast<Int_t>(counter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

Bool_t GetRBoolValue(const int counter)
{
    return static_cast<Bool_t>(counter % 2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

UInt_t GetRUIntValue(const int counter)
{
    return static_cast<UInt_t>(1234U) + static_cast<UInt_t>(counter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ULong64_t GetRULong64Value(const int counter)
{
    return static_cast<ULong64_t>(1234UL) + static_cast<ULong64_t>(counter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TString GetRTStringValue(const int counter)
{
    return TString("1 2 3 4 " + std::to_string(counter));
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<Float_t> GetRFloatVectorValue(const int counter)
{
    return std::vector<Float_t>{1.2f + static_cast<Float_t>(counter), -3.4f + static_cast<Float_t>(counter)};
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<Int_t> GetRIntVectorValue(const int counter)
{
    return std::vector<Int_t>{12 + static_cast<Int_t>(counter), -34 + static_cast<Int_t>(counter)};
}