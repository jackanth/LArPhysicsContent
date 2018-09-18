/**
 *  @file   larphysicscontent/LinkDef.h
 *
 *  @brief  Header file for linking custom ROOT classes.
 *
 *  $Log: $
 */

#ifndef LAR_LINK_DEF_H
#define LAR_LINK_DEF_H 1

#include "TString.h"
#include "Rtypes.h"
#include <vector>

#ifdef __MAKECINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class std::vector<Int_t>+;
#pragma link C++ class std::vector<Float_t>+;
#pragma link C++ class std::vector<UInt_t>+;
#pragma link C++ class std::vector<UInt64_t>+;
#pragma link C++ class std::vector<Bool_t>+;
#pragma link C++ class std::vector<TString>+;

#pragma link C++ class std::vector<std::vector<Int_t>>+;
#pragma link C++ class std::vector<std::vector<Float_t>>+;
#endif // #ifdef __MAKECINT__

#endif // #ifndef LAR_LINK_DEF_H