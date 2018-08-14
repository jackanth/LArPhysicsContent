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
#include <vector>

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class vector<int>+;
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<unsigned int>+;
#pragma link C++ class vector<unsigned long>+;
#pragma link C++ class vector<unsigned long long>+;
#pragma link C++ class vector<bool>+;
#pragma link C++ class vector<TString>+;

#pragma link C++ class vector<vector<int>>+;
#pragma link C++ class vector<vector<float>>+;
#endif

#endif // #ifndef LAR_LINK_DEF_H