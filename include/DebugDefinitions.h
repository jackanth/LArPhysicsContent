/**
 *  @file LArPhysicsContent/include/DebugDefinitions.h
 *
 *  @brief Temporary header file for debug-related definitions.
 * 
 *  $Log: $
 *
 */

#ifndef LAR_DEBUG_DEFINITIONS_H
#define LAR_DEBUG_DEFINITIONS_H 1

#include <iostream>

#define TEXT_NORMAL       "\033[0m"
#define TEXT_BOLD         "\033[1m"

#define TEXT_RED          "\033[0;31m"
#define TEXT_GREEN        "\033[0;32m"
#define TEXT_MAGENTA      "\033[0;35m"
#define TEXT_WHITE        "\033[0;37m"

#define TEXT_RED_BOLD     "\033[1;31m"
#define TEXT_GREEN_BOLD   "\033[1;32m"
#define TEXT_MAGENTA_BOLD "\033[1;35m"
#define TEXT_WHITE_BOLD   "\033[1;37m"

#define COUT(a) std::cout << a << TEXT_NORMAL << std::endl
#define CERR(a) std::cerr << TEXT_RED << a << TEXT_NORMAL << std::endl

#endif // #ifndef LAR_DEBUG_DEFINITIONS_H