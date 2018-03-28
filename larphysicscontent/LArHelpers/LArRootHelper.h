/**
 *  @file   larphysicscontent/LArHelpers/LArRootHelper.h
 *
 *  @brief  Header file for the lar root helper class.
 *
 *  $Log: $
 */
#ifndef LAR_ROOT_HELPER_H
#define LAR_ROOT_HELPER_H 1

#include "TNtuple.h"

namespace lar_physics_content
{

/**
 *  @brief LArRootHelper class.
 */
class LArRootHelper
{
public:
    /**
     *  @brief  Write a TNtuple object to a file
     *
     *  @param  pNtuple address of the TNtuple object
     *  @param  fileName the path to the new ROOT file
     *  @param  verboseMode whether to save in verbose mode
     */
    static void WriteNTuple(TNtuple *const pNtuple, const std::string &fileName, const bool verboseMode);

    /**
     *  @brief  Load a TNtuple object from a file
     *
     *  @param  filePath the path to the ROOT file
     *  @param  nTupleName the name of the TNtuple object to load
     *
     *  @return address of the loaded TNtuple object
     */
    static TNtuple * LoadNTupleFromFile(const std::string &filePath, const std::string &nTupleName);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_ROOT_HELPER_H
