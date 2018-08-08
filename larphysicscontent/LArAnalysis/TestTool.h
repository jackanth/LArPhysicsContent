/**
 *  @file   larphysicscontent/LArAnalysis/TestTool.h
 *
 *  @brief  Header file for the test tool class.
 *
 *  $Log: $
 */
#ifndef LAR_TEST_TOOL_H
#define LAR_TEST_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief  TestTool class
 */
class TestTool : public AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TestTool();

    bool Run(const Algorithm *const pAlgorithm);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_TEST_TOOL_H
