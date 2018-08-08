/**
 *  @file   larphysicscontent/LArAnalysis/TestTool.cc
 *
 *  @brief  Implementation of the test tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larphysicscontent/LArAnalysis/TestTool.h"

using namespace pandora;

namespace lar_physics_content
{
TestTool::TestTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TestTool::Run(const Algorithm *const pAlgorithm)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    (void)xmlHandle;
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_physics_content
