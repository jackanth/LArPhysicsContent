/**
 *  @file   larpandoracontent/LArPhysicsContent.cc
 *
 *  @brief  Factory implementations for physics content
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/Algorithm.h"
#include "Pandora/Pandora.h"

#include "AnalysisAlgorithm.h"
#include "AnalysisDataAlgorithm.h"
#include "WriteAnalysisParticlesAlgorithm.h"
#include "BirksHitSelectionTool.h"

#include "LArPhysicsContent.h"

#define LAR_ALGORITHM_LIST(d)                                                                                                   \
    d("LArAnalysis",                AnalysisAlgorithm)                                                                          \
    d("LArWriteAnalysisParticles",  WriteAnalysisParticlesAlgorithm)                                                            \
    d("LArAnalysisData",            AnalysisDataAlgorithm)

#define LAR_ALGORITHM_TOOL_LIST(d)                                                                                              \
    d("LArBirksHitSelection", BirksHitSelectionTool)

#define FACTORY Factory

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_physics_content
{

#define LAR_CONTENT_CREATE_ALGORITHM_FACTORY(a, b)                                                                              \
class b##FACTORY : public pandora::AlgorithmFactory                                                                             \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::Algorithm *CreateAlgorithm() const {return new b;};                                                                \
};

LAR_ALGORITHM_LIST(LAR_CONTENT_CREATE_ALGORITHM_FACTORY)

//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY(a, b)                                                                         \
class b##FACTORY : public pandora::AlgorithmToolFactory                                                                         \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::AlgorithmTool *CreateAlgorithmTool() const {return new b;};                                                        \
};

LAR_ALGORITHM_TOOL_LIST(LAR_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY)

} // namespace lar_physics_content

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_REGISTER_ALGORITHM(a, b)                                                                                     \
{                                                                                                                                \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmFactory(pandora, a, new lar_physics_content::b##FACTORY)); \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                              \
        return statusCode;                                                                                                       \
}

#define LAR_CONTENT_REGISTER_ALGORITHM_TOOL(a, b)                                                                                    \
{                                                                                                                                    \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmToolFactory(pandora, a, new lar_physics_content::b##FACTORY)); \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                                  \
        return statusCode;                                                                                                           \
}

pandora::StatusCode LArPhysicsContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_ALGORITHM_LIST(LAR_CONTENT_REGISTER_ALGORITHM);
    LAR_ALGORITHM_TOOL_LIST(LAR_CONTENT_REGISTER_ALGORITHM_TOOL);
    return pandora::STATUS_CODE_SUCCESS;
}