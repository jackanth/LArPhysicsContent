/**
 *  @file   larphysicscontent/LArAnalysis/LeeAnalysisNtupleTool.h
 *
 *  @brief  Header file for the LEE analysis ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_LEE_ANALYSIS_NTUPLE_TOOL_H
#define LAR_LEE_ANALYSIS_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

namespace lar_physics_content
{
/**
 *  @brief  LeeAnalysisNtupleTool class
 */
class LeeAnalysisNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    LeeAnalysisNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    LeeAnalysisNtupleTool(const LeeAnalysisNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    LeeAnalysisNtupleTool(LeeAnalysisNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    LeeAnalysisNtupleTool &operator=(const LeeAnalysisNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    LeeAnalysisNtupleTool &operator=(LeeAnalysisNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~LeeAnalysisNtupleTool() = default;

protected:
    std::vector<LArNtupleRecord> ProcessEvent(
        const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArInteractionValidationInfo> &spInteractionInfo) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget) override;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_physics_content

#endif // #ifndef LAR_LEE_ANALYSIS_NTUPLE_TOOL_H
