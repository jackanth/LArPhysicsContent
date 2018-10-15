/**
 *  @file   test/TestNtupleTool.cc
 *
 *  @brief  Implementation of the test ntuple tool class.
 *
 *  $Log: $
 */

#include "test/TestNtupleTool.h"
#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

namespace lar_physics_content
{
TestNtupleTool::TestNtupleTool() :
    NtupleVariableBaseTool(),
    m_eventCounter(0),
    m_neutrinoCounter(0),
    m_cosmicCounter(0),
    m_primaryCounter(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestNtupleTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return NtupleVariableBaseTool::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> TestNtupleTool::ProcessEvent(
    const pandora::PfoList &pfoList, const std::vector<std::shared_ptr<LArInteractionValidationInfo>> &eventValidationInfo)
{
    std::cout << "TestNtupleTool: Processing event (" << pfoList.size() << " PFOs, " << eventValidationInfo.size() << " interactions)" << std::endl;

    ++m_eventCounter;
    m_neutrinoCounter = 0;
    m_cosmicCounter   = 0;
    m_primaryCounter  = 0;
    return GetTestRecords(m_eventCounter - 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> TestNtupleTool::ProcessNeutrino(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
    const std::shared_ptr<LArInteractionValidationInfo> &spMcInteraction)
{
    std::cout << "TestNtupleTool: Processing neutrino (pPfo = " << pPfo << ", spMcInteraction = " << spMcInteraction << ", "
              << pfoList.size() << " PFOs)" << std::endl;

    return GetTestRecords(m_neutrinoCounter++);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> TestNtupleTool::ProcessPrimary(
    const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    std::cout << "TestNtupleTool: Processing primary (pPfo = " << pPfo << ", spMcTarget = " << spMcTarget << ", " << pfoList.size()
              << " PFOs)" << std::endl;

    return GetTestRecords(m_primaryCounter++);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> TestNtupleTool::ProcessCosmicRay(
    const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, const std::shared_ptr<LArMCTargetValidationInfo> &spMcTarget)
{
    std::cout << "TestNtupleTool: Processing cosmic ray (pPfo = " << pPfo << ", spMcTarget = " << spMcTarget << ", " << pfoList.size()
              << " PFOs)" << std::endl;

    return GetTestRecords(m_cosmicCounter++);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<LArNtupleRecord> TestNtupleTool::GetTestRecords(const int counter) const
{
    std::vector<LArNtupleRecord> records;

    records.emplace_back("RFloat", static_cast<LArNtupleRecord::RFloat>(-1.234) + static_cast<LArNtupleRecord::RFloat>(counter));
    records.emplace_back("RInt", static_cast<LArNtupleRecord::RInt>(-1234) + static_cast<LArNtupleRecord::RInt>(counter));

    records.emplace_back("RBool", static_cast<LArNtupleRecord::RBool>(counter % 2));
    records.emplace_back("RUInt", static_cast<LArNtupleRecord::RUInt>(1234U) + static_cast<LArNtupleRecord::RUInt>(counter));

    records.emplace_back("RULong64", static_cast<LArNtupleRecord::RULong64>(1234UL) + static_cast<LArNtupleRecord::RULong64>(counter));

    records.emplace_back("RTString", LArNtupleRecord::RTString("1 2 3 4 " + std::to_string(counter)));

    records.emplace_back("RFloatVector", LArNtupleRecord::RFloatVector{1.2f + static_cast<LArNtupleRecord::RFloat>(counter),
                                             -3.4f + static_cast<LArNtupleRecord::RFloat>(counter)});
    records.emplace_back("RIntVector",
        LArNtupleRecord::RIntVector{12 + static_cast<LArNtupleRecord::RInt>(counter), -34 + static_cast<LArNtupleRecord::RInt>(counter)});

    return records;
}

} // namespace lar_physics_content
