/**
 *  @file   test/TestNtupleTool.h
 *
 *  @brief  Header file for the test ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_TEST_NTUPLE_TOOL_H
#define LAR_TEST_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief  TestTool class
 */
class TestNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TestNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    TestNtupleTool(const TestNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    TestNtupleTool(TestNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    TestNtupleTool &operator=(const TestNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    TestNtupleTool &operator=(TestNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~TestNtupleTool() = default;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_eventCounter;    ///< The event counter
    int m_neutrinoCounter; ///< The neutrino counter
    int m_pfoCounter;      ///< The PFO counter
    int m_cosmicCounter;   ///< The cosmic counter
    int m_primaryCounter;  ///< The primary counter

    std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessParticle(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    /**
     *  @brief  Get the set of test records for a given counter
     *
     *  @param  counter the counter value
     *
     *  @return the test records
     */
    std::vector<LArNtupleRecord> GetTestRecords(const int counter) const;
};

} // namespace lar_physics_content

#endif // #ifndef LAR_TEST_NTUPLE_TOOL_H
