/**
 *  @file   larphysicscontent/LArAnalysis/CommonNtupleTool.h
 *
 *  @brief  Header file for the common ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_COMMON_NTUPLE_TOOL_H
#define LAR_COMMON_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArHelpers/LArAnalysisHelper.h"
#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

namespace lar_physics_content
{
/**
 *  @brief  CommonNtupleTool class
 */
class CommonNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CommonNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    CommonNtupleTool(const CommonNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    CommonNtupleTool(CommonNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    CommonNtupleTool &operator=(const CommonNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    CommonNtupleTool &operator=(CommonNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~CommonNtupleTool() = default;

protected:
    std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Produce records generic to every considered particle class
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceGenericPfoRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList) const;

    /**
     *  @brief  Produce records generic to every considered particle class except neutrinos
     *
     *  @param  pPfo optional address of the PFO
     *  @param  pfoList the list of all PFOs
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceNonNeutrinoPfoRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList) const;

    /**
     *  @brief  Get the fraction of fiducial 3D hits
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the fraction of hits
     */
    float GetFractionOfFiducialThreeDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get the direction at the vertex for a track
     *
     *  @param  pPfo address of the PFO
     *  @param  pVertex address of the vertex
     *
     *  @return the direction
     */
    pandora::CartesianVector GetTrackDirectionAtVertex(const pandora::ParticleFlowObject *const pPfo, const pandora::Vertex *const pVertex) const;

    /**
     *  @brief  Get the direction at the vertex for a shower
     *
     *  @param  pPfo address of the PFO
     *  @param  pVertex address of the vertex
     *
     *  @return the direction
     */
    pandora::CartesianVector GetShowerDirectionAtVertex(const pandora::ParticleFlowObject *const pPfo, const pandora::Vertex *const pVertex) const;

    /**
     *  @brief  Get the single vertex of a PFO, otherwise return nullptr
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the vertex
     */
    const pandora::Vertex *GetSingleVertex(const pandora::ParticleFlowObject *const pPfo) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector CommonNtupleTool::GetTrackDirectionAtVertex(
    const pandora::ParticleFlowObject *const pPfo, const pandora::Vertex *const pVertex) const
{
    if (const LArNtupleHelper::TrackFitSharedPtr &spTrackFit = this->GetTrackFit(pPfo))
        return LArAnalysisHelper::GetFittedDirectionAtPosition(*spTrackFit, pVertex->GetPosition());

    return pandora::CartesianVector(0.f, 0.f, 0.f);
}

} // namespace lar_physics_content

#endif // #ifndef LAR_COMMON_NTUPLE_TOOL_H
