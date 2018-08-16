/**
 *  @file   larphysicscontent/LArAnalysis/DefaultNtupleTool.h
 *
 *  @brief  Header file for the default ntuple tool class.
 *
 *  $Log: $
 */
#ifndef LAR_DEFAULT_NTUPLE_TOOL_H
#define LAR_DEFAULT_NTUPLE_TOOL_H 1

#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"

namespace lar_physics_content
{
/**
 *  @brief  DefaultNtupleTool class
 */
class DefaultNtupleTool : public NtupleVariableBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DefaultNtupleTool();

    /**
     *  @brief  Default copy constructor
     */
    DefaultNtupleTool(const DefaultNtupleTool &) = default;

    /**
     *  @brief  Default move constructor
     */
    DefaultNtupleTool(DefaultNtupleTool &&) = default;

    /**
     *  @brief  Default copy assignment operator
     */
    DefaultNtupleTool &operator=(const DefaultNtupleTool &) = default;

    /**
     *  @brief  Default move assignment operator
     */
    DefaultNtupleTool &operator=(DefaultNtupleTool &&) = default;

    /**
     *  @brief  Default destructor
     */
    ~DefaultNtupleTool() = default;

protected:
    std::vector<LArNtupleRecord> ProcessEvent(const pandora::PfoList &pfoList, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessCosmicRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

    std::vector<LArNtupleRecord> ProcessPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) override;

private:
    template <typename T>
    using PfoCache =
        std::unordered_map<const pandora::ParticleFlowObject *, std::decay_t<T>>; ///< Alias for a cache from PFO addresses to other objects

    pandora::CartesianVector m_fiducialCutLowMargins;  ///< The low-coordinate margins for the fiducial cut
    pandora::CartesianVector m_fiducialCutHighMargins; ///< The high-coordinate margins for the fiducial cut
    pandora::CartesianVector m_minFiducialCoordinates; ///< The minimum fiducial coordinates
    pandora::CartesianVector m_maxFiducialCoordinates; ///< The maximum fiducial coordinates

    mutable PfoCache<pandora::CaloHitList> m_cacheDownstreamThreeDHits; ///< The pfo cache of downstream 3D hits
    mutable PfoCache<pandora::CaloHitList> m_cacheDownstreamUHits;      ///< The pfo cache of downstream U hits
    mutable PfoCache<pandora::CaloHitList> m_cacheDownstreamVHits;      ///< The pfo cache of downstream V hits
    mutable PfoCache<pandora::CaloHitList> m_cacheDownstreamWHits;      ///< The pfo cache of downstream W hits

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Produce records generic to every considered particle class
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceGenericPfoRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList) const;

    /**
     *  @brief  Produce MC records generic to every considered particle class
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticle optional address of the MC particle
     *  @param  pMCParticleList optional list of all MC particles
     *
     *  @return the records
     */
    std::vector<LArNtupleRecord> ProduceGenericPfoMCRecords(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList,
        const pandora::MCParticle *const pMCParticle, const pandora::MCParticleList *const pMCParticleList) const;

    float GetFractionOfFiducialThreeDHits(const pandora::ParticleFlowObject *const pPfo) const;

    const pandora::CaloHitList &GetAllDownstreamThreeDHits(const pandora::ParticleFlowObject *const pPfo) const;

    pandora::CaloHitList GetAllDownstreamTwoDHits(const pandora::ParticleFlowObject *const pPfo) const;

    const pandora::CaloHitList &GetAllDownstreamUHits(const pandora::ParticleFlowObject *const pPfo) const;

    const pandora::CaloHitList &GetAllDownstreamVHits(const pandora::ParticleFlowObject *const pPfo) const;

    const pandora::CaloHitList &GetAllDownstreamWHits(const pandora::ParticleFlowObject *const pPfo) const;

    const pandora::CaloHitList &GetAllDownstreamHitsImpl(
        const pandora::ParticleFlowObject *const pPfo, const pandora::HitType hitType, PfoCache<pandora::CaloHitList> &cache) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &DefaultNtupleTool::GetAllDownstreamThreeDHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_3D, m_cacheDownstreamThreeDHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &DefaultNtupleTool::GetAllDownstreamUHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_VIEW_U, m_cacheDownstreamUHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &DefaultNtupleTool::GetAllDownstreamVHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_VIEW_V, m_cacheDownstreamVHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &DefaultNtupleTool::GetAllDownstreamWHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_VIEW_W, m_cacheDownstreamWHits);
}

} // namespace lar_physics_content

#endif // #ifndef LAR_DEFAULT_NTUPLE_TOOL_H
