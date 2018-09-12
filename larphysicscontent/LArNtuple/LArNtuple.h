/**
 *  @file   larphysicscontent/LArNtuple/LArNtuple.h
 *
 *  @brief  Header file for the lar ntuple class.
 *
 *  $Log: $
 */
#ifndef LAR_NTUPLE_H
#define LAR_NTUPLE_H 1

#include "larphysicscontent/LArNtuple/LArBranchPlaceholder.h"
#include "larphysicscontent/LArNtuple/LArNtupleRecord.h"
#include "larphysicscontent/LArNtuple/NtupleVariableBaseTool.h"
#include "larphysicscontent/LArObjects/LArRootRegistry.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "Objects/ParticleFlowObject.h"
#include "Pandora/Algorithm.h"

#include "TFile.h"
#include "TTree.h"

#include <any>
#include <deque>

namespace lar_physics_content
{

/**
 *  @briefi Forward declaration of the AnalysisNtupleAlgorithm class
 */
class AnalysisNtupleAlgorithm;

/**
 *  @brief  LArNtuple class
 */
class LArNtuple
{
public:
    using VectorRecordMap = std::multimap<std::string, LArNtupleRecord>; ///< Alias for a map from branch names to records for vector records

    /**
     * @brief  Deleted copy constructor
     */
    LArNtuple(const LArNtuple &) = delete;

    /**
     * @brief  Default move constructor
     */
    LArNtuple(LArNtuple &&) = default;

    /**
     * @brief  Deleted copy assignment operator
     */
    LArNtuple &operator=(const LArNtuple &) = delete;

    /**
     * @brief  Default move assignment operator
     */
    LArNtuple &operator=(LArNtuple &&) = default;

    /**
     * @brief  Virtual destructor
     */
    virtual ~LArNtuple();

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  filePath the TTree file path
     *  @param  treeName the TTree name
     *  @param  treeTitle the TTree title
     */
    LArNtuple(const std::string &filePath, const std::string &treeName, const std::string &treeTitle, const bool appendMode);

    /**
     *  @brief  Get an MC particle
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    virtual const pandora::MCParticle *GetMCParticle(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get an MC cosmic
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    virtual const pandora::MCParticle *GetMCCosmic(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get an MC primary
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    virtual const pandora::MCParticle *GetMCPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get an MC neutrino
     *
     *  @param  pPfo address of the PFO
     *  @param  pfoList the list of all PFOs
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    virtual const pandora::MCParticle *GetMCNeutrino(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    friend class AnalysisNtupleAlgorithm;
    friend class NtupleVariableBaseTool;

private:
    using MCParticleMapFn = std::function<const pandora::MCParticle *(const pandora::MCParticle *const)>; ///< Alias for an MC particle map function
    using BranchMap = std::unordered_map<std::string, LArBranchPlaceholder>; ///< Alias for a map from branch names to branch placeholders
    using Cache     = std::deque<std::any>; ///< Alias for a generic cache (deque to preserve pointers to elements)
    using VectorBranchTypeMap =
        std::unordered_map<LArNtupleHelper::VECTOR_BRANCH_TYPE, BranchMap>; ///< Alias for a map from vector branch types to their branch map

    template <typename T>
    using PfoCache =
        std::unordered_map<const pandora::ParticleFlowObject *, std::decay_t<T>>; ///< Alias for a cache from PFO addresses to other objects

    template <typename T>
    using RecordMapGetter = std::function<LArBranchPlaceholder::NtupleRecordMap<const std::decay_t<T> *>(const LArBranchPlaceholder &)>; ///< Alias for a record map getter function

    TFile *                                              m_pOutputTFile;           ///< The output TFile
    TTree *                                              m_pOutputTree;            ///< The output TTree
    BranchMap                                            m_scalarBranchMap;        ///< The scalar branch map
    BranchMap                                            m_vectorElementBranchMap; ///< The vector element branch map
    VectorBranchTypeMap                                  m_vectorBranchMaps; ///< The map from vector branch types to their branch maps
    bool                                                 m_addressesSet;     ///< Whether the addresses have been set
    bool                                                 m_ntupleEmpty;      ///< Whether the ntuple is empty
    bool                                                 m_areVectorElementsLocked;   ///< Whether scalar entries are locked
    unsigned int                                         m_trackSlidingFitWindow;     ///< The track sliding fit window size
    mutable Cache                                        m_cache;                     ///< The cache
    mutable PfoCache<const pandora::MCParticle *>        m_cacheMCParticles;          ///< The cached mappings from PFOs to MC particles
    mutable PfoCache<const pandora::MCParticle *>        m_cacheMCCosmics;            ///< The cached mappings from PFOs to MC cosmics
    mutable PfoCache<const pandora::MCParticle *>        m_cacheMCPrimaries;          ///< The cached mappings from PFOs to MC primaries
    mutable PfoCache<const pandora::MCParticle *>        m_cacheMCNeutrinos;          ///< The cached mappings from PFOs to MC neutrinos
    mutable PfoCache<pandora::CaloHitList>               m_cacheDownstreamThreeDHits; ///< The pfo cache of downstream 3D hits
    mutable PfoCache<pandora::CaloHitList>               m_cacheDownstreamUHits;      ///< The pfo cache of downstream U hits
    mutable PfoCache<pandora::CaloHitList>               m_cacheDownstreamVHits;      ///< The pfo cache of downstream V hits
    mutable PfoCache<pandora::CaloHitList>               m_cacheDownstreamWHits;      ///< The pfo cache of downstream W hits
    mutable PfoCache<pandora::PfoList>                   m_cacheDownstreamPfos;       ///< The pfo cache of downstream pfos
    mutable PfoCache<LArNtupleHelper::TrackFitSharedPtr> m_cacheTrackFits;            ///< The pfo cache of track fits
    std::shared_ptr<LArRootRegistry>                     m_spRegistry;                ///< The ROOT registry

    /**
     *  @brief  Add a scalar record to the cache
     *
     *  @param  record the record
     */
    void AddScalarRecord(const LArNtupleRecord &record);

    /**
     *  @brief  Add a vector record element to the cache
     *
     *  @param  record the record
     *  @param  type the vector type
     */
    void AddVectorRecordElement(const LArNtupleRecord &record, const LArNtupleHelper::VECTOR_BRANCH_TYPE type);

    /**
     *  @brief  Fill the vectors using the cached elements
     *
     *  @param  type the vector type
     */
    void FillVectors(const LArNtupleHelper::VECTOR_BRANCH_TYPE type);

    /**
     *  @brief  Push the vector records onto the stack
     *
     *  @param  type the vector type
     */
    void PushVectors(const LArNtupleHelper::VECTOR_BRANCH_TYPE type);

    /**
     *  @brief  Fill the TTree using the cache entries and reset the cache
     */
    void Fill();

    /**
     *  @brief  Reset the per-event ntuple state
     */
    void Reset();

    /**
     *  @brief  Get all 3D hits downstream of a PFO, including from the PFO itself
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamThreeDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all 2D hits downstream of a PFO, including from the PFO itself
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    pandora::CaloHitList GetAllDownstreamTwoDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all PFOs downstream of a PFO, including the PFO itself
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the downstream PFOs
     */
    const pandora::PfoList &GetAllDownstreamPfos(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all U hits downstream of a PFO, including from the PFO itself
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamUHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all V hits downstream of a PFO, including from the PFO itself
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamVHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Get all W hits downstream of a PFO, including from the PFO itself
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamWHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Create a vector of objects from a vector of shared pointers to their records
     *
     *  @param  records the record shared pointers
     *
     *  @return the vector
     */
    template <typename T>
    std::vector<std::decay_t<T>> CreateVector(const std::vector<LArBranchPlaceholder::NtupleRecordSPtr> &records) const;

    /**
     *  @brief  Cache an object an return a reference to the cached object
     *
     *  @param  obj the object to cache
     *
     *  @return reference to the cached object
     */
    template <typename T>
    std::decay_t<T> &CacheObject(T &&obj) const;

    /**
     *  @brief  Get the current vector branch map
     *
     *  @param  type the vector type
     *
     *  @return the current branch map
     */
    BranchMap &GetVectorBranchMap(const LArNtupleHelper::VECTOR_BRANCH_TYPE type);

    /**
     *  @brief  Add a branch or set the branch address in advance of the first fill
     *
     *  @param  branchName the branch name
     *  @param  obj the object whose address to use
     *  @param  splitMode whether to create the branch in split mode
     */
    template <typename TOBJECT>
    void AddBranch(const std::string &branchName, std::decay_t<TOBJECT> &obj, const bool splitMode) const;

    /**
     *  @brief  Push a scalar to a branch
     *
     *  @param  branchName the branch name
     *  @param  branchPlaceholder the branch placeholder
     *  @param  splitMode whether to create the branch in split mode
     */
    template <typename T>
    void PushScalarToBranch(const std::string &branchName, LArBranchPlaceholder &branchPlaceholder, const bool splitMode) const;

    /**
     *  @brief  Push a vector to a branch
     *
     *  @param  branchName the branch name
     *  @param  branchPlaceholder the branch placeholder
     *  @param  splitMode whether to create the branch in split mode
     */
    template <typename T>
    void PushVectorToBranch(const std::string &branchName, LArBranchPlaceholder &branchPlaceholder, const bool splitMode) const;

    /**
     *  @brief  Push an object to a branch (i.e. cache it in the right place and add the branch if required)
     *
     *  @param  branchName the branch name
     *  @param  branchPlaceholder the branch placeholder
     *  @param  object the object to push
     *  @param  splitMode whether to create the branch in split mode
     */
    template <typename TOBJ, typename T>
    void PushToBranch(const std::string &branchName, LArBranchPlaceholder &branchPlaceholder, T &&object, const bool splitMode) const;

    /**
     *  @brief  Set the vector branch addresses
     *
     *  @return the number of records added
     */
    std::size_t SetVectorBranchAddresses();

    /**
     *  @brief  Set the scalar branch addresses
     *
     *  @return the number of records added
     */
    std::size_t SetScalarBranchAddresses();

    /**
     *  @brief  Validate a record and add it to the ntuple
     *
     *  @param  branchMap the branch map to populate
     *  @param  record the record
     */
    void ValidateAndAddRecord(BranchMap &branchMap, const LArNtupleRecord &record);

    /**
     *  @brief  Instantiate the TFile object
     *
     *  @param  appendMode whether we are in append mode
     *  @param  filePath the file path
     */
    void InstantiateTFile(const bool appendMode, const std::string &filePath);

    /**
     *  @brief  Instantiate the TTree object
     *
     *  @param  appendMode whether we are in append mode
     *  @param  treeName the name of the tree
     *  @param  treeTitle the title of the tree
     *  @param  filePath the file path
     */
    void InstantiateTTree(const bool appendMode, const std::string &treeName, const std::string &treeTitle, const std::string &filePath);

    /**
     *  @brief  Get an MC particle (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    const pandora::MCParticle *GetMCParticleWrapper(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get an MC cosmic (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    const pandora::MCParticle *GetMCCosmicWrapper(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get an MC primary (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    const pandora::MCParticle *GetMCPrimaryWrapper(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get an MC neutrino (wrapper method)
     *
     *  @param  pPfo address of the PFO
     *  @param  pMCParticleList optional pointer to the MC particle list
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    const pandora::MCParticle *GetMCNeutrinoWrapper(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList);

    /**
     *  @brief  Get the MCParticle weight map for a set of CaloHits
     *
     *  @param  caloHitList the list of CaloHits
     *  @param  mapFn the MC particle mapping function
     *
     *  @return the MCParticle weight map
     */
    pandora::MCParticleWeightMap GetMCParticleWeightMap(const pandora::CaloHitList &caloHitList, const MCParticleMapFn &mapFn) const;

    /**
     *  @brief  Get all hits downstream of a PFO (implementation method)
     *
     *  @param  pPfo address of the PFO
     *  @param  hitType the hit type
     *  @param  the cache to use
     *
     *  @return the hits
     */
    const pandora::CaloHitList &GetAllDownstreamHitsImpl(
        const pandora::ParticleFlowObject *const pPfo, const pandora::HitType hitType, PfoCache<pandora::CaloHitList> &cache) const;

    /**
     *  @brief  Get an MC particle (implementation method)
     *
     *  @param  caloHitList the CaloHit list
     *  @param  mapFn the MC particle mapping function
     *
     *  @return address of the corresponding MCParticle, if one can be found
     */
    const pandora::MCParticle *GetMCParticleImpl(const pandora::CaloHitList &caloHitList, const MCParticleMapFn &mapFn) const;

    /**
     *  @brief  Get all 2D hits of a PFO
     *
     *  @param  pPfo address of the PFO
     *
     *  @return the hits
     */
    pandora::CaloHitList GetAllTwoDHits(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Wrapper for checking a cache before getting an object
     *
     *  @param  pPfo address of the PFO
     *  @param  cache the cache
     *  @parma  getter the getter function
     *
     *  @return the object
     */
    template <typename T>
    const std::decay_t<T> &CacheWrapper(const pandora::ParticleFlowObject *const pPfo, PfoCache<std::decay_t<T>> &cache,
        const std::function<std::decay_t<T>()> &getter) const;

    /**
     *  @brief  Retrieve a scalar record
     *
     *  @param  branchName the branch name
     *
     *  @return the record shared pointer
     */
    LArBranchPlaceholder::NtupleRecordSPtr GetScalarRecord(const std::string &branchName) const;

    /**
     *  @brief  Retrieve a vector record element
     *
     *  @param  type the vector branch type
     *  @param  branchName the branch name
     *  @param  pPfo address of the associated PFO
     *
     *  @return the record shared pointer
     */
    LArBranchPlaceholder::NtupleRecordSPtr GetVectorRecordElement(
        const LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Retrieve a vector record element
     *
     *  @param  type the vector branch type
     *  @param  branchName the branch name
     *  @param  pMCParticle address of the associated MC particle
     *
     *  @return the record shared pointer
     */
    LArBranchPlaceholder::NtupleRecordSPtr GetVectorRecordElement(
        const LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Retrieve a vector record element (implementation method)
     *
     *  @param  type the vector branch type
     *  @param  branchName the branch name
     *  @param  pParticle address of the associated particle
     *  @param  recordMapGetter the record map getter function
     *
     *  @return the record shared pointer
     */
    template <typename TPARTICLE>
    LArBranchPlaceholder::NtupleRecordSPtr GetVectorRecordElementImpl(const LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName,
        const std::decay_t<TPARTICLE> *const pParticle, const RecordMapGetter<std::decay_t<TPARTICLE>> &recordMapGetter) const;

    /**
     *  @brief  Retrieve a branch placeholder
     *
     *  @param  type the vector branch type
     *  @param  branchName the branch name
     *
     *  @return the branch placeholder
     */
    const LArBranchPlaceholder &GetBranchPlaceholder(const LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName) const;

    /**
     *  @brief  Retrieve a branch placeholder
     *
     *  @param  branchMap the branch map
     *  @param  branchName the branch name
     *
     *  @return the branch placeholder
     */
    const LArBranchPlaceholder &GetBranchPlaceholder(const BranchMap &branchMap, const std::string &branchName) const;

    /**
     *  @brief  Get a track fit, first trying to use the cache
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pPfo address of the PFO
     *
     *  @return shared pointer to the track fit, or nullptr if fit fails
     */
    const LArNtupleHelper::TrackFitSharedPtr &GetTrackFit(const pandora::Pandora &pandoraInstance, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Calculate a track fit
     *
     *  @param  pandoraInstance the instance of Pandora
     *  @param  pPfo address of the PFO
     *  @param  slidingFitWindow the sliding fit window
     *
     *  @return shared pointer to the track fit, or nullptr if fit fails
     */
    LArNtupleHelper::TrackFitSharedPtr CalculateTrackFit(
        const pandora::Pandora &pandoraInstance, const pandora::ParticleFlowObject *const pPfo, const unsigned int slidingFitWindow) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArNtuple::GetMCCosmic(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const)
{
    const pandora::CaloHitList caloHitList = this->GetAllDownstreamTwoDHits(pPfo);
    return this->GetMCParticleImpl(caloHitList, lar_content::LArMCParticleHelper::GetPrimaryMCParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArNtuple::GetMCPrimary(const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const)
{
    const pandora::CaloHitList caloHitList = this->GetAllDownstreamTwoDHits(pPfo);
    return this->GetMCParticleImpl(caloHitList, lar_content::LArMCParticleHelper::GetPrimaryMCParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArNtuple::GetAllDownstreamThreeDHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_3D, m_cacheDownstreamThreeDHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArNtuple::GetAllDownstreamUHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_VIEW_U, m_cacheDownstreamUHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArNtuple::GetAllDownstreamVHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_VIEW_V, m_cacheDownstreamVHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &LArNtuple::GetAllDownstreamWHits(const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetAllDownstreamHitsImpl(pPfo, pandora::TPC_VIEW_W, m_cacheDownstreamWHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
std::vector<std::decay_t<T>> LArNtuple::CreateVector(const std::vector<LArBranchPlaceholder::NtupleRecordSPtr> &records) const
{
    using T_D = std::decay_t<T>;
    std::vector<T_D> vector;

    for (const LArBranchPlaceholder::NtupleRecordSPtr &spRecord : records)
        vector.emplace_back(spRecord->Value<T_D>());

    return vector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline std::decay_t<T> &LArNtuple::CacheObject(T &&obj) const
{
    m_cache.emplace_back(obj);
    return std::any_cast<std::decay_t<T> &>(m_cache.back());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArNtuple::BranchMap &LArNtuple::GetVectorBranchMap(const LArNtupleHelper::VECTOR_BRANCH_TYPE type)
{
    return m_vectorBranchMaps.emplace(type, BranchMap()).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TOBJECT>
inline void LArNtuple::AddBranch(const std::string &branchName, std::decay_t<TOBJECT> &obj, const bool splitMode) const
{
    using TOBJECT_D = std::decay_t<TOBJECT>;

    if (m_ntupleEmpty)
        m_pOutputTree->Branch(branchName.c_str(), &obj, 32000, splitMode ? 99 : -1);

    else // we have loaded this non-empty TTree from a file and now need to tie up all the branches with new addresses
    {
        TBranch *pBranch = m_pOutputTree->GetBranch(branchName.c_str());

        if (!pBranch)
        {
            std::cerr << "LArNtuple: Could not append to existing TTree because no existing branch matched '" << branchName << "'" << std::endl;
            throw pandora::STATUS_CODE_NOT_FOUND;
        }

        // For any non-fundamental types, we need to cache a nullptr to the cached object and use a pointer to this cached nullptr as the
        // argument to SetAddress, then set the nullptr to point to the object
        if constexpr (std::is_fundamental_v<TOBJECT_D>)
            pBranch->SetAddress(&obj);

        else
        {
            TOBJECT_D *&pObj = this->CacheObject<TOBJECT_D *>(nullptr);
            pBranch->SetAddress(&pObj);
            pObj = &obj;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void LArNtuple::PushScalarToBranch(const std::string &branchName, LArBranchPlaceholder &branchPlaceholder, const bool splitMode) const
{
    using T_D = std::decay_t<T>;
    this->PushToBranch<T_D>(branchName, branchPlaceholder, branchPlaceholder.GetNtupleScalarRecord()->Value<T_D>(), splitMode);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void LArNtuple::PushVectorToBranch(const std::string &branchName, LArBranchPlaceholder &branchPlaceholder, const bool splitMode) const
{
    using T_D = std::decay_t<T>;
    this->PushToBranch<std::vector<T_D>>(branchName, branchPlaceholder, this->CreateVector<T_D>(branchPlaceholder.GetNtupleVectorRecord()), splitMode);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TOBJ, typename T>
void LArNtuple::PushToBranch(const std::string &branchName, LArBranchPlaceholder &branchPlaceholder, T &&object, const bool splitMode) const
{
    using TOBJ_D = std::decay_t<TOBJ>;

    if (m_addressesSet)
    {
        std::any *const pCacheElement = branchPlaceholder.CacheElement();

        if (!pCacheElement)
        {
            std::cerr << "LArNtuple: Cache element address has not been set" << std::endl;
            throw pandora::STATUS_CODE_FAILURE;
        }

        std::any_cast<TOBJ_D &>(*pCacheElement) = std::forward<T>(object);
    }

    else
    {
        TOBJ_D &cachedObject = this->CacheObject(std::forward<T>(object));
        branchPlaceholder.CacheElement(&m_cache.back());
        this->AddBranch<TOBJ_D>(branchName, cachedObject, splitMode);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArNtuple::GetMCParticleWrapper(
    const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList)
{
    return this->CacheWrapper<const pandora::MCParticle *>(
        pPfo, m_cacheMCParticles, [&]() { return this->GetMCParticle(pPfo, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArNtuple::GetMCCosmicWrapper(
    const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList)
{
    return this->CacheWrapper<const pandora::MCParticle *>(pPfo, m_cacheMCCosmics, [&]() { return this->GetMCCosmic(pPfo, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArNtuple::GetMCPrimaryWrapper(
    const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList)
{
    return this->CacheWrapper<const pandora::MCParticle *>(
        pPfo, m_cacheMCPrimaries, [&]() { return this->GetMCPrimary(pPfo, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArNtuple::GetMCNeutrinoWrapper(
    const pandora::ParticleFlowObject *const pPfo, const pandora::MCParticleList *const pMCParticleList)
{
    return this->CacheWrapper<const pandora::MCParticle *>(
        pPfo, m_cacheMCNeutrinos, [&]() { return this->GetMCNeutrino(pPfo, pMCParticleList); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::decay_t<T> &LArNtuple::CacheWrapper(
    const pandora::ParticleFlowObject *const pPfo, PfoCache<std::decay_t<T>> &cache, const std::function<std::decay_t<T>()> &getter) const
{
    // Check the cache first (can return nullptr)
    const auto findIter = cache.find(pPfo);

    if (findIter != cache.end())
        return findIter->second;

    // Not in cache; find it and cache it
    return cache.emplace(pPfo, getter()).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArBranchPlaceholder::NtupleRecordSPtr LArNtuple::GetVectorRecordElement(
    const LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const pandora::ParticleFlowObject *const pPfo) const
{
    return this->GetVectorRecordElementImpl<pandora::ParticleFlowObject>(
        type, branchName, pPfo, [](const LArBranchPlaceholder &branchPlaceholder) { return branchPlaceholder.GetPfoRecordMap(); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArBranchPlaceholder::NtupleRecordSPtr LArNtuple::GetVectorRecordElement(
    const LArNtupleHelper::VECTOR_BRANCH_TYPE type, const std::string &branchName, const pandora::MCParticle *const pMCParticle) const
{
    return this->GetVectorRecordElementImpl<pandora::MCParticle>(type, branchName, pMCParticle,
        [](const LArBranchPlaceholder &branchPlaceholder) { return branchPlaceholder.GetMCParticleRecordMap(); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TPARTICLE>
LArBranchPlaceholder::NtupleRecordSPtr LArNtuple::GetVectorRecordElementImpl(const LArNtupleHelper::VECTOR_BRANCH_TYPE type,
    const std::string &branchName, const std::decay_t<TPARTICLE> *const pParticle, const RecordMapGetter<std::decay_t<TPARTICLE>> &recordMapGetter) const
{
    const LArBranchPlaceholder &branchPlaceholder = this->GetBranchPlaceholder(type, branchName);
    const auto &                pfoRecordMap      = recordMapGetter(branchPlaceholder);
    const auto                  recordMapFindIter = pfoRecordMap.find(pParticle);

    if (recordMapFindIter == pfoRecordMap.end())
    {
        std::cerr << "LArNtuple: Could not find PFO amongst vector elements for branch '" << branchName << "'" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    const LArBranchPlaceholder::NtupleRecordSPtr spRecord = recordMapFindIter->second;

    if (!spRecord)
    {
        std::cerr << "LArNtuple: Found vector record element by name '" << branchName << "' but it had not yet been filled" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }

    return spRecord;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArNtupleHelper::TrackFitSharedPtr &LArNtuple::GetTrackFit(
    const pandora::Pandora &pandoraInstance, const pandora::ParticleFlowObject *const pPfo) const
{
    return this->CacheWrapper<LArNtupleHelper::TrackFitSharedPtr>(
        pPfo, m_cacheTrackFits, [&]() { return this->CalculateTrackFit(pandoraInstance, pPfo, m_trackSlidingFitWindow); });
}

} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_H