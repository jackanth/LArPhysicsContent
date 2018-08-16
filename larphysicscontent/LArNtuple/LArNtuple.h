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
     * @brief  Destructor
     */
    ~LArNtuple();

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
     *  @brief  Add a scalar record to the cache
     *
     *  @param  record The record
     */
    void AddScalarRecord(const LArNtupleRecord &record);

    /**
     *  @brief  Add a vector record element to the cache
     *
     *  @param  record The record
     */
    void AddVectorRecordElement(const LArNtupleRecord &record);

    /**
     *  @brief  Fill the vectors using the cached elements
     */
    void FillVectors();

    /**
     *  @brief  Push the vector records onto the stack
     */
    void PushVectors();

    /**
     *  @brief  Fill the TTree using the cache entries and reset the cache
     */
    void Fill();

    /**
     *  @brief  Reset the per-event ntuple state
     */
    void Reset();

    friend class AnalysisNtupleAlgorithm;

public:
    using BranchMap = std::unordered_map<std::string, LArBranchPlaceholder>; ///< Alias for a map from branch names to branch placeholders
    using Cache     = std::deque<std::any>; ///< Alias for a generic cache (deque to preserve pointers to elements)

    TFile *                m_pOutputTFile;            ///< The output TFile
    TTree *                m_pOutputTree;             ///< The output TTree
    BranchMap              m_scalarBranchMap;         ///< The scalar branch map
    BranchMap              m_vectorElementBranchMap;  ///< The vector element branch map
    std::vector<BranchMap> m_vectorBranchMaps;        ///< The vector branch map
    std::size_t            m_vectorBranchMapIndex;    ///< The vector branch map index
    bool                   m_addressesSet;            ///< Whether the addresses have been set
    bool                   m_ntupleEmpty;             ///< Whether the ntuple is empty
    bool                   m_areVectorElementsLocked; ///< Whether scalar entries are locked
    mutable Cache          m_cache;                   ///< The cache

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
     *  @return the current branch map
     */
    BranchMap &GetCurrentVectorBranchMap();

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
};

//------------------------------------------------------------------------------------------------------------------------------------------
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

template <typename TOBJECT>
inline void LArNtuple::AddBranch(const std::string &branchName, std::decay_t<TOBJECT> &obj, const bool splitMode) const
{
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
        if constexpr (std::is_fundamental_v<std::decay_t<TOBJECT>>)
            pBranch->SetAddress(&obj);

        else
        {
            std::decay_t<TOBJECT *> &pObj = this->CacheObject<std::decay_t<TOBJECT *>>(nullptr);
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
        if (branchPlaceholder.CacheIndex() >= m_cache.size())
        {
            std::cerr << "LArNtuple: Cache index was out of bounds" << std::endl;
            throw pandora::STATUS_CODE_FAILURE;
        }

        std::any_cast<TOBJ_D &>(m_cache.at(branchPlaceholder.CacheIndex())) = std::forward<T>(object);
    }

    else
    {
        TOBJ_D &cachedObject = this->CacheObject(std::forward<T>(object));
        branchPlaceholder.CacheIndex(m_cache.size() - 1UL);
        this->AddBranch<TOBJ_D>(branchName, cachedObject, splitMode);
    }
}

} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_H
