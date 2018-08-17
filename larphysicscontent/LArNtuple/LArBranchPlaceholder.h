/**
 *  @file   larphysicscontent/LArNtuple/LArBranchPlaceholder.h
 *
 *  @brief  Header file for the lar branch placeholder class.
 *
 *  $Log: $
 */
#ifndef LAR_BRANCH_PLACEHOLDER
#define LAR_BRANCH_PLACEHOLDER 1

#include "larphysicscontent/LArNtuple/LArNtupleRecord.h"
#include <memory>

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the LArNtuple class
 */
class LArNtuple;

/**
 *  @brief  Branch placeholder class
 */
class LArBranchPlaceholder
{
public:
    using NtupleRecordSPtr = std::shared_ptr<LArNtupleRecord>; ///< Alias for a shared pointer to an ntuple record

    /**
     * @brief  Default copy constructor
     */
    LArBranchPlaceholder(const LArBranchPlaceholder &) = default;

    /**
     * @brief  Default move constructor
     */
    LArBranchPlaceholder(LArBranchPlaceholder &&) = default;

    /**
     * @brief  Default copy assignment operator
     */
    LArBranchPlaceholder &operator=(const LArBranchPlaceholder &) = default;

    /**
     * @brief  Default move assignment operator
     */
    LArBranchPlaceholder &operator=(LArBranchPlaceholder &&) = default;

    /**
     * @brief  Default destructor
     */
    ~LArBranchPlaceholder() = default;

    /**
     *  @brief  Set the ntuple scalar record
     *
     *  @param  record the record to add
     */
    void SetNtupleScalarRecord(const LArNtupleRecord &record);

    /**
     *  @brief  Push the stored ntuple scalar record into the vector
     */
    void PushNtupleScalarRecord();

    /**
     *  @brief  Clear the ntuple records
     */
    void ClearNtupleRecords() noexcept;

    /**
     *  @brief  Get the vector of LArNtuple records
     *
     *  @return the vector of ntuple record shared pointers
     */
    const std::vector<NtupleRecordSPtr> &GetNtupleVectorRecord() const noexcept;

    /**
     *  @brief  Get the scalar record shared pointer if one exists
     *
     *  @return the record shared pointer, or nullptr if there is none
     */
    NtupleRecordSPtr GetNtupleScalarRecord() const;

protected:
    /**
     *  @brief  Constructor
     *
     *  @param  record the first LArNtupleRecord
     */
    explicit LArBranchPlaceholder(const LArNtupleRecord &record) noexcept;

    /**
     *  @brief  Get the value type
     *
     *  @return the value type
     */
    LArNtupleRecord::VALUE_TYPE ValueType() const noexcept;

    /**
     *  @brief  Get the cache index if it has been set
     *
     *  @return the cache index
     */
    std::size_t CacheIndex() const;

    /**
     *  @brief  Set the cache index
     *
     *  @param cacheIndex the cache index
     */
    void CacheIndex(const std::size_t cacheIndex) noexcept;

    /**
     *  @brief  Get whether the cache index has been set
     *
     *  @return whether the cache index has been set
     */
    bool CacheIndexSet() const noexcept;

    friend class LArNtuple;

private:
    NtupleRecordSPtr              m_spNtupleScalarRecord; ///< Shared pointer to the ntuple scalar record
    std::vector<NtupleRecordSPtr> m_ntupleVectorRecord;   ///< The vector record shared pointers
    LArNtupleRecord::VALUE_TYPE   m_valueType;            ///< The branch's value type
    std::size_t                   m_cacheIndex;           ///< The cache index
    bool                          m_cacheIndexSet;        ///< Whether the cache index has been set
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArBranchPlaceholder::ClearNtupleRecords() noexcept
{
    m_ntupleVectorRecord.clear();
    m_spNtupleScalarRecord = nullptr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::vector<LArBranchPlaceholder::NtupleRecordSPtr> &LArBranchPlaceholder::GetNtupleVectorRecord() const noexcept
{
    return m_ntupleVectorRecord;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArBranchPlaceholder::NtupleRecordSPtr LArBranchPlaceholder::GetNtupleScalarRecord() const
{
    return m_spNtupleScalarRecord;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArNtupleRecord::VALUE_TYPE LArBranchPlaceholder::ValueType() const noexcept
{
    return m_valueType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArBranchPlaceholder::CacheIndex(const std::size_t cacheIndex) noexcept
{
    m_cacheIndex    = cacheIndex;
    m_cacheIndexSet = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArBranchPlaceholder::CacheIndexSet() const noexcept
{
    return m_cacheIndexSet;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_BRANCH_PLACEHOLDER
