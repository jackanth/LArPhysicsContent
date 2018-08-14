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
     *  @brief  Push the ntuple scalar record
     *
     *  @param  record the record to add
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
     *  @brief  Assuming this is a scalar record (so the number of records must be 0 or 1), get the single record shared pointer if exists
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
    LArBranchPlaceholder(const LArNtupleRecord &record) noexcept;

    /**
     *  @brief  Get the value type
     *
     *  @return the value type
     */
    LArNtupleRecord::VALUE_TYPE ValueType() const noexcept;

    std::size_t CacheIndex() const;
    void CacheIndex(const std::size_t cacheIndex);

    friend class LArNtuple;

private:
    NtupleRecordSPtr m_spNtupleScalarRecord;            ///< Shared pointer to the ntuple scalar record
    std::vector<NtupleRecordSPtr> m_ntupleVectorRecord; ///< The vector record shared pointers
    LArNtupleRecord::VALUE_TYPE m_valueType;            ///< The branch's value type
    std::size_t m_cacheIndex;
    bool m_cacheIndexSet;
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

inline std::size_t LArBranchPlaceholder::CacheIndex() const
{
    if (m_cacheIndexSet)
        return m_cacheIndex;

    std::cerr << "LArBranchPlaceholder: ..." << std::endl;
    throw STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArBranchPlaceholder::CacheIndex(const std::size_t cacheIndex)
{
    m_cacheIndex = cacheIndex;
    m_cacheIndexSet = true;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_BRANCH_PLACEHOLDER
