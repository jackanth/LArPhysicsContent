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

#include <any>
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
    template <typename T>
    using NtupleRecordMap = std::unordered_map<std::decay_t<T>, NtupleRecordSPtr>; ///< Alias for a map to ntuple record shared pointers

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
     *  @return the cache element pointer
     */
    std::any *CacheElement() const noexcept;

    /**
     *  @brief  Set the cache element pointer
     *
     *  @param  pCacheElement address of the element
     */
    void CacheElement(std::any *const pCacheElement) noexcept;

    /**
     *  @brief  Get the PFO record map
     *
     *  @return the record map
     */
    const NtupleRecordMap<const pandora::ParticleFlowObject *> &GetPfoRecordMap() const noexcept;

    /**
     *  @brief  Get the MC particle record map
     *
     *  @return the record map
     */
    const NtupleRecordMap<const pandora::MCParticle *> &GetMCParticleRecordMap() const noexcept;

    friend class LArNtuple;

private:
    NtupleRecordSPtr                                     m_spNtupleScalarRecord; ///< Shared pointer to the ntuple scalar record
    std::vector<NtupleRecordSPtr>                        m_ntupleVectorRecord;   ///< The vector record shared pointers
    LArNtupleRecord::VALUE_TYPE                          m_valueType;            ///< The branch's value type
    std::any *                                           m_pCacheElement;        ///< The cache element pointer
    NtupleRecordMap<const pandora::ParticleFlowObject *> m_pfoRecordMap;         ///< The map from PFOs to record shared pointers
    NtupleRecordMap<const pandora::MCParticle *>         m_mcParticleRecordMap;  ///< The map from MCParticles to record shared pointers

    /**
     *  @brief  Populate the maps using the vector entries
     */
    void PopulateMaps();
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArBranchPlaceholder::ClearNtupleRecords() noexcept
{
    m_ntupleVectorRecord.clear();
    m_spNtupleScalarRecord = nullptr;
    m_pfoRecordMap.clear();
    m_mcParticleRecordMap.clear();
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

inline void LArBranchPlaceholder::CacheElement(std::any *const pCacheElement) noexcept
{
    m_pCacheElement = pCacheElement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::any *LArBranchPlaceholder::CacheElement() const noexcept
{
    return m_pCacheElement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArBranchPlaceholder::NtupleRecordMap<const pandora::ParticleFlowObject *> &LArBranchPlaceholder::GetPfoRecordMap() const noexcept
{
    return m_pfoRecordMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArBranchPlaceholder::NtupleRecordMap<const pandora::MCParticle *> &LArBranchPlaceholder::GetMCParticleRecordMap() const noexcept
{
    return m_mcParticleRecordMap;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_BRANCH_PLACEHOLDER
