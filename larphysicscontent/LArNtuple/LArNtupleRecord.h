/**
 *  @file   larphysicscontent/LArNtple/LArNtupleRecord.h
 *
 *  @brief  Header file for the lar ntuple record class.
 *
 *  $Log: $
 */
#ifndef LAR_NTUPLE_RECORD_H
#define LAR_NTUPLE_RECORD_H 1

#include "Pandora/AlgorithmHeaders.h"

#include "Rtypes.h"
#include "TString.h"

#include <cassert>
#include <variant>

namespace lar_physics_content
{

/**
 *  @brief  Forward declaration of the LArNtuple class
 */
class LArNtuple;

/**
 *  @brief  Forward declaration of the LArBranchPlaceholder class
 */
class LArBranchPlaceholder;

/**
 *  @brief  Forward declaration of the NtupleVariableBaseTool class
 */
class NtupleVariableBaseTool;

/**
 *  @brief  LArNtupleRecord class
 */
class LArNtupleRecord
{
public:
    /**
     *  @brief The different kinds of value type
     */
    enum class VALUE_TYPE : unsigned
    {
        R_FLOAT        = 0U, ///< The ROOT float type
        R_INT          = 1U, ///< The ROOT int type
        R_BOOL         = 2U, ///< The ROOT bool type
        R_UINT         = 3U, ///< The ROOT uint type
        R_ULONG64      = 4U, ///< The ROOT ulong64 type
        R_TSTRING      = 5U, ///< The ROOT TString type
        R_FLOAT_VECTOR = 6U, ///< A vector of ROOT float types
        R_INT_VECTOR   = 7U  ///< A vector of ROOT int types
        R_FLOAT_MATRIX = 8U, ///< A 2D matrix of ROOT float types
        R_INT_MATRIX   = 9U  ///< A 2D matrix of ROOT int types
    };

    using RFloat       = Float_t;              ///< Alias for a ROOT float type
    using RInt         = Int_t;                ///< Alias for a ROOT int type
    using RBool        = Bool_t;               ///< Alias for a ROOT bool type
    using RUInt        = UInt_t;               ///< Alias for a ROOT uint type
    using RULong64     = ULong64_t;            ///< Alias for a ROOT ulong64 type
    using RTString     = TString;              ///< Alias for a ROOT TString type
    using RFloatVector = std::vector<Float_t>; ///< Alias for a vector of ROOT float types
    using RIntVector   = std::vector<Int_t>;   ///< Alias for a vector of ROOT int types
    using RFloatMatrix = std::vector<std::vector<Float_t>>; ///< Alias for a 2D matrix of ROOT float types
    using RIntMatrix   = std::vector<std::vector<Int_t>>;   ///< Alias for a 2D matrix of ROOT int types

    /**
     *  @brief  Constructor
     *
     *  @param  branchName the branch name
     *  @param  value the value
     *  @param  writeToNtuple whether to write this record to the ntuple
     */
    template <typename TVALUE, typename = std::enable_if_t<std::is_same_v<std::decay_t<TVALUE>, RFloat> || std::is_same_v<std::decay_t<TVALUE>, RInt> ||
                                                           std::is_same_v<std::decay_t<TVALUE>, RBool> || std::is_same_v<std::decay_t<TVALUE>, RUInt> ||
                                                           std::is_same_v<std::decay_t<TVALUE>, RULong64> || std::is_same_v<std::decay_t<TVALUE>, RTString> ||
                                                           std::is_same_v<std::decay_t<TVALUE>, RFloatVector> || std::is_same_v<std::decay_t<TVALUE>, RIntVector> ||
                                                           std::is_same_v<std::decay_t<TVALUE>, RFloatMatrix> || std::is_same_v<std::decay_t<TVALUE>, RIntMatrix>>>
    LArNtupleRecord(std::string branchName, TVALUE &&value, const bool writeToNtuple = true) noexcept;

    /**
     * @brief  Default copy constructor
     */
    LArNtupleRecord(const LArNtupleRecord &) = default;

    /**
     * @brief  Default move constructor
     */
    LArNtupleRecord(LArNtupleRecord &&) = default;

    /**
     * @brief  Default copy assignment operator
     */
    LArNtupleRecord &operator=(const LArNtupleRecord &) = default;

    /**
     * @brief  Default move assignment operator
     */
    LArNtupleRecord &operator=(LArNtupleRecord &&) = default;

    /**
     * @brief  Default destructor
     */
    ~LArNtupleRecord() = default;

    /**
     *  @brief  Get the branch name
     *
     *  @return the branch name
     */
    const std::string &BranchName() const noexcept;

protected:
    /**
     *  @brief  Get the value type
     *
     *  @return the value type
     */
    VALUE_TYPE ValueType() const noexcept;

    /**
     *  @brief  Get the value
     *
     *  @return the value
     */
    template <typename TVALUE>
    const std::decay_t<TVALUE> &Value() const;

    /**
     *  @brief  Add a prefix to the branch name
     *
     *  @param  prefix the prefix
     */
    void AddBranchNamePrefix(const std::string &prefix);

    /**
     *  @brief  Get the PFO
     *
     *  @return address of the PFO if there is one
     */
    const pandora::ParticleFlowObject *GetPfo() const noexcept;

    /**
     *  @brief  Set the PFO
     *
     *  @param  pPfo address of the PFO
     */
    void SetPfo(const pandora::ParticleFlowObject *const pPfo) noexcept;

    /**
     *  @brief  Get the MC particle
     *
     *  @return address of the MC particle if there is one
     */
    const pandora::MCParticle *GetMCParticle() const noexcept;

    /**
     *  @brief  Set the MC particle
     *
     *  @param  pMCParticle address of the MC particle
     */
    void SetMCParticle(const pandora::MCParticle *const pMCParticle) noexcept;

    /**
     *  @brief  Get whether to write this record to the ntuple
     *
     *  @return whether to write this to the ntuple
     */
    bool WriteToNtuple() const noexcept;

    friend class LArNtuple;
    friend class LArBranchPlaceholder;
    friend class NtupleVariableBaseTool;

private:
    using VariantType = std::variant<RFloat, RInt, RBool, RUInt, RULong64, RTString, RFloatVector, RIntVector, RFloatMatrix, RIntMatrix>; ///< Alias for the variant type

    VALUE_TYPE                         m_valueType;     ///< The value type
    std::string                        m_branchName;    ///< The branch name
    VariantType                        m_value;         ///< The value
    const pandora::ParticleFlowObject *m_pPfo;          ///< Address of the associated PFO
    const pandora::MCParticle *        m_pMCParticle;   ///< Address of the associated MC particle
    bool                               m_writeToNtuple; ///< Whether to write this record to the ntuple
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TVALUE, typename>
LArNtupleRecord::LArNtupleRecord(std::string branchName, TVALUE &&value, const bool writeToNtuple) noexcept :
    m_valueType(VALUE_TYPE::R_FLOAT),
    m_branchName(std::move_if_noexcept(branchName)),
    m_value(value),
    m_pPfo(nullptr),
    m_pMCParticle(nullptr),
    m_writeToNtuple(writeToNtuple)
{
    using TVALUE_D = std::decay_t<TVALUE>;

    if constexpr (std::is_same_v<TVALUE_D, RFloat>)
        m_valueType = VALUE_TYPE::R_FLOAT;

    else if (std::is_same_v<TVALUE_D, RInt>)
        m_valueType = VALUE_TYPE::R_INT;

    else if (std::is_same_v<TVALUE_D, RBool>)
        m_valueType = VALUE_TYPE::R_BOOL;

    else if (std::is_same_v<TVALUE_D, RUInt>)
        m_valueType = VALUE_TYPE::R_UINT;

    else if (std::is_same_v<TVALUE_D, RULong64>)
        m_valueType = VALUE_TYPE::R_ULONG64;

    else if (std::is_same_v<TVALUE_D, RTString>)
        m_valueType = VALUE_TYPE::R_TSTRING;

    else if (std::is_same_v<TVALUE_D, RFloatVector>)
        m_valueType = VALUE_TYPE::R_FLOAT_VECTOR;

    else if (std::is_same_v<TVALUE_D, RIntVector>)
        m_valueType = VALUE_TYPE::R_INT_VECTOR;

    else if (std::is_same_v<TVALUE_D, RFloatMatrix>)
        m_valueType = VALUE_TYPE::R_FLOAT_MATRIX;

    else if (std::is_same_v<TVALUE_D, RIntMatrix>)
        m_valueType = VALUE_TYPE::R_INT_MATRIX;

    else // unreachable
        assert(false && "LArNtupleRecord: Unknown value type");
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &LArNtupleRecord::BranchName() const noexcept
{
    return m_branchName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArNtupleRecord::VALUE_TYPE LArNtupleRecord::ValueType() const noexcept
{
    return m_valueType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TVALUE>
const std::decay_t<TVALUE> &LArNtupleRecord::Value() const
{
    try
    {
        return std::get<std::decay_t<TVALUE>>(m_value);
    }

    catch (const std::bad_variant_access &)
    {
        std::cerr << "LArNtupleRecord: Invalid value type for branch '" << m_branchName << "'" << std::endl;
        throw pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    catch (...)
    {
        std::cerr << "LArNtupleRecord: Failed to get value for branch '" << m_branchName << "'" << std::endl;
        throw pandora::STATUS_CODE_FAILURE;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArNtupleRecord::AddBranchNamePrefix(const std::string &prefix)
{
    if (prefix.empty())
    {
        std::cerr << "LArNtupleRecord: Cannot add empty branch name prefix" << std::endl;
        throw pandora::STATUS_CODE_NOT_ALLOWED;
    }

    m_branchName = prefix + m_branchName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *LArNtupleRecord::GetPfo() const noexcept
{
    return m_pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArNtupleRecord::SetPfo(const pandora::ParticleFlowObject *const pPfo) noexcept
{
    m_pPfo = pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *LArNtupleRecord::GetMCParticle() const noexcept
{
    return m_pMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArNtupleRecord::SetMCParticle(const pandora::MCParticle *const pMCParticle) noexcept
{
    m_pMCParticle = pMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArNtupleRecord::WriteToNtuple() const noexcept
{
    return m_writeToNtuple;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_RECORD_H