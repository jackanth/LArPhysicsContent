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

using namespace pandora;

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
    enum class VALUE_TYPE
    {
        R_FLOAT,        ///< The ROOT float type
        R_INT,          ///< The ROOT int type
        R_BOOL,         ///< The ROOT bool type
        R_UINT,         ///< The ROOT uint type
        R_ULONG64,      ///< The ROOT ulong64 type
        R_TSTRING,      ///< The ROOT TString type
        R_FLOAT_VECTOR, ///< A vector of ROOT float types
        R_INT_VECTOR    ///< A vector of ROOT int types
    };

    using RFloat       = Float_t;              ///< Alias for a ROOT float type
    using RInt         = Int_t;                ///< Alias for a ROOT int type
    using RBool        = Bool_t;               ///< Alias for a ROOT bool type
    using RUInt        = UInt_t;               ///< Alias for a ROOT uint type
    using RULong64     = ULong64_t;            ///< Alias for a ROOT ulong64 type
    using RTString     = TString;              ///< Alias for a ROOT TString type
    using RFloatVector = std::vector<Float_t>; ///< Alias for a vector of ROOT float types
    using RIntVector   = std::vector<Int_t>;   ///< Alias for a vector of ROOT int types

    /**
     *  @brief  Constructor
     *
     *  @param  branchName The branch name
     *  @param  value The value
     */
    template <typename TVALUE, typename = std::enable_if_t<std::is_same_v<std::decay_t<TVALUE>, RFloat> || std::is_same_v<std::decay_t<TVALUE>, RInt> ||
                                                           std::is_same_v<std::decay_t<TVALUE>, RBool> || std::is_same_v<std::decay_t<TVALUE>, RUInt> ||
                                                           std::is_same_v<std::decay_t<TVALUE>, RULong64> || std::is_same_v<std::decay_t<TVALUE>, RTString> ||
                                                           std::is_same_v<std::decay_t<TVALUE>, RFloatVector> || std::is_same_v<std::decay_t<TVALUE>, RIntVector>>>
    LArNtupleRecord(std::string branchName, TVALUE &&value) noexcept;

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
     *  @return The branch name
     */
    const std::string &BranchName() const noexcept;

protected:
    /**
     *  @brief  Get the value type
     *
     *  @return The value type
     */
    VALUE_TYPE ValueType() const noexcept;

    /**
     *  @brief  Get the value
     *
     *  @return The value
     */
    template <typename TVALUE>
    std::decay_t<TVALUE> Value() const;

    /**
     *  @brief  Add a prefix to the branch name
     *
     *  @param  prefix the prefix
     */
    void AddBranchNamePrefix(const std::string &prefix);

    friend class LArNtuple;
    friend class LArBranchPlaceholder;
    friend class NtupleVariableBaseTool;

private:
    using VariantType = std::variant<RFloat, RInt, RBool, RUInt, RULong64, RTString, RFloatVector, RIntVector>; ///< Alias for the variant type

    VALUE_TYPE  m_valueType;  ///< The value type
    std::string m_branchName; ///< The branch name
    VariantType m_value;      ///< The value
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TVALUE, typename>
LArNtupleRecord::LArNtupleRecord(std::string branchName, TVALUE &&value) noexcept :
    m_valueType(VALUE_TYPE::R_FLOAT),
    m_branchName(std::move_if_noexcept(branchName)),
    m_value(value)
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
std::decay_t<TVALUE> LArNtupleRecord::Value() const
{
    try
    {
        return std::get<std::decay_t<TVALUE>>(m_value);
    }

    catch (const std::bad_variant_access &)
    {
        std::cerr << "LArNtupleRecord: Invalid value type" << std::endl;
        throw pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    catch (...)
    {
        std::cerr << "LArNtupleRecord: Failed to get value" << std::endl;
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

    m_branchName = prefix + '_' + m_branchName;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_NTUPLE_RECORD_H