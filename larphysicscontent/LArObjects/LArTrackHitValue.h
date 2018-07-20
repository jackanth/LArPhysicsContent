/**
 *  @file   larphysicscontent/LArObjects/LArTrackHitValue.h
 *
 *  @brief  Header file for the track hit value class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_HIT_VALUE_H
#define LAR_TRACK_HIT_VALUE_H 1

#include "Objects/CaloHit.h"

#include <vector>

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief LArTrackHitValue class
 */
class LArTrackHitValue
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pCaloHit address of the corresponding CaloHit
     *  @param  coordinate the projected coordinate along the track fit (cm)
     *  @param  threeDDistance the 3D distance this hit corresponds to (cm)
     *  @param  caloValue a calorimetric quantity of the hit, e.g. energy or integrated ADC
     */
    LArTrackHitValue(const CaloHit *const pCaloHit, const float coordinate, const float threeDDistance, const float caloValue) noexcept;

    /**
     *  @brief  Get the CaloHit
     *
     *  @return address of the CaloHit
     */
    const CaloHit *GetCaloHit() const noexcept;

    /**
     *  @brief  Get the coordinate
     *
     *  @return the coordinate
     */
    float Coordinate() const noexcept;

    /**
     *  @brief  Get the 3D distance this his corresponds to
     *
     *  @return the 3D distance
     */
    float ThreeDDistance() const noexcept;

    /**
     *  @brief  Get the calorimetric quantity
     *
     *  @return the calorimetric quantity
     */
    float CaloValue() const noexcept;

    /**
     *  @brief  Set the calorimetric quantity
     *
     *  @param  caloValue the calorimetric quantity
     */
    void CaloValue(const float caloValue) noexcept;

private:
    const CaloHit *m_pCaloHit; ///< Pointer to the corresponding CaloHit
    float m_coordinate;        ///< The coordinate along the track fit (cm)
    float m_threeDDistance;    ///< The 3D distance this hit corresponds to (cm)
    float m_caloValue;         ///< The hit's calorimetric quantity, e.g. energy or integrated ADC
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArTrackHitValue::LArTrackHitValue(const CaloHit *const pCaloHit, const float coordinate, const float threeDDistance, const float caloValue) noexcept
    : m_pCaloHit(pCaloHit), m_coordinate(coordinate), m_threeDDistance(threeDDistance), m_caloValue(caloValue)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CaloHit *LArTrackHitValue::GetCaloHit() const noexcept
{
    return m_pCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArTrackHitValue::Coordinate() const noexcept
{
    return m_coordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArTrackHitValue::ThreeDDistance() const noexcept
{
    return m_threeDDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArTrackHitValue::CaloValue() const noexcept
{
    return m_caloValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArTrackHitValue::CaloValue(const float caloValue) noexcept
{
    m_caloValue = caloValue;
}

} // namespace lar_physics_content

#endif // #ifndef LAR_TRACK_HIT_VALUE_H
