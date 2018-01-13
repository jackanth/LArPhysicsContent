/**
 *  @file LArPhysicsContent/include/LArTrackHitEnergy.h
 *
 *  @brief Header file for the track hit energy class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_HIT_ENERGY_H
#define LAR_TRACK_HIT_ENERGY_H 1

#include "Objects/CaloHit.h"

#include <vector>

using namespace pandora;

namespace lar_physics_content
{
/**
 *  @brief LArTrackHitEnergy class.
 * 
 */
class LArTrackHitEnergy
{
public:
    using Vector = std::vector<LArTrackHitEnergy>; ///< Alias for a vector of LArTrackHitEnergy objects.

    /**
     *  @brief Constructor.
     * 
     *  @param pCaloHit Pointer to the corresponding CaloHit.
     *  @param coordinate The projected coordinate along the track fit (cm).
     *  @param uncorrectedEnergy The hit energy before recombination correction (MeV).
     *  @param correctedEnergy The hit energy after recombination correction (MeV).
     * 
     */
    LArTrackHitEnergy(const CaloHit *const pCaloHit, const float coordinate, const float uncorrectedEnergy,
                   const float correctedEnergy) noexcept;
    
    /**
     *  @brief Get the CaloHit.
     * 
     *  @return Pointer to the CaloHit.
     * 
     */
    const CaloHit * GetCaloHit() const noexcept;
    
    /**
     *  @brief Get the coordinate.
     * 
     *  @return The coordinate.
     * 
     */
    float Coordinate() const noexcept;
    
    /**
     *  @brief Get the uncorrected energy.
     * 
     *  @return The uncorrected energy.
     * 
     */
    float UncorrectedEnergy() const noexcept;
    
    /**
     *  @brief Get the corrected energy.
     * 
     *  @return The corrected energy.
     * 
     */
    float CorrectedEnergy() const noexcept;
    
    /**
     *  @brief Get the apply-correction boolean.
     * 
     *  @return The apply-correction boolean.
     * 
     */
    bool ApplyCorrection() const noexcept;
    
    /**
     *  @brief Set the apply-correction boolean.
     * 
     *  @param applyCorrection The apply-correction boolean.
     * 
     */
    void ApplyCorrection(const bool applyCorrection) noexcept;
    
    private:
    const CaloHit * m_pCaloHit;          ///< Pointer to the corresponding CaloHit.
    float           m_coordinate;        ///< The coordinate along the track fit (cm).
    float           m_uncorrectedEnergy; ///< The hit energy, not corrected for recombination (MeV).
    float           m_correctedEnergy;   ///< The hit energy, corrected for recombination (MeV).
    bool            m_applyCorrection;   ///< Whether to use the recombination-corrected energy.
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArTrackHitEnergy::LArTrackHitEnergy(const CaloHit *const pCaloHit, const float coordinate, const float uncorrectedEnergy,
                                     const float correctedEnergy) noexcept :
    m_pCaloHit{pCaloHit},
    m_coordinate{coordinate},
    m_uncorrectedEnergy{uncorrectedEnergy},
    m_correctedEnergy{correctedEnergy},
    m_applyCorrection{true}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CaloHit * LArTrackHitEnergy::GetCaloHit() const noexcept
{
    return this->m_pCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArTrackHitEnergy::Coordinate() const noexcept
{
    return this->m_coordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArTrackHitEnergy::UncorrectedEnergy() const noexcept
{
    return this->m_uncorrectedEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArTrackHitEnergy::CorrectedEnergy() const noexcept
{
    return this->m_correctedEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArTrackHitEnergy::ApplyCorrection() const noexcept
{
    return this->m_applyCorrection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArTrackHitEnergy::ApplyCorrection(const bool applyCorrection) noexcept
{
    this->m_applyCorrection = applyCorrection;
}
} // namespace lar_physics_content

#endif // #ifndef LAR_TRACK_HIT_ENERGY_H
