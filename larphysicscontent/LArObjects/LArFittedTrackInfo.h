/**
 *  @file   larphysicscontent/LArObjects/LArFittedTrackInfo.h
 *
 *  @brief  Header file for the lar fitted track info class.
 *
 *  $Log: $
 */
#ifndef LAR_FITTED_TRACK_INFO_H
#define LAR_FITTED_TRACK_INFO_H 1

#include "larphysicscontent/LArObjects/LArTrackHitValue.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace pandora;
using namespace lar_content;

namespace lar_physics_content
{

/**
 *  @brief LArFittedTrackInfo class
 */
class LArFittedTrackInfo
{
    public:
    using TrackHitValueVector  = std::vector<LArTrackHitValue>; ///< Alias for a vector of LArTrackHitValues
    /**
     *  @brief  Constructor
     *
     *  @param  pPfo address of the PFO
     *  @param  hitChargeVector the hit charge vector
     *  @param  fit the 3D fit
     *  @param  range the range
     */
    LArFittedTrackInfo(const ParticleFlowObject *const pPfo, TrackHitValueVector hitChargeVector, ThreeDSlidingFitResult fit, const float range) noexcept;

    /**
     *  @brief  Get the address of the PFO
     *
     *  @return the PFO address
     */
    const ParticleFlowObject * Pfo() const noexcept;

    /**
     *  @brief  Get the hit charge vector
     *
     *  @return the hit charge vector
     */
    const TrackHitValueVector & HitChargeVector() const noexcept;

    /**
     *  @brief  Get the 3D fit
     *
     *  @return the 3D fit
     */
    const ThreeDSlidingFitResult & Fit() const noexcept;

    /**
     *  @brief  Get the range
     *
     *  @return the range
     */
    float Range() const noexcept;

private:
    const ParticleFlowObject    *m_pPfo;               ///< Address of the corresponding PFO
    TrackHitValueVector          m_hitChargeVector;    ///< The sorted vector of hit charges
    ThreeDSlidingFitResult       m_fit;                ///< The fit result
    float                        m_range;              ///< The range
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArFittedTrackInfo::LArFittedTrackInfo(const ParticleFlowObject *const pPfo, TrackHitValueVector hitChargeVector, ThreeDSlidingFitResult fit,
    const float range) noexcept :
    m_pPfo(pPfo),
    m_hitChargeVector(std::move_if_noexcept(hitChargeVector)),
    m_fit(std::move_if_noexcept(fit)),
    m_range(range)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ParticleFlowObject * LArFittedTrackInfo::Pfo() const noexcept
{
    return m_pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArFittedTrackInfo::TrackHitValueVector & LArFittedTrackInfo::HitChargeVector() const noexcept
{
    return m_hitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDSlidingFitResult & LArFittedTrackInfo::Fit() const noexcept
{
    return m_fit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArFittedTrackInfo::Range() const noexcept
{
    return m_range;
}
} // namespace lar_physics_content

#endif // #ifndef LAR_FITTED_TRACK_INFO_H
