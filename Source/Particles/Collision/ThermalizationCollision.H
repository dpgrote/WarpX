/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_THERMALIZATIONCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_THERMALIZATIONCOLLISION_H_

#include "CollisionBase.H"
#include "Particles/MultiParticleContainer_fwd.H"

#include <AMReX_REAL.H>

#include <iosfwd>

using namespace amrex::literals;

namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::ParticleReal LambertW(amrex::ParticleReal x) {
        // Calculate the Lambert W on the n=-1 branch
        // using Fritsch's iterationmethod.
        // Starting with wn=-1.1 puts the result on the appropriate branch.
        // See https://arxiv.org/pdf/1003.1628.pdf
        // 6 iterations should be good enough
        amrex::ParticleReal wn = -1.1_rt;
        for (int i=0 ; i < 6 ; i++) {
            amrex::ParticleReal zn = std::log(x/wn) - wn;
            amrex::ParticleReal qn = 2._rt*(1._rt + wn)*(1._rt + wn + 2._rt*zn/3._rt);
            amrex::ParticleReal en = (zn/(1._rt + wn))*((qn - zn)/(qn - 2._rt*zn));
            wn = wn*(1._rt + en);
        }
        return wn;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::ParticleReal invertv3gaussian(amrex::ParticleReal rr) {
        // This returns a random number with the v**3*exp(-v**2) distribution.
        if (rr == 0._rt) {
            return 0._rt;
        } else {
            amrex::ParticleReal x = (rr - 1._rt)/std::exp(1._rt);
            return std::sqrt(-LambertW(x) - 1._rt);
        }
    }
}

class ThermalizationCollision
    : public CollisionBase
{
public:

    ThermalizationCollision (std::string collision_name);

    /** Perform the collisions
     *
     * @param cur_time Current time
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc) override;

private:

    amrex::Real m_thermal_temperature;
    amrex::Real m_mean_free_path;
    bool m_use_thermal_v;
    bool m_use_full_velocity;

};

#endif // WARPX_PARTICLES_COLLISION_THERMALIZATIONCOLLISION_H_
