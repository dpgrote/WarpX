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
#include "Particles/WarpXParticleContainer_fwd.H"

#include <AMReX_REAL.H>

#include <AMReX_BaseFwd.H>

#include <string>

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
