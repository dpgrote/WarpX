/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ThermalizationCollision.H"
#include "Particles/Collision/CollisionBase.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_Extension.H>
/* #include <AMReX_GpuLaunchFunctsC.H> */
#include <AMReX_GpuQualifiers.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_Vector.H>

#include <stdlib.h>
#include <array>
#include <cmath>
#include <string>

ThermalizationCollision::ThermalizationCollision (std::string const collision_name)
    : CollisionBase(collision_name)
{

    amrex::ParmParse pp(collision_name);

    getWithParser(pp, "thermal_temperature", m_thermal_temperature);
    getWithParser(pp, "mean_free_path", m_mean_free_path);

    m_use_thermal_v = true;
    m_use_full_velocity = true;
    pp.query("use_thermal_v", m_use_thermal_v);
    pp.query("use_full_velocity", m_use_full_velocity);

}

void
ThermalizationCollision::doCollisions (amrex::Real cur_time, MultiParticleContainer* mypc)
{
    const amrex::Real dt = WarpX::GetInstance().getdt(0);
    if ( int(std::floor(cur_time/dt)) % m_ndt != 0 ) return;

    if (m_mean_free_path == 0._rt) {
        return;
    }

    for (auto& species_name : m_species_names) {

        auto& species = mypc->GetParticleContainerFromName(species_name);

        amrex::Real const mean_free_path = m_mean_free_path;
        amrex::Real const thermal_velocity = std::sqrt(m_thermal_temperature*PhysConst::q_e/species.getMass());
        int const ndt = m_ndt;

        // Loop over refinement levels
        for (int lev = 0; lev <= species.finestLevel(); ++lev){

            amrex::Real const dt_lev = WarpX::GetInstance().getdt(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (WarpXParIter pti(species, lev); pti.isValid(); ++pti)
            {
                auto& attribs = pti.GetAttribs();
                amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
                amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
                long np = pti.numParticles();

                if (m_use_thermal_v) {

                    // Use the thermal velocity to calculate the probability
                    // Velocity is replaced by a normal distribution
                    const bool use_full_velocity = m_use_full_velocity;
                    amrex::ParallelForRNG( np,
                        [=] AMREX_GPU_DEVICE (int ip, amrex::RandomEngine const& engine) noexcept
                        {
                            // Ignore the gamma factor
                            if (thermal_velocity*dt_lev*ndt/mean_free_path > amrex::Random(engine)) {
                                ux[ip] = amrex::RandomNormal(0._rt, thermal_velocity, engine);
                                if (use_full_velocity) {
                                    uy[ip] = amrex::RandomNormal(0._rt, thermal_velocity, engine);
                                    uz[ip] = amrex::RandomNormal(0._rt, thermal_velocity, engine);
                                }
                            }
                        }
                    );

                } else {

                    // Use the particle velocity to calculate the probability
                    // Velocity is replaced by a v*exp(-v**2/2*thermal_velocity**2) distribution
                    if (m_use_full_velocity) {

                        amrex::ParallelForRNG( np,
                            [=] AMREX_GPU_DEVICE (int ip, amrex::RandomEngine const& engine) noexcept
                            {
                                // Ignore the gamma factor since gamma ~ 1
                                amrex::ParticleReal v = std::sqrt(ux[ip]*ux[ip] + uy[ip]*uy[ip] + uz[ip]*uz[ip]);
                                if (v*dt_lev*ndt/mean_free_path > amrex::Random(engine)) {
                                    amrex::ParticleReal const rr = amrex::Random(engine);
                                    amrex::ParticleReal const vv = invertv3gaussian(rr);
                                    amrex::ParticleReal const vnew = std::sqrt(2._rt)*thermal_velocity*vv;
                                    amrex::ParticleReal const phi = 2._rt*MathConst::pi*amrex::Random(engine);
                                    amrex::ParticleReal const theta = std::acos(1._rt - 2._rt*amrex::Random(engine));
                                    ux[ip] = vnew*std::cos(phi)*std::sin(theta);
                                    uy[ip] = vnew*std::sin(phi)*std::sin(theta);
                                    uz[ip] = vnew*std::cos(theta);
                                }
                            }
                        );

                    } else {

                        amrex::ParallelForRNG( np,
                            [=] AMREX_GPU_DEVICE (int ip, amrex::RandomEngine const& engine) noexcept
                            {
                                // Ignore the gamma factor since gamma ~ 1
                                amrex::ParticleReal v = std::abs(ux[ip]);
                                if (v*dt_lev*ndt/mean_free_path > amrex::Random(engine)) {
                                    // (1 - rand) is used since rand includes 0, avoiding 1/0.
                                    // rand does not include 1 so (1-rand) will never be 0.
                                    amrex::ParticleReal const uxrand = 1._rt - amrex::Random(engine);
                                    amrex::ParticleReal const vnew = thermal_velocity*std::sqrt(2._rt*std::log(1._rt/uxrand));
                                    amrex::ParticleReal const ss = (amrex::Random(engine) < 0.5_rt ? -1._rt : +1._rt);
                                    ux[ip] = ss*vnew;
                                }
                            }
                        );
                    }

                }
            }
        }
    }
}

