/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ThermalizationCollision.H"
#include "WarpX.H"

using namespace amrex::literals;

ThermalizationCollision::ThermalizationCollision (std::string const collision_name)
    : CollisionBase(collision_name)
{

    amrex::ParmParse pp(collision_name);

    getWithParser(pp, "thermal_temperature", m_thermal_temperature);
    getWithParser(pp, "mean_free_path", m_mean_free_path);

    m_use_thermal_v = true;
    pp.query("use_thermal_v", m_use_thermal_v);

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

            amrex::Real const dt = WarpX::GetInstance().getdt(lev);

            for (WarpXParIter pti(species, lev); pti.isValid(); ++pti)
            {
                auto& attribs = pti.GetAttribs();
                amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
                long np = pti.numParticles();

                if (m_use_thermal_v) {

                    // Use the thermal velocity to calculate the probability
                    // Velocity is replaced by a normal distribution
                    amrex::ParallelForRNG( np,
                        [=] AMREX_GPU_DEVICE (int ip, amrex::RandomEngine const& engine) noexcept
                        {
                            // Ignore the gamma factor
                            if (thermal_velocity*dt*ndt/mean_free_path > amrex::Random(engine)) {
                                ux[ip] = amrex::RandomNormal(0._rt, thermal_velocity, engine);
                            }
                        }
                    );

                } else {

                    // Use the particle velocity to calculate the probability
                    // Velocity is replaced by a v*exp(-v**2/2*thermal_velocity**2) distribution
                    amrex::ParallelForRNG( np,
                        [=] AMREX_GPU_DEVICE (int ip, amrex::RandomEngine const& engine) noexcept
                        {
                            // Ignore the gamma factor
                            if (std::abs(ux[ip])*dt*ndt/mean_free_path > amrex::Random(engine)) {
                                amrex::ParticleReal const ss = (amrex::Random(engine) < 0.5_rt ? -1._rt : +1._rt);
                                // (1 - rand) is used since rand includes 0, avoiding 1/0.
                                // rand does not include 1 so (1-rand) will never be 0.
                                amrex::ParticleReal const uxrand = 1._rt - amrex::Random(engine);
                                ux[ip] = ss*thermal_velocity*std::sqrt(2._rt*std::log(1._rt/uxrand));
                            }
                        }
                    );

                }
            }
        }
    }
}

