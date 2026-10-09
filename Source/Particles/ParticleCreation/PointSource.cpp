/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Particles/PhysicalParticleContainer.H"

#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/profiler/ProfilerWrapper.H>

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
    [[nodiscard]] amrex::Vector<amrex::ParticleReal>
    build_point_source_cdf (
        long nparticles,
        amrex::ParticleReal width,
        amrex::ParticleReal height,
        amrex::ParticleReal p)
    {
        using namespace amrex::literals;

        constexpr long min_samples = 1024;
        constexpr long max_samples = 1000000;
        const long nsamples = std::max(min_samples, std::min(10L*nparticles, max_samples));

        amrex::Vector<amrex::ParticleReal> cdf(nsamples+1, 0._prt);
        const amrex::ParticleReal dx = width/static_cast<amrex::ParticleReal>(nsamples);

        auto const pdf = [=] (amrex::ParticleReal x)
        {
            const amrex::ParticleReal theta = std::atan2(x, height);
            return std::pow(std::cos(theta), p);
        };

        amrex::ParticleReal prev_f = pdf(-0.5_prt*width);
        for (long i = 1; i <= nsamples; ++i) {
            const amrex::ParticleReal x = -0.5_prt*width + static_cast<amrex::ParticleReal>(i)*dx;
            const amrex::ParticleReal f = pdf(x);
            cdf[i] = cdf[i-1] + 0.5_prt*(prev_f + f)*dx;
            prev_f = f;
        }

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(cdf.back() > 0._prt,
            "Point source CDF normalization must be strictly positive.");

        for (auto& val : cdf) {
            val /= cdf.back();
        }
        cdf.back() = 1._prt;
        return cdf;
    }

    [[nodiscard]] amrex::ParticleReal sample_aperture_position (
        amrex::ParticleReal u,
        amrex::ParticleReal width,
        amrex::Vector<amrex::ParticleReal> const& cdf)
    {
        using namespace amrex::literals;

        auto const it = std::upper_bound(cdf.begin(), cdf.end(), u);
        if (it == cdf.begin()) {
            return -0.5_prt*width;
        }
        if (it == cdf.end()) {
            return 0.5_prt*width;
        }

        const long upper = static_cast<long>(std::distance(cdf.begin(), it));
        const long lower = upper - 1;
        const amrex::ParticleReal x0 = -0.5_prt*width + width*static_cast<amrex::ParticleReal>(lower)
            / static_cast<amrex::ParticleReal>(cdf.size()-1);
        const amrex::ParticleReal x1 = -0.5_prt*width + width*static_cast<amrex::ParticleReal>(upper)
            / static_cast<amrex::ParticleReal>(cdf.size()-1);
        const amrex::ParticleReal c0 = cdf[lower];
        const amrex::ParticleReal c1 = cdf[upper];

        if (c1 <= c0) {
            return x0;
        }

        const amrex::ParticleReal alpha = (u - c0)/(c1 - c0);
        return x0 + alpha*(x1 - x0);
    }
}

void
PhysicalParticleContainer::AddPointSource (PlasmaInjector const& plasma_injector)
{
    using namespace amrex::literals;

    ABLASTR_PROFILE("PhysicalParticleContainer::AddPointSource()");

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_1D_Z) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    amrex::ignore_unused(plasma_injector);
    WARPX_ABORT_WITH_MESSAGE(
        "point_source injection is currently supported only in Cartesian XZ and 3D geometries.");
#else
    const long nparticles = plasma_injector.point_source_nparticles;
    if (nparticles == 0) {
        return;
    }

    const int nprocs = amrex::ParallelDescriptor::NProcs();
    const int myproc = amrex::ParallelDescriptor::MyProc();
    const long base_nparticles = nparticles / nprocs;
    const long remainder = nparticles % nprocs;
    const long local_nparticles = base_nparticles + (myproc < remainder ? 1 : 0);

    if (WarpX::gamma_boost > 1._prt) {
        WARPX_ABORT_WITH_MESSAGE(
            "point_source injection is not yet implemented for boosted-frame simulations.");
    }

    const amrex::ParticleReal width = plasma_injector.point_source_width;
    const amrex::ParticleReal height = plasma_injector.point_source_height;
    const amrex::ParticleReal aperture_z = plasma_injector.point_source_aperture_z;
    const amrex::ParticleReal p = plasma_injector.point_source_p;
    const amrex::ParticleReal vdrift = plasma_injector.point_source_vdrift;
    const amrex::ParticleReal vparallelrms = plasma_injector.point_source_vparallelrms;
    const amrex::ParticleReal vperprms = plasma_injector.point_source_vperprms;
    const amrex::ParticleReal initial_vparallelrms = plasma_injector.point_source_initial_vparallelrms;
    const amrex::ParticleReal initial_vperprms = plasma_injector.point_source_initial_vperprms;
    const amrex::ParticleReal taucycle = plasma_injector.point_source_taucycle;
    const amrex::ParticleReal weight = plasma_injector.point_source_weight;

    amrex::Gpu::HostVector<amrex::ParticleReal> particle_x;
    amrex::Gpu::HostVector<amrex::ParticleReal> particle_y;
    amrex::Gpu::HostVector<amrex::ParticleReal> particle_z;
    amrex::Gpu::HostVector<amrex::ParticleReal> particle_ux;
    amrex::Gpu::HostVector<amrex::ParticleReal> particle_uy;
    amrex::Gpu::HostVector<amrex::ParticleReal> particle_uz;
    amrex::Gpu::HostVector<amrex::ParticleReal> particle_w;

    particle_x.reserve(local_nparticles);
    particle_y.reserve(local_nparticles);
    particle_z.reserve(local_nparticles);
    particle_ux.reserve(local_nparticles);
    particle_uy.reserve(local_nparticles);
    particle_uz.reserve(local_nparticles);
    particle_w.reserve(local_nparticles);

    const amrex::Vector<amrex::ParticleReal> cdf = build_point_source_cdf(nparticles, width, height, p);

    for (long i = 0; i < local_nparticles; ++i) {
        const amrex::ParticleReal t = amrex::Random()*taucycle;
        const amrex::ParticleReal aperture_x = sample_aperture_position(amrex::Random(), width, cdf);
        const amrex::ParticleReal theta = std::atan2(aperture_x, height);
        const amrex::ParticleReal position_vparallel = vdrift + amrex::RandomNormal(0._prt, vparallelrms);
        const amrex::ParticleReal position_vperp = amrex::RandomNormal(0._prt, vperprms);
        const amrex::ParticleReal initial_vparallel = vdrift + amrex::RandomNormal(0._prt, initial_vparallelrms);
        const amrex::ParticleReal initial_vperp = amrex::RandomNormal(0._prt, initial_vperprms);

        const amrex::ParticleReal position_vx = position_vparallel*std::sin(theta) +
            position_vperp*std::cos(theta);
        const amrex::ParticleReal position_vz = position_vparallel*std::cos(theta) -
            position_vperp*std::sin(theta);
        const amrex::ParticleReal vx = initial_vparallel*std::sin(theta) +
            initial_vperp*std::cos(theta);
        const amrex::ParticleReal vz = initial_vparallel*std::cos(theta) -
            initial_vperp*std::sin(theta);
        const amrex::ParticleReal position_v2 = position_vx*position_vx + position_vz*position_vz;
        const amrex::ParticleReal v2 = vx*vx + vz*vz;

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(position_v2 < PhysConst::c2,
            "Point source placement speed must remain below the speed of light.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v2 < PhysConst::c2,
            "Point source particle speed must remain below the speed of light.");

        const amrex::ParticleReal gamma = 1._prt/std::sqrt(1._prt - v2/PhysConst::c2);
        const amrex::ParticleReal x = aperture_x + position_vx*t;
        const amrex::ParticleReal z = aperture_z + position_vz*t;
        const amrex::ParticleReal ux = gamma*vx;
        const amrex::ParticleReal uz = gamma*vz;

#if defined(WARPX_DIM_3D)
        CheckAndAddParticle(x, 0._prt, z, ux, 0._prt, uz, weight,
                            particle_x, particle_y, particle_z,
                            particle_ux, particle_uy, particle_uz,
                            particle_w);
#elif defined(WARPX_DIM_XZ)
        CheckAndAddParticle(x, 0._prt, z, ux, 0._prt, uz, weight,
                            particle_x, particle_y, particle_z,
                            particle_ux, particle_uy, particle_uz,
                            particle_w);
#endif
    }

    auto const np = static_cast<long>(particle_z.size());
    const amrex::Vector<amrex::ParticleReal> xp(particle_x.data(), particle_x.data() + np);
    const amrex::Vector<amrex::ParticleReal> yp(particle_y.data(), particle_y.data() + np);
    const amrex::Vector<amrex::ParticleReal> zp(particle_z.data(), particle_z.data() + np);
    const amrex::Vector<amrex::ParticleReal> uxp(particle_ux.data(), particle_ux.data() + np);
    const amrex::Vector<amrex::ParticleReal> uyp(particle_uy.data(), particle_uy.data() + np);
    const amrex::Vector<amrex::ParticleReal> uzp(particle_uz.data(), particle_uz.data() + np);

    amrex::Vector<amrex::Vector<amrex::ParticleReal>> attr;
    const amrex::Vector<amrex::ParticleReal> wp(particle_w.data(), particle_w.data() + np);
    attr.push_back(wp);

    const amrex::Vector<amrex::Vector<int>> attr_int;

    AddNParticles(0, np, xp, yp, zp, uxp, uyp, uzp,
                  1, attr, 0, attr_int, 1);
#endif
}
