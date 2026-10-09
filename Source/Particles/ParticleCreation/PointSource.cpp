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

using namespace amrex;
using namespace amrex::literals;

namespace
{
    [[nodiscard]] amrex::Vector<amrex::Real>
    build_point_source_cdf (long nparticles, amrex::Real width, amrex::Real height, amrex::Real p)
    {
        constexpr long min_samples = 1024;
        constexpr long max_samples = 1000000;
        const long nsamples = std::max(min_samples, std::min(10L*nparticles, max_samples));

        amrex::Vector<amrex::Real> cdf(nsamples+1, 0._rt);
        const amrex::Real dx = width/static_cast<amrex::Real>(nsamples);

        auto const pdf = [=] (amrex::Real x)
        {
            const amrex::Real theta = std::atan2(x, height);
            return std::pow(std::cos(theta), p);
        };

        amrex::Real prev_f = pdf(-0.5_rt*width);
        for (long i = 1; i <= nsamples; ++i) {
            const amrex::Real x = -0.5_rt*width + static_cast<amrex::Real>(i)*dx;
            const amrex::Real f = pdf(x);
            cdf[i] = cdf[i-1] + 0.5_rt*(prev_f + f)*dx;
            prev_f = f;
        }

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(cdf.back() > 0._rt,
            "Point source CDF normalization must be strictly positive.");

        for (auto& val : cdf) {
            val /= cdf.back();
        }
        cdf.back() = 1._rt;
        return cdf;
    }

    [[nodiscard]] amrex::Real sample_aperture_position (
        amrex::Real u,
        amrex::Real width,
        amrex::Vector<amrex::Real> const& cdf)
    {
        auto const it = std::upper_bound(cdf.begin(), cdf.end(), u);
        if (it == cdf.begin()) {
            return -0.5_rt*width;
        }
        if (it == cdf.end()) {
            return 0.5_rt*width;
        }

        const long upper = static_cast<long>(std::distance(cdf.begin(), it));
        const long lower = upper - 1;
        const amrex::Real x0 = -0.5_rt*width + width*static_cast<amrex::Real>(lower)
            / static_cast<amrex::Real>(cdf.size()-1);
        const amrex::Real x1 = -0.5_rt*width + width*static_cast<amrex::Real>(upper)
            / static_cast<amrex::Real>(cdf.size()-1);
        const amrex::Real c0 = cdf[lower];
        const amrex::Real c1 = cdf[upper];

        if (c1 <= c0) {
            return x0;
        }

        const amrex::Real alpha = (u - c0)/(c1 - c0);
        return x0 + alpha*(x1 - x0);
    }
}

void
PhysicalParticleContainer::AddPointSource (PlasmaInjector const& plasma_injector)
{
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

    if (WarpX::gamma_boost > 1._rt) {
        WARPX_ABORT_WITH_MESSAGE(
            "point_source injection is not yet implemented for boosted-frame simulations.");
    }

    const amrex::Real width = plasma_injector.point_source_width;
    const amrex::Real height = plasma_injector.point_source_height;
    const amrex::Real aperture_z = plasma_injector.point_source_aperture_z;
    const amrex::Real p = plasma_injector.point_source_p;
    const amrex::Real vdrift = plasma_injector.point_source_vdrift;
    const amrex::Real vparallelrms = plasma_injector.point_source_vparallelrms;
    const amrex::Real vperprms = plasma_injector.point_source_vperprms;
    const amrex::Real taucycle = plasma_injector.point_source_taucycle;
    const amrex::ParticleReal weight = plasma_injector.point_source_weight;

    amrex::Gpu::HostVector<ParticleReal> particle_x;
    amrex::Gpu::HostVector<ParticleReal> particle_y;
    amrex::Gpu::HostVector<ParticleReal> particle_z;
    amrex::Gpu::HostVector<ParticleReal> particle_ux;
    amrex::Gpu::HostVector<ParticleReal> particle_uy;
    amrex::Gpu::HostVector<ParticleReal> particle_uz;
    amrex::Gpu::HostVector<ParticleReal> particle_w;

    if (ParallelDescriptor::IOProcessor()) {
        particle_x.reserve(nparticles);
        particle_y.reserve(nparticles);
        particle_z.reserve(nparticles);
        particle_ux.reserve(nparticles);
        particle_uy.reserve(nparticles);
        particle_uz.reserve(nparticles);
        particle_w.reserve(nparticles);

        const amrex::Vector<amrex::Real> cdf = build_point_source_cdf(nparticles, width, height, p);

        for (long i = 0; i < nparticles; ++i) {
            const amrex::Real t = amrex::Random()*taucycle;
            const amrex::Real aperture_x = sample_aperture_position(amrex::Random(), width, cdf);
            const amrex::Real theta = std::atan2(aperture_x, height);
            const amrex::Real vparallel = vdrift + amrex::RandomNormal(0._rt, vparallelrms);
            const amrex::Real vperp = amrex::RandomNormal(0._rt, vperprms);

            const amrex::Real vx = vparallel*std::sin(theta) + vperp*std::cos(theta);
            const amrex::Real vz = vparallel*std::cos(theta) - vperp*std::sin(theta);
            const amrex::Real v2 = vx*vx + vz*vz;

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v2 < PhysConst::c2,
                "Point source particle speed must remain below the speed of light.");

            const amrex::Real gamma = 1._rt/std::sqrt(1._rt - v2/PhysConst::c2);
            const amrex::ParticleReal x = aperture_x + vx*t;
            const amrex::ParticleReal z = aperture_z + vz*t;
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
    }

    auto const np = static_cast<long>(particle_z.size());
    const amrex::Vector<ParticleReal> xp(particle_x.data(), particle_x.data() + np);
    const amrex::Vector<ParticleReal> yp(particle_y.data(), particle_y.data() + np);
    const amrex::Vector<ParticleReal> zp(particle_z.data(), particle_z.data() + np);
    const amrex::Vector<ParticleReal> uxp(particle_ux.data(), particle_ux.data() + np);
    const amrex::Vector<ParticleReal> uyp(particle_uy.data(), particle_uy.data() + np);
    const amrex::Vector<ParticleReal> uzp(particle_uz.data(), particle_uz.data() + np);

    amrex::Vector<amrex::Vector<ParticleReal>> attr;
    const amrex::Vector<ParticleReal> wp(particle_w.data(), particle_w.data() + np);
    attr.push_back(wp);

    const amrex::Vector<amrex::Vector<int>> attr_int;

    AddNParticles(0, np, xp, yp, zp, uxp, uyp, uzp,
                  1, attr, 0, attr_int, 1);
#endif
}
