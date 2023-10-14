/* Copyright 2022 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "Diagnostics/MultiDiagnostics.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "Evolve/WarpXDtType.H"
#include "Evolve/WarpXPushType.H"
#ifdef WARPX_USE_PSATD
#   ifdef WARPX_DIM_RZ
#       include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#   else
#       include "FieldSolver/SpectralSolver/SpectralSolver.H"
#   endif
#endif
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Python/WarpX_py.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"

#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <vector>

using namespace amrex;
using ablastr::utils::SignalHandling;

void
WarpX::EvolveImplicitPicardInit (const int lev)
{

    if (lev == 0) {
        // Add space to save the positions and velocities at the start of the time steps
        for (auto const& pc : *mypc) {
#if (AMREX_SPACEDIM >= 2)
            pc->AddRealComp("x_n");
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
            pc->AddRealComp("y_n");
#endif
            pc->AddRealComp("z_n");
            pc->AddRealComp("ux_n");
            pc->AddRealComp("uy_n");
            pc->AddRealComp("uz_n");
        }
    }

    // Initialize MultiFabs to hold the E and B fields at the start of the time steps
    // Only one refinement level is supported
    const int nlevs_max = maxLevel() + 1;
    Efield_n.resize(nlevs_max);
    Bfield_n.resize(nlevs_max);
    Efield_save.resize(nlevs_max);
    Bfield_save.resize(nlevs_max);

    // Strange, the WarpX::DistributionMap(0) is not consistent with Ex_fp.DistributionMap()???

    // Note that the *_fp will be the n+theta and n+1 time level
    AllocInitMultiFabFromModel(Efield_n[lev][0], *Efield_fp[0][0], lev, "Efield_n[0]");
    AllocInitMultiFabFromModel(Efield_n[lev][1], *Efield_fp[0][1], lev, "Efield_n[1]");
    AllocInitMultiFabFromModel(Efield_n[lev][2], *Efield_fp[0][2], lev, "Efield_n[2]");
    AllocInitMultiFabFromModel(Bfield_n[lev][0], *Bfield_fp[0][0], lev, "Bfield_n[0]");
    AllocInitMultiFabFromModel(Bfield_n[lev][1], *Bfield_fp[0][1], lev, "Bfield_n[1]");
    AllocInitMultiFabFromModel(Bfield_n[lev][2], *Bfield_fp[0][2], lev, "Bfield_n[2]");

    AllocInitMultiFabFromModel(Efield_save[lev][0], *Efield_fp[0][0], lev, "Efield_save[0]");
    AllocInitMultiFabFromModel(Efield_save[lev][1], *Efield_fp[0][1], lev, "Efield_save[1]");
    AllocInitMultiFabFromModel(Efield_save[lev][2], *Efield_fp[0][2], lev, "Efield_save[2]");
    AllocInitMultiFabFromModel(Bfield_save[lev][0], *Bfield_fp[0][0], lev, "Bfield_save[0]");
    AllocInitMultiFabFromModel(Bfield_save[lev][1], *Bfield_fp[0][1], lev, "Bfield_save[1]");
    AllocInitMultiFabFromModel(Bfield_save[lev][2], *Bfield_fp[0][2], lev, "Bfield_save[2]");

}

void
WarpX::EvolveImplicitPicard (int numsteps)
{
    WARPX_PROFILE_REGION("WarpX::EvolveImplicitPicard()");
    WARPX_PROFILE("WarpX::EvolveImplicitPicard()");

    Real cur_time = t_new[0];

    int numsteps_max;
    if (numsteps < 0) {  // Note that the default argument is numsteps = -1
        numsteps_max = max_step;
    } else {
        numsteps_max = istep[0] + numsteps;
    }

    bool early_params_checked = false; // check typos in inputs after step 1

    static Real evolve_time = 0;

    const int step_begin = istep[0];
    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {
        WARPX_PROFILE("WarpX::EvolveImplicitPicard::step");
        Real evolve_time_beg_step = amrex::second();

        SignalHandling::CheckSignals();

        multi_diags->NewIteration();

        // Start loop on time steps
        if (verbose) {
            amrex::Print() << "STEP " << step+1 << " starts ...\n";
        }
        ExecutePythonCallback("beforestep");

        amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(0);
        if (cost) {
            if (step > 0 && load_balance_intervals.contains(step+1))
            {
                LoadBalance();

                // Reset the costs to 0
                ResetCosts();
            }
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                cost = WarpX::getCosts(lev);
                if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
                {
                    // Perform running average of the costs
                    // (Giving more importance to most recent costs; only needed
                    // for timers update, heuristic load balance considers the
                    // instantaneous costs)
                    for (int i : cost->IndexArray())
                    {
                        (*cost)[i] *= (1._rt - 2._rt/load_balance_intervals.localPeriod(step+1));
                    }
                }
            }
        }

        // We have B^{n} and E^{n}.
        // Particles have p^{n} and x^{n}.

        // E and B are up-to-date inside the domain only
        /* FillBoundaryE(guard_cells.ng_alloc_EB); */
        /* FillBoundaryB(guard_cells.ng_alloc_EB); */

        // Save the values at the start of the time step,
        // copying particle data to x_n etc.
        for (auto const& pc : *mypc) {
            SaveParticlesAtStepStart (*pc, 0);
        }

        // Save the fields at the start of the step
        amrex::MultiFab::Copy(*Efield_n[0][0], *Efield_fp[0][0], 0, 0, ncomps, Efield_fp[0][0]->nGrowVect());
        amrex::MultiFab::Copy(*Efield_n[0][1], *Efield_fp[0][1], 0, 0, ncomps, Efield_fp[0][1]->nGrowVect());
        amrex::MultiFab::Copy(*Efield_n[0][2], *Efield_fp[0][2], 0, 0, ncomps, Efield_fp[0][2]->nGrowVect());
        amrex::MultiFab::Copy(*Bfield_n[0][0], *Bfield_fp[0][0], 0, 0, ncomps, Bfield_fp[0][0]->nGrowVect());
        amrex::MultiFab::Copy(*Bfield_n[0][1], *Bfield_fp[0][1], 0, 0, ncomps, Bfield_fp[0][1]->nGrowVect());
        amrex::MultiFab::Copy(*Bfield_n[0][2], *Bfield_fp[0][2], 0, 0, ncomps, Bfield_fp[0][2]->nGrowVect());

        // Start the iterations
        amrex::Real deltaE = 1._rt;
        amrex::Real deltaB = 1._rt;
        int iteration_count = 0;
        while (iteration_count < max_picard_iterations &&
               (deltaE > picard_iteration_tolerance || deltaB > picard_iteration_tolerance)) {
            iteration_count++;

            // Advance the particle positions by 1/2 dt,
            // particle velocities by dt, then take average of old and new v,
            // deposit currents, giving J at n+1/2
            // This uses Efield_fp and Bfield_fp, the field at n+1/2 from the previous iteration.
            bool skip_deposition = false;
            PushType push_type = PushType::Implicit;
            PushParticlesandDepose(cur_time, skip_deposition, push_type);

            SyncCurrentAndRho();

            if (picard_iteration_tolerance > 0. || iteration_count == max_picard_iterations) {
                // Save the E at n+1/2 from the previous iteration so that the change
                // in this iteration can be calculated
                amrex::MultiFab::Copy(*Efield_save[0][0], *Efield_fp[0][0], 0, 0, ncomps, 0);
                amrex::MultiFab::Copy(*Efield_save[0][1], *Efield_fp[0][1], 0, 0, ncomps, 0);
                amrex::MultiFab::Copy(*Efield_save[0][2], *Efield_fp[0][2], 0, 0, ncomps, 0);
            }

            // Copy Efield_n into Efield_fp since EvolveE updates Efield_fp in place
            amrex::MultiFab::Copy(*Efield_fp[0][0], *Efield_n[0][0], 0, 0, ncomps, Efield_n[0][0]->nGrowVect());
            amrex::MultiFab::Copy(*Efield_fp[0][1], *Efield_n[0][1], 0, 0, ncomps, Efield_n[0][1]->nGrowVect());
            amrex::MultiFab::Copy(*Efield_fp[0][2], *Efield_n[0][2], 0, 0, ncomps, Efield_n[0][2]->nGrowVect());

            // Updates Efield_fp so it holds the new E at n+1/2
            EvolveE(0.5_rt*dt[0]);
            // WarpX::sync_nodal_points is used to avoid instability
            FillBoundaryE(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
            ApplyEfieldBoundary(0, PatchType::fine);

            if (picard_iteration_tolerance > 0. || iteration_count == max_picard_iterations) {
                // Save the B at n+1/2 from the previous iteration so that the change
                // in this iteration can be calculated
                amrex::MultiFab::Copy(*Bfield_save[0][0], *Bfield_fp[0][0], 0, 0, ncomps, 0);
                amrex::MultiFab::Copy(*Bfield_save[0][1], *Bfield_fp[0][1], 0, 0, ncomps, 0);
                amrex::MultiFab::Copy(*Bfield_save[0][2], *Bfield_fp[0][2], 0, 0, ncomps, 0);
            }

            // Copy Bfield_n into Bfield_fp since EvolveB updates Bfield_fp in place
            amrex::MultiFab::Copy(*Bfield_fp[0][0], *Bfield_n[0][0], 0, 0, ncomps, Bfield_n[0][0]->nGrowVect());
            amrex::MultiFab::Copy(*Bfield_fp[0][1], *Bfield_n[0][1], 0, 0, ncomps, Bfield_n[0][1]->nGrowVect());
            amrex::MultiFab::Copy(*Bfield_fp[0][2], *Bfield_n[0][2], 0, 0, ncomps, Bfield_n[0][2]->nGrowVect());

            // This updates Bfield_fp so it holds the new B at n+1/2
            EvolveB(0.5_rt*dt[0], DtType::Full);
            // WarpX::sync_nodal_points is used to avoid instability
            FillBoundaryB(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
            ApplyBfieldBoundary(0, PatchType::fine, DtType::Full);

            // The B field update needs
            if (num_mirrors>0){
                applyMirrors(cur_time);
                // E : guard cells are NOT up-to-date from the mirrors
                // B : guard cells are NOT up-to-date from the mirrors
            }

            if (picard_iteration_tolerance > 0. || iteration_count == max_picard_iterations) {
                // Calculate the change in E and B from this iteration
                // deltaE = abs(Enew - Eold)/max(abs(Enew))
                Efield_save[0][0]->minus(*Efield_fp[0][0], 0, ncomps, 0);
                Efield_save[0][1]->minus(*Efield_fp[0][1], 0, ncomps, 0);
                Efield_save[0][2]->minus(*Efield_fp[0][2], 0, ncomps, 0);
                amrex::Real maxE0 = std::max(1._rt, Efield_fp[0][0]->norm0(0, 0));
                amrex::Real maxE1 = std::max(1._rt, Efield_fp[0][1]->norm0(0, 0));
                amrex::Real maxE2 = std::max(1._rt, Efield_fp[0][2]->norm0(0, 0));
                amrex::Real deltaE0 = Efield_save[0][0]->norm0(0, 0)/maxE0;
                amrex::Real deltaE1 = Efield_save[0][1]->norm0(0, 0)/maxE1;
                amrex::Real deltaE2 = Efield_save[0][2]->norm0(0, 0)/maxE2;
                /* amrex::Print() << "deltaE " << iteration_count << " " << deltaE0 << " " << deltaE1 << " " << deltaE2 << "\n"; */
                deltaE = std::max(std::max(deltaE0, deltaE1), deltaE2);
                Bfield_save[0][0]->minus(*Bfield_fp[0][0], 0, ncomps, 0);
                Bfield_save[0][1]->minus(*Bfield_fp[0][1], 0, ncomps, 0);
                Bfield_save[0][2]->minus(*Bfield_fp[0][2], 0, ncomps, 0);
                amrex::Real maxB0 = std::max(1._rt, Bfield_fp[0][0]->norm0(0, 0));
                amrex::Real maxB1 = std::max(1._rt, Bfield_fp[0][1]->norm0(0, 0));
                amrex::Real maxB2 = std::max(1._rt, Bfield_fp[0][2]->norm0(0, 0));
                amrex::Real deltaB0 = Bfield_save[0][0]->norm0(0, 0)/maxB0;
                amrex::Real deltaB1 = Bfield_save[0][1]->norm0(0, 0)/maxB1;
                amrex::Real deltaB2 = Bfield_save[0][2]->norm0(0, 0)/maxB2;
                /* amrex::Print() << "deltaB " << iteration_count << " " << deltaB0 << " " << deltaB1 << " " << deltaB2 << "\n"; */
                deltaB = std::max(std::max(deltaB0, deltaB1), deltaB2);
                amrex::Print() << "Max delta " << iteration_count << " " << deltaE << " " << deltaB << "\n";
            }

            // Now, the particle positions and velocities and the Efield_fp and Bfield_fp hold
            // the new values at n+1/2
        }

        amrex::Print() << "Picard iterations = " << iteration_count << ", Eerror =  " << deltaE << ", Berror =  " << deltaB << "\n";

        // Advance particles to step n+1
        for (auto const& pc : *mypc) {
            FinishImplicitParticleUpdate(*pc, 0);
        }

        // Advance fields to step n+1
        FinishImplicitFieldUpdate(Efield_fp, Efield_n);
        FinishImplicitFieldUpdate(Bfield_fp, Bfield_n);
        // WarpX::sync_nodal_points is used to avoid instability
        FillBoundaryE(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);
        FillBoundaryB(guard_cells.ng_alloc_EB, WarpX::sync_nodal_points);

        // Run multi-physics modules:
        // ionization, Coulomb collisions, QED
        doFieldIonization();
        ExecutePythonCallback("beforecollisions");
        mypc->doCollisions( cur_time, dt[0] );
        ExecutePythonCallback("aftercollisions");
#ifdef WARPX_QED
        doQEDEvents();
        mypc->doQEDSchwinger();
#endif

        for (int lev = 0; lev <= max_level; ++lev) {
            ++istep[lev];
        }

        cur_time += dt[0];

        ShiftGalileanBoundary();

        // sync up time
        for (int i = 0; i <= max_level; ++i) {
            t_new[i] = cur_time;
        }
        multi_diags->FilterComputePackFlush( step, false, true );

        bool move_j = is_synchronized;
        // If is_synchronized we need to shift j too so that next step we can evolve E by dt/2.
        // We might need to move j because we are going to make a plotfile.
        MoveWindow(step+1, move_j);

        mypc->ContinuousFluxInjection(cur_time, dt[0]);

        mypc->ApplyBoundaryConditions();

        // interact the particles with EB walls (if present)
#ifdef AMREX_USE_EB
        mypc->ScrapeParticles(amrex::GetVecOfConstPtrs(m_distance_to_eb));
#endif

        m_particle_boundary_buffer->gatherParticles(*mypc, amrex::GetVecOfConstPtrs(m_distance_to_eb));

        // Implicit solver: particles can move by an arbitrary number of cells
        mypc->Redistribute();

        if (sort_intervals.contains(step+1)) {
            if (verbose) {
                amrex::Print() << Utils::TextMsg::Info("re-sorting particles");
            }
            mypc->SortParticlesByBin(sort_bin_size);
        }

        // afterstep callback runs with the updated global time. It is included
        // in the evolve timing.
        ExecutePythonCallback("afterstep");

        /// reduced diags
        if (reduced_diags->m_plot_rd != 0)
        {
            reduced_diags->LoadBalance();
            reduced_diags->ComputeDiags(step);
            reduced_diags->WriteToFile(step);
        }
        multi_diags->FilterComputePackFlush( step );

        // execute afterdiagnostic callbacks
        ExecutePythonCallback("afterdiagnostics");

        // inputs: unused parameters (e.g. typos) check after step 1 has finished
        if (!early_params_checked) {
            amrex::Print() << "\n"; // better: conditional \n based on return value
            amrex::ParmParse().QueryUnusedInputs();

            //Print the warning list right after the first step.
            amrex::Print() <<
                ablastr::warn_manager::GetWMInstance().PrintGlobalWarnings("FIRST STEP");
            early_params_checked = true;
        }

        // create ending time stamp for calculating elapsed time each iteration
        Real evolve_time_end_step = amrex::second();
        evolve_time += evolve_time_end_step - evolve_time_beg_step;

        HandleSignals();

        if (verbose) {
            amrex::Print()<< "STEP " << step+1 << " ends." << " TIME = " << cur_time
                        << " DT = " << dt[0] << "\n";
            amrex::Print()<< "Evolve time = " << evolve_time
                      << " s; This step = " << evolve_time_end_step-evolve_time_beg_step
                      << " s; Avg. per step = " << evolve_time/(step-step_begin+1) << " s\n\n";
        }

        if (cur_time >= stop_time - 1.e-3*dt[0] || SignalHandling::TestAndResetActionRequestFlag(SignalHandling::SIGNAL_REQUESTS_BREAK)) {
            break;
        }

        // End loop on time steps
    }
    // This if statement is needed for PICMI, which allows the Evolve routine to be
    // called multiple times, otherwise diagnostics will be done at every call,
    // regardless of the diagnostic period parameter provided in the inputs.
    if (istep[0] == max_step || (stop_time - 1.e-3*dt[0] <= cur_time && cur_time < stop_time + dt[0])) {
        multi_diags->FilterComputePackFlushLastTimestep( istep[0] );
    }
}

void
WarpX::SaveParticlesAtStepStart (WarpXParticleContainer& pc, const int lev)
{

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {

    auto particle_comps = pc.getParticleComps();

    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {

        const auto getPosition = GetParticlePosition(pti);

        auto& attribs = pti.GetAttribs();
        ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
        ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
        ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

#if (AMREX_SPACEDIM >= 2)
        ParticleReal* x_n = pti.GetAttribs(particle_comps["x_n"]).dataPtr();
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
        ParticleReal* y_n = pti.GetAttribs(particle_comps["y_n"]).dataPtr();
#endif
        ParticleReal* z_n = pti.GetAttribs(particle_comps["z_n"]).dataPtr();
        ParticleReal* ux_n = pti.GetAttribs(particle_comps["ux_n"]).dataPtr();
        ParticleReal* uy_n = pti.GetAttribs(particle_comps["uy_n"]).dataPtr();
        ParticleReal* uz_n = pti.GetAttribs(particle_comps["uz_n"]).dataPtr();

        const long np = pti.numParticles();

        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            amrex::ParticleReal xp, yp, zp;
            getPosition(ip, xp, yp, zp);

#if (AMREX_SPACEDIM >= 2)
            x_n[ip] = xp;
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
            y_n[ip] = yp;
#endif
            z_n[ip] = zp;

            ux_n[ip] = ux[ip];
            uy_n[ip] = uy[ip];
            uz_n[ip] = uz[ip];

        });

    }
    }
}

void
WarpX::FinishImplicitParticleUpdate (WarpXParticleContainer& pc, const int lev)
{

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {

    auto particle_comps = pc.getParticleComps();

    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {

        const auto getPosition = GetParticlePosition(pti);
        const auto setPosition = SetParticlePosition(pti);

        auto& attribs = pti.GetAttribs();
        ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
        ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
        ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

#if (AMREX_SPACEDIM >= 2)
        ParticleReal* x_n = pti.GetAttribs(particle_comps["x_n"]).dataPtr();
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
        ParticleReal* y_n = pti.GetAttribs(particle_comps["y_n"]).dataPtr();
#endif
        ParticleReal* z_n = pti.GetAttribs(particle_comps["z_n"]).dataPtr();
        ParticleReal* ux_n = pti.GetAttribs(particle_comps["ux_n"]).dataPtr();
        ParticleReal* uy_n = pti.GetAttribs(particle_comps["uy_n"]).dataPtr();
        ParticleReal* uz_n = pti.GetAttribs(particle_comps["uz_n"]).dataPtr();

        const long np = pti.numParticles();

        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            amrex::ParticleReal xp, yp, zp;
            getPosition(ip, xp, yp, zp);

#if (AMREX_SPACEDIM >= 2)
            xp = 2._rt*xp - x_n[ip];
#endif
#if defined(WARPX_DIM_3D) || defined(WARPX_DIM_RZ)
            yp = 2._rt*yp - y_n[ip];
#endif
            zp = 2._rt*zp - z_n[ip];

            ux[ip] = 2._rt*ux[ip] - ux_n[ip];
            uy[ip] = 2._rt*uy[ip] - uy_n[ip];
            uz[ip] = 2._rt*uz[ip] - uz_n[ip];

            setPosition(ip, xp, yp, zp);
        });

    }
    }
}

void
WarpX::FinishImplicitFieldUpdate(amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& Field_fp,
                                 amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& Field_n)
{

    for (int lev = 0; lev <= finest_level; ++lev) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
       for ( amrex::MFIter mfi(*Field_fp[lev][0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {

            amrex::Array4<amrex::Real> const& Fx = Field_fp[lev][0]->array(mfi);
            amrex::Array4<amrex::Real> const& Fy = Field_fp[lev][1]->array(mfi);
            amrex::Array4<amrex::Real> const& Fz = Field_fp[lev][2]->array(mfi);

            amrex::Array4<amrex::Real> const& Fx_n = Field_n[lev][0]->array(mfi);
            amrex::Array4<amrex::Real> const& Fy_n = Field_n[lev][1]->array(mfi);
            amrex::Array4<amrex::Real> const& Fz_n = Field_n[lev][2]->array(mfi);

            amrex::Box tbx = mfi.tilebox(Field_fp[lev][0]->ixType().toIntVect());
            amrex::Box tby = mfi.tilebox(Field_fp[lev][1]->ixType().toIntVect());
            amrex::Box tbz = mfi.tilebox(Field_fp[lev][2]->ixType().toIntVect());

            /* tbx.grow(Field_fp[lev][0]->nGrowVect()); */
            /* tby.grow(Field_fp[lev][1]->nGrowVect()); */
            /* tbz.grow(Field_fp[lev][2]->nGrowVect()); */

            amrex::ParallelFor(
            tbx, ncomps, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Fx(i,j,k,n) = 2._rt*Fx(i,j,k,n) - Fx_n(i,j,k,n);
            },
            tby, ncomps, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Fy(i,j,k,n) = 2._rt*Fy(i,j,k,n) - Fy_n(i,j,k,n);
            },
            tbz, ncomps, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Fz(i,j,k,n) = 2._rt*Fz(i,j,k,n) - Fz_n(i,j,k,n);
            });
        }
    }
}
