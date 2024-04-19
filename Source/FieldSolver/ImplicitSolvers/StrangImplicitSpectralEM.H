/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef STRANG_IMPLICIT_SPECTRALEM_H_
#define STRANG_IMPLICIT_SPECTRALEM_H_

#include "FieldSolver/ImplicitSolvers/WarpXSolverVec.H"
#include "NonlinearSolvers/NonlinearSolverLibrary.H"

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "ImplicitSolver.H"

/** @file
 *  Implicit spectral electromagnetic time solver class. This is a fully implicit
 *  algorithm where both the fields and particles are treated implicitly.
 *
 *  The time stencil is
 *    Advance (Eg^n, Bg^n) -> (Eg^{n+1/2}, Bg^{n+1/2}) source free // E transverse
 *    Iterate:
 *      Eg^{n+1} = Eg^n + c^2*dt*( - mu0*Jg^{n+1/2} ) // E longitudinal
 *      xp^{n+1} = xp^n + dt*up^{n+1/2}/(0.5*(gammap^n + gammap^{n+1}))
 *      up^{n+1} = up^n + dt*qp/mp*(Ep^{n+1/2} + up^{n+1/2}/gammap^{n+1/2} x Bp^{n+1/2})
 *    Advance (Eg^n+1/2, Bg^n+1/2) -> (Eg^{n+1}, Bg^{n+1}) source free // E transverse
 *
 *  The algorithm is exactly energy conserving for.
 *  The algorithm is numerially stable for any time step.
 *  I.e., the CFL condition for light waves does not
 *  have to be satisifed and the time step is not limited by the plasma period. However, how
 *  efficiently the algorithm can use large time steps depends strongly on the nonlinear solver.
 *  Furthermore, the time step should always be such that particles do not travel outside the
 *  ghost region of the box they live in, which is an MPI-related limitation. The time step
 *  is always limited by the need to resolve the appropriate physics.
 *
 */

class StrangImplicitSpectralEM : public ImplicitSolver
{
public:

    StrangImplicitSpectralEM() = default;

    ~StrangImplicitSpectralEM() override {
        // If the function below is not called, then the following
        // error message occurs after the simulation finishes:
        // malloc_consolidate(): unaligned fastbin chunk detected
        WarpXSolverVec::clearDotMask();
    }

    // Prohibit Move and Copy operations
    StrangImplicitSpectralEM(const StrangImplicitSpectralEM&) = delete;
    StrangImplicitSpectralEM& operator=(const StrangImplicitSpectralEM&) = delete;
    StrangImplicitSpectralEM(StrangImplicitSpectralEM&&) = delete;
    StrangImplicitSpectralEM& operator=(StrangImplicitSpectralEM&&) = delete;

    void Define ( WarpX*  a_WarpX ) override;

    [[nodiscard]] bool IsDefined () const override { return m_is_defined; }

    void PrintParameters () const override;

    void Initialize () override {;}

    void GetParticleSolverParams ( int&  a_max_particle_iterations,
                                   amrex::ParticleReal&  a_particle_tolerance ) override
    {
        a_max_particle_iterations = m_max_particle_iterations;
        a_particle_tolerance = m_particle_tolerance;
    }

    void OneStep ( amrex::Real  a_time,
                   amrex::Real  a_dt,
                   int          a_step ) override;

    void PreRHSOp ( const WarpXSolverVec&  a_E,
                          amrex::Real      a_time,
                          amrex::Real      a_dt,
                          int              a_nl_iter,
                          bool             a_from_jacobian ) override;

    void ComputeRHS ( WarpXSolverVec&  a_Erhs,
                const WarpXSolverVec&  a_E,
                      amrex::Real      a_time,
                      amrex::Real      a_dt ) override;

private:

    bool m_is_defined = false;

    WarpX* m_WarpX;

    /**
     * \brief tolerance used by the Picard method used to obtain a self-consistent
     *  update of the particle positions and velocities for given E and B on the grid
     */
    amrex::ParticleReal m_particle_tolerance = 1.0e-10;

    /**
     * \brief maximum iterations for the Picard method used to obtain a self-consistent
     *  update of the particle positions and velocities for given E and B on the grid
     */
    int m_max_particle_iterations = 20;

    NonlinearSolverType m_nlsolver_type = NonlinearSolverType::Picard;
    std::unique_ptr<NonlinearSolver<WarpXSolverVec,StrangImplicitSpectralEM>> m_nlsolver;

    WarpXSolverVec m_E, m_Eold;

    void UpdateWarpXState (const WarpXSolverVec& a_E,
                           amrex::Real a_time);

    void FinishFieldUpdate (amrex::Real a_new_time);

};

#endif
