/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "StrangImplicitSpectralEM.H"
#include "WarpX.H"

void StrangImplicitSpectralEM::Define ( WarpX* const a_WarpX )
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_defined,
        "StrangImplicitSpectralEM object is already defined!");

    // Retain a pointer back to main WarpX class
    m_WarpX = a_WarpX;

    // Define E vectors
    m_E.Define( m_WarpX->getEfield_fp_vec() );
    m_Eold.Define( m_WarpX->getEfield_fp_vec() );

    // Need to define the WarpXSolverVec owned dot_mask to do dot
    // product correctly for linear and nonlinear solvers
    amrex::Vector<amrex::Geometry> const & Geom = m_WarpX->Geom();
    m_E.SetDotMask(Geom);

    // Parse implicit solver parameters
    amrex::ParmParse const pp("implicit_evolve");

    std::string nlsolver_type_str;
    pp.query("nonlinear_solver", nlsolver_type_str);
    if (nlsolver_type_str=="picard") {
        m_nlsolver_type = NonlinearSolverType::Picard;
        m_max_particle_iterations = 1;
        m_particle_tolerance = 0.0;
    }
    else if (nlsolver_type_str=="newton") {
        m_nlsolver_type = NonlinearSolverType::Newton;
        pp.query("max_particle_iterations", m_max_particle_iterations);
        pp.query("particle_tolerance", m_particle_tolerance);
    }
    else {
        WARPX_ABORT_WITH_MESSAGE(
            "invalid nonlinear_solver specified. Valid options are picard and newton.");
    }

    // Define the nonlinear solver
    if (m_nlsolver_type == NonlinearSolverType::Picard) {
        m_nlsolver = std::make_unique<PicardSolver<WarpXSolverVec,StrangImplicitSpectralEM>>();
        m_nlsolver->Define(m_E, this);
    }
    else if (m_nlsolver_type == NonlinearSolverType::Newton) {
        m_nlsolver = std::make_unique<NewtonSolver<WarpXSolverVec,StrangImplicitSpectralEM>>();
        m_nlsolver->Define(m_E, this);
    }

    m_is_defined = true;
}

void StrangImplicitSpectralEM::PrintParameters () const
{
    if (!m_WarpX->Verbose()) { return; }
    amrex::Print() << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "----------- IMPLICIT SPECTRAL EM SOLVER PARAMETERS --------" << std::endl;
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << "max particle iterations:    " << m_max_particle_iterations << std::endl;
    amrex::Print() << "particle tolerance:         " << m_particle_tolerance << std::endl;
    if (m_nlsolver_type==NonlinearSolverType::Picard) {
        amrex::Print() << "Nonlinear solver type:      Picard" << std::endl;
    }
    else if (m_nlsolver_type==NonlinearSolverType::Newton) {
        amrex::Print() << "Nonlinear solver type:      Newton" << std::endl;
    }
    m_nlsolver->PrintParams();
    amrex::Print() << "-----------------------------------------------------------" << std::endl;
    amrex::Print() << std::endl;
}

void StrangImplicitSpectralEM::OneStep ( amrex::Real a_time,
                                         amrex::Real a_dt,
                                         int a_step )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_step);

    // Fields have E^{n} and B^{n}
    // Particles have p^{n} and x^{n}.

    // Save the values at the start of the time step,
    m_WarpX->SaveParticlesAtImplicitStepStart();

    // Advance the fields to time n+1/2 source free
    m_WarpX->SpectralSourceFreeFieldAdvance();

    // Save the fields at the start of the step
    m_Eold.Copy( m_WarpX->getEfield_fp_vec() );
    m_E = m_Eold; // initial guess for E

    amrex::Real const half_time = a_time + 0.5_rt*a_dt;

    // Solve nonlinear system for E at t_{n+1/2}
    // Particles will be advanced to t_{n+1/2}
    m_nlsolver->Solve( m_E, m_Eold, half_time, a_dt );

    // Update WarpX owned Efield_fp and Bfield_fp to t_{n+1/2}
    UpdateWarpXState( m_E, half_time );

    // Update field boundary probes prior to updating fields to t_{n+1}
    //UpdateBoundaryProbes( a_dt )

    // Advance particles from time n+1/2 to time n+1
    m_WarpX->FinishImplicitParticleUpdate();

    // Advance E and B fields from time n+1/2 to time n+1
    amrex::Real const new_time = a_time + a_dt;
    FinishFieldUpdate( new_time );

    // Advance the fields to time n+1 source free
    m_WarpX->SpectralSourceFreeFieldAdvance();

}

void StrangImplicitSpectralEM::PreRHSOp ( WarpXSolverVec const & a_E,
                                          amrex::Real a_time,
                                          amrex::Real a_dt,
                                          int a_nl_iter,
                                          bool a_from_jacobian )
{
    amrex::ignore_unused(a_E);

    // update derived variable B and then update WarpX owned Efield_fp and Bfield_fp
    UpdateWarpXState( a_E, a_time );

    // Advance the particle positions by 1/2 dt,
    // particle velocities by dt, then take average of old and new v,
    // deposit currents, giving J at n+1/2 used in ComputeRHSE below
    m_WarpX->PreRHSOp( a_time, a_dt, a_nl_iter, a_from_jacobian );

}

void StrangImplicitSpectralEM::ComputeRHS ( WarpXSolverVec& a_Erhs,
                                            WarpXSolverVec const & a_E,
                                            amrex::Real a_time,
                                            amrex::Real a_dt )
{

    using namespace amrex::literals;
    amrex::ignore_unused(a_E, a_time);

    /* m_WarpX->ComputeRHSE_OnlyJ(0.5_rt*a_dt, a_Erhs); */

    amrex::Real constexpr coeff = PhysConst::c * PhysConst::c * PhysConst::mu0;

    // For Strang split implicit PSATD, the RHS = -dt*mu*c**2*J
    a_Erhs.Copy(m_WarpX->getcurrent_fp_vec());
    a_Erhs.scale(-coeff * 0.5_rt*a_dt);

}

void StrangImplicitSpectralEM::UpdateWarpXState (WarpXSolverVec const & a_E,
                                                 amrex::Real a_time)
{

    // Update Efield_fp owned by WarpX
    m_WarpX->SetElectricFieldAndApplyBCs( a_E );

    if (WarpX::num_mirrors > 0){
        m_WarpX->applyMirrors(a_time);
        // E : guard cells are NOT up-to-date from the mirrors
        // B : guard cells are NOT up-to-date from the mirrors
    }

}

void StrangImplicitSpectralEM::FinishFieldUpdate ( amrex::Real a_new_time )
{
    using namespace amrex::literals;
    amrex::ignore_unused(a_new_time);

    // Eg^{n+1} = 2*E_g^{n+1/2} - E_g^n

    amrex::Real const c0 = 1._rt/0.5_rt;
    amrex::Real const c1 = 1._rt - c0;
    m_E.linComb( c0, m_E, c1, m_Eold );
    m_WarpX->SetElectricFieldAndApplyBCs( m_E );

    if (WarpX::num_mirrors>0){
        m_WarpX->applyMirrors(a_new_time);
        // E : guard cells are NOT up-to-date from the mirrors
        // B : guard cells are NOT up-to-date from the mirrors
    }

}
