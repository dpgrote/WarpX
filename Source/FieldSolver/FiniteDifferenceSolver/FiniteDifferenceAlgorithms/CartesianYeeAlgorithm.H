/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_YEE_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_YEE_H_

#include "FieldAccessorFunctors.H"
#include "Utils/WarpXConst.H"
#include "Utils/Parser/ParserUtils.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <array>
#include <cmath>


/**
 * This struct contains only static functions to initialize the stencil coefficients
 * and to compute finite-difference derivatives for the Cartesian Yee algorithm.
 */
struct CartesianYeeAlgorithm {

    static void InitializeStencilCoefficients (
        std::array<amrex::Real,3>& cell_size,
        amrex::Vector<amrex::Real>& stencil_coefs_x,
        amrex::Vector<amrex::Real>& stencil_coefs_y,
        amrex::Vector<amrex::Real>& stencil_coefs_z ) {

        using namespace amrex;
        // Store the inverse cell size along each direction in the coefficients
        stencil_coefs_x.resize(1);
        stencil_coefs_x[0] = 1._rt/cell_size[0];
        stencil_coefs_y.resize(1);
        stencil_coefs_y[0] = 1._rt/cell_size[1];
        stencil_coefs_z.resize(1);
        stencil_coefs_z[0] = 1._rt/cell_size[2];
    }

    /**
     * Compute the maximum timestep, for which the scheme remains stable
     * (Courant-Friedrichs-Levy limit) */
    static amrex::Real ComputeMaxDt ( amrex::Real const * const dx ) {
        using namespace amrex::literals;
        amrex::Real const delta_t  = 1._rt / ( std::sqrt( AMREX_D_TERM(
                                           1._rt / (dx[0]*dx[0]),
                                         + 1._rt / (dx[1]*dx[1]),
                                         + 1._rt / (dx[2]*dx[2])
                                     ) ) * PhysConst::c );
        return delta_t;
    }

    /**
     * \brief Returns maximum number of guard cells required by the field-solve
     */
    static amrex::IntVect GetMaxGuardCell () {
        // The yee solver requires one guard cell in each dimension
        int ngx = 1, ngy = 1, ngz = 1;
        const amrex::ParmParse pp_grid("grid");
        utils::parser::queryWithParser(pp_grid, "nx_guard", ngx);
        utils::parser::queryWithParser(pp_grid, "ny_guard", ngy);
        utils::parser::queryWithParser(pp_grid, "nz_guard", ngz);

        return (amrex::IntVect(AMREX_D_DECL(ngx,ngy,ngz)));
    }

    /**
     * Perform derivative along x on a cell-centered grid, from a nodal field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const /*n_coefs_x*/,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;
#if (defined WARPX_DIM_1D_Z)
        amrex::ignore_unused(F, coefs_x, i, j, k, ncomp);
        return 0._rt; // 1D Cartesian: derivative along x is 0
#else
        amrex::Real const inv_dx = coefs_x[0];
        return inv_dx*( F(i+1,j,k,ncomp) - F(i,j,k,ncomp) );
#endif
    }

    /**
     * Perform derivative along x on a nodal grid, from a cell-centered field `F`*/
    template< typename T_Field>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDx (
        T_Field const& F,
        amrex::Real const * const coefs_x, int const /*n_coefs_x*/,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;
#if (defined WARPX_DIM_1D_Z)
        amrex::ignore_unused(F, coefs_x, i, j, k, ncomp);
        return 0._rt; // 1D Cartesian: derivative along x is 0
#else
        amrex::Real const inv_dx = coefs_x[0];
        const int ii = (i < F.end.x ? i : F.begin.x);
        const int im1 = (i > F.begin.x ? i - 1 : F.end.x - 1);
        return inv_dx*( F(ii,j,k,ncomp) - F(im1,j,k,ncomp) );
#endif
    }

    /**
     * Perform derivative along y on a cell-centered grid, from a nodal field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;
#if defined WARPX_DIM_3D
        Real const inv_dy = coefs_y[0];
        amrex::ignore_unused(n_coefs_y);
        return inv_dy*( F(i,j+1,k,ncomp) - F(i,j,k,ncomp) );
#elif (defined WARPX_DIM_XZ || WARPX_DIM_1D_Z)
        amrex::ignore_unused(F, coefs_y, n_coefs_y,
            i, j, k, ncomp);
        return 0._rt; // 1D and 2D Cartesian: derivative along y is 0
#else
        amrex::ignore_unused(F, coefs_y, n_coefs_y,
            i, j, k, ncomp);
#endif
    }

    /**
     * Perform derivative along y on a nodal grid, from a cell-centered field `F`*/
    template< typename T_Field>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDy (
        T_Field const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;
#if defined WARPX_DIM_3D
        Real const inv_dy = coefs_y[0];
        amrex::ignore_unused(n_coefs_y);
        return inv_dy*( F(i,j,k,ncomp) - F(i,j-1,k,ncomp) );
#elif (defined WARPX_DIM_XZ || WARPX_DIM_1D_Z)
        amrex::ignore_unused(F, coefs_y, n_coefs_y,
            i, j, k, ncomp);
        return 0._rt; // 1D and 2D Cartesian: derivative along y is 0
#else
        amrex::ignore_unused(F, coefs_y, n_coefs_y,
            i, j, k, ncomp);
#endif
    }

    /**
     * Perform derivative along z on a cell-centered grid, from a nodal field `F`*/
   AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const /*n_coefs_z*/,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;
        Real const inv_dz = coefs_z[0];
#if defined WARPX_DIM_3D
        return inv_dz*( F(i,j,k+1,ncomp) - F(i,j,k,ncomp) );
#elif (defined WARPX_DIM_XZ)
        return inv_dz*( F(i,j+1,k,ncomp) - F(i,j,k,ncomp) );
#elif (defined WARPX_DIM_1D_Z)
        return inv_dz*( F(i+1,j,k,ncomp) - F(i,j,k,ncomp) );
#endif
    }

    /**
     * Perform derivative along z on a nodal grid, from a cell-centered field `F`*/
    template< typename T_Field>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDz (
        T_Field const& F,
        amrex::Real const * const coefs_z, int const /*n_coefs_z*/,
        int const i, int const j, int const k, int const ncomp=0 ) {

        using namespace amrex;
        Real const inv_dz = coefs_z[0];
#if defined WARPX_DIM_3D
        return inv_dz*( F(i,j,k,ncomp) - F(i,j,k-1,ncomp) );
#elif (defined WARPX_DIM_XZ)
        const int jj = (j < F.end.y ? j : F.begin.y);
        const int jm1 = (j > F.begin.y ? j - 1 : F.end.y - 1);
        return inv_dz*( F(i,jj,k,ncomp) - F(i,jm1,k,ncomp) );
#elif (defined WARPX_DIM_1D_Z)
        return inv_dz*( F(i,j,k,ncomp) - F(i-1,j,k,ncomp) );
#endif
    }

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_YEE_H_
