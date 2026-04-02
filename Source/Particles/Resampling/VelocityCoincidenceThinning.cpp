/* Copyright 2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Particles/ShapeFactors.H"
#include "VelocityCoincidenceThinning.H"
#include "WarpX.H"


VelocityCoincidenceThinning::VelocityCoincidenceThinning (const std::string& species_name)
{
    using namespace amrex::literals;

    const amrex::ParmParse pp_species_name(species_name);

    utils::parser::queryWithParser(
        pp_species_name, "resampling_min_ppc", m_min_ppc
    );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_min_ppc >= 1,
        "Resampling min_ppc should be greater than or equal to 1"
    );

    amrex::ParticleReal target_weight = 0;
    if (utils::parser::queryWithParser(
        pp_species_name, "resampling_algorithm_target_weight", target_weight
    )) {
        // factor of 2 since each cluster is reduced to 2 particles
        m_cluster_weight = target_weight * 2.0_prt;
    }

    std::string velocity_grid_type_str = "spherical";
    pp_species_name.query(
        "resampling_algorithm_velocity_grid_type", velocity_grid_type_str
    );
    if (velocity_grid_type_str == "spherical") {
        m_velocity_grid_type = VelocityGridType::Spherical;
        utils::parser::getWithParser(
            pp_species_name, "resampling_algorithm_delta_ur", m_delta_ur
        );
        utils::parser::getWithParser(
            pp_species_name, "resampling_algorithm_n_theta", m_ntheta
        );
        utils::parser::getWithParser(
            pp_species_name, "resampling_algorithm_n_phi", m_nphi
        );
    }
    else if (velocity_grid_type_str == "cartesian") {
        m_velocity_grid_type = VelocityGridType::Cartesian;
        utils::parser::getArrWithParser(
            pp_species_name, "resampling_algorithm_delta_u", m_delta_u
        );
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_delta_u.size() == 3,
            "resampling_algorithm_delta_u must have three components.");
    }
    else {
        WARPX_ABORT_WITH_MESSAGE("Unkown velocity grid type.");
    }
}

void VelocityCoincidenceThinning::operator() (
    const amrex::Geometry& geom_lev, WarpXParIter& pti,
    const int lev, WarpXParticleContainer * const pc) const
{
    using namespace amrex::literals;

    int shape = WarpX::nox;

    auto& ptile = pc->ParticlesAt(lev, pti);
    const auto n_parts_in_tile = pti.numParticles();
    auto& soa = ptile.GetStructOfArrays();
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_3D)
    auto * const AMREX_RESTRICT x = soa.GetRealData(PIdx::x).data();
#elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    auto * const AMREX_RESTRICT x = soa.GetRealData(PIdx::x).data(); // rename to PIdx::r after PR #4667
#endif
#if defined(WARPX_DIM_3D)
    auto * const AMREX_RESTRICT y = soa.GetRealData(PIdx::y).data();
#endif
#if defined(WARPX_ZINDEX)
    auto * const AMREX_RESTRICT z = soa.GetRealData(PIdx::z).data();
#endif
    auto * const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
    auto * const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
    auto * const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();
    auto * const AMREX_RESTRICT w = soa.GetRealData(PIdx::w).data();
    auto * const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();

    amrex::Geometry geom_sub;
    amrex::Box cbx;
    if (shape % 2 == 0) {
        // Create a new geometry with cell sizes one half smaller
        const amrex::Box& dom = geom_lev.Domain();
        const int coord = geom_lev.Coord();

        // Physical bounds stay the same
        amrex::RealBox rb(geom_lev.ProbLo(), geom_lev.ProbHi());

        // Refine index-space domain by factor 2 in every direction
        amrex::Box fine_dom = amrex::refine(dom, 2);
        cbx = fine_dom;

        geom_sub = amrex::Geometry(fine_dom, rb, coord, geom_lev.isPeriodic());
    } else {
        geom_sub = geom_lev;
        cbx = pti.tilebox(amrex::IntVect::TheZeroVector());
    }
    auto const dxi = geom_lev.InvCellSizeArray();
    auto const dx = geom_lev.CellSizeArray();
    amrex::Box box = pti.validbox();
    const WarpX& warpx = WarpX::GetInstance();
    const amrex::IntVect& ng_rho = warpx.get_ng_depos_rho();
    box.grow(ng_rho);
    const amrex::XDim3 xyzmin = WarpX::LowerCorner(box, lev, 0._rt);
    Compute_shape_factor<2> const compute_shape_factor;

    // Using this function means that we must loop over the cells in the ParallelFor.
    auto bins = ParticleUtils::findParticlesInEachCell(geom_sub, cbx, ptile);

    const auto n_cells = static_cast<int>(bins.numBins());
    auto *const indices = bins.permutationPtr();
    auto *const cell_offsets = bins.offsetsPtr();

    const auto min_ppc = m_min_ppc;
    const auto cluster_weight = m_cluster_weight;
    const auto mass = pc->getMass();

    // check if species mass > 0
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        mass > 0,
        "VelocityCoincidenceThinning does not yet work for massless particles."
    );

    // create a GPU vector to hold the momentum cluster index for each particle
    amrex::Gpu::DeviceVector<int> momentum_bin_number(n_parts_in_tile);
    auto* momentum_bin_number_data = momentum_bin_number.data();

    // create a GPU vector to hold the index sorting for the momentum bins
    amrex::Gpu::DeviceVector<int> sorted_indices(n_parts_in_tile);
    auto* sorted_indices_data = sorted_indices.data();

    constexpr auto c2 = PhysConst::c * PhysConst::c;

    auto velocityBinCalculator = VelocityBinCalculator();
    velocityBinCalculator.velocity_grid_type = m_velocity_grid_type;
    if (m_velocity_grid_type == VelocityGridType::Spherical) {
        velocityBinCalculator.dur = m_delta_ur;
        velocityBinCalculator.n1 = m_ntheta;
        velocityBinCalculator.n2 = m_nphi;
        velocityBinCalculator.dutheta = 2.0_prt * MathConst::pi / m_ntheta;
        velocityBinCalculator.duphi = MathConst::pi / m_nphi;
    }
    else if (m_velocity_grid_type == VelocityGridType::Cartesian) {
        velocityBinCalculator.dux = m_delta_u[0];
        velocityBinCalculator.duy = m_delta_u[1];
        velocityBinCalculator.duz = m_delta_u[2];

        // get the minimum and maximum velocities to determine the velocity space
        // grid boundaries
        {
            using ReduceOpsT = amrex::TypeMultiplier<amrex::ReduceOps,
                                                     amrex::ReduceOpMin[3],
                                                     amrex::ReduceOpMax[2]>;
            using ReduceDataT = amrex::TypeMultiplier<amrex::ReduceData, amrex::ParticleReal[5]>;
            ReduceOpsT reduce_op;
            ReduceDataT reduce_data(reduce_op);
            using ReduceTuple = typename ReduceDataT::Type;
            reduce_op.eval(n_parts_in_tile, reduce_data, [=] AMREX_GPU_DEVICE(int i) -> ReduceTuple {
                return {ux[i], uy[i], uz[i], ux[i], uy[i]};
            });
            auto hv = reduce_data.value(reduce_op);
            velocityBinCalculator.ux_min = amrex::get<0>(hv);
            velocityBinCalculator.uy_min = amrex::get<1>(hv);
            velocityBinCalculator.uz_min = amrex::get<2>(hv);
            velocityBinCalculator.ux_max = amrex::get<3>(hv);
            velocityBinCalculator.uy_max = amrex::get<4>(hv);
        }

        velocityBinCalculator.n1 = static_cast<int>(
            std::ceil((velocityBinCalculator.ux_max - velocityBinCalculator.ux_min) / m_delta_u[0])
        );
        velocityBinCalculator.n2 = static_cast<int>(
            std::ceil((velocityBinCalculator.uy_max - velocityBinCalculator.uy_min) / m_delta_u[1])
        );
    }
    auto heapSort = HeapSort();

    // Loop over cells
    amrex::ParallelForRNG( n_cells,
        [=] AMREX_GPU_DEVICE (int i_cell, amrex::RandomEngine const& engine) noexcept
        {
            // The particles that are in the cell `i_cell` are
            // given by the `indices[cell_start:cell_stop]`
            const auto cell_start = static_cast<int>(cell_offsets[i_cell]);
            const auto cell_stop  = static_cast<int>(cell_offsets[i_cell+1]);
            const auto cell_numparts = cell_stop - cell_start;

            // do nothing for cells with less particles than min_ppc
            // (this intentionally includes skipping empty cells, too)
            if (cell_numparts < min_ppc) {
                return;
            }

            // Loop over particles and label them with the appropriate momentum bin
            // number. Also assign initial ordering to the sorted_indices array.
            velocityBinCalculator(
                ux, uy, uz, indices, momentum_bin_number_data, sorted_indices_data,
                cell_start, cell_stop
            );

            // sort indices based on comparing values in momentum_bin_number
            heapSort(sorted_indices_data, momentum_bin_number_data, cell_start, cell_numparts);

            // initialize variables used to hold cluster totals
            int particles_in_bin = 0;
            amrex::ParticleReal total_weight = 0._prt, total_energy = 0._prt;
#if !defined(WARPX_DIM_1D_Z)
            amrex::ParticleReal cluster_x = 0._prt;
#endif
#if defined(WARPX_DIM_3D)
            amrex::ParticleReal cluster_y = 0._prt;
#endif
#if defined(WARPX_ZINDEX)
            amrex::ParticleReal cluster_z = 0._prt;
#endif
            amrex::ParticleReal cluster_ux = 0._prt, cluster_uy = 0._prt, cluster_uz = 0._prt;

            // When shape == 2, save the deposition weights.
            // The location of the two particles will depend on these in order to conserve charge
            amrex::ParticleReal ss[3];
#if !defined(WARPX_DIM_1D_Z)
            amrex::ParticleReal cluster_wx[3] = {0._prt};
#endif
#if defined(WARPX_DIM_3D)
            amrex::ParticleReal cluster_wy[3] = {0._prt};
#endif
#if defined(WARPX_ZINDEX)
            amrex::ParticleReal cluster_wz[3] = {0._prt};
#endif

            // Finally, loop through the particles in the cell and merge
            // ones in the same momentum bin
            for (int i = cell_start; i < cell_stop; ++i)
            {
                particles_in_bin += 1;
                const auto part_idx = indices[sorted_indices_data[i]];

#if !defined(WARPX_DIM_1D_Z)
                cluster_x += w[part_idx]*x[part_idx];
#endif
#if defined(WARPX_DIM_3D)
                cluster_y += w[part_idx]*y[part_idx];
#endif
#if defined(WARPX_ZINDEX)
                cluster_z += w[part_idx]*z[part_idx];
#endif
                cluster_ux += w[part_idx]*ux[part_idx];
                cluster_uy += w[part_idx]*uy[part_idx];
                cluster_uz += w[part_idx]*uz[part_idx];
                total_weight += w[part_idx];
                total_energy += w[part_idx] * Algorithms::KineticEnergy(
                    ux[part_idx], uy[part_idx], uz[part_idx], mass
                );

                if (shape == 2) {
#if !defined(WARPX_DIM_1D_Z)
                    double const wx = (x[part_idx] - xyzmin.x)*dxi[0];
                    compute_shape_factor(ss, wx);
                    for (int iw = 0 ; iw < 3 ; iw++) {
                        cluster_wx[iw] += ss[iw]*w[part_idx];
                    }
#endif
#if defined(WARPX_DIM_3D)
                    double const wy = (y[part_idx] - xyzmin.y)*dxi[1];
                    compute_shape_factor(ss, wy);
                    for (int iw = 0 ; iw < 3 ; iw++) {
                        cluster_wy[iw] += ss[iw]*w[part_idx];
                    }
#endif
#if defined(WARPX_ZINDEX)
                    double const wz = (z[part_idx] - xyzmin.z)*dxi[WARPX_ZINDEX];
                    compute_shape_factor(ss, wz);
                    for (int iw = 0 ; iw < 3 ; iw++) {
                        cluster_wz[iw] += ss[iw]*w[part_idx];
                    }
#endif

                }

                // check if this is the last particle in the current momentum bin,
                // or if the next particle would push the current cluster weight
                // to exceed the maximum specified cluster weight
                if (
                    (i == cell_stop - 1)
                    || (momentum_bin_number_data[sorted_indices_data[i]] != momentum_bin_number_data[sorted_indices_data[i + 1]])
                    || (total_weight + w[indices[sorted_indices_data[i+1]]] > cluster_weight)
                ) {

                    // check if the bin has more than 2 particles in it
                    if ( particles_in_bin > 2 && total_weight > std::numeric_limits<amrex::ParticleReal>::min() ){

                        // get average quantities for the bin
#if !defined(WARPX_DIM_1D_Z)
                        cluster_x /= total_weight;
#endif
#if defined(WARPX_DIM_3D)
                        cluster_y /= total_weight;
#endif
#if defined(WARPX_ZINDEX)
                        cluster_z /= total_weight;
#endif
                        cluster_ux /= total_weight;
                        cluster_uy /= total_weight;
                        cluster_uz /= total_weight;

                        // perform merging of momentum bin particles
                        auto u_perp2 = cluster_ux*cluster_ux + cluster_uy*cluster_uy;
                        auto u_perp = std::sqrt(u_perp2);
                        auto cluster_u_mag2 = u_perp2 + cluster_uz*cluster_uz;
                        auto cluster_u_mag = std::sqrt(cluster_u_mag2);

                        // calculate required velocity magnitude to achieve
                        // energy conservation
                        auto v_mag2 = total_energy / total_weight * (
                            (total_energy / total_weight + 2._prt * mass * c2 )
                            / (mass * mass * c2)
                        );
                        auto v_perp = (v_mag2 > cluster_u_mag2) ? std::sqrt(v_mag2 - cluster_u_mag2) : 0_prt;

                        // choose random angle for new velocity vector
                        auto phi = amrex::Random(engine) * MathConst::pi;

                        // set new velocity components based on chosen phi
                        auto vx = v_perp * std::cos(phi);
                        auto vy = v_perp * std::sin(phi);

                        // calculate rotation angles to parallel coord. frame
                        auto cos_theta = (cluster_u_mag > 0._prt) ? cluster_uz / cluster_u_mag : 0._prt;
                        auto sin_theta = (cluster_u_mag > 0._prt) ? u_perp / cluster_u_mag : 0._prt;
                        auto cos_phi = (u_perp > 0._prt) ? cluster_ux / u_perp : 0._prt;
                        auto sin_phi = (u_perp > 0._prt) ? cluster_uy / u_perp : 0._prt;

                        // rotate new velocity vector to labframe
                        auto ux_new = (
                            vx * cos_theta * cos_phi - vy * sin_phi
                            + cluster_u_mag * sin_theta * cos_phi
                        );
                        auto uy_new = (
                            vx * cos_theta * sin_phi + vy * cos_phi
                            + cluster_u_mag * sin_theta * sin_phi
                        );
                        auto uz_new = -vx * sin_theta + cluster_u_mag * cos_theta;

                        // set the last two particles' attributes according to
                        // the bin's aggregate values
                        const auto part_idx2 = indices[sorted_indices_data[i - 1]];

                        if (shape == 2) {

                            // This generates the two particle's positions so that the deposited charge of the
                            // two particles is the same as that of the particles to be merged, thus conserving charge.
                            // The positions x1 and x2 are given by the solution of these two equations:
                            // W0 = 0.5*(0.5 - x1)**2*w1 + 0.5*(0.5 - x2)**2*w2
                            // W2 = 0.5*(0.5 + x1)**2*w1 + 0.5*(0.5 + x2)**2*w2
                            // x1 and x2 are scaled relative to the grid cell.
                            // W0 and W2 are the weights summed over the particles being merged.
                            // w1 and w2 are the weights of the two new particles.
                            // This solution assumes that the particles in the lower and upper halves of the grid
                            // cell are handled separately as each set of particles contributes to only three of the
                            // four terms with 2nd order weighting. In the lower half, the x1 and x2 should be between
                            // 0 and 1/2. In the upper half, the x1 and x2 are relative to the upper grid cell and
                            // should be in the range -1/2 to 0.
                            // In most cases, the solution (x1, x2) is good. However, if the distribution of merged
                            // particles is lopsided (with more particles towards one side than the other), one of
                            // the resulting x1 or x2 can land outside the half grid cell. In this case, the
                            // w1 and w2 are adjusted so that both lie within the half grid cell.
                            // In all cases, w1 is increased and w2 decreased - there are two pairs of solutions
                            // for the positions, and the pair is chosen so that increasing w1 will shift the position
                            // back into the half grid cell. The minimum w1 is calculated needed to shift the particle that
                            // is out of bounds to the boundary. The maximum w1 is also calculated that would shift the
                            // particle within the bounds to the other boundary, above which the particle would go out
                            // of bounds. The w1 is then set to the average of the minimum and the maximum.
                            // In 1D, there will always be a solution.
                            // In 2D and 3D, there may not be a solution that satisfies the constraint in all dimensions -
                            // in those cases, the merging is skipped.
                            auto adjust_weight_limits = [=] AMREX_GPU_DEVICE(int dir, amrex::ParticleReal gmin,
                                                                             amrex::ParticleReal cluster_avg,
                                                                             amrex::ParticleReal const * cluster_w,
                                                                             amrex::ParticleReal & w1_min,
                                                                             amrex::ParticleReal & w2_min)
                                                                             noexcept {
                                const amrex::ParticleReal W0 = cluster_w[0];
                                const amrex::ParticleReal W2 = cluster_w[2];
                                const amrex::ParticleReal A = W2 - W0;
                                const amrex::ParticleReal B = W0 + W2 - 0.25_prt*total_weight;

                                // ix will be used to determine which half of the cell the particles are in
                                const int ix = static_cast<int>((cluster_avg - gmin)*dxi[dir]*2._prt);
                                const amrex::ParticleReal xg_avg = (cluster_avg - gmin)*dxi[dir] - (ix % 2 == 1 ? 0.5_prt : 0._prt);
                                const bool x_lopsided_high = (xg_avg >= 0.25_prt);
                                // The lower and upper bounds of the half grid cell
                                const amrex::ParticleReal xlower = (ix % 2 == 0 ? 0._prt : -0.5_prt);
                                const amrex::ParticleReal xhigher = (ix % 2 == 0 ? 0.5_prt : 0._prt);
                                const amrex::ParticleReal numer_high = xhigher*total_weight - A;
                                const amrex::ParticleReal numer_low = A - xlower*total_weight;
                                const amrex::ParticleReal denom = (total_weight*B - A*A);
                                amrex::ParticleReal w1_min_x, w2_min_x;
                                if (x_lopsided_high) {
                                    w1_min_x = denom*total_weight/(denom + numer_high*numer_high);
                                    w2_min_x = denom*total_weight/(denom + numer_low*numer_low);
                                } else {
                                    w1_min_x = denom*total_weight/(denom + numer_low*numer_low);
                                    w2_min_x = denom*total_weight/(denom + numer_high*numer_high);
                                }
                                w1_min = std::max(w1_min, w1_min_x);
                                w2_min = std::max(w2_min, w2_min_x);
                            };

                            auto shape2_positions = [=] AMREX_GPU_DEVICE(int dir, amrex::ParticleReal gmin,
                                                                         amrex::ParticleReal cluster_avg,
                                                                         amrex::ParticleReal const * cluster_w,
                                                                         amrex::ParticleReal & w1,
                                                                         amrex::ParticleReal & w2,
                                                                         amrex::ParticleReal & pos1,
                                                                         amrex::ParticleReal & pos2) noexcept {

                                const int ix = static_cast<int>((cluster_avg - gmin)*dxi[dir]*2._prt);
                                const amrex::ParticleReal xg_avg = (cluster_avg - gmin)*dxi[dir] - (ix % 2 == 1 ? 0.5_prt : 0._prt);
                                const bool x_lopsided_high = (xg_avg >= 0.25_prt);
                                const amrex::ParticleReal W0 = cluster_w[0];
                                const amrex::ParticleReal W2 = cluster_w[2];
                                const amrex::ParticleReal A = W2 - W0;
                                const amrex::ParticleReal B = W0 + W2 - 0.25_prt*total_weight;

                                amrex::ParticleReal x1, x2;
                                if (x_lopsided_high) {
                                    x1 = (A + std::sqrt(w2/w1*(total_weight*B - A*A)))/total_weight;
                                    x2 = (A - std::sqrt(w1/w2*(total_weight*B - A*A)))/total_weight;
                                } else {
                                    x1 = (A - std::sqrt(w2/w1*(total_weight*B - A*A)))/total_weight;
                                    x2 = (A + std::sqrt(w1/w2*(total_weight*B - A*A)))/total_weight;
                                }
                                const int i_sub = static_cast<int>((cluster_avg - gmin)*dxi[dir] + 0.5_prt);
                                pos1 = gmin + (i_sub + x1)*dx[dir];
                                pos2 = gmin + (i_sub + x2)*dx[dir];
                            };

                            amrex::ParticleReal w1_min = 0._prt;
                            amrex::ParticleReal w2_min = 0._prt;

#if !defined(WARPX_DIM_1D_Z)
                            adjust_weight_limits(0, xyzmin.x, cluster_x, cluster_wx, w1_min, w2_min);
#endif
#if defined(WARPX_DIM_3D)
                            adjust_weight_limits(1, xyzmin.y, cluster_y, cluster_wy, w1_min, w2_min);
#endif
#if defined(WARPX_ZINDEX)
                            adjust_weight_limits(WARPX_ZINDEX, xyzmin.z, cluster_z, cluster_wz, w1_min, w2_min);
#endif

                            amrex::ParticleReal w1_max = total_weight - w2_min;

                            if (w1_min < w1_max) {
                                // Only do the merge if weight constraints is satified in all dimensions

                                amrex::ParticleReal w1, w2;

                                if (w1_min < total_weight*0.5_prt && w1_max > total_weight*0.5_prt) {
                                    // If possible, use equally weights particles
                                    w1 = total_weight*0.5_prt;
                                    w2 = total_weight*0.5_prt;
                                } else {
                                    w1 = 0.5_prt*(w1_min + w1_max);
                                    w2 = total_weight - w1;
                                }

#if !defined(WARPX_DIM_1D_Z)
                                shape2_positions(0, xyzmin.x, cluster_x, cluster_wx, w1, w2, x[part_idx], x[part_idx2]);
#endif
#if defined(WARPX_DIM_3D)
                                shape2_positions(1, xyzmin.y, cluster_y, cluster_wy, w1, w2, y[part_idx], y[part_idx2]);
#endif
#if defined(WARPX_ZINDEX)
                                shape2_positions(WARPX_ZINDEX, xyzmin.z, cluster_z, cluster_wz, w1, w2, z[part_idx], z[part_idx2]);
#endif

                                w[part_idx] = w1;
                                w[part_idx2] = w2;

                                ux[part_idx] = ux_new;
                                uy[part_idx] = uy_new;
                                uz[part_idx] = uz_new;
                                ux[part_idx2] = 2._prt * cluster_ux - ux_new;
                                uy[part_idx2] = 2._prt * cluster_uy - uy_new;
                                uz[part_idx2] = 2._prt * cluster_uz - uz_new;

                                // set ids of merged particles so they will be removed
                                for (int j = 2; j < particles_in_bin; ++j){
                                    idcpu[indices[sorted_indices_data[i - j]]] = amrex::ParticleIdCpus::Invalid;
                                }
                            }

                        } else {

                            w[part_idx] = total_weight / 2._prt;
                            w[part_idx2] = total_weight / 2._prt;

#if !defined(WARPX_DIM_1D_Z)
                            x[part_idx] = cluster_x;
                            x[part_idx2] = cluster_x;
#endif
#if defined(WARPX_DIM_3D)
                            y[part_idx] = cluster_y;
                            y[part_idx2] = cluster_y;
#endif
#if defined(WARPX_ZINDEX)
                            z[part_idx] = cluster_z;
                            z[part_idx2] = cluster_z;
#endif

                            ux[part_idx] = ux_new;
                            uy[part_idx] = uy_new;
                            uz[part_idx] = uz_new;
                            ux[part_idx2] = 2._prt * cluster_ux - ux_new;
                            uy[part_idx2] = 2._prt * cluster_uy - uy_new;
                            uz[part_idx2] = 2._prt * cluster_uz - uz_new;

                            // set ids of merged particles so they will be removed
                            for (int j = 2; j < particles_in_bin; ++j){
                                idcpu[indices[sorted_indices_data[i - j]]] = amrex::ParticleIdCpus::Invalid;
                            }

                        }
                    }

                    // restart the tallies
                    particles_in_bin = 0;
                    total_weight = 0._prt;
                    total_energy = 0._prt;
#if !defined(WARPX_DIM_1D_Z)
                    cluster_x = 0_prt;
#endif
#if defined(WARPX_DIM_3D)
                    cluster_y = 0_prt;
#endif
#if defined(WARPX_ZINDEX)
                    cluster_z = 0_prt;
#endif
                    cluster_ux = 0_prt;
                    cluster_uy = 0_prt;
                    cluster_uz = 0_prt;

                    if (shape == 2) {
                        for (int iw=0 ; iw < 3 ; iw++) {
#if !defined(WARPX_DIM_1D_Z)
                            cluster_wx[iw] = 0._prt;
#endif
#if defined(WARPX_DIM_3D)
                            cluster_wy[iw] = 0._prt;
#endif
#if defined(WARPX_ZINDEX)
                            cluster_wz[iw] = 0._prt;
#endif
                        }
                    }

                }
            }
        }
    );
}
