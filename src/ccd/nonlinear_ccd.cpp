#pragma once

#include <ipc/ccd/ccd.hpp>

#include <functional>

namespace ipc {

bool edge_edge_nonlinear_ccd(
    const std::function<Vector3I(const Interval&)>& ea0, // ea0([t0, t1])
    const std::function<Vector3I(const Interval&)>& ea1, // ea1([t0, t1])
    const std::function<Vector3I(const Interval&)>& eb0, // eb0([t0, t1])
    const std::function<Vector3I(const Interval&)>& eb1, // eb1([t0, t1])
    double& toi,
    const double tmax,
    const double min_distance,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    int dim = 3; // TODO
    assert(dim == 3);
    assert(min_distance >= 0);

    const long ea0i = bodyA.edges(edgeA_id, 0);
    const long ea1i = bodyA.edges(edgeA_id, 1);
    const long eb0i = bodyB.edges(edgeB_id, 0);
    const long eb1i = bodyB.edges(edgeB_id, 1);

    const interval t0 = _interval(0);
    const interval t1 = _interval(tmax);
    const interval t = _interval(0, tmax);

    double distance_t0 =
        sqrt(edge_edge_distance(ea0(t0), ea1(t0), eb0(t0), eb1(t0)));
    if (distance_t0 <= minimum_separation_distance) {
        spdlog::warn(
            "initial distance in edge-edge CCD is less than MS={:g}!",
            minimum_separation_distance);
        toi = 0;
        return true;
    }
    assert(distance_t0 > minimum_separation_distance);

#ifdef TIME_CCD_QUERIES
    igl::Timer timer;
    timer.start();
#endif

    bool is_impacting = false;
    PoseD poseA_ti0 = poseA_t0, poseB_ti0 = poseB_t0;
    double ti0 = 0;
    std::stack<double> ts;

    // Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = FIXED_NUM_PIECES; i > 0; i--) {
        ts.push(i / double(FIXED_NUM_PIECES) * earliest_toi);
    }
    int num_subdivisions = FIXED_NUM_PIECES;
#else
    ts.push(earliest_toi);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        double ti1 = ts.top();

        PoseD poseA_ti1 = PoseD::interpolate(poseA_t0, poseA_t1, ti1);
        PoseD poseB_ti1 = PoseD::interpolate(poseB_t0, poseB_t1, ti1);

        double distance_ti0 = sqrt(edge_edge_distance(
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 0)),
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 1)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 0)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 1))));

#ifdef USE_DECREASING_DISTANCE_CHECK
        if (distance_ti0 < DECREASING_DISTANCE_FACTOR * distance_t0
            && ti0 >= DECREASING_DISTANCE_MIN_TIME) {
            toi = ti0;
            is_impacting = true;
            break;
        }
#endif

        double min_distance = 0;
#ifndef USE_FIXED_PIECES
        Interval ti(0, 1);

        PoseI poseIA_ti0 = poseA_ti0.cast<Interval>();
        PoseI poseIA_ti1 = poseA_ti1.cast<Interval>();
        PoseI poseIB_ti0 = poseB_ti0.cast<Interval>();
        PoseI poseIB_ti1 = poseB_ti1.cast<Interval>();

        PoseI poseIA = PoseI::interpolate(poseIA_ti0, poseIA_ti1, ti);
        PoseI poseIB = PoseI::interpolate(poseIB_ti0, poseIB_ti1, ti);

        MatrixMax3I RA = poseIA.construct_rotation_matrix();
        MatrixMax3I RB = poseIB.construct_rotation_matrix();

        double min_ea_distance = 0;
        Vector3I ea0_ti0 = bodyA.world_vertex(poseA_ti0, ea0i).cast<Interval>();
        Vector3I ea0_ti1 = bodyA.world_vertex(poseA_ti1, ea0i).cast<Interval>();
        Vector3I ea0 = bodyA.world_vertex(RA, poseIA.position, ea0i);
        Interval d = (ea0 - ((ea0_ti1 - ea0_ti0) * ti + ea0_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_ea_distance = std::max(min_ea_distance, d.upper());

        Vector3I ea1_ti0 = bodyA.world_vertex(poseA_ti0, ea1i).cast<Interval>();
        Vector3I ea1_ti1 = bodyA.world_vertex(poseA_ti1, ea1i).cast<Interval>();
        Vector3I ea1 = bodyA.world_vertex(RA, poseIA.position, ea1i);
        d = (ea1 - ((ea1_ti1 - ea1_ti0) * ti + ea1_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_ea_distance = std::max(min_ea_distance, d.upper());

        double min_eb_distance = 0;
        Vector3I eb0_ti0 = bodyB.world_vertex(poseB_ti0, eb0i).cast<Interval>();
        Vector3I eb0_ti1 = bodyB.world_vertex(poseB_ti1, eb0i).cast<Interval>();
        Vector3I eb0 = bodyB.world_vertex(RB, poseIB.position, eb0i);
        d = (eb0 - ((eb0_ti1 - eb0_ti0) * ti + eb0_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_eb_distance = std::max(min_eb_distance, d.upper());

        Vector3I eb1_ti0 = bodyB.world_vertex(poseB_ti0, eb1i).cast<Interval>();
        Vector3I eb1_ti1 = bodyB.world_vertex(poseB_ti1, eb1i).cast<Interval>();
        Vector3I eb1 = bodyB.world_vertex(RB, poseIB.position, eb1i);
        d = (eb1 - ((eb1_ti1 - eb1_ti0) * ti + eb1_ti0)).norm();
        assert(abs(d.lower()) < 1e-12); // The endpoints are part of both curves
        min_eb_distance = std::max(min_eb_distance, d.upper());

        min_distance = min_ea_distance + min_eb_distance;

        if (min_distance >= TRAJECTORY_DISTANCE_FACTOR * distance_ti0
            && (num_subdivisions < MAX_NUM_SUBDIVISIONS || ti0 == 0)) {
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif
        min_distance += minimum_separation_distance;
        // spdlog::critical("min_distance={:g}", min_distance);

        double output_tolerance;
        // 0: normal ccd method which only checks t = [0,1]
        // 1: ccd with max_itr and t=[0, t_max]
        const int CCD_TYPE = 1;
        is_impacting = inclusion_ccd::edgeEdgeCCD_double(
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 0)),
            bodyA.world_vertex(poseA_ti0, bodyA.edges(edgeA_id, 1)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 0)),
            bodyB.world_vertex(poseB_ti0, bodyB.edges(edgeB_id, 1)),
            bodyA.world_vertex(poseA_ti1, bodyA.edges(edgeA_id, 0)),
            bodyA.world_vertex(poseA_ti1, bodyA.edges(edgeA_id, 1)),
            bodyB.world_vertex(poseB_ti1, bodyB.edges(edgeB_id, 0)),
            bodyB.world_vertex(poseB_ti1, bodyB.edges(edgeB_id, 1)),
            { { -1, -1, -1 } },        // rounding error
            min_distance,              // minimum separation distance
            toi,                       // time of impact
            LINEAR_CCD_TOL,            // delta
            1.0,                       // Maximum time to check
            LINEAR_CCD_MAX_ITERATIONS, // Maximum number of iterations
            output_tolerance,          // delta_actual
            TIGHT_INCLUSION_CCD_TYPE,  //
            TIGHT_INCLUSION_NO_ZERO_TOI);

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            if (toi == 0) {
                // This is impossible because distance_t0 > MS_DIST
                ts.push((ti1 + ti0) / 2);
                num_subdivisions++;
                // spdlog::warn(
                //     "failure=\"Edge-edge MSCCD says toi=0, but "
                //     "distance_t0={0:g}\" failsafe=\"spliting [{1:g}, "
                //     "{2:g}] into ([{1:g}, {3:g}], [{3:g}, {2:g}])\"",
                //     distance_t0, ti0, ti1, ts.top());
                continue;
            }
            break;
        }

        ts.pop();
        ti0 = ti1;
        poseA_ti0 = poseA_ti1;
        poseB_ti0 = poseB_ti1;
    }
    // spdlog::trace("ee_ccd_num_subdivision={:d}", num_subdivisions);

#ifdef TIME_CCD_QUERIES
    timer.stop();
    fmt::print("EE {:.16g}\n", timer.getElapsedTime());
#endif

    // This time of impact is very dangerous for convergence
    assert(!is_impacting || toi > 0);
    return is_impacting;
}

} // namespace ipc
