#include "nonlinear_ccd.hpp"

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <tight_inclusion/ccd.hpp>

#include <stack>

// #define USE_FIXED_PIECES

namespace ipc {

#ifdef USE_FIXED_PIECES
static constexpr size_t FIXED_NUM_PIECES = 100l;
#else
static constexpr size_t MAX_NUM_SUBDIVISIONS = 1000l;
#endif

// ============================================================================

#ifdef IPC_TOOLKIT_WITH_FILIB
double IntervalNonlinearTrajectory::max_distance_from_linear(
    const double t0, const double t1) const
{
    // Estimate max t ∈ [0, 1] ‖ p((t1 - t0) * t + t0) - lerp(p(t0), p(t1), t) ‖
    const VectorMax3I p_t0 = (*this)(t0).cast<filib::Interval>();
    const VectorMax3I p_t1 = (*this)(t1).cast<filib::Interval>();

    double max_d = 0;

    constexpr double n = 100.0;
    double ti0 = t0;
    for (int i = 1; i <= n; i++) {
        const double ti1 = i / n * (t1 - t0) + t0;

        const VectorMax3I p = (*this)(filib::Interval(ti0, ti1));
        const filib::Interval d = norm(
            p - ((p_t1 - p_t0) * filib::Interval(i / n, (i + 1) / n) + p_t0));

        max_d = std::max(max_d, d.SUP);
        ti0 = ti1;
    }

    return max_d;
}
#endif

// ============================================================================

bool conservative_piecewise_linear_ccd(
    const std::function<double(const double)>& distance,
    const std::function<double(const double, const double)>&
        max_distance_from_linear,
    const std::function<bool(
        const double /*ti0*/,
        const double /*ti1*/,
        const double /*min_distance*/,
        const bool /*no_zero_toi*/,
        double& /*toi*/)>& linear_ccd,
    double& toi,
    const double tmax,
    const double min_sep_distance,
    const double conservative_rescaling)
{
    const double distance_t0 = distance(0);
    if (check_initial_distance(distance_t0, min_sep_distance, toi)) {
        return true;
    }
    assert(distance_t0 > min_sep_distance);

    double ti0 = 0;
    std::stack<double> ts;

// Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = FIXED_NUM_PIECES; i > 0; i--) {
        ts.push(i / double(FIXED_NUM_PIECES) * tmax);
    }
    int num_subdivisions = FIXED_NUM_PIECES;
#else
    ts.push(tmax);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        const double ti1 = ts.top();

        const double distance_ti0 = distance(ti0);

        // If distance has decreased by a factor and the toi is not near zero,
        // then we can call this a collision.
        if (distance_ti0 < (1 - conservative_rescaling) * distance_t0
            && ti0 >= CCD_SMALL_TOI) {
            toi = ti0;
            logger().trace(
                "Distance small enough distance_ti0={:g}; toi={:g}",
                distance_ti0, toi);
            return true;
        }

        double min_distance = max_distance_from_linear(ti0, ti1);

#ifndef USE_FIXED_PIECES
        // Check if the minimum distance is too large and we need to subdivide
        // (Large distances cause the slow CCD)
        if ((min_distance
             >= std::min((1 - conservative_rescaling) * distance_ti0, 0.01))
            && (num_subdivisions < MAX_NUM_SUBDIVISIONS || ti0 == 0)) {
            logger().trace(
                "Subdividing at ti=[{:g}, {:g}] min_distance={:g} distance_ti0={:g}",
                ti0, ti1, min_distance, distance_ti0);
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif

        min_distance += min_sep_distance;

        const bool is_impacting =
            linear_ccd(ti0, ti1, min_distance, /*no_zero_toi=*/ti0 == 0, toi);

        logger().trace(
            "Evaluated at ti=[{:g}, {:g}] min_distance={:g} distance_ti0={:g}; result={}{}",
            ti0, ti1, min_distance, distance_ti0, is_impacting,
            is_impacting ? fmt::format(" toi={:g}", (ti1 - ti0) * toi + ti0)
                         : "");

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            if (toi == 0) {
                // This is impossible because distance_t0 > min_sep_distance
                ts.push((ti1 + ti0) / 2);
                num_subdivisions++;
                continue;
            }
            return true;
        }

        ts.pop();
        ti0 = ti1;
    }

    return false;
}

// ============================================================================

bool point_point_nonlinear_ccd(
    const NonlinearTrajectory& p0,
    const NonlinearTrajectory& p1,
    double& toi,
    const double tmax,
    const double min_distance,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(point_point_distance(p0(t), p1(t)));
        },
        [&](const double t0, const double t1) {
            return std::max(
                p0.max_distance_from_linear(t0, t1),
                p1.max_distance_from_linear(t0, t1));
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool no_zero_toi, double& _toi) {
            double output_tolerance;
            return ticcd::edgeEdgeCCD(
                to_3D(p0(ti0)), to_3D(p0(ti0)), to_3D(p1(ti0)), to_3D(p1(ti0)),
                to_3D(p0(ti1)), to_3D(p0(ti1)), to_3D(p1(ti1)), to_3D(p1(ti1)),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                _min_distance,                // minimum separation distance
                _toi,                         // time of impact
                tolerance,                    // delta
                1.0,                          // maximum time to check
                max_iterations,               // maximum number of iterations
                output_tolerance,             // delta_actual
                no_zero_toi);                 // no zero toi
        },
        toi, tmax, min_distance, conservative_rescaling);
}

bool point_edge_nonlinear_ccd(
    const NonlinearTrajectory& p,
    const NonlinearTrajectory& e0,
    const NonlinearTrajectory& e1,
    double& toi,
    const double tmax,
    const double min_distance,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(point_edge_distance(p(t), e0(t), e1(t)));
        },
        [&](const double t0, const double t1) {
            return p.max_distance_from_linear(t0, t1)
                + std::max(
                       e0.max_distance_from_linear(t0, t1),
                       e1.max_distance_from_linear(t0, t1));
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool no_zero_toi, double& _toi) {
            double output_tolerance;
            return ticcd::edgeEdgeCCD(
                to_3D(p(ti0)), to_3D(p(ti0)), to_3D(e0(ti0)), to_3D(e1(ti0)),
                to_3D(p(ti1)), to_3D(p(ti1)), to_3D(e0(ti1)), to_3D(e1(ti1)),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                _min_distance,                // minimum separation distance
                _toi,                         // time of impact
                tolerance,                    // delta
                1.0,                          // maximum time to check
                max_iterations,               // maximum number of iterations
                output_tolerance,             // delta_actual
                no_zero_toi);                 // no zero toi
        },
        toi, tmax, min_distance, conservative_rescaling);
}

bool edge_edge_nonlinear_ccd(
    const NonlinearTrajectory& ea0,
    const NonlinearTrajectory& ea1,
    const NonlinearTrajectory& eb0,
    const NonlinearTrajectory& eb1,
    double& toi,
    const double tmax,
    const double min_sep_distance,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(edge_edge_distance(ea0(t), ea1(t), eb0(t), eb1(t)));
        },
        [&](const double t0, const double t1) {
            return std::max(
                       ea0.max_distance_from_linear(t0, t1),
                       ea1.max_distance_from_linear(t0, t1))
                + std::max(
                       eb0.max_distance_from_linear(t0, t1),
                       eb1.max_distance_from_linear(t0, t1));
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool no_zero_toi, double& _toi) {
            double output_tolerance;
            return ticcd::edgeEdgeCCD(
                ea0(ti0), ea1(ti0), eb0(ti0), eb1(ti0), //
                ea0(ti1), ea1(ti1), eb0(ti1), eb1(ti1),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                _min_distance,                // minimum separation distance
                _toi,                         // time of impact
                tolerance,                    // delta
                1.0,                          // maximum time to check
                max_iterations,               // maximum number of iterations
                output_tolerance,             // delta_actual
                no_zero_toi);                 // no zero toi
        },
        toi, tmax, min_sep_distance, conservative_rescaling);
}

bool point_triangle_nonlinear_ccd(
    const NonlinearTrajectory& p,
    const NonlinearTrajectory& t0,
    const NonlinearTrajectory& t1,
    const NonlinearTrajectory& t2,
    double& toi,
    const double tmax,
    const double min_distance,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(point_triangle_distance(p(t), t0(t), t1(t), t2(t)));
        },
        [&](const double ti0, const double ti1) {
            return p.max_distance_from_linear(ti0, ti1)
                + std::max({ t0.max_distance_from_linear(ti0, ti1),
                             t1.max_distance_from_linear(ti0, ti1),
                             t2.max_distance_from_linear(ti0, ti1) });
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool no_zero_toi, double& _toi) {
            double output_tolerance;
            return ticcd::vertexFaceCCD(
                p(ti0), t0(ti0), t1(ti0), t2(ti0), //
                p(ti1), t0(ti1), t1(ti1), t2(ti1),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                _min_distance,                // minimum separation distance
                _toi,                         // time of impact
                tolerance,                    // delta
                1.0,                          // maximum time to check
                max_iterations,               // maximum number of iterations
                output_tolerance,             // delta_actual
                no_zero_toi);                 // no zero toi
        },
        toi, tmax, min_distance, conservative_rescaling);
}

} // namespace ipc
