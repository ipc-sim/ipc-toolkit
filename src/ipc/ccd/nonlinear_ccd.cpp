#include "nonlinear_ccd.hpp"

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <tight_inclusion/ccd.hpp>

#include <stack>

// #define USE_FIXED_PIECES

namespace ipc {

namespace {
    inline Eigen::Vector3d to_3D(const VectorMax3d& v)
    {
        assert(v.size() == 2 || v.size() == 3);
        return v.size() == 2 ? Eigen::Vector3d(v.x(), v.y(), 0) : v.head<3>();
    }

    inline bool check_initial_distance(
        const double initial_distance, const double min_distance, double& toi)
    {
        if (initial_distance > min_distance) {
            return false;
        }

        logger().warn(
            "Initial distance {} ≤ d_min={}, returning toi=0!",
            initial_distance, min_distance);

        toi = 0; // Initially touching

        return true;
    }
} // namespace

// ============================================================================

#ifdef USE_FIXED_PIECES
static constexpr size_t FIXED_NUM_PIECES = 100;
#else
static constexpr size_t MAX_NUM_SUBDIVISIONS = 1'000l;
#endif

// Tolerance for small time of impact which triggers further refinement
static constexpr double SMALL_TOI = 1e-6;

// ============================================================================

double
NonlinearTrajectory::max_distance_from_linear(const filib::Interval& t) const
{
    // Estimate max t ∈ [0, 1] ‖ p((t1 - t0) * t + t0) - lerp(p(t0), p(t1), t) ‖
    const VectorMax3I p_t0 = (*this)(t.INF).cast<filib::Interval>();
    const VectorMax3I p_t1 = (*this)(t.SUP).cast<filib::Interval>();

    double max_d = 0;

    constexpr double n = 100.0;
    double ti0 = t.INF;
    for (int i = 1; i <= n; i++) {
        const double ti1 = i / n * (t.SUP - t.INF) + t.INF;

        const VectorMax3I p = (*this)(filib::Interval(ti0, ti1));
        const filib::Interval d = norm(
            p - ((p_t1 - p_t0) * filib::Interval(i / n, (i + 1) / n) + p_t0));

        max_d = std::max(max_d, d.SUP);
        ti0 = ti1;
    }

    return max_d;
}

// ============================================================================

bool conservative_piecewise_linear_ccd(
    const std::function<double(double)>& distance,
    const std::function<double(const filib::Interval&)>&
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
    const filib::Interval t(0, tmax);

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
        const filib::Interval ti(ti0, ti1);

        const double distance_ti0 = distance(ti0);

        // If distance has decreased by a factor and the toi is not near zero,
        // then we can call this a collision.
        if (distance_ti0 < (1 - conservative_rescaling) * distance_t0
            && ti0 >= SMALL_TOI) {
            logger().trace(
                "Distance small enough distance_ti0={:g}; toi={:g}",
                distance_ti0, toi);
            toi = ti0;
            return true;
        }

        double min_distance = max_distance_from_linear(ti);

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

        double output_tolerance;
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
        [&](const filib::Interval& t) {
            return std::max(
                p0.max_distance_from_linear(t), p1.max_distance_from_linear(t));
        },
        [&](const double ti0, const double ti1, const double min_distance,
            const bool no_zero_toi, double& toi) {
            double output_tolerance;
            return ticcd::edgeEdgeCCD(
                to_3D(p0(ti0)), to_3D(p0(ti0)), to_3D(p1(ti0)), to_3D(p1(ti0)),
                to_3D(p0(ti1)), to_3D(p0(ti1)), to_3D(p1(ti1)), to_3D(p1(ti1)),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                min_distance,                 // minimum separation distance
                toi,                          // time of impact
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
        [&](const filib::Interval& t) {
            return p.max_distance_from_linear(t)
                + std::max(
                       e0.max_distance_from_linear(t),
                       e1.max_distance_from_linear(t));
        },
        [&](const double ti0, const double ti1, const double min_distance,
            const bool no_zero_toi, double& toi) {
            double output_tolerance;
            return ticcd::edgeEdgeCCD(
                to_3D(p(ti0)), to_3D(p(ti0)), to_3D(e0(ti0)), to_3D(e1(ti0)),
                to_3D(p(ti1)), to_3D(p(ti1)), to_3D(e0(ti1)), to_3D(e1(ti1)),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                min_distance,                 // minimum separation distance
                toi,                          // time of impact
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
        [&](const filib::Interval& t) {
            return std::max(
                       ea0.max_distance_from_linear(t),
                       ea1.max_distance_from_linear(t))
                + std::max(
                       eb0.max_distance_from_linear(t),
                       eb1.max_distance_from_linear(t));
        },
        [&](const double ti0, const double ti1, const double min_distance,
            const bool no_zero_toi, double& toi) {
            double output_tolerance;
            return ticcd::edgeEdgeCCD(
                ea0(ti0), ea1(ti0), eb0(ti0), eb1(ti0), //
                ea0(ti1), ea1(ti1), eb0(ti1), eb1(ti1),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                min_distance,                 // minimum separation distance
                toi,                          // time of impact
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
        [&](const filib::Interval& t) {
            return p.max_distance_from_linear(t)
                + std::max({ t0.max_distance_from_linear(t),
                             t1.max_distance_from_linear(t),
                             t2.max_distance_from_linear(t) });
        },
        [&](const double ti0, const double ti1, const double min_distance,
            const bool no_zero_toi, double& toi) {
            double output_tolerance;
            return ticcd::vertexFaceCCD(
                p(ti0), t0(ti0), t1(ti0), t2(ti0), //
                p(ti1), t0(ti1), t1(ti1), t2(ti1),
                Eigen::Array3d::Constant(-1), // rounding error (auto)
                min_distance,                 // minimum separation distance
                toi,                          // time of impact
                tolerance,                    // delta
                1.0,                          // maximum time to check
                max_iterations,               // maximum number of iterations
                output_tolerance,             // delta_actual
                no_zero_toi);                 // no zero toi
        },
        toi, tmax, min_distance, conservative_rescaling);
}

} // namespace ipc
