#include "ccd.hpp"

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
#include <ipc/ccd/inexact_point_edge.hpp>
#include <CTCD.h>
#else
#include <tight_inclusion/ccd.hpp>
#endif

#include <algorithm> // std::min/max
#include <array>

namespace ipc {

#ifndef IPC_TOOLKIT_WITH_INEXACT_CCD
/// Scale the distance tolerance to be at most this fraction of the initial
/// distance.
static constexpr double INITIAL_DISTANCE_TOLERANCE_SCALE = 0.5;
#endif

/// Special value for max_iterations to run tight inclusion without a maximum
/// number of iterations.
static constexpr long TIGHT_INCLUSION_UNLIMITED_ITERATIONS = -1;

bool ccd_strategy(
    const std::function<bool(
        long /*max_iterations*/,
        double /*min_distance*/,
        bool /*no_zero_toi*/,
        double& /*toi*/)>& ccd,
    const long max_iterations,
    const double min_distance,
    const double initial_distance,
    const double conservative_rescaling,
    double& toi)
{
    if (check_initial_distance(initial_distance, min_distance, toi)) {
        return true;
    }

    double min_effective_distance =
        (1.0 - conservative_rescaling) * (initial_distance - min_distance);
#ifndef IPC_TOOLKIT_WITH_INEXACT_CCD
    // Tight Inclusion performs better when the minimum separation is small
    min_effective_distance = std::min(min_effective_distance, 1e-4);
#endif
    min_effective_distance += min_distance;

    assert(min_effective_distance < initial_distance);

    // Do not use no_zero_toi because the minimum distance is arbitrary and can
    // be removed if the query is challenging (i.e., produces small ToI).
    bool is_impacting =
        ccd(max_iterations, min_effective_distance, /*no_zero_toi=*/false, toi);

    // #ifndef IPC_TOOLKIT_WITH_INEXACT_CCD
    //     // Tight inclusion will have higher accuracy and better performance
    //     // if we shrink the minimum distance. The value 1e-10 is arbitrary.
    //     while (is_impacting && toi < CCD_SMALL_TOI && min_distance > 1e-10) {
    //         min_distance /= 10;
    //         is_impacting =
    //             ccd(max_iterations, min_distance, /*no_zero_toi=*/false,
    //             toi);
    //     }
    // #endif

    if (is_impacting && toi < CCD_SMALL_TOI) {
        is_impacting = ccd(
            /*max_iterations=*/TIGHT_INCLUSION_UNLIMITED_ITERATIONS,
            /*min_distance=*/min_distance, /*no_zero_toi=*/true, toi);

        if (is_impacting) {
            toi *= conservative_rescaling;
            assert(toi != 0);
        }
    }

    return is_impacting;
}

bool point_point_ccd_3D(
    const Eigen::Vector3d& p0_t0,
    const Eigen::Vector3d& p1_t0,
    const Eigen::Vector3d& p0_t1,
    const Eigen::Vector3d& p1_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance = sqrt(point_point_distance(p0_t0, p1_t0));

    if (p0_t0 == p0_t1 && p1_t0 == p1_t1) { // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

#ifndef IPC_TOOLKIT_WITH_INEXACT_CCD
    const double adjusted_tolerance = std::min(
        INITIAL_DISTANCE_TOLERANCE_SCALE * initial_distance, tolerance);
#endif

    const auto ccd = [&](long _max_iterations, double _min_distance,
                         bool no_zero_toi, double& _toi) -> bool {
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
        return CTCD::vertexVertexCTCD(
            p0_t0, p1_t0, p0_t1, p1_t1, min_distance, toi);
#else
        double output_tolerance;
        // NOTE: Use degenerate edge-edge
        return ticcd::edgeEdgeCCD(
            p0_t0, p0_t0, p1_t0, p1_t0, p0_t1, p0_t1, p1_t1, p1_t1,
            Eigen::Array3d::Constant(-1), // rounding error (auto)
            _min_distance,                // minimum separation distance
            _toi,                         // time of impact
            adjusted_tolerance,           // delta
            tmax,                         // maximum time to check
            _max_iterations,              // maximum number of iterations
            output_tolerance,             // delta_actual
            no_zero_toi);
#endif
    };

    return ccd_strategy(
        ccd, max_iterations, min_distance, initial_distance,
        conservative_rescaling, toi);
}

bool point_point_ccd(
    const VectorMax3d& p0_t0,
    const VectorMax3d& p1_t0,
    const VectorMax3d& p0_t1,
    const VectorMax3d& p1_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    assert(p0_t0.size() == p1_t0.size());
    assert(p0_t0.size() == p0_t1.size());
    assert(p0_t0.size() == p1_t1.size());
    return point_point_ccd_3D(
        to_3D(p0_t0), to_3D(p1_t0), to_3D(p0_t1), to_3D(p1_t1), toi,
        min_distance, tmax, tolerance, max_iterations, conservative_rescaling);
}

bool point_edge_ccd_3D(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& e0_t0,
    const Eigen::Vector3d& e1_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& e0_t1,
    const Eigen::Vector3d& e1_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance =
        sqrt(point_edge_distance(p_t0, e0_t0, e1_t0));

    if (p_t0 == p_t1 && e0_t0 == e0_t1 && e1_t0 == e1_t1) { // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

#ifndef IPC_TOOLKIT_WITH_INEXACT_CCD
    const double adjusted_tolerance = std::min(
        INITIAL_DISTANCE_TOLERANCE_SCALE * initial_distance, tolerance);
#endif

    const auto ccd = [&](long _max_iterations, double _min_distance,
                         bool no_zero_toi, double& _toi) -> bool {
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
        return CTCD::vertexEdgeCTCD(
            p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, min_distance, toi);
#else
        double output_tolerance = tolerance;
        // NOTE: Use degenerate edge-edge
        bool is_impacting = ticcd::edgeEdgeCCD(
            p_t0, p_t0, e0_t0, e1_t0, p_t1, p_t1, e0_t1, e1_t1,
            Eigen::Array3d::Constant(-1), // rounding error (auto)
            _min_distance,                // minimum separation distance
            _toi,                         // time of impact
            adjusted_tolerance,           // delta
            tmax,                         // maximum time to check
            _max_iterations,              // maximum number of iterations
            output_tolerance,             // delta_actual
            no_zero_toi);
        if (adjusted_tolerance < output_tolerance && toi < CCD_SMALL_TOI) {
            logger().trace(
                "ticcd::edgeEdgeCCD exceeded iteration limit (min_dist={:g} "
                "max_iterations={:d} input_tol={:g} output_tol={:g} toi={:g})",
                min_distance, max_iterations, adjusted_tolerance,
                output_tolerance, toi);
        }
        return is_impacting;
#endif
    };

    return ccd_strategy(
        ccd, max_iterations, min_distance, initial_distance,
        conservative_rescaling, toi);
}

bool point_edge_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    assert(p_t1.size() == p_t0.size());
    assert(e0_t0.size() == p_t0.size() && e1_t0.size() == p_t0.size());
    assert(e0_t1.size() == p_t0.size() && e1_t1.size() == p_t0.size());
    return point_edge_ccd_3D(
        to_3D(p_t0), to_3D(e0_t0), to_3D(e1_t0), to_3D(p_t1), to_3D(e0_t1),
        to_3D(e1_t1), toi, min_distance, tmax, tolerance, max_iterations,
        conservative_rescaling);
}

bool edge_edge_ccd(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance =
        sqrt(edge_edge_distance(ea0_t0, ea1_t0, eb0_t0, eb1_t0));

    if (ea0_t0 == ea0_t1 && ea1_t0 == ea1_t1 && eb0_t0 == eb0_t1
        && eb1_t0 == eb1_t1) { // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

#ifndef IPC_TOOLKIT_WITH_INEXACT_CCD
    const double adjusted_tolerance = std::min(
        INITIAL_DISTANCE_TOLERANCE_SCALE * initial_distance, tolerance);
#endif

    const auto ccd = [&](long _max_iterations, double _min_distance,
                         bool no_zero_toi, double& _toi) -> bool {
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
        return CTCD::edgeEdgeCTCD(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
            min_distance, toi);
#else
        double output_tolerance;
        bool is_impacting = ticcd::edgeEdgeCCD(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
            Eigen::Array3d::Constant(-1), // rounding error (auto)
            _min_distance,                // minimum separation distance
            _toi,                         // time of impact
            adjusted_tolerance,           // delta
            tmax,                         // maximum time to check
            _max_iterations,              // maximum number of iterations
            output_tolerance,             // delta_actual
            no_zero_toi);
        if (adjusted_tolerance < output_tolerance && toi < CCD_SMALL_TOI) {
            logger().trace(
                "ticcd::edgeEdgeCCD exceeded iteration limit (min_dist={:g} "
                "max_iterations={:d} input_tol={:g} output_tol={:g} toi={:g})",
                min_distance, max_iterations, adjusted_tolerance,
                output_tolerance, toi);
        }
        return is_impacting;
#endif
    };

    return ccd_strategy(
        ccd, max_iterations, min_distance, initial_distance,
        conservative_rescaling, toi);
}

bool point_triangle_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& t0_t0,
    const Eigen::Vector3d& t1_t0,
    const Eigen::Vector3d& t2_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& t0_t1,
    const Eigen::Vector3d& t1_t1,
    const Eigen::Vector3d& t2_t1,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling)
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance =
        sqrt(point_triangle_distance(p_t0, t0_t0, t1_t0, t2_t0));

    if (p_t0 == p_t1 && t0_t0 == t0_t1 && t1_t0 == t1_t1 && t2_t0 == t2_t1) {
        // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

#ifndef IPC_TOOLKIT_WITH_INEXACT_CCD
    const double adjusted_tolerance = std::min(
        INITIAL_DISTANCE_TOLERANCE_SCALE * initial_distance, tolerance);
#endif

    const auto ccd = [&](long _max_iterations, double _min_distance,
                         bool no_zero_toi, double& _toi) -> bool {
#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
        return CTCD::vertexFaceCTCD(
            p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, //
            min_distance, toi);
#else
        double output_tolerance;
        bool is_impacting = ticcd::vertexFaceCCD(
            p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1,
            Eigen::Array3d::Constant(-1), // rounding error (auto)
            _min_distance,                // minimum separation distance
            _toi,                         // time of impact
            adjusted_tolerance,           // delta
            tmax,                         // maximum time to check
            _max_iterations,              // maximum number of iterations
            output_tolerance,             // delta_actual
            no_zero_toi);
        if (adjusted_tolerance < output_tolerance && toi < CCD_SMALL_TOI) {
            logger().trace(
                "ticcd::vertexFaceCCD exceeded iteration limit (min_dist={:g} "
                "max_iterations={:d} input_tol={:g} output_tol={:g} toi={:g})",
                min_distance, max_iterations, adjusted_tolerance,
                output_tolerance, toi);
        }
        return is_impacting;
#endif
    };

    return ccd_strategy(
        ccd, max_iterations, min_distance, initial_distance,
        conservative_rescaling, toi);
}

// -----------------------------------------------------------------------------

bool check_initial_distance(
    const double initial_distance, const double min_distance, double& toi)
{
    if (initial_distance > min_distance) {
        return false;
    }

    logger().warn(
        "Initial distance {} ≤ d_min={}, returning toi=0!", initial_distance,
        min_distance);

    toi = 0; // Initially touching

    return true;
}

} // namespace ipc
