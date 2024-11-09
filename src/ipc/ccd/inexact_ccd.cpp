#include "inexact_ccd.hpp"

#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD

#include <ipc/ccd/check_initial_distance.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <CTCD.h>

namespace ipc {

InexactCCD::InexactCCD(const double conservative_rescaling)
    : conservative_rescaling(conservative_rescaling)
{
}

bool InexactCCD::ccd_strategy(
    const std::function<bool(double /*min_distance*/, double& /*toi*/)>& ccd,
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
    min_effective_distance += min_distance;

    assert(min_effective_distance < initial_distance);

    // Do not use no_zero_toi because the minimum distance is arbitrary and can
    // be removed if the query is challenging (i.e., produces small ToI).
    bool is_impacting = ccd(min_effective_distance, toi);

    if (is_impacting && toi < SMALL_TOI) {
        is_impacting = ccd(/*min_distance=*/min_distance, toi);

        if (is_impacting) {
            toi *= conservative_rescaling;
            assert(toi != 0);
        }
    }

    return is_impacting;
}

bool InexactCCD::point_point_ccd_3D(
    const Eigen::Vector3d& p0_t0,
    const Eigen::Vector3d& p1_t0,
    const Eigen::Vector3d& p0_t1,
    const Eigen::Vector3d& p1_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance = sqrt(point_point_distance(p0_t0, p1_t0));

    if (p0_t0 == p0_t1 && p1_t0 == p1_t1) { // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

    const auto ccd = [&](double _min_distance, double& _toi) -> bool {
        return CTCD::vertexVertexCTCD(
            p0_t0, p1_t0, p0_t1, p1_t1, _min_distance, _toi);
    };

    return ccd_strategy(
        ccd, min_distance, initial_distance, conservative_rescaling, toi);
}

bool InexactCCD::point_point_ccd(
    const VectorMax3d& p0_t0,
    const VectorMax3d& p1_t0,
    const VectorMax3d& p0_t1,
    const VectorMax3d& p1_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    assert(p0_t0.size() == p1_t0.size());
    assert(p0_t0.size() == p0_t1.size());
    assert(p0_t0.size() == p1_t1.size());
    return point_point_ccd_3D(
        to_3D(p0_t0), to_3D(p1_t0), to_3D(p0_t1), to_3D(p1_t1), toi,
        min_distance, tmax);
}

bool InexactCCD::point_edge_ccd_3D(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& e0_t0,
    const Eigen::Vector3d& e1_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& e0_t1,
    const Eigen::Vector3d& e1_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance =
        sqrt(point_edge_distance(p_t0, e0_t0, e1_t0));

    if (p_t0 == p_t1 && e0_t0 == e0_t1 && e1_t0 == e1_t1) { // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

    const auto ccd = [&](double _min_distance, double& _toi) -> bool {
        return CTCD::vertexEdgeCTCD(
            p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, _min_distance, _toi);
    };

    return ccd_strategy(
        ccd, min_distance, initial_distance, conservative_rescaling, toi);
}

bool InexactCCD::point_edge_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    assert(p_t1.size() == p_t0.size());
    assert(e0_t0.size() == p_t0.size() && e1_t0.size() == p_t0.size());
    assert(e0_t1.size() == p_t0.size() && e1_t1.size() == p_t0.size());
    return point_edge_ccd_3D(
        to_3D(p_t0), to_3D(e0_t0), to_3D(e1_t0), to_3D(p_t1), to_3D(e0_t1),
        to_3D(e1_t1), toi, min_distance, tmax);
}

bool InexactCCD::edge_edge_ccd(
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
    const double tmax) const
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance =
        sqrt(edge_edge_distance(ea0_t0, ea1_t0, eb0_t0, eb1_t0));

    if (ea0_t0 == ea0_t1 && ea1_t0 == ea1_t1 && eb0_t0 == eb0_t1
        && eb1_t0 == eb1_t1) { // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

    const auto ccd = [&](double _min_distance, double& _toi) -> bool {
        return CTCD::edgeEdgeCTCD(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
            _min_distance, _toi);
    };

    return ccd_strategy(
        ccd, min_distance, initial_distance, conservative_rescaling, toi);
}

bool InexactCCD::point_triangle_ccd(
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
    const double tmax) const
{
    assert(tmax >= 0 && tmax <= 1.0);

    const double initial_distance =
        sqrt(point_triangle_distance(p_t0, t0_t0, t1_t0, t2_t0));

    if (p_t0 == p_t1 && t0_t0 == t0_t1 && t1_t0 == t1_t1 && t2_t0 == t2_t1) {
        // No motion
        return check_initial_distance(initial_distance, min_distance, toi);
    }

    const auto ccd = [&](double _min_distance, double& _toi) -> bool {
        return CTCD::vertexFaceCTCD(
            p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, _min_distance,
            _toi);
    };

    return ccd_strategy(
        ccd, min_distance, initial_distance, conservative_rescaling, toi);
}

} // namespace ipc

#endif