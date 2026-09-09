#include "distance_type_exact.hpp"

#include <ipc/utils/logger.hpp>

#include <Eigen/Geometry>

#ifdef IPC_TOOLKIT_WITH_GEOGRAM
#include "fp_filters.h"

#include <geogram/numerics/exact_geometry.h>
// geogram 1.10 no longer pulls this in transitively via exact_geometry.h.
#include <geogram/numerics/predicates.h>
#endif

namespace ipc {

#ifdef IPC_TOOLKIT_WITH_GEOGRAM
using ExReal = GEO::expansion_nt; // exact scalar type
using ExVec3 = GEO::vec3E;        // exact vector type

inline void init_pck()
{
    struct PckInit {
        PckInit() { GEO::PCK::initialize(); }
    };
    static PckInit _;
}

inline ExVec3 make_exact(Eigen::ConstRef<VectorMax3d> v)
{
    ExReal x { v.x() };
    ExReal y { v.y() };
    ExReal z { v.size() < 3 ? 0 : v.z() }; // compatibility with 2D vectors
    return ExVec3(std::move(x), std::move(y), std::move(z));
}

int dot3_3d(
    Eigen::ConstRef<VectorMax3d> _p0,
    Eigen::ConstRef<VectorMax3d> _p1,
    Eigen::ConstRef<VectorMax3d> _p2)
{
    // Evaluates the sign of dot(p1-p0, p2-p0)
    const int s = dot3_3d_filter(_p0.data(), _p1.data(), _p2.data());
    if (s != FPG_UNCERTAIN_VALUE) {
        return s;
    }
    logger().trace("dot3_3d filter uncertain - fallback to exact arithmetic");
    const ExVec3 p0 = make_exact(_p0);
    const ExVec3 p1 = make_exact(_p1);
    const ExVec3 p2 = make_exact(_p2);
    const ExReal ss = dot(p1 - p0, p2 - p0);
    return (ss > 0) ? 1 : ((ss < 0) ? -1 : 0);
}

int dot3_2d(
    Eigen::ConstRef<VectorMax3d> _p0,
    Eigen::ConstRef<VectorMax3d> _p1,
    Eigen::ConstRef<VectorMax3d> _p2)
{
    // Evaluates the sign of dot(p1-p0, p2-p0)
    const int s = dot3_2d_filter(_p0.data(), _p1.data(), _p2.data());
    if (s != FPG_UNCERTAIN_VALUE) {
        return s;
    }
    logger().trace("dot3_2d filter uncertain - fallback to exact arithmetic");
    const ExVec3 p0 = make_exact(_p0);
    const ExVec3 p1 = make_exact(_p1);
    const ExVec3 p2 = make_exact(_p2);
    const ExReal ss = dot(p1 - p0, p2 - p0);
    return (ss > 0) ? 1 : ((ss < 0) ? -1 : 0);
}

int cross_dot_cross_1(
    Eigen::ConstRef<VectorMax3d> _p0,
    Eigen::ConstRef<VectorMax3d> _p1,
    Eigen::ConstRef<VectorMax3d> _p2,
    Eigen::ConstRef<VectorMax3d> _p3)
{
    /*
    Evaluates the sign of dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p0)) =
    = dot(p1-p0, p3-p0) * dot(p1-p0, p2-p0) - dot(p1-p0, p1-p0) * dot(p2-p0,
    p3-p0)
    */
    const int s = cross_dot_cross_1_3d_filter(
        _p0.data(), _p1.data(), _p2.data(), _p3.data());
    if (s != FPG_UNCERTAIN_VALUE) {
        return s;
    }
    logger().trace(
        "cross_dot_cross_1 filter uncertain - fallback to exact arithmetic");
    const ExVec3 p0 = make_exact(_p0);
    const ExVec3 p1 = make_exact(_p1);
    const ExVec3 p2 = make_exact(_p2);
    const ExVec3 p3 = make_exact(_p3);
    const ExReal ss = dot(cross(p1 - p0, p2 - p0), cross(p3 - p0, p1 - p0));
    return (ss > 0) ? 1 : ((ss < 0) ? -1 : 0);
}

int cross_dot_cross_2(
    Eigen::ConstRef<VectorMax3d> _p0,
    Eigen::ConstRef<VectorMax3d> _p1,
    Eigen::ConstRef<VectorMax3d> _p2,
    Eigen::ConstRef<VectorMax3d> _p3)
{
    /*
    Evaluates the sign of dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p2)) =
    = dot(p1-p0, p3-p0) * dot(p1-p2, p2-p0) - dot(p1-p0, p1-p2) * dot(p2-p0,
    p3-p0)
    */
    const int s = cross_dot_cross_2_3d_filter(
        _p0.data(), _p1.data(), _p2.data(), _p3.data());
    if (s != FPG_UNCERTAIN_VALUE) {
        return s;
    }
    logger().trace(
        "cross_dot_cross_1 filter uncertain - fallback to exact arithmetic");
    const ExVec3 p0 = make_exact(_p0);
    const ExVec3 p1 = make_exact(_p1);
    const ExVec3 p2 = make_exact(_p2);
    const ExVec3 p3 = make_exact(_p3);
    const ExReal ss = dot(cross(p1 - p0, p2 - p0), cross(p3 - p0, p1 - p2));
    return (ss > 0) ? 1 : ((ss < 0) ? -1 : 0);
}

static PointEdgeDistanceType point_edge_distance_type_predicate(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    init_pck();
    assert(p.size() == e0.size() && p.size() == e1.size());
    if (p.size() == 2) {
        if (dot3_2d(e0, p, e1) <= 0) {
            return PointEdgeDistanceType::P_E0;
        } else if (dot3_2d(e1, p, e0) <= 0) {
            return PointEdgeDistanceType::P_E1;
        } else {
            return PointEdgeDistanceType::P_E;
        }
    } else {
        if (dot3_3d(e0, p, e1) <= 0) {
            return PointEdgeDistanceType::P_E0;
        } else if (dot3_3d(e1, p, e0) <= 0) {
            return PointEdgeDistanceType::P_E1;
        } else {
            return PointEdgeDistanceType::P_E;
        }
    }
}
#endif // IPC_TOOLKIT_WITH_GEOGRAM

PointEdgeDistanceType point_edge_distance_type_exact(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
#ifdef IPC_TOOLKIT_WITH_GEOGRAM
    if (DistanceTypeConfig::instance().use_standard()) {
        return point_edge_distance_type(p, e0, e1);
    }
    return point_edge_distance_type_predicate(p, e0, e1);
#else
    return point_edge_distance_type(p, e0, e1);
#endif
}

#ifdef IPC_TOOLKIT_WITH_GEOGRAM
static PointTriangleDistanceType point_triangle_distance_type_predicate(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
    init_pck();
    const int dot01 = dot3_3d(t0, p, t1);
    const int dot02 = dot3_3d(t0, p, t2);
    if (dot01 <= 0 && dot02 <= 0) {
        return PointTriangleDistanceType::P_T0;
    }
    const int dot12 = dot3_3d(t1, p, t2);
    const int dot10 = dot3_3d(t1, p, t0);
    if (dot12 <= 0 && dot10 <= 0) {
        return PointTriangleDistanceType::P_T1;
    }
    const int dot20 = dot3_3d(t2, p, t0);
    const int dot21 = dot3_3d(t2, p, t1);
    if (dot20 <= 0 && dot21 <= 0) {
        return PointTriangleDistanceType::P_T2;
    }

    if (cross_dot_cross_1(t0, t1, t2, p) >= 0 && dot01 > 0 && dot10 > 0) {
        return PointTriangleDistanceType::P_E0;
    }
    if (cross_dot_cross_1(t1, t2, t0, p) >= 0 && dot12 > 0 && dot21 > 0) {
        return PointTriangleDistanceType::P_E1;
    }
    if (cross_dot_cross_1(t2, t0, t1, p) >= 0 && dot20 > 0 && dot02 > 0) {
        return PointTriangleDistanceType::P_E2;
    }

    return PointTriangleDistanceType::P_T;
}
#endif // IPC_TOOLKIT_WITH_GEOGRAM

PointTriangleDistanceType point_triangle_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
#ifdef IPC_TOOLKIT_WITH_GEOGRAM
    if (DistanceTypeConfig::instance().use_standard()) {
        return point_triangle_distance_type(p, t0, t1, t2);
    }
    return point_triangle_distance_type_predicate(p, t0, t1, t2);
#else
    return point_triangle_distance_type(p, t0, t1, t2);
#endif
}

bool is_almost_parallel_edge_edge(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d u = ea1 - ea0;
    const Eigen::Vector3d v = eb1 - eb0;
    const double cross_norm_sqr = u.cross(v).squaredNorm();
    const double a = u.squaredNorm();
    const double c = v.squaredNorm();
    // Relative sin² test: parallel when sin²(θ) < PARALLEL_THRESHOLD.
    // Scaling by a*c (rather than max(1, a*c)) keeps this scale-invariant.
    return cross_norm_sqr < a * c * PARALLEL_THRESHOLD;
}

bool is_parallel_edge_edge(
    Eigen::ConstRef<Eigen::Vector3d> _ea0,
    Eigen::ConstRef<Eigen::Vector3d> _ea1,
    Eigen::ConstRef<Eigen::Vector3d> _eb0,
    Eigen::ConstRef<Eigen::Vector3d> _eb1)
{
#ifdef IPC_TOOLKIT_WITH_GEOGRAM
    if constexpr (PARALLEL_THRESHOLD == 0.0) {
        init_pck();
        // TODO use a zero filter?
        const int s = cross_null_3d_filter(
            _ea0.data(), _ea1.data(), _eb0.data(), _eb1.data());
        if (s != FPG_UNCERTAIN_VALUE) {
            return false;
        }
        const ExVec3 ea0 = make_exact(_ea0);
        const ExVec3 ea1 = make_exact(_ea1);
        const ExVec3 eb0 = make_exact(_eb0);
        const ExVec3 eb1 = make_exact(_eb1);
        const ExReal cross_norm_sqr = cross(ea1 - ea0, eb1 - eb0).length2();
        return cross_norm_sqr == 0;
    } else {
        return is_almost_parallel_edge_edge(_ea0, _ea1, _eb0, _eb1);
    }
#else
    // Without geogram the exact test is unavailable; PARALLEL_THRESHOLD must
    // be non-zero for the thresholded test to be meaningful.
    static_assert(
        PARALLEL_THRESHOLD != 0.0,
        "PARALLEL_THRESHOLD == 0 requires the exact predicates (geogram).");
    return is_almost_parallel_edge_edge(_ea0, _ea1, _eb0, _eb1);
#endif
}

#ifdef IPC_TOOLKIT_WITH_GEOGRAM
static EdgeEdgeDistanceType edge_edge_distance_type_predicate(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    init_pck();

    const PointEdgeDistanceType dt_ea0 =
        point_edge_distance_type_exact(ea0, eb0, eb1);
    const PointEdgeDistanceType dt_ea1 =
        point_edge_distance_type_exact(ea1, eb0, eb1);

    if (dt_ea0 == PointEdgeDistanceType::P_E0 && dot3_3d(ea0, eb0, ea1) <= 0) {
        return EdgeEdgeDistanceType::EA0_EB0;
    }
    if (dt_ea0 == PointEdgeDistanceType::P_E1 && dot3_3d(ea0, eb1, ea1) <= 0) {
        return EdgeEdgeDistanceType::EA0_EB1;
    }
    if (dt_ea1 == PointEdgeDistanceType::P_E0 && dot3_3d(ea1, eb0, ea0) <= 0) {
        return EdgeEdgeDistanceType::EA1_EB0;
    }
    if (dt_ea1 == PointEdgeDistanceType::P_E1 && dot3_3d(ea1, eb1, ea0) <= 0) {
        return EdgeEdgeDistanceType::EA1_EB1;
    }

    const PointEdgeDistanceType dt_eb0 =
        point_edge_distance_type_exact(eb0, ea0, ea1);
    const PointEdgeDistanceType dt_eb1 =
        point_edge_distance_type_exact(eb1, ea0, ea1);

    if (dt_eb0 == PointEdgeDistanceType::P_E
        && cross_dot_cross_2(eb0, ea0, ea1, eb1) >= 0) {
        return EdgeEdgeDistanceType::EA_EB0;
    }
    if (dt_eb1 == PointEdgeDistanceType::P_E
        && cross_dot_cross_2(eb1, ea0, ea1, eb0) >= 0) {
        return EdgeEdgeDistanceType::EA_EB1;
    }
    if (dt_ea0 == PointEdgeDistanceType::P_E
        && cross_dot_cross_2(ea0, eb0, eb1, ea1) >= 0) {
        return EdgeEdgeDistanceType::EA0_EB;
    }
    if (dt_ea1 == PointEdgeDistanceType::P_E
        && cross_dot_cross_2(ea1, eb0, eb1, ea0) >= 0) {
        return EdgeEdgeDistanceType::EA1_EB;
    }

    return EdgeEdgeDistanceType::EA_EB;
}
#endif // IPC_TOOLKIT_WITH_GEOGRAM

EdgeEdgeDistanceType edge_edge_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
#ifdef IPC_TOOLKIT_WITH_GEOGRAM
    if (DistanceTypeConfig::instance().use_standard()) {
        return edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }
    return edge_edge_distance_type_predicate(ea0, ea1, eb0, eb1);
#else
    return edge_edge_distance_type(ea0, ea1, eb0, eb1);
#endif
}

} // namespace ipc
