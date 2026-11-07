#include "distance_type_exact.hpp"

#include <ipc/utils/logger.hpp>

#include <Eigen/Geometry>

// fp_filters.h must come first: the PSM #defines FPG_UNCERTAIN_VALUE (to the
// same value), which would otherwise clobber the constant's declaration.
// clang-format off
#include "fp_filters.h"
#include "exact_vector.hpp"
// clang-format on

namespace ipc {

namespace {

    // The predicates below take raw coordinate pointers and live in this
    // translation unit, so a classification's whole chain of filtered
    // predicates inlines into its public entry point. The entry point converts
    // its arguments, checks DistanceTypeConfig, and initializes the exact
    // arithmetic once, instead of once per nested point-edge test. Only the
    // exact fallbacks, taken when a filter is inconclusive, stay out of line.

    using ExReal = ExactReal;
    using ExVec3 = ExactVector3;

    void init_predicates() { init_exact_arithmetic(); }

    /// @brief Exact copy of a 2D or 3D point (z = 0 in 2D).
    ExVec3 make_exact(const double* v, const int dim)
    {
        return ExVec3::from(v, dim);
    }

    int exact_sign(const ExReal& s) { return int(s.sign()); }

    // -- Exact fallbacks ------------------------------------------------------

    /// @brief Sign of dot(p1-p0, p2-p0) in exact arithmetic.
    int dot3_exact(
        const double* _p0, const double* _p1, const double* _p2, const int dim)
    {
        logger().trace("dot3 filter uncertain - fallback to exact arithmetic");
        const ExVec3 p0 = make_exact(_p0, dim);
        const ExVec3 p1 = make_exact(_p1, dim);
        const ExVec3 p2 = make_exact(_p2, dim);
        return exact_sign((p1 - p0).dot(p2 - p0));
    }

    /// @brief Sign of dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p0)) in exact
    /// arithmetic.
    int cross_dot_cross_1_exact(
        const double* _p0,
        const double* _p1,
        const double* _p2,
        const double* _p3)
    {
        logger().trace(
            "cross_dot_cross_1 filter uncertain - fallback to exact arithmetic");
        const ExVec3 p0 = make_exact(_p0, 3);
        const ExVec3 p1 = make_exact(_p1, 3);
        const ExVec3 p2 = make_exact(_p2, 3);
        const ExVec3 p3 = make_exact(_p3, 3);
        return exact_sign(
            (p1 - p0).cross(p2 - p0).dot((p3 - p0).cross(p1 - p0)));
    }

    /// @brief Sign of dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p2)) in exact
    /// arithmetic.
    int cross_dot_cross_2_exact(
        const double* _p0,
        const double* _p1,
        const double* _p2,
        const double* _p3)
    {
        logger().trace(
            "cross_dot_cross_2 filter uncertain - fallback to exact arithmetic");
        const ExVec3 p0 = make_exact(_p0, 3);
        const ExVec3 p1 = make_exact(_p1, 3);
        const ExVec3 p2 = make_exact(_p2, 3);
        const ExVec3 p3 = make_exact(_p3, 3);
        return exact_sign(
            (p1 - p0).cross(p2 - p0).dot((p3 - p0).cross(p1 - p2)));
    }

    // -- Filtered predicates --------------------------------------------------

    /// @brief Sign of dot(p1-p0, p2-p0).
    template <int dim>
    inline int dot3(const double* p0, const double* p1, const double* p2)
    {
        int s;
        if constexpr (dim == 2) {
            s = dot3_2d_filter(p0, p1, p2);
        } else {
            s = dot3_3d_filter(p0, p1, p2);
        }
        return s != FPG_UNCERTAIN_VALUE ? s : dot3_exact(p0, p1, p2, dim);
    }

    /// @brief Sign of dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p0))
    /// = dot(p1-p0, p3-p0) * dot(p1-p0, p2-p0)
    ///   - dot(p1-p0, p1-p0) * dot(p2-p0, p3-p0).
    inline int cross_dot_cross_1(
        const double* p0, const double* p1, const double* p2, const double* p3)
    {
        const int s = cross_dot_cross_1_3d_filter(p0, p1, p2, p3);
        return s != FPG_UNCERTAIN_VALUE
            ? s
            : cross_dot_cross_1_exact(p0, p1, p2, p3);
    }

    /// @brief Sign of dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p2))
    /// = dot(p1-p0, p3-p0) * dot(p1-p2, p2-p0)
    ///   - dot(p1-p0, p1-p2) * dot(p2-p0, p3-p0).
    inline int cross_dot_cross_2(
        const double* p0, const double* p1, const double* p2, const double* p3)
    {
        const int s = cross_dot_cross_2_3d_filter(p0, p1, p2, p3);
        return s != FPG_UNCERTAIN_VALUE
            ? s
            : cross_dot_cross_2_exact(p0, p1, p2, p3);
    }

    // -- Classifiers ----------------------------------------------------------

    template <int dim>
    inline PointEdgeDistanceType
    point_edge_type(const double* p, const double* e0, const double* e1)
    {
        if (dot3<dim>(e0, p, e1) <= 0) {
            return PointEdgeDistanceType::P_E0;
        } else if (dot3<dim>(e1, p, e0) <= 0) {
            return PointEdgeDistanceType::P_E1;
        } else {
            return PointEdgeDistanceType::P_E;
        }
    }

    inline PointTriangleDistanceType point_triangle_type(
        const double* p, const double* t0, const double* t1, const double* t2)
    {
        const int dot01 = dot3<3>(t0, p, t1);
        const int dot02 = dot3<3>(t0, p, t2);
        if (dot01 <= 0 && dot02 <= 0) {
            return PointTriangleDistanceType::P_T0;
        }
        const int dot12 = dot3<3>(t1, p, t2);
        const int dot10 = dot3<3>(t1, p, t0);
        if (dot12 <= 0 && dot10 <= 0) {
            return PointTriangleDistanceType::P_T1;
        }
        const int dot20 = dot3<3>(t2, p, t0);
        const int dot21 = dot3<3>(t2, p, t1);
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

    inline EdgeEdgeDistanceType edge_edge_type(
        const double* ea0,
        const double* ea1,
        const double* eb0,
        const double* eb1)
    {
        using PE = PointEdgeDistanceType;

        const PE dt_ea0 = point_edge_type<3>(ea0, eb0, eb1);
        const PE dt_ea1 = point_edge_type<3>(ea1, eb0, eb1);

        if (dt_ea0 == PE::P_E0 && dot3<3>(ea0, eb0, ea1) <= 0) {
            return EdgeEdgeDistanceType::EA0_EB0;
        }
        if (dt_ea0 == PE::P_E1 && dot3<3>(ea0, eb1, ea1) <= 0) {
            return EdgeEdgeDistanceType::EA0_EB1;
        }
        if (dt_ea1 == PE::P_E0 && dot3<3>(ea1, eb0, ea0) <= 0) {
            return EdgeEdgeDistanceType::EA1_EB0;
        }
        if (dt_ea1 == PE::P_E1 && dot3<3>(ea1, eb1, ea0) <= 0) {
            return EdgeEdgeDistanceType::EA1_EB1;
        }

        const PE dt_eb0 = point_edge_type<3>(eb0, ea0, ea1);
        const PE dt_eb1 = point_edge_type<3>(eb1, ea0, ea1);

        if (dt_eb0 == PE::P_E && cross_dot_cross_2(eb0, ea0, ea1, eb1) >= 0) {
            return EdgeEdgeDistanceType::EA_EB0;
        }
        if (dt_eb1 == PE::P_E && cross_dot_cross_2(eb1, ea0, ea1, eb0) >= 0) {
            return EdgeEdgeDistanceType::EA_EB1;
        }
        if (dt_ea0 == PE::P_E && cross_dot_cross_2(ea0, eb0, eb1, ea1) >= 0) {
            return EdgeEdgeDistanceType::EA0_EB;
        }
        if (dt_ea1 == PE::P_E && cross_dot_cross_2(ea1, eb0, eb1, ea0) >= 0) {
            return EdgeEdgeDistanceType::EA1_EB;
        }

        return EdgeEdgeDistanceType::EA_EB;
    }

} // namespace

PointEdgeDistanceType point_edge_distance_type_exact(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    if (!DistanceTypeConfig::instance().use_standard()) {
        assert(p.size() == e0.size() && p.size() == e1.size());
        init_predicates();
        return p.size() == 2
            ? point_edge_type<2>(p.data(), e0.data(), e1.data())
            : point_edge_type<3>(p.data(), e0.data(), e1.data());
    }
    return point_edge_distance_type(p, e0, e1);
}

PointTriangleDistanceType point_triangle_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2)
{
    if (!DistanceTypeConfig::instance().use_standard()) {
        init_predicates();
        return point_triangle_type(p.data(), t0.data(), t1.data(), t2.data());
    }
    return point_triangle_distance_type(p, t0, t1, t2);
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
    if constexpr (PARALLEL_THRESHOLD == 0.0) {
        init_predicates();
        // TODO use a zero filter?
        const int s = cross_null_3d_filter(
            _ea0.data(), _ea1.data(), _eb0.data(), _eb1.data());
        if (s != FPG_UNCERTAIN_VALUE) {
            return false;
        }
        const ExVec3 ea0 = make_exact(_ea0.data(), 3);
        const ExVec3 ea1 = make_exact(_ea1.data(), 3);
        const ExVec3 eb0 = make_exact(_eb0.data(), 3);
        const ExVec3 eb1 = make_exact(_eb1.data(), 3);
        return exact_sign((ea1 - ea0).cross(eb1 - eb0).squared_norm()) == 0;
    } else {
        return is_almost_parallel_edge_edge(_ea0, _ea1, _eb0, _eb1);
    }
}

EdgeEdgeDistanceType edge_edge_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    if (!DistanceTypeConfig::instance().use_standard()) {
        init_predicates();
        return edge_edge_type(ea0.data(), ea1.data(), eb0.data(), eb1.data());
    }
    return edge_edge_distance_type(ea0, ea1, eb0, eb1);
}

} // namespace ipc
