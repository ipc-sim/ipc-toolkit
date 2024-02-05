#pragma once
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/distance/distance_type.hpp>

#include <ipc/smooth_contact/primitives/point3.hpp>
#include <ipc/smooth_contact/primitives/face.hpp>
#include <ipc/smooth_contact/primitives/edge.hpp>

namespace ipc {
template <typename PrimitiveA, typename PrimitiveB>
struct PrimitiveDistType { };

template <> struct PrimitiveDistType<Point3, Point3> {
    using type = PointPointDistanceType;
};

template <> struct PrimitiveDistType<Edge3, Point3> {
    using type = PointEdgeDistanceType;
};

template <> struct PrimitiveDistType<Face, Point3> {
    using type = PointTriangleDistanceType;
};

template <> struct PrimitiveDistType<Edge3, Edge3> {
    using type = EdgeEdgeDistanceType;
};

template <typename PrimitiveA, typename PrimitiveB> class PrimitiveDistance {
    static_assert(
        PrimitiveA::dim == PrimitiveB::dim,
        "Primitives must have the same dimension");
    constexpr static int dim = PrimitiveA::dim;
    constexpr static int n_core_dofs =
        PrimitiveA::n_core_points * PrimitiveA::dim
        + PrimitiveB::n_core_points * PrimitiveB::dim;

public:
    static typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type
    compute_distance_type(const Vector<double, n_core_dofs>& x);

    static double compute_distance(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const long& a,
        const long& b,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);

    // points from primitiveA to primitiveB
    static Vector<double, dim> compute_closest_direction(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const long& a,
        const long& b,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
};

template <typename PrimitiveA, typename PrimitiveB, typename T>
class PrimitiveDistanceTemplate {
    static_assert(
        PrimitiveA::dim == PrimitiveB::dim,
        "Primitives must have the same dimension");
    constexpr static int dim = PrimitiveA::dim;
    constexpr static int n_core_dofs =
        PrimitiveA::n_core_points * PrimitiveA::dim
        + PrimitiveB::n_core_points * PrimitiveB::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype);
    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr);
};

} // namespace ipc