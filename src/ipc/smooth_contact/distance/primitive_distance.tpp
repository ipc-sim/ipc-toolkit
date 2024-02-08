#pragma once

#include "primitive_distance.hpp"

namespace ipc {

template <typename T> class PrimitiveDistanceTemplate<Face, Point3, T> {
    static_assert(
        Face::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Face::dim;
    constexpr static int n_core_dofs =
        Face::n_core_points * Face::dim + Point3::n_core_points * Point3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Face, Point3>::type dtype)
    {
        return point_triangle_closest_point_direction<T>(
            x.tail(3) /* point */, x.head(3), x.segment(3, 3),
            x.segment(6, 3) /* face */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return point_face_mollifier<T>(
            x.tail(3) /* point */, x.head(3), x.segment(3, 3),
            x.segment(6, 3) /* face */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Edge3, T> {
    static_assert(
        Edge3::dim == Edge3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Edge3::dim;
    constexpr static int n_core_dofs =
        Edge3::n_core_points * Edge3::dim + Edge3::n_core_points * Edge3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Edge3>::type dtype)
    {
        return edge_edge_closest_point_direction<T>(
            x.head(3) /* edge 0 */, x.segment(3, 3) /* edge 0 */,
            x.segment(6, 3) /* edge 1 */, x.tail(3) /* edge 1 */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        std::array<HEAVISIDE_TYPE, 4> types;
        types.fill(HEAVISIDE_TYPE::VARIANT);
        return edge_edge_mollifier<T>(
            x.head(3) /* edge 0 */, x.segment(3, 3) /* edge 0 */,
            x.segment(6, 3) /* edge 1 */, x.tail(3) /* edge 1 */, types,
            dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Point3, T> {
    static_assert(
        Edge3::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Edge3::dim;
    constexpr static int n_core_dofs =
        Edge3::n_core_points * Edge3::dim + Point3::n_core_points * Point3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Point3>::type dtype)
    {
        return PointEdgeDistance<T, dim>::point_edge_closest_point_direction(
            x.tail(3) /* point */, x.head(3) /* edge */,
            x.segment(3, 3) /* edge */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return point_edge_mollifier<T>(
            x.tail(3) /* point */, x.segment(3, 3) /* edge */,
            x.head(3) /* edge */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Point3, Point3, T> {
    static_assert(
        Point3::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Point3::dim;
    constexpr static int n_core_dofs = Point3::n_core_points * Point3::dim
        + Point3::n_core_points * Point3::dim;

public:
    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Point3, Point3>::type dtype)
    {
        return x.tail(3) - x.head(3);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return T(1.);
    }
};

}
