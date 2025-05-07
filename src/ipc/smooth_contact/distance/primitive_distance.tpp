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
    static T compute_distance(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Face, Point3>::type dtype)
    {
        return point_triangle_sqr_distance<T>(
            x.template tail<3>() /* point */, x.template head<3>(), x.template segment<3>(3),
            x.template segment<3>(6) /* face */, dtype);
    }

    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Face, Point3>::type dtype)
    {
        return point_triangle_closest_point_direction<T>(
            x.template tail<3>() /* point */, x.template head<3>(), x.template segment<3>(3),
            x.template segment<3>(6) /* face */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return point_face_mollifier<T>(
            x.template tail<3>() /* point */, x.template head<3>(), x.template segment<3>(3),
            x.template segment<3>(6) /* face */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Edge3, T> {
    static_assert(
        Edge3::dim == Edge3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Edge3::dim;
    constexpr static int n_core_dofs =
        Edge3::n_core_points * Edge3::dim + Edge3::n_core_points * Edge3::dim;

public:
    static T compute_distance(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Edge3>::type dtype)
    {
        return edge_edge_sqr_distance<T>(
            x.template head<3>() /* edge 0 */, x.template segment<3>(3) /* edge 0 */,
            x.template segment<3>(6) /* edge 1 */, x.template tail<3>() /* edge 1 */, dtype);
    }

    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Edge3>::type dtype)
    {
        return edge_edge_closest_point_direction<T>(
            x.template head<3>() /* edge 0 */, x.template segment<3>(3) /* edge 0 */,
            x.template segment<3>(6) /* edge 1 */, x.template tail<3>() /* edge 1 */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        std::array<HEAVISIDE_TYPE, 4> types;
        types.fill(HEAVISIDE_TYPE::VARIANT);
        return edge_edge_mollifier<T>(
            x.template head<3>() /* edge 0 */, x.template segment<3>(3) /* edge 0 */,
            x.template segment<3>(6) /* edge 1 */, x.template tail<3>() /* edge 1 */, types,
            dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge2, Point2, T> {
    static_assert(
        Edge2::dim == Point2::dim, "Primitives must have the same dimension");
    constexpr static int dim = Point2::dim;
    constexpr static int n_core_dofs = Edge2::n_core_points * Edge2::dim
        + Point2::n_core_points * Point2::dim;

public:
    static T compute_distance(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge2, Point2>::type dtype)
    {
        return PointEdgeDistance<T, dim>::point_edge_sqr_distance(
            x.template tail<2>() /* point */, x.template head<2>() /* edge */,
            x.template segment<2>(2) /* edge */, dtype);
    }

    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge2, Point2>::type dtype)
    {
        return PointEdgeDistance<T, dim>::point_edge_closest_point_direction(
            x.template tail<2>() /* point */, x.template head<2>() /* edge */,
            x.template segment<2>(2) /* edge */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return point_edge_mollifier<T, 2>(
            x.template tail<2>() /* point */, x.template segment<2>(2) /* edge */,
            x.template head<2>() /* edge */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Point3, T> {
    static_assert(
        Edge3::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Edge3::dim;
    constexpr static int n_core_dofs =
        Edge3::n_core_points * Edge3::dim + Point3::n_core_points * Point3::dim;

public:
    static T compute_distance(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Point3>::type dtype)
    {
        return PointEdgeDistance<T, dim>::point_edge_sqr_distance(
            x.template tail<3>() /* point */, x.template head<3>() /* edge */,
            x.template segment<3>(3) /* edge */, dtype);
    }

    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Edge3, Point3>::type dtype)
    {
        return PointEdgeDistance<T, dim>::point_edge_closest_point_direction(
            x.template tail<3>() /* point */, x.template head<3>() /* edge */,
            x.template segment<3>(3) /* edge */, dtype);
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return point_edge_mollifier<T, 3>(
            x.template tail<3>() /* point */, x.template segment<3>(3) /* edge */,
            x.template head<3>() /* edge */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Point2, Point2, T> {
    static_assert(
        Point2::dim == Point2::dim, "Primitives must have the same dimension");
    constexpr static int dim = Point2::dim;
    constexpr static int n_core_dofs = Point2::n_core_points * Point2::dim
        + Point2::n_core_points * Point2::dim;

public:
    static T compute_distance(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Point2, Point2>::type dtype)
    {
        return (x.template tail<2>() - x.template head<2>()).squaredNorm();
    }

    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Point2, Point2>::type dtype)
    {
        return x.template tail<2>() - x.template head<2>();
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return T(1.);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Point3, Point3, T> {
    static_assert(
        Point3::dim == Point3::dim, "Primitives must have the same dimension");
    constexpr static int dim = Point3::dim;
    constexpr static int n_core_dofs = Point3::n_core_points * Point3::dim
        + Point3::n_core_points * Point3::dim;

public:
    static T compute_distance(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Point3, Point3>::type dtype)
    {
        return (x.template tail<3>() - x.template head<3>()).squaredNorm();
    }

    static Vector<T, dim> compute_closest_direction(
        const Vector<T, n_core_dofs>& x,
        typename PrimitiveDistType<Point3, Point3>::type dtype)
    {
        return x.template tail<3>() - x.template head<3>();
    }

    static T mollifier(const Vector<T, n_core_dofs>& x, const T& dist_sqr)
    {
        return T(1.);
    }
};

}
