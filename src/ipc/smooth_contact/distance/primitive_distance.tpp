#pragma once

#include "primitive_distance.hpp"

namespace ipc {

template <typename T> class PrimitiveDistanceTemplate<Face, Point3, T> {
    static_assert(
        Face::DIM == Point3::DIM, "Primitives must have the same dimension");
    constexpr static int DIM = Face::DIM;
    constexpr static int N_CORE_DOFS =
        Face::N_CORE_POINTS * Face::DIM + Point3::N_CORE_POINTS * Point3::DIM;

public:
    static T compute_distance(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Face, Point3>::type dtype)
    {
        return point_triangle_sqr_distance<T>(
            x.template tail<3>() /* point */, x.template head<3>(),
            x.template segment<3>(3), x.template segment<3>(6) /* face */,
            dtype);
    }

    static Vector<T, DIM> compute_closest_direction(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Face, Point3>::type dtype)
    {
        return point_triangle_closest_point_direction<T>(
            x.template tail<3>() /* point */, x.template head<3>(),
            x.template segment<3>(3), x.template segment<3>(6) /* face */,
            dtype);
    }

    static T mollifier(const Vector<T, N_CORE_DOFS>& x, const T& dist_sqr)
    {
        return point_face_mollifier<T>(
            x.template tail<3>() /* point */, x.template head<3>(),
            x.template segment<3>(3), x.template segment<3>(6) /* face */,
            dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Edge3, T> {
    static_assert(
        Edge3::DIM == Edge3::DIM, "Primitives must have the same dimension");
    constexpr static int DIM = Edge3::DIM;
    constexpr static int N_CORE_DOFS =
        Edge3::N_CORE_POINTS * Edge3::DIM + Edge3::N_CORE_POINTS * Edge3::DIM;

public:
    static T compute_distance(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Edge3, Edge3>::type dtype)
    {
        return edge_edge_sqr_distance<T>(
            x.template head<3>() /* edge 0 */,
            x.template segment<3>(3) /* edge 0 */,
            x.template segment<3>(6) /* edge 1 */,
            x.template tail<3>() /* edge 1 */, dtype);
    }

    static Vector<T, DIM> compute_closest_direction(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Edge3, Edge3>::type dtype)
    {
        return edge_edge_closest_point_direction<T>(
            x.template head<3>() /* edge 0 */,
            x.template segment<3>(3) /* edge 0 */,
            x.template segment<3>(6) /* edge 1 */,
            x.template tail<3>() /* edge 1 */, dtype);
    }

    static T mollifier(const Vector<T, N_CORE_DOFS>& x, const T& dist_sqr)
    {
        std::array<HeavisideType, 4> types;
        types.fill(HeavisideType::VARIANT);
        return edge_edge_mollifier<T>(
            x.template head<3>() /* edge 0 */,
            x.template segment<3>(3) /* edge 0 */,
            x.template segment<3>(6) /* edge 1 */,
            x.template tail<3>() /* edge 1 */, types, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge2, Point2, T> {
    static_assert(
        Edge2::DIM == Point2::DIM, "Primitives must have the same dimension");
    constexpr static int DIM = Point2::DIM;
    constexpr static int N_CORE_DOFS =
        Edge2::N_CORE_POINTS * Edge2::DIM + Point2::N_CORE_POINTS * Point2::DIM;

public:
    static T compute_distance(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Edge2, Point2>::type dtype)
    {
        return PointEdgeDistance<T, DIM>::point_edge_sqr_distance(
            x.template tail<2>() /* point */, x.template head<2>() /* edge */,
            x.template segment<2>(2) /* edge */, dtype);
    }

    static Vector<T, DIM> compute_closest_direction(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Edge2, Point2>::type dtype)
    {
        return PointEdgeDistance<T, DIM>::point_edge_closest_point_direction(
            x.template tail<2>() /* point */, x.template head<2>() /* edge */,
            x.template segment<2>(2) /* edge */, dtype);
    }

    static T mollifier(const Vector<T, N_CORE_DOFS>& x, const T& dist_sqr)
    {
        return point_edge_mollifier<T, 2>(
            x.template tail<2>() /* point */,
            x.template segment<2>(2) /* edge */,
            x.template head<2>() /* edge */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Edge3, Point3, T> {
    static_assert(
        Edge3::DIM == Point3::DIM, "Primitives must have the same dimension");
    constexpr static int DIM = Edge3::DIM;
    constexpr static int N_CORE_DOFS =
        Edge3::N_CORE_POINTS * Edge3::DIM + Point3::N_CORE_POINTS * Point3::DIM;

public:
    static T compute_distance(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Edge3, Point3>::type dtype)
    {
        return PointEdgeDistance<T, DIM>::point_edge_sqr_distance(
            x.template tail<3>() /* point */, x.template head<3>() /* edge */,
            x.template segment<3>(3) /* edge */, dtype);
    }

    static Vector<T, DIM> compute_closest_direction(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Edge3, Point3>::type dtype)
    {
        return PointEdgeDistance<T, DIM>::point_edge_closest_point_direction(
            x.template tail<3>() /* point */, x.template head<3>() /* edge */,
            x.template segment<3>(3) /* edge */, dtype);
    }

    static T mollifier(const Vector<T, N_CORE_DOFS>& x, const T& dist_sqr)
    {
        return point_edge_mollifier<T, 3>(
            x.template tail<3>() /* point */,
            x.template segment<3>(3) /* edge */,
            x.template head<3>() /* edge */, dist_sqr);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Point2, Point2, T> {
    static_assert(
        Point2::DIM == Point2::DIM, "Primitives must have the same dimension");
    constexpr static int DIM = Point2::DIM;
    constexpr static int N_CORE_DOFS = Point2::N_CORE_POINTS * Point2::DIM
        + Point2::N_CORE_POINTS * Point2::DIM;

public:
    static T compute_distance(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Point2, Point2>::type dtype)
    {
        return (x.template tail<2>() - x.template head<2>()).squaredNorm();
    }

    static Vector<T, DIM> compute_closest_direction(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Point2, Point2>::type dtype)
    {
        return x.template tail<2>() - x.template head<2>();
    }

    static T mollifier(const Vector<T, N_CORE_DOFS>& x, const T& dist_sqr)
    {
        return T(1.);
    }
};

template <typename T> class PrimitiveDistanceTemplate<Point3, Point3, T> {
    static_assert(
        Point3::DIM == Point3::DIM, "Primitives must have the same dimension");
    constexpr static int DIM = Point3::DIM;
    constexpr static int N_CORE_DOFS = Point3::N_CORE_POINTS * Point3::DIM
        + Point3::N_CORE_POINTS * Point3::DIM;

public:
    static T compute_distance(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Point3, Point3>::type dtype)
    {
        return (x.template tail<3>() - x.template head<3>()).squaredNorm();
    }

    static Vector<T, DIM> compute_closest_direction(
        const Vector<T, N_CORE_DOFS>& x,
        typename PrimitiveDistType<Point3, Point3>::type dtype)
    {
        return x.template tail<3>() - x.template head<3>();
    }

    static T mollifier(const Vector<T, N_CORE_DOFS>& x, const T& dist_sqr)
    {
        return T(1.);
    }
};

} // namespace ipc
