#include "primitive_distance.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <ipc/smooth_contact/distance/edge_edge.hpp>
#include <ipc/smooth_contact/distance/point_edge.hpp>
#include <ipc/smooth_contact/distance/point_face.hpp>

namespace ipc {

// template <typename PrimitiveA, typename PrimitiveB>
// std::tuple<
//     Vector<double, PrimitiveDistance<PrimitiveA, PrimitiveB>::dim>, 
//     Eigen::Matrix<double, PrimitiveDistance<PrimitiveA, PrimitiveB>::dim, PrimitiveDistance<PrimitiveA, PrimitiveB>::n_core_dofs>>
// PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction_gradient(
//     const Vector<double, PrimitiveDistance<PrimitiveA, PrimitiveB>::n_core_dofs>& x,
//     typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
// {
//     DiffScalarBase::setVariableCount(n_core_dofs);
//     using T = ADGrad<n_core_dofs>;
//     const Vector<T, n_core_dofs> X = slice_positions<T, n_core_dofs, 1>(x);
//     const Vector<T, dim> d = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(X, dtype);

//     Vector<double, dim> out;
//     Eigen::Matrix<double, dim, n_core_dofs> J = Eigen::Matrix<double, dim, n_core_dofs>::Zero();
//     for (int i = 0; i < dim; i++)
//     {
//         out(i) = d(i).getValue();
//         J.row(i) = d.getGradient();
//     }
//     return std::make_tuple(out, J);
// }

// template <typename PrimitiveA, typename PrimitiveB>
// std::tuple<
//     Vector<double, PrimitiveDistance<PrimitiveA, PrimitiveB>::dim>, 
//     Eigen::Matrix<double, PrimitiveDistance<PrimitiveA, PrimitiveB>::dim, PrimitiveDistance<PrimitiveA, PrimitiveB>::n_core_dofs>,
//     std::array<Eigen::Matrix<double, PrimitiveDistance<PrimitiveA, PrimitiveB>::n_core_dofs, PrimitiveDistance<PrimitiveA, PrimitiveB>::n_core_dofs>, PrimitiveDistance<PrimitiveA, PrimitiveB>::dim>>
// PrimitiveDistance<PrimitiveA, PrimitiveB>::compute_closest_direction_hessian(
//     const Vector<double, PrimitiveDistance<PrimitiveA, PrimitiveB>::n_core_dofs>& x,
//     typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type dtype)
// {
//     DiffScalarBase::setVariableCount(n_core_dofs);
//     using T = ADHessian<n_core_dofs>;
//     const Vector<T, n_core_dofs> X = slice_positions<T, n_core_dofs, 1>(x);
//     const Vector<T, dim> d = PrimitiveDistanceTemplate<PrimitiveA, PrimitiveB, T>::compute_closest_direction(X, dtype);

//     Vector<double, dim> out;
//     Eigen::Matrix<double, dim, n_core_dofs> J = Eigen::Matrix<double, dim, n_core_dofs>::Zero();
//     std::array<Eigen::Matrix<double, n_core_dofs, n_core_dofs>, dim> H;
//     for (int i = 0; i < dim; i++)
//     {
//         out(i) = d(i).getValue();
//         J.row(i) = d.getGradient();
//         H[i] = d.getHessian();
//     }
//     return std::make_tuple(out, J, H);
// }

template <>
typename PrimitiveDistType<Face, Point3>::type
PrimitiveDistance<Face, Point3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return point_triangle_distance_type(
        x.tail(3) /* point */, x.head(3), x.segment(3, 3),
        x.segment(6, 3) /* face */);
}

template <>
typename PrimitiveDistType<Edge3, Edge3>::type
PrimitiveDistance<Edge3, Edge3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return edge_edge_distance_type(
        x.head(3) /* edge 0 */, x.segment(3, 3) /* edge 0 */,
        x.segment(6, 3) /* edge 1 */, x.tail(3) /* edge 1 */);
}

template <>
typename PrimitiveDistType<Edge3, Point3>::type
PrimitiveDistance<Edge3, Point3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return point_edge_distance_type(
        x.tail(3) /* point */, x.head(3) /* edge */,
        x.segment(3, 3) /* edge */);
}

template <>
typename PrimitiveDistType<Point3, Point3>::type
PrimitiveDistance<Point3, Point3>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return PointPointDistanceType::AUTO;
}

template <>
double PrimitiveDistance<Face, Point3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Face, Point3>::type dtype)
{
    return point_triangle_distance(
        V.row(b), V.row(mesh.faces()(a, 0)), V.row(mesh.faces()(a, 1)),
        V.row(mesh.faces()(a, 2)), dtype);
}

template <>
double PrimitiveDistance<Edge3, Edge3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    return edge_edge_distance(
        V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)),
        V.row(mesh.edges()(b, 0)), V.row(mesh.edges()(b, 1)), dtype);
}

template <>
double PrimitiveDistance<Edge3, Point3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Point3>::type dtype)
{
    return point_edge_distance(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
double PrimitiveDistance<Point3, Point3>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    return point_point_distance(V.row(a), V.row(b));
}

template <>
Vector<double, PrimitiveDistance<Face, Point3>::dim>
PrimitiveDistance<Face, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Face, Point3>::type dtype)
{
    return point_triangle_closest_point_direction<double>(
        V.row(b), V.row(mesh.faces()(a, 0)), V.row(mesh.faces()(a, 1)),
        V.row(mesh.faces()(a, 2)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Edge3, Edge3>::dim>
PrimitiveDistance<Edge3, Edge3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    return edge_edge_closest_point_direction<double>(
        V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)),
        V.row(mesh.edges()(b, 0)), V.row(mesh.edges()(b, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Edge3, Point3>::dim>
PrimitiveDistance<Edge3, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge3, Point3>::type dtype)
{
    return PointEdgeDistance<double, 3>::point_edge_closest_point_direction(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Point3, Point3>::dim>
PrimitiveDistance<Point3, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    return V.row(b) - V.row(a);
}

#ifndef DERIVATIVES_WITH_AUTODIFF

template <>
std::tuple<Vector<double, PrimitiveDistance<Edge3, Edge3>::dim>, Eigen::Matrix<double, PrimitiveDistance<Edge3, Edge3>::dim, PrimitiveDistance<Edge3, Edge3>::n_core_dofs>,
    std::array<Eigen::Matrix<double, PrimitiveDistance<Edge3, Edge3>::n_core_dofs, PrimitiveDistance<Edge3, Edge3>::n_core_dofs>, PrimitiveDistance<Edge3, Edge3>::dim>>
PrimitiveDistance<Edge3, Edge3>::compute_closest_direction_hessian(
    const Vector<double, n_core_dofs>& x,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    assert(dtype == EdgeEdgeDistanceType::EA_EB);
    return line_line_closest_point_direction_hessian(
        x.head<3>() /* edge 0 */, x.segment<3>(3) /* edge 0 */,
        x.segment<3>(6) /* edge 1 */, x.tail<3>() /* edge 1 */);
}

template <>
std::tuple<Vector<double, PrimitiveDistance<Point3, Point3>::dim>, Eigen::Matrix<double, PrimitiveDistance<Point3, Point3>::dim, PrimitiveDistance<Point3, Point3>::n_core_dofs>,
    std::array<Eigen::Matrix<double, PrimitiveDistance<Point3, Point3>::n_core_dofs, PrimitiveDistance<Point3, Point3>::n_core_dofs>, PrimitiveDistance<Point3, Point3>::dim>>
PrimitiveDistance<Point3, Point3>::compute_closest_direction_hessian(
    const Vector<double, n_core_dofs>& x,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    Vector3d out = x.tail(3) - x.head(3);
    Eigen::Matrix<double, 3, 6> J = Eigen::Matrix<double, 3, 6>::Zero();
    J.leftCols<3>().diagonal().array() = -1;
    J.rightCols<3>().diagonal().array() = 1;
    std::array<Eigen::Matrix<double, 6, 6>, 3> H;
    for (auto &h : H)
        h.setZero();
    return std::make_tuple(out, J, H);
}

#endif

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

template class PrimitiveDistance<Edge3, Point3>;
template class PrimitiveDistance<Edge3, Edge3>;
template class PrimitiveDistance<Point3, Point3>;
template class PrimitiveDistance<Face, Point3>;

template class PrimitiveDistanceTemplate<Edge3, Point3, double>;
template class PrimitiveDistanceTemplate<Edge3, Edge3, double>;
template class PrimitiveDistanceTemplate<Point3, Point3, double>;
template class PrimitiveDistanceTemplate<Face, Point3, double>;

template class PrimitiveDistanceTemplate<Edge3, Point3, ADGrad<9>>;
template class PrimitiveDistanceTemplate<Edge3, Edge3, ADGrad<12>>;
template class PrimitiveDistanceTemplate<Point3, Point3, ADGrad<6>>;
template class PrimitiveDistanceTemplate<Face, Point3, ADGrad<12>>;

template class PrimitiveDistanceTemplate<Edge3, Point3, ADHessian<9>>;
template class PrimitiveDistanceTemplate<Edge3, Edge3, ADHessian<12>>;
template class PrimitiveDistanceTemplate<Point3, Point3, ADHessian<6>>;
template class PrimitiveDistanceTemplate<Face, Point3, ADHessian<12>>;
} // namespace ipc
