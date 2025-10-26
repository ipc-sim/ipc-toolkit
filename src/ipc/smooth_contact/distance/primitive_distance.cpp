#include "primitive_distance.hpp"

#include <ipc/config.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/smooth_contact/distance/edge_edge.hpp>
#include <ipc/smooth_contact/distance/point_edge.hpp>
#include <ipc/smooth_contact/distance/point_face.hpp>

namespace ipc {

template <>
typename PrimitiveDistType<Face, Point3>::type
PrimitiveDistance<Face, Point3>::compute_distance_type(
    const Vector<double, N_CORE_DOFS>& x)
{
    return point_triangle_distance_type(
        x.tail(3) /* point */, x.head(3), x.segment(3, 3),
        x.segment(6, 3) /* face */);
}

template <>
typename PrimitiveDistType<Edge3, Edge3>::type
PrimitiveDistance<Edge3, Edge3>::compute_distance_type(
    const Vector<double, N_CORE_DOFS>& x)
{
    return edge_edge_distance_type(
        x.head(3) /* edge 0 */, x.segment(3, 3) /* edge 0 */,
        x.segment(6, 3) /* edge 1 */, x.tail(3) /* edge 1 */);
}

template <>
typename PrimitiveDistType<Edge3, Point3>::type
PrimitiveDistance<Edge3, Point3>::compute_distance_type(
    const Vector<double, N_CORE_DOFS>& x)
{
    return point_edge_distance_type(
        x.tail(3) /* point */, x.head(3) /* edge */,
        x.segment(3, 3) /* edge */);
}

template <>
typename PrimitiveDistType<Edge2, Point2>::type
PrimitiveDistance<Edge2, Point2>::compute_distance_type(
    const Vector<double, N_CORE_DOFS>& x)
{
    return PointEdgeDistanceType::AUTO;
}

template <>
typename PrimitiveDistType<Point2, Point2>::type
PrimitiveDistance<Point2, Point2>::compute_distance_type(
    const Vector<double, N_CORE_DOFS>& x)
{
    return PointPointDistanceType::P_P;
}

template <>
typename PrimitiveDistType<Point3, Point3>::type
PrimitiveDistance<Point3, Point3>::compute_distance_type(
    const Vector<double, N_CORE_DOFS>& x)
{
    return PointPointDistanceType::P_P;
}

template <>
Vector<double, PrimitiveDistance<Face, Point3>::DIM>
PrimitiveDistance<Face, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const index_t a,
    const index_t b,
    typename PrimitiveDistType<Face, Point3>::type dtype)
{
    return point_triangle_closest_point_direction<double>(
        V.row(b), V.row(mesh.faces()(a, 0)), V.row(mesh.faces()(a, 1)),
        V.row(mesh.faces()(a, 2)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Edge3, Edge3>::DIM>
PrimitiveDistance<Edge3, Edge3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const index_t a,
    const index_t b,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    return edge_edge_closest_point_direction<double>(
        V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)),
        V.row(mesh.edges()(b, 0)), V.row(mesh.edges()(b, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Edge3, Point3>::DIM>
PrimitiveDistance<Edge3, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const index_t a,
    const index_t b,
    typename PrimitiveDistType<Edge3, Point3>::type dtype)
{
    return PointEdgeDistance<double, 3>::point_edge_closest_point_direction(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Edge2, Point2>::DIM>
PrimitiveDistance<Edge2, Point2>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const index_t a,
    const index_t b,
    typename PrimitiveDistType<Edge2, Point2>::type dtype)
{
    return PointEdgeDistance<double, 2>::point_edge_closest_point_direction(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Point2, Point2>::DIM>
PrimitiveDistance<Point2, Point2>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const index_t a,
    const index_t b,
    typename PrimitiveDistType<Point2, Point2>::type dtype)
{
    return V.row(b) - V.row(a);
}

template <>
Vector<double, PrimitiveDistance<Point3, Point3>::DIM>
PrimitiveDistance<Point3, Point3>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const index_t a,
    const index_t b,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    return V.row(b) - V.row(a);
}

#ifndef IPC_TOOLKIT_DEBUG_AUTODIFF

template <>
std::tuple<
    Vector<double, PrimitiveDistance<Edge3, Edge3>::DIM>,
    Eigen::Matrix<
        double,
        PrimitiveDistance<Edge3, Edge3>::DIM,
        PrimitiveDistance<Edge3, Edge3>::N_CORE_DOFS>,
    std::array<
        Eigen::Matrix<
            double,
            PrimitiveDistance<Edge3, Edge3>::N_CORE_DOFS,
            PrimitiveDistance<Edge3, Edge3>::N_CORE_DOFS>,
        PrimitiveDistance<Edge3, Edge3>::DIM>>
PrimitiveDistance<Edge3, Edge3>::compute_closest_direction_hessian(
    const Vector<double, N_CORE_DOFS>& x,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    assert(dtype == EdgeEdgeDistanceType::EA_EB);
    return line_line_closest_point_direction_hessian(
        x.head<3>() /* edge 0 */, x.segment<3>(3) /* edge 0 */,
        x.segment<3>(6) /* edge 1 */, x.tail<3>() /* edge 1 */);
}

template <>
std::tuple<
    Vector<double, PrimitiveDistance<Point2, Point2>::DIM>,
    Eigen::Matrix<
        double,
        PrimitiveDistance<Point2, Point2>::DIM,
        PrimitiveDistance<Point2, Point2>::N_CORE_DOFS>,
    std::array<
        Eigen::Matrix<
            double,
            PrimitiveDistance<Point2, Point2>::N_CORE_DOFS,
            PrimitiveDistance<Point2, Point2>::N_CORE_DOFS>,
        PrimitiveDistance<Point2, Point2>::DIM>>
PrimitiveDistance<Point2, Point2>::compute_closest_direction_hessian(
    const Vector<double, N_CORE_DOFS>& x,
    typename PrimitiveDistType<Point2, Point2>::type dtype)
{
    Eigen::Vector2d out = x.tail(2) - x.head(2);
    Eigen::Matrix<double, 2, 4> J = Eigen::Matrix<double, 2, 4>::Zero();
    J.leftCols<2>().diagonal().array() = -1;
    J.rightCols<2>().diagonal().array() = 1;
    std::array<Eigen::Matrix<double, 4, 4>, 2> H;
    for (auto& h : H) {
        h.setZero();
    }
    return std::make_tuple(out, J, H);
}

template <>
std::tuple<
    Vector<double, PrimitiveDistance<Point3, Point3>::DIM>,
    Eigen::Matrix<
        double,
        PrimitiveDistance<Point3, Point3>::DIM,
        PrimitiveDistance<Point3, Point3>::N_CORE_DOFS>,
    std::array<
        Eigen::Matrix<
            double,
            PrimitiveDistance<Point3, Point3>::N_CORE_DOFS,
            PrimitiveDistance<Point3, Point3>::N_CORE_DOFS>,
        PrimitiveDistance<Point3, Point3>::DIM>>
PrimitiveDistance<Point3, Point3>::compute_closest_direction_hessian(
    const Vector<double, N_CORE_DOFS>& x,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    Eigen::Vector3d out = x.tail(3) - x.head(3);
    Eigen::Matrix<double, 3, 6> J = Eigen::Matrix<double, 3, 6>::Zero();
    J.leftCols<3>().diagonal().array() = -1;
    J.rightCols<3>().diagonal().array() = 1;
    std::array<Eigen::Matrix<double, 6, 6>, 3> H;
    for (auto& h : H) {
        h.setZero();
    }
    return std::make_tuple(out, J, H);
}

template <>
GradType<PrimitiveDistance<Edge3, Edge3>::N_CORE_DOFS + 1>
PrimitiveDistance<Edge3, Edge3>::compute_mollifier_gradient(
    const Vector<double, N_CORE_DOFS>& x, const double dist_sqr)
{
    const auto otypes = edge_edge_mollifier_type(
        x.head(3), x.segment(3, 3), x.segment(6, 3), x.tail(3), dist_sqr);

    return edge_edge_mollifier_gradient(
        x.head(3), x.segment(3, 3), x.segment(6, 3), x.tail(3), otypes,
        dist_sqr);
}

template <>
HessianType<PrimitiveDistance<Edge3, Edge3>::N_CORE_DOFS + 1>
PrimitiveDistance<Edge3, Edge3>::compute_mollifier_hessian(
    const Vector<double, N_CORE_DOFS>& x, const double dist_sqr)
{
    const auto otypes = edge_edge_mollifier_type(
        x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), dist_sqr);

    return edge_edge_mollifier_hessian(
        x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), otypes,
        dist_sqr);
}

template <>
GradType<PrimitiveDistance<Face, Point3>::N_CORE_DOFS + 1>
PrimitiveDistance<Face, Point3>::compute_mollifier_gradient(
    const Vector<double, N_CORE_DOFS>& x, const double dist_sqr)
{
    const auto [val, grad] = point_face_mollifier_gradient(
        x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dist_sqr);
    Vector<int, 13> indices;
    indices << 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 12;
    return std::make_tuple(val, grad(indices));
}

template <>
HessianType<PrimitiveDistance<Face, Point3>::N_CORE_DOFS + 1>
PrimitiveDistance<Face, Point3>::compute_mollifier_hessian(
    const Vector<double, N_CORE_DOFS>& x, const double dist_sqr)
{
    const auto [val, grad, hess] = point_face_mollifier_hessian(
        x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dist_sqr);
    Vector<int, 13> indices;
    indices << 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 12;
    return std::make_tuple(val, grad(indices), hess(indices, indices));
}

#endif

template class PrimitiveDistance<Edge2, Point2>;
template class PrimitiveDistance<Point2, Point2>;

template class PrimitiveDistance<Edge3, Point3>;
template class PrimitiveDistance<Edge3, Edge3>;
template class PrimitiveDistance<Point3, Point3>;
template class PrimitiveDistance<Face, Point3>;
} // namespace ipc
