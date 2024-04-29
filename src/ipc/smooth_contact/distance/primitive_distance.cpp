#include "primitive_distance.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <ipc/smooth_contact/distance/edge_edge.hpp>
#include <ipc/smooth_contact/distance/point_edge.hpp>
#include <ipc/smooth_contact/distance/point_face.hpp>

namespace ipc {

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
typename PrimitiveDistType<Edge2, Point2>::type
PrimitiveDistance<Edge2, Point2>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return PointEdgeDistanceType::AUTO;
}

template <>
typename PrimitiveDistType<Point2, Point2>::type
PrimitiveDistance<Point2, Point2>::compute_distance_type(
    const Vector<double, n_core_dofs>& x)
{
    return PointPointDistanceType::AUTO;
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
GradType<PrimitiveDistance<Face, Point3>::n_core_dofs> 
PrimitiveDistance<Face, Point3>::compute_distance_gradient(
    const Vector<double, PrimitiveDistance<Face, Point3>::n_core_dofs>& x,
    typename PrimitiveDistType<Face, Point3>::type dtype)
{
    const double dist = point_triangle_distance(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dtype);
    Vector<double, n_core_dofs> dist_grad = point_triangle_distance_gradient(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dtype);
    dist_grad = dist_grad({3,4,5,6,7,8,9,10,11,0,1,2}).eval();
    return {dist, dist_grad};
}

template <>
HessianType<PrimitiveDistance<Face, Point3>::n_core_dofs> 
PrimitiveDistance<Face, Point3>::compute_distance_hessian(
    const Vector<double, PrimitiveDistance<Face, Point3>::n_core_dofs>& x,
    typename PrimitiveDistType<Face, Point3>::type dtype)
{
    const double dist = point_triangle_distance(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dtype);
    Vector<double, n_core_dofs> dist_grad = point_triangle_distance_gradient(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dtype);
    Eigen::Matrix<double, n_core_dofs, n_core_dofs> dist_hess = point_triangle_distance_hessian(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dtype);
    const std::vector<int> ind{3,4,5,6,7,8,9,10,11,0,1,2};
    dist_grad = dist_grad(ind).eval();
    dist_hess = dist_hess(ind, ind).eval();
    return {dist, dist_grad, dist_hess};
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
GradType<PrimitiveDistance<Edge3, Edge3>::n_core_dofs> 
PrimitiveDistance<Edge3, Edge3>::compute_distance_gradient(
    const Vector<double, PrimitiveDistance<Edge3, Edge3>::n_core_dofs>& x,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    return {
        edge_edge_distance(
            x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), dtype),
        edge_edge_distance_gradient(
            x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), dtype)};
}

template <>
HessianType<PrimitiveDistance<Edge3, Edge3>::n_core_dofs> 
PrimitiveDistance<Edge3, Edge3>::compute_distance_hessian(
    const Vector<double, PrimitiveDistance<Edge3, Edge3>::n_core_dofs>& x,
    typename PrimitiveDistType<Edge3, Edge3>::type dtype)
{
    return {
        edge_edge_distance(
            x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), dtype),
        edge_edge_distance_gradient(
            x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), dtype),
        edge_edge_distance_hessian(
            x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), dtype)};
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
GradType<PrimitiveDistance<Edge3, Point3>::n_core_dofs> 
PrimitiveDistance<Edge3, Point3>::compute_distance_gradient(
    const Vector<double, PrimitiveDistance<Edge3, Point3>::n_core_dofs>& x,
    typename PrimitiveDistType<Edge3, Point3>::type dtype)
{
    const double dist = point_edge_distance(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), dtype);
    Vector<double, n_core_dofs> dist_grad = point_edge_distance_gradient(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), dtype);
    dist_grad = dist_grad({3,4,5,6,7,8,0,1,2}).eval();
    return {dist, dist_grad};
}

template <>
HessianType<PrimitiveDistance<Edge3, Point3>::n_core_dofs> 
PrimitiveDistance<Edge3, Point3>::compute_distance_hessian(
    const Vector<double, PrimitiveDistance<Edge3, Point3>::n_core_dofs>& x,
    typename PrimitiveDistType<Edge3, Point3>::type dtype)
{
    const double dist = point_edge_distance(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), dtype);
    Vector<double, n_core_dofs> dist_grad = point_edge_distance_gradient(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), dtype);
    Eigen::Matrix<double, n_core_dofs, n_core_dofs> dist_hess = point_edge_distance_hessian(
            x.tail<3>(), x.head<3>(), x.segment<3>(3), dtype);
    const std::vector<int> ind{3,4,5,6,7,8,0,1,2};
    dist_grad = dist_grad(ind).eval();
    dist_hess = dist_hess(ind, ind).eval();
    return {dist, dist_grad, dist_hess};
}

template <>
GradType<PrimitiveDistance<Point3, Point3>::n_core_dofs> 
PrimitiveDistance<Point3, Point3>::compute_distance_gradient(
    const Vector<double, PrimitiveDistance<Point3, Point3>::n_core_dofs>& x,
    typename PrimitiveDistType<Point3, Point3>::type dtype)
{
    const double dist = (x.head<3>() - x.tail<3>()).squaredNorm();
    Vector<double, n_core_dofs> dist_grad;
    dist_grad << 2 * (x.head<3>() - x.tail<3>()), -2 * (x.head<3>() - x.tail<3>());
    return {dist, dist_grad};
}

template <>
double PrimitiveDistance<Edge2, Point2>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge2, Point2>::type dtype)
{
    return point_edge_distance(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
double PrimitiveDistance<Point2, Point2>::compute_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Point2, Point2>::type dtype)
{
    return point_point_distance(V.row(a), V.row(b));
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
Vector<double, PrimitiveDistance<Edge2, Point2>::dim>
PrimitiveDistance<Edge2, Point2>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Edge2, Point2>::type dtype)
{
    return PointEdgeDistance<double, 2>::point_edge_closest_point_direction(
        V.row(b), V.row(mesh.edges()(a, 0)), V.row(mesh.edges()(a, 1)), dtype);
}

template <>
Vector<double, PrimitiveDistance<Point2, Point2>::dim>
PrimitiveDistance<Point2, Point2>::compute_closest_direction(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const long& a,
    const long& b,
    typename PrimitiveDistType<Point2, Point2>::type dtype)
{
    return V.row(b) - V.row(a);
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
std::tuple<Vector<double, PrimitiveDistance<Point2, Point2>::dim>, Eigen::Matrix<double, PrimitiveDistance<Point2, Point2>::dim, PrimitiveDistance<Point2, Point2>::n_core_dofs>,
    std::array<Eigen::Matrix<double, PrimitiveDistance<Point2, Point2>::n_core_dofs, PrimitiveDistance<Point2, Point2>::n_core_dofs>, PrimitiveDistance<Point2, Point2>::dim>>
PrimitiveDistance<Point2, Point2>::compute_closest_direction_hessian(
    const Vector<double, n_core_dofs>& x,
    typename PrimitiveDistType<Point2, Point2>::type dtype)
{
    Vector2d out = x.tail(2) - x.head(2);
    Eigen::Matrix<double, 2, 4> J = Eigen::Matrix<double, 2, 4>::Zero();
    J.leftCols<2>().diagonal().array() = -1;
    J.rightCols<2>().diagonal().array() = 1;
    std::array<Eigen::Matrix<double, 4, 4>, 2> H;
    for (auto &h : H)
        h.setZero();
    return std::make_tuple(out, J, H);
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

template <>
GradType<PrimitiveDistance<Edge3, Edge3>::n_core_dofs + 1>
PrimitiveDistance<Edge3, Edge3>::compute_mollifier_gradient(
    const Vector<double, n_core_dofs>& x, const double dist_sqr)
{
    const auto otypes = edge_edge_mollifier_type(
        x.head(3), x.segment(3, 3), x.segment(6, 3),
        x.tail(3), dist_sqr);

    return edge_edge_mollifier_gradient(
        x.head(3), x.segment(3, 3), x.segment(6, 3),
        x.tail(3), otypes, dist_sqr);
}

template <>
HessianType<PrimitiveDistance<Edge3, Edge3>::n_core_dofs + 1>
PrimitiveDistance<Edge3, Edge3>::compute_mollifier_hessian(
    const Vector<double, n_core_dofs>& x, const double dist_sqr)
{
    const auto otypes = edge_edge_mollifier_type(
        x.head<3>(), x.segment<3>(3), x.segment<3>(6),
        x.tail<3>(), dist_sqr);

    return edge_edge_mollifier_hessian(
        x.head<3>(), x.segment<3>(3), x.segment<3>(6),
        x.tail<3>(), otypes, dist_sqr);
}

template <>
GradType<PrimitiveDistance<Face, Point3>::n_core_dofs + 1>
PrimitiveDistance<Face, Point3>::compute_mollifier_gradient(
    const Vector<double, n_core_dofs>& x, const double dist_sqr)
{
    const auto [val, grad] = point_face_mollifier_gradient(
        x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dist_sqr);
    Vector<int, 13> indices;
    indices << 3,4,5,6,7,8,9,10,11,0,1,2,12;
    return std::make_tuple(val, grad(indices));
}

template <>
HessianType<PrimitiveDistance<Face, Point3>::n_core_dofs + 1>
PrimitiveDistance<Face, Point3>::compute_mollifier_hessian(
    const Vector<double, n_core_dofs>& x, const double dist_sqr)
{
    const auto [val, grad, hess] = point_face_mollifier_hessian(
        x.tail<3>(), x.head<3>(), x.segment<3>(3), x.segment<3>(6), dist_sqr);
    Vector<int, 13> indices;
    indices << 3,4,5,6,7,8,9,10,11,0,1,2,12;
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
