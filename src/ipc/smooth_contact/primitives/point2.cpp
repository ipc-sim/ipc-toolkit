#include "point2.hpp"

namespace ipc {

Point2::Point2(const long &id, 
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const VectorMax3d& d,
    const double &alpha,
    const double &beta)
: Primitive(id, alpha, beta)
{
    _vert_ids.resize(3);
    _vert_ids[0] = id;

    if (mesh.vertex_edge_adjacencies()[id].size() != 2)
        logger().error("Invalid number of vertex neighbor in 2D! {} should be 2.", mesh.vertex_edge_adjacencies()[id].size());
    for (long i : mesh.vertex_edge_adjacencies()[id])
    {
        if (mesh.edges()(i, 0) == id)
            _vert_ids[1] = mesh.edges()(i, 1);
        else if (mesh.edges()(i, 1) == id)
            _vert_ids[2] = mesh.edges()(i, 0);
        else
            logger().error("Wrong edge-vertex adjacency!");
    }

    is_active_ = smooth_point2_term_type(vertices.row(id), d, vertices.row(_vert_ids[1]), vertices.row(_vert_ids[2]), _alpha, _beta);
}

int Point2::n_vertices() const
{
    return 3;
}

double Point2::potential(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
{
    return smooth_point2_term<double>(x.segment<dim>(0), d, x.segment<dim>(dim), x.segment<dim>(dim), _alpha, _beta);
}
Vector<double, -1, Point2::max_size+Point2::dim> Point2::grad(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
{
    DiffScalarBase::setVariableCount(4*dim);
    using T = ADGrad<4*dim>;
    Vector<double, 4*dim> tmp;
    tmp << d, x;
    Eigen::Matrix<T, 4, dim> X = slice_positions<T, 4, dim>(tmp);
    return smooth_point2_term<T>(X.row(1), X.row(0), X.row(2), X.row(3), _alpha, _beta).getGradient();
}
MatrixMax<double, Point2::max_size+Point2::dim, Point2::max_size+Point2::dim> Point2::hessian(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
{
    DiffScalarBase::setVariableCount(4*dim);
    using T = ADHessian<4*dim>;
    Vector<double, 4*dim> tmp;
    tmp << d, x;
    Eigen::Matrix<T, 4, dim> X = slice_positions<T, 4, dim>(tmp);
    return smooth_point2_term<T>(X.row(1), X.row(0), X.row(2), X.row(3), _alpha, _beta).getHessian();
}

template <class scalar>
scalar smooth_point2_term(
    const Eigen::Ref<const Vector2<scalar>>& v,
    const Eigen::Ref<const Vector2<scalar>>& direc,
    const Eigen::Ref<const Vector2<scalar>>& e0,
    const Eigen::Ref<const Vector2<scalar>>& e1,
    const double &alpha,
    const double &beta)
{
    const Vector2<scalar> dn = -direc.normalized();
    const Vector2<scalar> t0 = (e0 - v).normalized(), t1 = (v - e1).normalized();

    const scalar tangent_term = Math<scalar>::smooth_heaviside(dn.dot(t0), alpha, beta) *
                        Math<scalar>::smooth_heaviside(-dn.dot(t1), alpha, beta);

    const scalar tmp = Math<scalar>::smooth_heaviside(-Math<scalar>::cross2(dn, t0), alpha, beta) + 
                         Math<scalar>::smooth_heaviside(-Math<scalar>::cross2(dn, t1), alpha, beta);
    const scalar normal_term = Math<scalar>::smooth_heaviside(tmp - 1., alpha, 0);

    return tangent_term * normal_term * ((e0 - v).norm() + (e1 - v).norm()) / 2.;
}

template double smooth_point2_term(
    const Eigen::Ref<const Vector2d>& v,
    const Eigen::Ref<const Vector2d>& direc,
    const Eigen::Ref<const Vector2d>& e0,
    const Eigen::Ref<const Vector2d>& e1,
    const double &alpha,
    const double &beta);
template ADGrad<12> smooth_point2_term(
    const Eigen::Ref<const Vector2<ADGrad<12>>>& v,
    const Eigen::Ref<const Vector2<ADGrad<12>>>& direc,
    const Eigen::Ref<const Vector2<ADGrad<12>>>& e0,
    const Eigen::Ref<const Vector2<ADGrad<12>>>& e1,
    const double &alpha,
    const double &beta);
template ADHessian<12> smooth_point2_term(
    const Eigen::Ref<const Vector2<ADHessian<12>>>& v,
    const Eigen::Ref<const Vector2<ADHessian<12>>>& direc,
    const Eigen::Ref<const Vector2<ADHessian<12>>>& e0,
    const Eigen::Ref<const Vector2<ADHessian<12>>>& e1,
    const double &alpha,
    const double &beta);
template ADGrad<10> smooth_point2_term(
    const Eigen::Ref<const Vector2<ADGrad<10>>>& v,
    const Eigen::Ref<const Vector2<ADGrad<10>>>& direc,
    const Eigen::Ref<const Vector2<ADGrad<10>>>& e0,
    const Eigen::Ref<const Vector2<ADGrad<10>>>& e1,
    const double &alpha,
    const double &beta);
template ADHessian<10> smooth_point2_term(
    const Eigen::Ref<const Vector2<ADHessian<10>>>& v,
    const Eigen::Ref<const Vector2<ADHessian<10>>>& direc,
    const Eigen::Ref<const Vector2<ADHessian<10>>>& e0,
    const Eigen::Ref<const Vector2<ADHessian<10>>>& e1,
    const double &alpha,
    const double &beta);

bool smooth_point2_term_type(
    const Eigen::Ref<const Vector2d>& v,
    const Eigen::Ref<const Vector2d>& direc,
    const Eigen::Ref<const Vector2d>& e0,
    const Eigen::Ref<const Vector2d>& e1,
    const double &alpha,
    const double &beta)
{
    const Vector2d dn = -direc.normalized();
    const Vector2d t0 = (e0 - v).normalized(), t1 = (v - e1).normalized();

    if (dn.dot(t0) <= -alpha || -dn.dot(t1) <= -alpha)
        return false;

    const double tmp = Math<double>::smooth_heaviside(-Math<double>::cross2(dn, t0), alpha, beta) + 
                         Math<double>::smooth_heaviside(-Math<double>::cross2(dn, t1), alpha, beta);
    if (tmp <= 1. - alpha)
        return false;

    return true;
}
}