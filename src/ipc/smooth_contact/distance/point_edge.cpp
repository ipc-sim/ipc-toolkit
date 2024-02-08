#include "point_edge.hpp"

namespace ipc {
template <typename scalar, int dim>
scalar PointEdgeDistance<scalar, dim>::point_point_sqr_distance(
    const Eigen::Ref<const Vector<scalar, dim>>& a,
    const Eigen::Ref<const Vector<scalar, dim>>& b)
{
    return (a - b).squaredNorm();
}

template <typename scalar, int dim>
scalar PointEdgeDistance<scalar, dim>::point_line_sqr_distance(
    const Eigen::Ref<const Vector<scalar, dim>>& p,
    const Eigen::Ref<const Vector<scalar, dim>>& e0,
    const Eigen::Ref<const Vector<scalar, dim>>& e1)
{
    if constexpr (dim == 2)
        return Math<scalar>::sqr(Math<scalar>::cross2(e0 - p, e1 - p))
            / (e1 - e0).squaredNorm();
    else
        return (e0 - p).cross(e1 - p).squaredNorm() / (e1 - e0).squaredNorm();
}

template <typename scalar, int dim>
scalar PointEdgeDistance<scalar, dim>::point_edge_sqr_distance(
    const Eigen::Ref<const Vector<scalar, dim>>& p,
    const Eigen::Ref<const Vector<scalar, dim>>& e0,
    const Eigen::Ref<const Vector<scalar, dim>>& e1,
    const PointEdgeDistanceType dtype)
{
    switch (dtype) {
    case PointEdgeDistanceType::P_E:
        return point_line_sqr_distance(p, e0, e1);
    case PointEdgeDistanceType::P_E0:
        return point_point_sqr_distance(p, e0);
    case PointEdgeDistanceType::P_E1:
        return point_point_sqr_distance(p, e1);
    case PointEdgeDistanceType::AUTO:
    default:
        const Vector<scalar, dim> t = e1 - e0;
        const Vector<scalar, dim> pos = p - e0;
        const scalar s = pos.dot(t) / t.squaredNorm();
        return (pos - Math<scalar>::L_ns(s) * t).squaredNorm();
    }
}

// template <int size>
// ADHessian<size> point_edge_sqr_distance(
//     const Eigen::Ref<const Vector<ADHessian<size>, dim>>& p,
//     const Eigen::Ref<const Vector<ADHessian<size>, dim>>& e0,
//     const Eigen::Ref<const Vector<ADHessian<size>, dim>>& e1,
//     const PointEdgeDistanceType dtype)
// {
//     std::cout << "This is partial specialization\n";
//     switch (dtype)
//     {
//     case PointEdgeDistanceType::P_E:
//         return point_line_sqr_distance<ADHessian<size>>(p, e0, e1);
//     case PointEdgeDistanceType::P_E0:
//         return point_point_sqr_distance<ADHessian<size>>(p, e0);
//     case PointEdgeDistanceType::P_E1:
//         return point_point_sqr_distance<ADHessian<size>>(p, e1);
//     case PointEdgeDistanceType::AUTO:
//     default:
//         const Vector<ADHessian<size>> t = e1 - e0;
//         const Vector<ADHessian<size>> pos = p - e0;
//         const ADHessian<size> s = pos.dot(t) / t.squaredNorm();
//         return (pos - Math<ADHessian<size>>::L_ns(s) * t).squaredNorm();
//     }
// }

template <typename scalar, int dim>
Vector<scalar, dim>
PointEdgeDistance<scalar, dim>::point_line_closest_point_direction(
    const Eigen::Ref<const Vector<scalar, dim>>& p,
    const Eigen::Ref<const Vector<scalar, dim>>& e0,
    const Eigen::Ref<const Vector<scalar, dim>>& e1)
{
    const Vector<scalar, dim> d = p - e0;
    const Vector<scalar, dim> t = e1 - e0;
    return d - (d.dot(t) / t.squaredNorm()) * t;
}

template <typename scalar, int dim>
Vector<scalar, dim>
PointEdgeDistance<scalar, dim>::point_edge_closest_point_direction(
    const Eigen::Ref<const Vector<scalar, dim>>& p,
    const Eigen::Ref<const Vector<scalar, dim>>& e0,
    const Eigen::Ref<const Vector<scalar, dim>>& e1,
    const PointEdgeDistanceType& dtype)
{
    switch (dtype) {
    case PointEdgeDistanceType::P_E:
        return point_line_closest_point_direction(p, e0, e1);
    case PointEdgeDistanceType::P_E0:
        return p - e0;
    case PointEdgeDistanceType::P_E1:
        return p - e1;
    case PointEdgeDistanceType::AUTO:
    default:
        Vector<scalar, dim> t = e1 - e0;
        const Vector<scalar, dim> pos = p - e0;
        const scalar s = pos.dot(t) / t.squaredNorm();
        return pos - Math<scalar>::L_ns(s) * t;
    }
}

template <int dim>
std::tuple<Vector<double, dim>, Eigen::Matrix<double, dim, dim * dim>>
PointEdgeDistanceDerivatives<dim>::point_line_closest_point_direction_grad(
    const Eigen::Ref<const Vector<double, dim>>& p,
    const Eigen::Ref<const Vector<double, dim>>& e0,
    const Eigen::Ref<const Vector<double, dim>>& e1)
{
    using T = ADGrad<dim * dim>;
    DiffScalarBase::setVariableCount(dim * dim);
    const Vector<T, dim> pT = slice_positions<T, 1, dim>(p);
    const Vector<T, dim> e0T = slice_positions<T, 1, dim>(e0, dim);
    const Vector<T, dim> e1T = slice_positions<T, 1, dim>(e1, 2 * dim);

    const Vector<T, dim> out =
        PointEdgeDistance<T, dim>::point_line_closest_point_direction(
            pT, e0T, e1T);
    Vector<double, dim> val;
    Eigen::Matrix<double, dim, dim * dim> grad;
    for (int i = 0; i < dim; i++) {
        val(i) = out(i).getValue();
        grad.row(i) = out(i).getGradient();
    }

    return std::make_tuple(val, grad);
}

template <int dim>
std::tuple<
    Vector<double, dim>,
    Eigen::Matrix<double, dim, dim * dim>,
    std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim>>
PointEdgeDistanceDerivatives<dim>::point_line_closest_point_direction_hessian(
    const Eigen::Ref<const Vector<double, dim>>& p,
    const Eigen::Ref<const Vector<double, dim>>& e0,
    const Eigen::Ref<const Vector<double, dim>>& e1)
{
    using T = ADHessian<dim * dim>;
    DiffScalarBase::setVariableCount(dim * dim);
    Vector<T, dim> pT = slice_positions<T, 1, dim>(p);
    Vector<T, dim> e0T = slice_positions<T, 1, dim>(e0, dim);
    Vector<T, dim> e1T = slice_positions<T, 1, dim>(e1, 2 * dim);

    Vector<T, dim> out =
        PointEdgeDistance<T, dim>::point_line_closest_point_direction(
            pT, e0T, e1T);
    Vector<double, dim> val;
    Eigen::Matrix<double, dim, dim * dim> grad;
    std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim> hess;
    for (int i = 0; i < dim; i++) {
        val(i) = out(i).getValue();
        grad.row(i) = out(i).getGradient();
        hess[i] = out(i).getHessian();
    }

    return std::make_tuple(val, grad, hess);
}

template <int dim>
Eigen::Matrix<double, dim, dim * dim>
PointEdgeDistanceDerivatives<dim>::point_edge_closest_point_direction_grad(
    const Eigen::Ref<const Vector<double, dim>>& p,
    const Eigen::Ref<const Vector<double, dim>>& e0,
    const Eigen::Ref<const Vector<double, dim>>& e1,
    const PointEdgeDistanceType& dtype)
{
    using T = ADGrad<dim * dim>;
    DiffScalarBase::setVariableCount(dim * dim);
    Vector<T, dim> pT = slice_positions<T, 1, dim>(p);
    Vector<T, dim> e0T = slice_positions<T, 1, dim>(e0, dim);
    Vector<T, dim> e1T = slice_positions<T, 1, dim>(e1, 2 * dim);

    Vector<T, dim> out =
        PointEdgeDistance<T, dim>::point_edge_closest_point_direction(
            pT, e0T, e1T);
    Eigen::Matrix<double, dim, dim * dim> grad;
    for (int i = 0; i < dim; i++)
        grad.row(i) = out(i).getGradient();

    return grad;
}

template <int dim>
std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim>
PointEdgeDistanceDerivatives<dim>::point_edge_closest_point_direction_hessian(
    const Eigen::Ref<const Vector<double, dim>>& p,
    const Eigen::Ref<const Vector<double, dim>>& e0,
    const Eigen::Ref<const Vector<double, dim>>& e1,
    const PointEdgeDistanceType& dtype)
{
    using T = ADHessian<dim * dim>;
    DiffScalarBase::setVariableCount(dim * dim);
    Vector<T, dim> pT = slice_positions<T, 1, dim>(p);
    Vector<T, dim> e0T = slice_positions<T, 1, dim>(e0, dim);
    Vector<T, dim> e1T = slice_positions<T, 1, dim>(e1, 2 * dim);

    Vector<T, dim> out =
        PointEdgeDistance<T, dim>::point_line_closest_point_direction(
            pT, e0T, e1T);
    std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim> hess;
    for (int i = 0; i < dim; i++)
        hess[i] = out(i).getHessian();

    return hess;
}

template class PointEdgeDistance<double, 2>;
template class PointEdgeDistance<double, 3>;

template class PointEdgeDistance<ADGrad<4>, 2>;
template class PointEdgeDistance<ADHessian<4>, 2>;

template class PointEdgeDistance<ADGrad<10>, 2>;
template class PointEdgeDistance<ADHessian<10>, 2>;

template class PointEdgeDistance<ADGrad<9>, 3>;
template class PointEdgeDistance<ADHessian<9>, 3>;

template class PointEdgeDistance<ADGrad<12>, 3>;
template class PointEdgeDistance<ADHessian<12>, 3>;

// #ifdef DERIVATIVES_WITH_AUTODIFF
template class PointEdgeDistance<ADGrad<13>, 3>;
template class PointEdgeDistance<ADHessian<13>, 3>;

template class PointEdgeDistance<ADGrad<15>, 3>;
template class PointEdgeDistance<ADHessian<15>, 3>;
// #endif

template class PointEdgeDistanceDerivatives<2>;
template class PointEdgeDistanceDerivatives<3>;
} // namespace ipc
