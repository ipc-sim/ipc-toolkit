#include "edge3.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
// d is a vector from any point on the edge to the point outside of the edge
Edge3::Edge3(
    const long& id,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const VectorMax3d& d,
    const ParameterType& param)
    : Primitive(id, param)
{
    auto ids = mesh.find_edge_adjacent_vertices(id);
    _vert_ids = std::vector<long>(ids.begin(), ids.begin() + ids.size());
    is_active_ = smooth_edge3_term_type(
        d.normalized(), vertices.row(_vert_ids[0]), vertices.row(_vert_ids[1]),
        vertices.row(_vert_ids[2]), vertices.row(_vert_ids[3]), _param,
        otypes);
}

int Edge3::n_vertices() const { return n_edge_neighbors_3d; }

double Edge3::potential(
    const Eigen::Ref<const Vector3d>& d,
    const Eigen::Ref<const Vector12d>& x) const
{
#ifdef DERIVATIVES_WITH_AUTODIFF
    return smooth_edge3_term_template<double>(
        d.normalized(), x.head<3>(), x.segment<3>(3), x.segment<3>(6),
        x.tail<3>(), _param, otypes);
#else
    return smooth_edge3_term(
        d, x.head<3>(), x.segment<3>(3), x.segment<3>(6),
        x.tail<3>(), _param, otypes);
#endif
}
Vector15d Edge3::grad(
    const Eigen::Ref<const Vector3d>& d,
    const Eigen::Ref<const Vector12d>& x) const
{
#ifdef DERIVATIVES_WITH_AUTODIFF
    Vector15d tmp;
    tmp << d, x;
    DiffScalarBase::setVariableCount(15);
    using T = ADGrad<15>;
    auto X = slice_positions<T, 5, 3>(tmp);
    return smooth_edge3_term_template<T>(
               X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3),
               X.row(4), _param, otypes)
        .getGradient();
#else
    return std::get<1>(smooth_edge3_term_gradient(
        d, x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), _param, otypes));
#endif
}
Matrix15d Edge3::hessian(
    const Eigen::Ref<const Vector3d>& d,
    const Eigen::Ref<const Vector12d>& x) const
{
#ifdef DERIVATIVES_WITH_AUTODIFF
    Vector15d tmp;
    tmp << d, x;
    DiffScalarBase::setVariableCount(15);
    using T = ADHessian<15>;
    auto X = slice_positions<T, 5, 3>(tmp);
    return smooth_edge3_term_template<T>(
               X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3),
               X.row(4), _param, otypes)
        .getHessian();
#else
    return std::get<2>(smooth_edge3_term_hessian(
        d, x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), _param, otypes));
#endif
}

bool smooth_edge3_term_type(
    const Eigen::Ref<const Vector3d>& dn,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const ParameterType& param,
    ORIENTATION_TYPES& otypes)
{
    otypes.set_size(2);

    const Vector3d t0 =
        PointEdgeDistance<double, 3>::point_line_closest_point_direction(
            f0, e0, e1);
    const Vector3d t1 =
        PointEdgeDistance<double, 3>::point_line_closest_point_direction(
            f1, e0, e1);
    otypes.tangent_type(0) =
        otypes.compute_type(-dn.dot(t0) / t0.norm(), param.alpha_t, param.beta_t);
    otypes.tangent_type(1) =
        otypes.compute_type(-dn.dot(t1) / t1.norm(), param.alpha_t, param.beta_t);
    if (otypes.tangent_type(0) == HEAVISIDE_TYPE::ZERO
        || otypes.tangent_type(1) == HEAVISIDE_TYPE::ZERO)
        return false;

    const Vector3d n0 = (e0 - f0).cross(e1 - f0);
    const Vector3d n1 = -(e0 - f1).cross(e1 - f1);
    const double tmp0 = dn.dot(n0) / n0.norm();
    const double tmp1 = dn.dot(n1) / n1.norm();
    otypes.normal_type(0) = otypes.compute_type(tmp0, param.alpha_n, param.beta_n);
    otypes.normal_type(1) = otypes.compute_type(tmp1, param.alpha_n, param.beta_n);
    const double sum = Math<double>::smooth_heaviside(tmp0, param.alpha_n, param.beta_n)
        + Math<double>::smooth_heaviside(tmp1, param.alpha_n, param.beta_n);
    if (sum <= 1 - param.alpha_n)
        return false;
    else if (sum >= 1) {
        otypes.normal_type(0) = HEAVISIDE_TYPE::ONE;
        otypes.normal_type(1) = HEAVISIDE_TYPE::ONE;
    }

    return true;
}

double smooth_edge3_normal_term(
    const Eigen::Ref<const Vector3d>& dn,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes)
{
    if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE
        || otypes.normal_type(1) == HEAVISIDE_TYPE::ONE)
        return 1.;

    return Math<double>::smooth_heaviside(
        negative_orientation_penalty(e0 - f0, e1 - f0, dn, alpha, beta)
            + negative_orientation_penalty(e0 - f1, e1 - f1, -dn, alpha, beta)
            - 1,
        alpha, 0);
}

std::tuple<double, Vector15d> smooth_edge3_normal_term_gradient(
    const Eigen::Ref<const Vector3d>& dn,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes)
{
    double val = 1.;
    Vector15d gradient = Vector15d::Zero();

    if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE
        || otypes.normal_type(1) == HEAVISIDE_TYPE::ONE)
        return std::make_tuple(val, gradient);

    val = 0.;

    if (otypes.normal_type(0) == HEAVISIDE_TYPE::VARIANT) {
        const auto [y, dy] = negative_orientation_penalty_grad(
            e0 - f0, e1 - f0, dn, alpha, beta);

        val += y;

        // dn
        gradient.head<3>() += dy.tail<3>();
        // e0, e1
        gradient.segment<6>(3) += dy.head<6>();
        // f0
        gradient.segment<3>(9) -= dy.segment<3>(0) + dy.segment<3>(3);
    }

    if (otypes.normal_type(1) == HEAVISIDE_TYPE::VARIANT) {
        const auto [y, dy] = negative_orientation_penalty_grad(
            e0 - f1, e1 - f1, -dn, alpha, beta);

        val += y;

        // dn
        gradient.head<3>() -= dy.tail<3>();
        // e0, e1
        gradient.segment<6>(3) += dy.head<6>();
        // f1
        gradient.segment<3>(12) -= dy.segment<3>(0) + dy.segment<3>(3);
    }

    gradient *= Math<double>::smooth_heaviside_grad(val - 1, alpha, 0);
    val = Math<double>::smooth_heaviside(val - 1, alpha, 0);
    return std::make_tuple(val, gradient);
}

std::tuple<double, Vector15d, Matrix15d> smooth_edge3_normal_term_hessian(
    const Eigen::Ref<const Vector3d>& dn,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes)
{
    double value = 1.;
    Vector15d gradient = Vector15d::Zero();
    Matrix15d hessian = Matrix15d::Zero();

    if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE
        || otypes.normal_type(1) == HEAVISIDE_TYPE::ONE)
        return std::make_tuple(value, gradient, hessian);

    value = 0.;

    for (int d : { 0, 1 }) {
        const Eigen::Ref<const Vector3d> f = (d == 0) ? f0 : f1;
        const int sign = (d == 0) ? 1 : -1;
        const int f_id = (d == 0) ? 9 : 12;
        if (otypes.normal_type(d) == HEAVISIDE_TYPE::VARIANT) {
            auto [y, dy, ddy] = negative_orientation_penalty_hess(
                e0 - f, e1 - f, sign * dn, alpha, beta);
            dy.tail<3>() *= sign;
            ddy.bottomRows<3>() *= sign;
            ddy.rightCols<3>() *= sign;

            value += y;

            // dn
            gradient.head<3>() += dy.tail<3>();
            // e0, e1
            gradient.segment<6>(3) += dy.head<6>();
            // f
            gradient.segment<3>(f_id) -= dy.segment<3>(0) + dy.segment<3>(3);

            // dn
            hessian.block<3, 3>(0, 0) += ddy.block<3, 3>(6, 6);
            // e0, e1
            hessian.block<6, 6>(3, 3) += ddy.block<6, 6>(0, 0);
            // f
            hessian.block<3, 3>(f_id, f_id) += ddy.block<3, 3>(0, 0)
                + ddy.block<3, 3>(3, 3) + ddy.block<3, 3>(0, 3)
                + ddy.block<3, 3>(3, 0);
            // dn & e0, e1
            hessian.block<3, 6>(0, 3) += ddy.block<3, 6>(6, 0);
            hessian.block<6, 3>(3, 0) += ddy.block<6, 3>(0, 6);
            // dn & f
            hessian.block<3, 3>(0, f_id) -=
                ddy.block<3, 3>(6, 0) + ddy.block<3, 3>(6, 3);
            hessian.block<3, 3>(f_id, 0) -=
                ddy.block<3, 3>(0, 6) + ddy.block<3, 3>(3, 6);
            // f & e0, e1
            hessian.block<6, 3>(3, f_id) -=
                ddy.block<6, 3>(0, 0) + ddy.block<6, 3>(0, 3);
            hessian.block<3, 6>(f_id, 3) -=
                ddy.block<3, 6>(0, 0) + ddy.block<3, 6>(3, 0);
        }
    }

    const double hess_val =
        Math<double>::smooth_heaviside_hess(value - 1, alpha, 0);
    const double grad_val =
        Math<double>::smooth_heaviside_grad(value - 1, alpha, 0);

    hessian = gradient * hess_val * gradient.transpose() + grad_val * hessian;
    gradient = gradient * grad_val;
    
    // verify with autodiff
    // {
    //     DiffScalarBase::setVariableCount(15);
    //     using T = ADHessian<15>;
    //     Vector15d tmp;
    //     tmp << dn, e0, e1, f0, f1;
    //     auto X = slice_positions<T, 15, 1>(tmp);
    //     const Vector3<T> n0 = (X.segment<3>(3) - X.segment<3>(9)).cross(X.segment<3>(6) - X.segment<3>(9));
    //     const Vector3<T> n1 = -(X.segment<3>(3) - X.tail<3>()).cross(X.segment<3>(6) - X.tail<3>());
    //     T normal_term = Math<T>::smooth_heaviside(
    //         (Math<T>::smooth_heaviside(X.head<3>().dot(n0) / n0.norm(), alpha, beta)
    //         + Math<T>::smooth_heaviside(X.head<3>().dot(n1) / n1.norm(), alpha, beta)
    //         - 1), alpha, 0);
        
    //     if ((normal_term.getHessian() - hessian).norm() > 1e-6 * hessian.norm()) {
    //         std::cout << "gradient: \n" << normal_term.getGradient().transpose() << "\n vs \n" << gradient.transpose() << std::endl;
    //         std::cout << "hessian mismatch: \n" << normal_term.getHessian() << "\n vs \n" << hessian << std::endl;
    //     }
    // }

    return std::make_tuple(
        Math<double>::smooth_heaviside(value - 1, alpha, 0),
        gradient,
        hessian);
}

double smooth_edge3_tangent_term(
    const Eigen::Ref<const Vector3d>& dn,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes)
{
    double tangent_term = 1.;
    if (otypes.tangent_type(0) != HEAVISIDE_TYPE::ONE) {
        const Vector3d t0 =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f0, e0, e1);
        tangent_term *= opposite_direction_penalty(t0, -dn, alpha, beta);
    }
    if (otypes.tangent_type(1) != HEAVISIDE_TYPE::ONE) {
        const Vector3d t1 =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f1, e0, e1);
        tangent_term *= opposite_direction_penalty(t1, -dn, alpha, beta);
    }

    return tangent_term;
}

std::tuple<double, Vector15d> smooth_edge3_tangent_term_gradient(
    const Eigen::Ref<const Vector3d>& dn,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes)
{
    Vector2d vals;
    vals << 1., 1.;
    std::array<Vector<double, 15>, 2> grads;
    for (auto& g : grads)
        g.setZero();

    for (int d : { 0, 1 }) {
        const Eigen::Ref<const Vector3d> f = d == 0 ? f0 : f1;
        if (otypes.tangent_type(d) != HEAVISIDE_TYPE::ONE) {
            const auto [t, g] = PointEdgeDistanceDerivatives<
                3>::point_line_closest_point_direction_grad(f, e0, e1);
            const auto [tmp_val, tmp_grad] =
                opposite_direction_penalty_grad(t, -dn, alpha, beta);

            vals[d] = tmp_val;

            Vector<double, 12> gradient_tmp;
            gradient_tmp << -tmp_grad.tail<3>(),
                g.transpose() * tmp_grad.head<3>();

            Vector<int, 12> indices;
            indices << 0, 1, 2, 9, 10, 11, 3, 4, 5, 6, 7, 8;
            if (d == 1)
                indices.segment<3>(3).array() += 3;

            grads[d](indices) = gradient_tmp;
        }
    }

    return std::make_tuple(
        vals.prod(), vals[0] * grads[1] + vals[1] * grads[0]);
}

std::tuple<double, Vector15d, Matrix15d> smooth_edge3_tangent_term_hessian(
    const Eigen::Ref<const Vector3d>& dn,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const double alpha,
    const double beta,
    const ORIENTATION_TYPES& otypes)
{
    Vector2d vals;
    vals << 1., 1.;
    std::array<Vector<double, 15>, 2> grads;
    std::array<Eigen::Matrix<double, 15, 15>, 2> hesses;
    for (auto& g : grads)
        g.setZero();
    for (auto& h : hesses)
        h.setZero();

    for (int d : { 0, 1 }) {
        const Eigen::Ref<const Vector3d> f = d == 0 ? f0 : f1;
        if (otypes.tangent_type(d) != HEAVISIDE_TYPE::ONE) {
            const auto [t, g, h] = PointEdgeDistanceDerivatives<
                3>::point_line_closest_point_direction_hessian(f, e0, e1);
            auto [tmp_val, tmp_grad, tmp_hess] =
                opposite_direction_penalty_hess(t, -dn, alpha, beta);
            tmp_grad.tail<3>() *= -1;
            tmp_hess.rightCols<3>() *= -1;
            tmp_hess.bottomRows<3>() *= -1;

            vals[d] = tmp_val;

            Vector<double, 12> gradient_tmp;
            gradient_tmp << tmp_grad.tail<3>(),
                g.transpose() * tmp_grad.head<3>();

            Eigen::Matrix<double, 12, 12> hessian_tmp;
            // dn
            hessian_tmp.block<3, 3>(0, 0) = tmp_hess.block<3, 3>(3, 3);
            // f, e0, e1
            hessian_tmp.block<9, 9>(3, 3) =
                g.transpose() * tmp_hess.block<3, 3>(0, 0) * g
                + tmp_grad(0) * h[0] + tmp_grad(1) * h[1] + tmp_grad(2) * h[2];

            // mixed 2nd derivatives
            hessian_tmp.block<3, 9>(0, 3) = tmp_hess.block<3, 3>(3, 0) * g;
            hessian_tmp.block<9, 3>(3, 0) =
                g.transpose() * tmp_hess.block<3, 3>(0, 3);

            Vector<int, 12> indices;
            indices << 0, 1, 2, 9, 10, 11, 3, 4, 5, 6, 7, 8;
            if (d == 1)
                indices.segment<3>(3).array() += 3;

            hesses[d](indices, indices) = hessian_tmp;
            grads[d](indices) = gradient_tmp;
        }
    }

    return std::make_tuple(
        vals.prod(), grads[0] * vals[1] + grads[1] * vals[0],
        hesses[0] * vals[1] + hesses[1] * vals[0]
            + grads[0] * grads[1].transpose()
            + grads[1] * grads[0].transpose());
}

double smooth_edge3_term(
    const Eigen::Ref<const Vector3d>& direc,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const ParameterType& param,
    const ORIENTATION_TYPES& otypes)
{
    const Vector3d dn = direc.normalized();
    double tangent_term =
        smooth_edge3_tangent_term(dn, e0, e1, f0, f1, param.alpha_t, param.beta_t, otypes);
    double normal_term =
        smooth_edge3_normal_term(dn, e0, e1, f0, f1, param.alpha_n, param.beta_n, otypes);

    return (e1 - e0).squaredNorm() * tangent_term * normal_term;
}

std::tuple<double, Vector15d> smooth_edge3_term_gradient(
    const Eigen::Ref<const Vector3d>& direc,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const ParameterType& param,
    const ORIENTATION_TYPES& otypes)
{
    assert(otypes.size() == 2);

    const auto [dn, dn_grad] = normalize_vector_grad(direc);

    auto [t_term, t_grad] = smooth_edge3_tangent_term_gradient(
        dn, e0, e1, f0, f1, param.alpha_t, param.beta_t, otypes);
    auto [n_term, n_grad] = smooth_edge3_normal_term_gradient(dn, e0, e1, f0,
    f1, param.alpha_n, param.beta_n, otypes);

    t_grad.head<3>() = dn_grad * t_grad.head<3>();
    n_grad.head<3>() = dn_grad * n_grad.head<3>();

    const double weight = (e1 - e0).squaredNorm();
    const double val = weight * t_term * n_term;

    Vector15d gradient =
        (weight * t_term) * n_grad + (weight * n_term) * t_grad;
    gradient.segment<3>(3) += (2 * t_term * n_term) * (e0 - e1);
    gradient.segment<3>(6) += (2 * t_term * n_term) * (e1 - e0);

    return std::make_tuple(val, gradient);
}

std::tuple<double, Vector15d, Matrix15d> smooth_edge3_term_hessian(
    const Eigen::Ref<const Vector3d>& direc,
    const Eigen::Ref<const Vector3d>& e0,
    const Eigen::Ref<const Vector3d>& e1,
    const Eigen::Ref<const Vector3d>& f0,
    const Eigen::Ref<const Vector3d>& f1,
    const ParameterType& param,
    const ORIENTATION_TYPES& otypes)
{
    assert(otypes.size() == 2);

    const auto [dn, dn_grad, dn_hess] = normalize_vector_hess(direc);

    auto [t_term, t_grad, t_hess] = smooth_edge3_tangent_term_hessian(
        dn, e0, e1, f0, f1, param.alpha_t, param.beta_t, otypes);
    auto [n_term, n_grad, n_hess] = smooth_edge3_normal_term_hessian(dn, e0,
    e1, f0, f1, param.alpha_n, param.beta_n, otypes);

    t_hess.topRows<3>() = dn_grad * t_hess.topRows<3>();
    t_hess.leftCols<3>() = t_hess.leftCols<3>() * dn_grad;
    t_hess.topLeftCorner<3, 3>() += dn_hess[0] * t_grad(0)
        + dn_hess[1] * t_grad(1) + dn_hess[2] * t_grad(2);
    t_grad.head<3>() = dn_grad * t_grad.head<3>();

    n_hess.topRows<3>() = dn_grad * n_hess.topRows<3>();
    n_hess.leftCols<3>() = n_hess.leftCols<3>() * dn_grad;
    n_hess.topLeftCorner<3, 3>() += dn_hess[0] * n_grad(0)
        + dn_hess[1] * n_grad(1) + dn_hess[2] * n_grad(2);
    n_grad.head<3>() = dn_grad * n_grad.head<3>();

    const double weight = (e1 - e0).squaredNorm();

    Vector15d gradient = t_term * n_grad + n_term * t_grad;

    Matrix15d hessian = (weight * t_term) * n_hess + (weight * n_term) * t_hess
        + (weight * t_grad) * n_grad.transpose()
        + (weight * n_grad) * t_grad.transpose();

    const double tmp_val = 2 * t_term * n_term;
    hessian.block<6, 6>(3, 3).diagonal().array() += tmp_val;
    hessian.block<3, 3>(3, 6).diagonal().array() -= tmp_val;
    hessian.block<3, 3>(6, 3).diagonal().array() -= tmp_val;

    Vector6d weight_grad;
    weight_grad << 2 * (e0 - e1), 2 * (e1 - e0);
    hessian.middleCols<6>(3) += gradient * weight_grad.transpose();
    hessian.middleRows<6>(3) += weight_grad * gradient.transpose();

    gradient *= weight;
    gradient.segment<3>(3) += tmp_val * (e0 - e1);
    gradient.segment<3>(6) += tmp_val * (e1 - e0);

    return std::make_tuple(weight * t_term * n_term, gradient, hessian);
}

template <typename scalar>
scalar smooth_edge3_term_template(
    const Eigen::Ref<const Vector3<scalar>>& dn,
    const Eigen::Ref<const Vector3<scalar>>& e0,
    const Eigen::Ref<const Vector3<scalar>>& e1,
    const Eigen::Ref<const Vector3<scalar>>& f0,
    const Eigen::Ref<const Vector3<scalar>>& f1,
    const ParameterType& param,
    const ORIENTATION_TYPES& otypes)
{
    scalar tangent_term = scalar(1.);
    if (otypes.tangent_type(0) != HEAVISIDE_TYPE::ONE) {
        const Vector3<scalar> t0 =
            PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                f0, e0, e1);
        tangent_term = tangent_term
            * Math<scalar>::smooth_heaviside(
                           -dn.dot(t0) / t0.norm(), param.alpha_t, param.beta_t);
    }
    if (otypes.tangent_type(1) != HEAVISIDE_TYPE::ONE) {
        const Vector3<scalar> t1 =
            PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                f1, e0, e1);
        tangent_term = tangent_term
            * Math<scalar>::smooth_heaviside(
                           -dn.dot(t1) / t1.norm(), param.alpha_t, param.beta_t);
    }

    scalar normal_term = scalar(1.);
    if (otypes.normal_type(0) != HEAVISIDE_TYPE::ONE
        && otypes.normal_type(1) != HEAVISIDE_TYPE::ONE) {
        const Vector3<scalar> n0 = (e0 - f0).cross(e1 - f0);
        const Vector3<scalar> n1 = -(e0 - f1).cross(e1 - f1);
        normal_term = Math<scalar>::smooth_heaviside(
            (Math<scalar>::smooth_heaviside(dn.dot(n0) / n0.norm(), param.alpha_n, param.beta_n)
             + Math<scalar>::smooth_heaviside(
                 dn.dot(n1) / n1.norm(), param.alpha_n, param.beta_n)
             - 1),
            param.alpha_n, 0);
    }

    return (e1 - e0).squaredNorm() * tangent_term * normal_term;
}

#ifdef DERIVATIVES_WITH_AUTODIFF
template double smooth_edge3_term_template(
    const Eigen::Ref<const Vector3<double>>& dn,
    const Eigen::Ref<const Vector3<double>>& e0,
    const Eigen::Ref<const Vector3<double>>& e1,
    const Eigen::Ref<const Vector3<double>>& f0,
    const Eigen::Ref<const Vector3<double>>& f1,
    const ParameterType& param,
    const ORIENTATION_TYPES& otypes);
template ADGrad<15> smooth_edge3_term_template(
    const Eigen::Ref<const Vector3<ADGrad<15>>>& dn,
    const Eigen::Ref<const Vector3<ADGrad<15>>>& e0,
    const Eigen::Ref<const Vector3<ADGrad<15>>>& e1,
    const Eigen::Ref<const Vector3<ADGrad<15>>>& f0,
    const Eigen::Ref<const Vector3<ADGrad<15>>>& f1,
    const ParameterType& param,
    const ORIENTATION_TYPES& otypes);
template ADHessian<15> smooth_edge3_term_template(
    const Eigen::Ref<const Vector3<ADHessian<15>>>& dn,
    const Eigen::Ref<const Vector3<ADHessian<15>>>& e0,
    const Eigen::Ref<const Vector3<ADHessian<15>>>& e1,
    const Eigen::Ref<const Vector3<ADHessian<15>>>& f0,
    const Eigen::Ref<const Vector3<ADHessian<15>>>& f1,
    const ParameterType& param,
    const ORIENTATION_TYPES& otypes);
#endif
} // namespace ipc
