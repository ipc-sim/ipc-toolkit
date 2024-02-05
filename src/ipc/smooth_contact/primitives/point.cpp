#include "point.hpp"
#include <ipc/utils/AutodiffTypes.hpp>

namespace ipc {
    namespace {
        template <int size, typename T = ADGrad<size * 3>>
        T evaluate(const Vector<double, size * 3> &input, const double &alpha, const double &beta, const ORIENTATION_TYPES &otypes)
        {
            DiffScalarBase::setVariableCount(size * 3);
            auto X = slice_positions<T, size, 3>(input);
            return smooth_point3_term<T, size - 2>(X.row(1), X.row(0), X.bottomRows(size - 2), alpha, beta, otypes);
        }

    }

    Point3::Point3(const long &id, 
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const VectorMax3d& d,
        const double &alpha,
        const double &beta)
    : Primitive(id, alpha, beta)
    {
        auto neighbor_ids = mesh.find_vertex_adjacent_vertices(id);
        n_neighbors = neighbor_ids.size();
        _vert_ids.reserve(1 + n_neighbors);
        _vert_ids.push_back(id);
        _vert_ids.insert( _vert_ids.end(), neighbor_ids.begin(), neighbor_ids.end() );

        if (_vert_ids.size() > n_vert_neighbors_3d)
            logger().error("Too many neighbors for point3 primitive! {} > {}", _vert_ids.size(), n_vert_neighbors_3d);

        Eigen::Vector3d v = vertices.row(id);
        MatrixMax<double, n_vert_neighbors_3d, 3> neighbors(n_neighbors, dim);
        int k = 0;
        for (long i : neighbor_ids)
            neighbors.row(k++) = vertices.row(i);

        is_active_ = smooth_point3_term_type(v, d, neighbors, _alpha, _beta, otypes);
    }

    int Point3::n_vertices() const
    {
        return n_neighbors + 1;
    }
    
    double Point3::potential(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
    {
        const Eigen::Matrix<double, -1, dim> X = slice_positions<double, -1, dim>(x);
        return smooth_point3_term<double, -1>(X.row(0), d, X.bottomRows(n_neighbors), _alpha, _beta, otypes);
    }
    Vector<double, -1, Point3::max_size+Point3::dim> Point3::grad(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
    {
        const Eigen::Matrix<double, -1, dim> X = slice_positions<double, -1, dim>(x);
        const auto [val, grad] = smooth_point3_term_gradient(d, X.row(0), X.bottomRows(n_neighbors), _alpha, _beta, otypes);
        return grad;
    }
    MatrixMax<double, Point3::max_size+Point3::dim, Point3::max_size+Point3::dim> Point3::hessian(const Vector<double, dim> &d, const Vector<double, -1, max_size> &x) const
    {
        const auto X = slice_positions<double, -1, dim>(x);
        const auto [val, grad, hess] = smooth_point3_term_hessian(d, X.row(0), X.bottomRows(n_neighbors), _alpha, _beta, otypes);
        return hess;
    }

std::tuple<double, Eigen::VectorXd> smooth_point3_term_gradient(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& neighbors,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    const int n_neighbors = neighbors.rows();
    const int n_dofs = (n_neighbors + 2) * 3;
    const int n_neighbor_dofs = n_neighbors * 3;
    assert(n_neighbors > 2);
    assert(otypes.size() == n_neighbors);

    const Eigen::Matrix<double, -1, 3, Eigen::RowMajor> tangents = neighbors.rowwise() - v;
    const Eigen::VectorXd tangents_vec = Eigen::Map<const Eigen::VectorXd>(tangents.data(), tangents.size());
    auto [dn, dn_grad] = normalize_vector_grad(direc);
    dn *= -1;
    dn_grad *= -1;

    const auto [tangent_term, tangent_grad] = smooth_point3_term_tangent_gradient(dn, tangents, alpha, beta, otypes);
    const auto [normal_term, normal_grad] = smooth_point3_term_normal_gradient(dn, tangents, alpha, beta, otypes);

    double val = tangent_term * normal_term;

    // gradient wrt. [dn, tangents]
    Eigen::VectorXd grad_tmp = tangent_grad * normal_term + normal_grad * tangent_term;

    const double weight = tangents.squaredNorm() / 3.;
    grad_tmp *= weight;
    grad_tmp.tail(n_neighbor_dofs) += (2. / 3. * val) * tangents_vec;

    val *= weight;

    // gradient wrt. [direc, v, neighbors]
    Eigen::VectorXd grad(n_dofs);
    grad.head<3>() = dn_grad * grad_tmp.head<3>();
    for (int d = 0; d < 3; d++)
        grad(d + 3) = -grad_tmp(Eigen::seqN(d + 3, n_neighbors, 3)).sum();
    grad.tail(n_neighbor_dofs) = grad_tmp.tail(n_neighbor_dofs);

    return std::tuple(val, grad);
}

std::tuple<double, Eigen::VectorXd, Eigen::MatrixXd> smooth_point3_term_hessian(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& neighbors,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    const int n_neighbors = neighbors.rows();
    const int n_dofs = (n_neighbors + 2) * 3;
    const int n_neighbor_dofs = n_neighbors * 3;
    assert(n_neighbors > 2);
    assert(otypes.size() == n_neighbors);

    const Eigen::Matrix<double, -1, 3, Eigen::RowMajor> tangents = neighbors.rowwise() - v;
    const Eigen::VectorXd tangents_vec = Eigen::Map<const Eigen::VectorXd>(tangents.data(), tangents.size());
    auto [dn, dn_grad, dn_hess] = normalize_vector_hess(direc);
    dn *= -1;
    dn_grad *= -1;
    for (auto &mat : dn_hess)
        mat *= -1;

    const auto [tangent_term, tangent_grad, tangent_hess] = smooth_point3_term_tangent_hessian(dn, tangents, alpha, beta, otypes);
    const auto [normal_term, normal_grad, normal_hess] = smooth_point3_term_normal_hessian(dn, tangents, alpha, beta, otypes);

    double val = tangent_term * normal_term;

    // gradient wrt. [dn, tangents]
    Eigen::VectorXd grad_tmp = tangent_grad * normal_term + normal_grad * tangent_term;

    // hessian wrt. [dn, tangents]
    Eigen::MatrixXd hess_tmp = tangent_hess * normal_term + normal_hess * tangent_term + 
                            tangent_grad * normal_grad.transpose() + normal_grad * tangent_grad.transpose();

    const double weight = tangents.squaredNorm() / 3.;
    hess_tmp *= weight;
    hess_tmp.bottomRightCorner(n_neighbor_dofs, n_neighbor_dofs).diagonal().array() += 2. / 3. * val;
    Eigen::VectorXd weight_grad(3 + n_neighbors * 3);
    weight_grad << Eigen::Vector3d::Zero(), 2. / 3. * tangents_vec;
    hess_tmp += weight_grad * grad_tmp.transpose() + grad_tmp * weight_grad.transpose();

    grad_tmp *= weight;
    grad_tmp.tail(n_neighbor_dofs) += (2. / 3. * val) * tangents_vec;

    val *= weight;

    // gradient wrt. [direc, v, neighbors]
    Eigen::VectorXd grad(n_dofs);
    grad.head<3>() = dn_grad * grad_tmp.head<3>();
    for (int d = 0; d < 3; d++)
        grad(d + 3) = -grad_tmp(Eigen::seqN(d + 3, n_neighbors, 3)).sum();
    grad.tail(n_neighbor_dofs) = grad_tmp.tail(n_neighbor_dofs);

    // hessian wrt. [direc, v, neighbors]
    Eigen::MatrixXd hess;
    hess.setZero(n_dofs, n_dofs);

    hess.topLeftCorner<3, 3>() = dn_grad * hess_tmp.topLeftCorner<3, 3>() * dn_grad.transpose() + 
                                dn_hess[0] * grad_tmp(0) + dn_hess[1] * grad_tmp(1) + dn_hess[2] * grad_tmp(2);
    hess.bottomRightCorner(n_neighbor_dofs, n_neighbor_dofs) = hess_tmp.bottomRightCorner(n_neighbor_dofs, n_neighbor_dofs);

    hess.block(0, 6, 3, n_neighbor_dofs) = dn_grad * hess_tmp.topRightCorner(3, n_neighbor_dofs);
    hess.block(6, 0, n_neighbor_dofs, 3) = hess_tmp.bottomLeftCorner(n_neighbor_dofs, 3) * dn_grad;

    for (int k = 0, id = 3; k < n_neighbors; k++, id += 3)
    {
        hess.block<3, 3>(0, 3) -= dn_grad * hess_tmp.block<3, 3>(0, id);
        hess.block<3, 3>(3, 0) -= hess_tmp.block<3, 3>(id, 0) * dn_grad;
        for (int l = 0, jd = 3; l < n_neighbors; l++, jd += 3)
        {
            hess.block<3, 3>(3, 3) += hess_tmp.block<3, 3>(id, jd);
            hess.block<3, 3>(3, 3 + id) -= hess_tmp.block<3, 3>(jd, id);
            hess.block<3, 3>(3 + id, 3) -= hess_tmp.block<3, 3>(id, jd);
        }
    }

    return std::tuple(val, grad, hess);
}

double smooth_point3_term_tangent(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    double tangent_term = 1.;
    for (int a = 0; a < tangents.rows(); a++)
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::VARIANT)
            tangent_term *= opposite_direction_penalty(tangents.row(a), direc, alpha, beta);

    return tangent_term;
}

std::tuple<double, Eigen::VectorXd> smooth_point3_term_tangent_gradient(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    const int nn = tangents.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Ones(nn, 1);
    Eigen::VectorXd acc_val_1 = Eigen::VectorXd::Ones(nn, 1);
    std::vector<Vector6d> tmp_grad(nn);
    for (int a = 0; a < nn; a++)
    {
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::VARIANT)
        {
            std::tie(values(a), tmp_grad[a]) = opposite_direction_penalty_grad(tangents.row(a), direc, alpha, beta);
            for (int b = 0; b < nn; b++)
                if (otypes.tangent_type(b) == HEAVISIDE_TYPE::VARIANT && b != a)
                    acc_val_1(b) *= values(a);
        }
    }

    Eigen::VectorXd tangent_grad = Eigen::VectorXd::Zero((nn + 1) * 3);
    for (int a = 0; a < nn; a++)
    {
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::VARIANT)
        {
            const int id = (a + 1) * 3;
            tangent_grad.segment<3>(id) = tmp_grad[a].head<3>() * acc_val_1(a);
            tangent_grad.segment<3>(0) += tmp_grad[a].tail<3>() * acc_val_1(a);
        }
    }

    return std::make_tuple(values.prod(), tangent_grad);
}

std::tuple<double, Eigen::VectorXd, Eigen::MatrixXd> smooth_point3_term_tangent_hessian(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    const int nn = tangents.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Ones(nn, 1);
    Eigen::VectorXd acc_val_1 = Eigen::VectorXd::Ones(nn, 1);
    Eigen::MatrixXd acc_val_2 = Eigen::MatrixXd::Ones(nn, nn);
    std::vector<Vector6d> tmp_grad(nn);
    std::vector<Matrix6d> tmp_hess(nn);
    for (int a = 0; a < nn; a++)
    {
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::VARIANT)
        {
            std::tie(values(a), tmp_grad[a], tmp_hess[a]) = opposite_direction_penalty_hess(tangents.row(a), direc, alpha, beta);
            for (int b = 0; b < nn; b++)
            {
                if (b != a)
                {
                    acc_val_1(b) *= values(a);
                
                    for (int c = 0; c < nn; c++)
                        if (a != c)
                            acc_val_2(b, c) *= values(a);
                }
            }
        }
    }

    Eigen::VectorXd tangent_grad = Eigen::VectorXd::Zero((nn + 1) * 3);
    Eigen::MatrixXd tangent_hess = Eigen::MatrixXd::Zero(tangent_grad.size(), tangent_grad.size());
    Eigen::Vector3d tmp;
    for (int a = 0; a < nn; a++)
    {
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::VARIANT)
        {
            const int id = (a + 1) * 3;
            tangent_grad.segment<3>(id) = tmp_grad[a].head<3>() * acc_val_1(a);
            tangent_grad.segment<3>(0) += tmp_grad[a].tail<3>() * acc_val_1(a);

            tmp.setZero();
            for (int b = 0; b < nn; b++)
                if (otypes.tangent_type(b) == HEAVISIDE_TYPE::VARIANT && b != a)
                {
                    tmp += tmp_grad[b].tail<3>() * acc_val_2(a, b);
                    tangent_hess.block<3, 3>(id, (b + 1) * 3) = tmp_grad[a].head<3>() * tmp_grad[b].head<3>().transpose() * acc_val_2(a, b);
                    tangent_hess.block<3, 3>((b + 1) * 3, id) = tmp_grad[b].head<3>() * tmp_grad[a].head<3>().transpose() * acc_val_2(a, b);
                }

            tangent_hess.block<3, 3>(id, id) = tmp_hess[a].block<3, 3>(0, 0) * acc_val_1(a);
            tangent_hess.block<3, 3>(0, 0)  += tmp_hess[a].block<3, 3>(3, 3) * acc_val_1(a);
            tangent_hess.block<3, 3>(0, 0)  += tmp_grad[a].tail<3>() * tmp.transpose();

            tangent_hess.block<3, 3>(id, 0) += tmp_hess[a].block<3, 3>(0, 3) * acc_val_1(a);
            tangent_hess.block<3, 3>(0, id) += tmp_hess[a].block<3, 3>(3, 0) * acc_val_1(a);
            tangent_hess.block<3, 3>(id, 0) += tmp_grad[a].head<3>() * tmp.transpose();
            tangent_hess.block<3, 3>(0, id) += tmp * tmp_grad[a].head<3>().transpose();
        }
    }

    return std::make_tuple(values.prod(), tangent_grad, tangent_hess);
}

double smooth_point3_term_normal(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE)
        return 1.;
    
    double normal_term = 0.;
    for (int a = 0; a < tangents.rows(); a++)
    {
        const Eigen::Ref<const RowVector3<double>> t_prev = tangents.row((a + tangents.rows() - 1) % tangents.rows());
        const Eigen::Ref<const RowVector3<double>> t = tangents.row(a);

        if (otypes.normal_type(a) == HEAVISIDE_TYPE::VARIANT)
            normal_term += negative_orientation_penalty(t_prev, t, -direc, alpha, beta);
    }

    return Math<double>::smooth_heaviside(normal_term - 1, beta + alpha, 0);
}

std::tuple<double, Eigen::VectorXd> smooth_point3_term_normal_gradient(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    Eigen::VectorXd grad = Eigen::VectorXd::Zero((tangents.rows() + 1) * 3);
    if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE)
        return std::make_tuple(1., grad);
    
    double normal_term = 0.;
    for (int a = 0; a < tangents.rows(); a++)
    {
        const int a_prev = (a + tangents.rows() - 1) % tangents.rows();
        const Eigen::Ref<const RowVector3<double>> t_prev = tangents.row(a_prev);
        const Eigen::Ref<const RowVector3<double>> t = tangents.row(a);

        if (otypes.normal_type(a) == HEAVISIDE_TYPE::VARIANT)
        {
            const int id_prev = (a_prev + 1) * 3;
            const int id = (a + 1) * 3;
            auto [y, dy, ddy] = negative_orientation_penalty_hess(t_prev, t, -direc, alpha, beta);
            
            normal_term += y;
            
            grad.segment<3>(id_prev) += dy.head(3);
            grad.segment<3>(id) += dy.segment(3, 3);
            grad.head<3>() -= dy.tail(3);
        }
    }

    const double val = Math<double>::smooth_heaviside(normal_term - 1, beta + alpha, 0);
    const double grad_val = Math<double>::smooth_heaviside_grad(normal_term - 1, beta + alpha, 0);

    return std::make_tuple(val, grad * grad_val);
}

std::tuple<double, Eigen::VectorXd, Eigen::MatrixXd> smooth_point3_term_normal_hessian(
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Ref<const Eigen::Matrix<double, -1, 3>>& tangents,
    const double &alpha,
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    Eigen::VectorXd grad = Eigen::VectorXd::Zero((tangents.rows() + 1) * 3);
    Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(grad.size(), grad.size());
    if (otypes.normal_type(0) == HEAVISIDE_TYPE::ONE)
        return std::make_tuple(1., grad, hess);
    
    double normal_term = 0.;
    for (int a = 0; a < tangents.rows(); a++)
    {
        const int a_prev = (a + tangents.rows() - 1) % tangents.rows();
        const Eigen::Ref<const RowVector3<double>> t_prev = tangents.row(a_prev);
        const Eigen::Ref<const RowVector3<double>> t = tangents.row(a);

        if (otypes.normal_type(a) == HEAVISIDE_TYPE::VARIANT)
        {
            const int id_prev = (a_prev + 1) * 3;
            const int id = (a + 1) * 3;
            auto [y, dy, ddy] = negative_orientation_penalty_hess(t_prev, t, -direc, alpha, beta);
            
            normal_term += y;
            
            grad.segment<3>(id_prev) += dy.head(3);
            grad.segment<3>(id) += dy.segment(3, 3);
            grad.head<3>() -= dy.tail(3);
            
            hess.block<3, 3>(id_prev, id_prev) += ddy.block<3, 3>(0, 0);
            hess.block<3, 3>(id_prev, id) += ddy.block<3, 3>(0, 3);
            hess.block<3, 3>(id_prev, 0) += ddy.block<3, 3>(0, 6);

            hess.block<3, 3>(id, id_prev) += ddy.block<3, 3>(3, 0);
            hess.block<3, 3>(id, id) += ddy.block<3, 3>(3, 3);
            hess.block<3, 3>(id, 0) += ddy.block<3, 3>(3, 6);

            hess.block<3, 3>(0, id_prev) += ddy.block<3, 3>(6, 0);
            hess.block<3, 3>(0, id) += ddy.block<3, 3>(6, 3);
            hess.block<3, 3>(0, 0) -= ddy.block<3, 3>(6, 6);
        }
    }

    const double val = Math<double>::smooth_heaviside(normal_term - 1, beta + alpha, 0);
    const double grad_val = Math<double>::smooth_heaviside_grad(normal_term - 1, beta + alpha, 0);
    const double hess_val = Math<double>::smooth_heaviside_hess(normal_term - 1, beta + alpha, 0);

    return std::make_tuple(val, grad * grad_val, grad * hess_val * grad.transpose() + grad_val * hess);
}

bool smooth_point2_term_type(
    const Eigen::Ref<const Vector2<double>>& v,
    const Eigen::Ref<const Vector2<double>>& direc,
    const Eigen::Ref<const Vector2<double>>& e0,
    const Eigen::Ref<const Vector2<double>>& e1,
    const double &alpha,
    const double &beta)
{
    const Vector2<double> t0 = (e0 - v).normalized(), t1 = (v - e1).normalized();

    if (direc.dot(t0) <= -alpha || -direc.dot(t1) <= -alpha)
        return false;

    const double tmp = Math<double>::smooth_heaviside(-Math<double>::cross2(direc, t0), alpha, beta) + 
                         Math<double>::smooth_heaviside(-Math<double>::cross2(direc, t1), alpha, beta);
    if (tmp <= 1. - alpha)
        return false;

    return true;
}

bool smooth_point3_term_type(
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Matrix<double, -1, 3> &neighbors,
    const double &alpha,
    const double &beta,
    ORIENTATION_TYPES &otypes)
{
    RowVector3<double> t, t_prev;
    assert(neighbors.rows() > 2);
    otypes.set_size(neighbors.rows());

    const RowVector3<double> dn = direc.normalized();
    double normal_term = 0;
    t_prev = neighbors.row(neighbors.rows()-1) - v;
    for (int a = 0; a < neighbors.rows(); a++)
    {
        t = neighbors.row(a) - v;
        otypes.tangent_type(a) = otypes.compute_type(-dn.dot(t) / t.norm(), alpha, beta);
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::ZERO)
            return false;

        const double tmp = dn.dot(t_prev.cross(t).normalized());
        otypes.normal_type(a) = otypes.compute_type(tmp, alpha, beta);
        normal_term += Math<double>::smooth_heaviside(tmp, alpha, beta);

        std::swap(t, t_prev);
    }

    if (normal_term >= 1)
    {
        for (int a = 0; a < neighbors.rows(); a++)
            otypes.normal_type(a) = HEAVISIDE_TYPE::ONE;
    }

    return normal_term > 1 - alpha;
}
}
