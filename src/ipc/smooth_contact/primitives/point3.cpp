#include "point3.hpp"

#include <ipc/utils/autodiff_types.hpp>

namespace ipc {

Point3::Point3(
    const index_t id,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const VectorMax3d& d,
    const SmoothContactParameters& params)
    : Primitive(id, params)
{
    orientable =
        mesh.is_orient_vertex(id) && mesh.vertices_to_faces()[id].size() > 0;

    // Build index mapping from all vertices to one-ring neighbors
    {
        global_to_local_vids[id] = 0;
        local_to_global_vids.push_back(id);

        int i = 0;
        edges.setZero(mesh.vertices_to_edges()[id].size(), 2);
        for (int eid : mesh.vertices_to_edges()[id]) {
            i++;
            index_t neighbor_id = mesh.edges()(eid, 0) == id
                ? mesh.edges()(eid, 1)
                : mesh.edges()(eid, 0);
            assert(
                global_to_local_vids.find(neighbor_id)
                == global_to_local_vids.end());
            global_to_local_vids[neighbor_id] = i;
            edges.row(i - 1) = Eigen::RowVector2i(0, i);
            local_to_global_vids.push_back(neighbor_id);
        }
    }

    // If the vertex is on some faces, record the local indices of the face
    // vertices
    if (mesh.vertices_to_faces()[id].size() > 0) {
        int i = 0;
        faces.setZero(mesh.vertices_to_faces()[id].size(), 3);
        for (auto f : mesh.vertices_to_faces()[id]) {
            int lv = 0;
            for (; lv < 3; lv++) {
                if (mesh.faces()(f, lv) == id) {
                    faces.row(i) << 0,
                        global_to_local_vids[mesh.faces()(f, (lv + 1) % 3)],
                        global_to_local_vids[mesh.faces()(f, (lv + 2) % 3)];
                    break;
                }
            }
            assert(lv < 3);
            i++;
        }
    }

    n_neighbors = local_to_global_vids.size() - 1;
    m_vertex_ids = local_to_global_vids;

    if (m_vertex_ids.size() > N_VERT_NEIGHBORS_3D)
        logger().error(
            "Too many neighbors for point3 primitive! {} > {}! Increase N_VERT_NEIGHBORS_3D in common.hpp",
            m_vertex_ids.size(), N_VERT_NEIGHBORS_3D);

    m_is_active =
        smooth_point3_term_type(vertices(local_to_global_vids, Eigen::all), d);
}

int Point3::n_vertices() const { return local_to_global_vids.size(); }

double Point3::potential(
    const Vector<double, DIM>& d, const Vector<double, -1, MAX_SIZE>& x) const
{
    const Eigen::Matrix<double, -1, DIM> X =
        slice_positions<double, -1, DIM>(x);
    return smooth_point3_term<double, -1>(X, d);
}

Vector<double, -1, Point3::MAX_SIZE + Point3::DIM> Point3::grad(
    const Vector<double, DIM>& d, const Vector<double, -1, MAX_SIZE>& x) const
{
#ifdef DERIVATIVES_WITH_AUTODIFF
    using T = ADGrad<-1>;
    Eigen::VectorXd tmp(x.size() + d.size());
    tmp << d, x;
    DiffScalarBase::setVariableCount(tmp.size());
    const Eigen::Matrix<T, -1, DIM> X = slice_positions<T, -1, DIM>(tmp);
    return smooth_point3_term<T, -1>(X.bottomRows(X.rows() - 1), X.row(0))
        .getGradient();
#else
    const Eigen::Matrix<double, -1, DIM> X =
        slice_positions<double, -1, DIM>(x);
    const auto [val, grad] = smooth_point3_term_gradient(d, X, params);
    return grad;
#endif
}

MatrixMax<
    double,
    Point3::MAX_SIZE + Point3::DIM,
    Point3::MAX_SIZE + Point3::DIM>
Point3::hessian(
    const Vector<double, DIM>& d, const Vector<double, -1, MAX_SIZE>& x) const
{
#ifdef DERIVATIVES_WITH_AUTODIFF
    using T = ADHessian<-1>;
    Eigen::VectorXd tmp(x.size() + d.size());
    tmp << d, x;
    DiffScalarBase::setVariableCount(tmp.size());
    const Eigen::Matrix<T, -1, DIM> X = slice_positions<T, -1, DIM>(tmp);
    return smooth_point3_term<T, -1>(X.bottomRows(X.rows() - 1), X.row(0))
        .getHessian();
#else
    const auto X = slice_positions<double, -1, DIM>(x);
    const auto [val, grad, hess] = smooth_point3_term_hessian(d, X, params);
    return hess;
#endif
}

GradType<-1> Point3::smooth_point3_term_tangent_gradient(
    Eigen::ConstRef<RowVector3<double>> direc,
    Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
    const double& alpha,
    const double& beta) const
{
    const int nn = tangents.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Ones(nn, 1);
    Eigen::VectorXd acc_val_1 = Eigen::VectorXd::Ones(nn, 1);
    std::vector<Vector6d> tmp_grad(nn);
    for (int a = 0; a < nn; a++) {
        if (otypes.tangent_type(a) == HeavisideType::VARIANT) {
            std::tie(values(a), tmp_grad[a]) = opposite_direction_penalty_grad(
                tangents.row(a), direc, alpha, beta);
            for (int b = 0; b < nn; b++)
                if (otypes.tangent_type(b) == HeavisideType::VARIANT && b != a)
                    acc_val_1(b) *= values(a);
        }
    }

    Eigen::VectorXd tangent_grad = Eigen::VectorXd::Zero((nn + 1) * 3);
    for (int a = 0; a < nn; a++) {
        if (otypes.tangent_type(a) == HeavisideType::VARIANT) {
            const int id = (a + 1) * 3;
            tangent_grad.segment<3>(id) = tmp_grad[a].head<3>() * acc_val_1(a);
            tangent_grad.segment<3>(0) += tmp_grad[a].tail<3>() * acc_val_1(a);
        }
    }

    return std::make_tuple(values.prod(), tangent_grad);
}

HessianType<-1> Point3::smooth_point3_term_tangent_hessian(
    Eigen::ConstRef<RowVector3<double>> direc,
    Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
    const double& alpha,
    const double& beta) const
{
    const int nn = tangents.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Ones(nn, 1);
    Eigen::VectorXd acc_val_1 = Eigen::VectorXd::Ones(nn, 1);
    Eigen::MatrixXd acc_val_2 = Eigen::MatrixXd::Ones(nn, nn);
    std::vector<Vector6d> tmp_grad(nn);
    std::vector<Matrix6d> tmp_hess(nn);
    for (int a = 0; a < nn; a++) {
        if (otypes.tangent_type(a) == HeavisideType::VARIANT) {
            std::tie(values(a), tmp_grad[a], tmp_hess[a]) =
                opposite_direction_penalty_hess(
                    tangents.row(a), direc, alpha, beta);
            for (int b = 0; b < nn; b++) {
                if (b != a) {
                    acc_val_1(b) *= values(a);

                    for (int c = 0; c < nn; c++)
                        if (a != c)
                            acc_val_2(b, c) *= values(a);
                }
            }
        }
    }

    Eigen::VectorXd tangent_grad = Eigen::VectorXd::Zero((nn + 1) * 3);
    Eigen::MatrixXd tangent_hess =
        Eigen::MatrixXd::Zero(tangent_grad.size(), tangent_grad.size());
    Eigen::Vector3d tmp;
    for (int a = 0; a < nn; a++) {
        if (otypes.tangent_type(a) == HeavisideType::VARIANT) {
            const int id = (a + 1) * 3;
            tangent_grad.segment<3>(id) = tmp_grad[a].head<3>() * acc_val_1(a);
            tangent_grad.segment<3>(0) += tmp_grad[a].tail<3>() * acc_val_1(a);

            tmp.setZero();
            for (int b = 0; b < nn; b++)
                if (otypes.tangent_type(b) == HeavisideType::VARIANT
                    && b != a) {
                    tmp += tmp_grad[b].tail<3>() * acc_val_2(a, b);
                    tangent_hess.block<3, 3>(id, (b + 1) * 3) =
                        tmp_grad[a].head<3>()
                        * tmp_grad[b].head<3>().transpose() * acc_val_2(a, b);
                    tangent_hess.block<3, 3>((b + 1) * 3, id) =
                        tmp_grad[b].head<3>()
                        * tmp_grad[a].head<3>().transpose() * acc_val_2(a, b);
                }

            tangent_hess.block<3, 3>(id, id) =
                tmp_hess[a].block<3, 3>(0, 0) * acc_val_1(a);
            tangent_hess.block<3, 3>(0, 0) +=
                tmp_hess[a].block<3, 3>(3, 3) * acc_val_1(a);
            tangent_hess.block<3, 3>(0, 0) +=
                tmp_grad[a].tail<3>() * tmp.transpose();

            tangent_hess.block<3, 3>(id, 0) +=
                tmp_hess[a].block<3, 3>(0, 3) * acc_val_1(a);
            tangent_hess.block<3, 3>(0, id) +=
                tmp_hess[a].block<3, 3>(3, 0) * acc_val_1(a);
            tangent_hess.block<3, 3>(id, 0) +=
                tmp_grad[a].head<3>() * tmp.transpose();
            tangent_hess.block<3, 3>(0, id) +=
                tmp * tmp_grad[a].head<3>().transpose();
        }
    }

    return std::make_tuple(values.prod(), tangent_grad, tangent_hess);
}

GradType<-1> Point3::smooth_point3_term_normal_gradient(
    Eigen::ConstRef<RowVector3<double>> direc,
    Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
    const double& alpha,
    const double& beta) const
{
    Eigen::VectorXd grad =
        Eigen::VectorXd::Zero(tangents.size() + direc.size());
    if (!orientable || otypes.normal_type(0) == HeavisideType::ONE)
        return std::make_tuple(1., grad);

    double normal_term = 0.;
    for (int a = 0; a < faces.rows(); a++) {
        const Eigen::Ref<const RowVector3<double>> t1 =
            tangents.row(faces(a, 1) - 1);
        const Eigen::Ref<const RowVector3<double>> t2 =
            tangents.row(faces(a, 2) - 1);

        if (otypes.normal_type(a) == HeavisideType::VARIANT) {
            const int id1 = faces(a, 1) * 3;
            const int id2 = faces(a, 2) * 3;
            const auto [y, dy] =
                negative_orientation_penalty_grad(t1, t2, -direc, alpha, beta);

            normal_term += y;

            grad.segment<3>(id1) += dy.head<3>();
            grad.segment<3>(id2) += dy.segment<3>(3);
            grad.head<3>() -= dy.tail<3>();
        }
    }

    // autodiff
    // TODO: replace with efficient code
    // double normal_term = 0;
    // {
    //     using T = ADGrad<-1>;
    //     Eigen::VectorXd tmp(direc.size() + tangents.size());
    //     tmp.head<3>() = direc;
    //     for (int i = 0; i < tangents.rows(); i++)
    //         tmp.segment<3>(3 * i + 3) = tangents.row(i);

    //     DiffScalarBase::setVariableCount(tmp.size());
    //     const Eigen::Matrix<T, -1, DIM> X = slice_positions<T, -1, DIM>(tmp);

    //     T normal_term_ad(0.);
    //     for (int a = 0; a < faces.rows(); a++) {
    //         if (otypes.normal_type(a) == HeavisideType::VARIANT)
    //             normal_term_ad =
    //                 normal_term_ad
    //                 + Math<T>::smooth_heaviside(
    //                     -X.row(0).dot(X.row(faces(a, 1))
    //                                       .cross(X.row(faces(a, 2)))
    //                                       .normalized()),
    //                     params.alpha_n, params.beta_n);
    //     }

    //     normal_term = normal_term_ad.getValue();
    //     grad = normal_term_ad.getGradient();
    // }

    const double val = Math<double>::smooth_heaviside(normal_term - 1, 1., 0);
    const double grad_val =
        Math<double>::smooth_heaviside_grad(normal_term - 1, 1., 0);

    return std::make_tuple(val, grad * grad_val);
}

HessianType<-1> Point3::smooth_point3_term_normal_hessian(
    Eigen::ConstRef<RowVector3<double>> direc,
    Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> tangents,
    const double& alpha,
    const double& beta) const
{
    Eigen::VectorXd grad = Eigen::VectorXd::Zero((tangents.rows() + 1) * 3);
    Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(grad.size(), grad.size());
    if (!orientable || otypes.normal_type(0) == HeavisideType::ONE)
        return std::make_tuple(1., grad, hess);

    // double normal_term = 0.;
    // for (int a = 0; a < faces.rows(); a++) {
    //     const Eigen::Ref<const RowVector3<double>> t1 = tangents.row(faces(a,
    //     1) - 1); const Eigen::Ref<const RowVector3<double>> t2 =
    //     tangents.row(faces(a, 2) - 1);

    //     if (otypes.normal_type(a) == HeavisideType::VARIANT)
    //     {
    //         const int id1 = faces(a, 1) * 3;
    //         const int id2 = faces(a, 2) * 3;
    //         auto [y, dy, ddy] = negative_orientation_penalty_hess(
    //             t1, t2, -direc, alpha, beta);

    //         normal_term += y;

    //         grad.segment<3>(id1) += dy.head(3);
    //         grad.segment<3>(id2) += dy.segment<3>(3);
    //         grad.head<3>() -= dy.tail<3>();

    //         hess.block<3, 3>(id1, id1) += ddy.block<3, 3>(0, 0);
    //         hess.block<3, 3>(id1, id2) += ddy.block<3, 3>(0, 3);
    //         hess.block<3, 3>(id1, 0) += ddy.block<3, 3>(0, 6);

    //         hess.block<3, 3>(id2, id1) += ddy.block<3, 3>(3, 0);
    //         hess.block<3, 3>(id2, id2) += ddy.block<3, 3>(3, 3);
    //         hess.block<3, 3>(id2, 0) += ddy.block<3, 3>(3, 6);

    //         hess.block<3, 3>(0, id1) += ddy.block<3, 3>(6, 0);
    //         hess.block<3, 3>(0, id2) += ddy.block<3, 3>(6, 3);
    //         hess.block<3, 3>(0, 0) -= ddy.block<3, 3>(6, 6);
    //     }
    // }

    // autodiff
    // TODO: replace with efficient code
    double normal_term = 0;
    {
        using T = ADHessian<-1>;
        Eigen::VectorXd tmp(direc.size() + tangents.size());
        tmp.head<3>() = direc;
        for (int i = 0; i < tangents.rows(); i++)
            tmp.segment<3>(3 * i + 3) = tangents.row(i);

        DiffScalarBase::setVariableCount(tmp.size());
        const Eigen::Matrix<T, -1, DIM> X = slice_positions<T, -1, DIM>(tmp);

        T normal_term_ad(0.);
        for (int a = 0; a < faces.rows(); a++) {
            if (otypes.normal_type(a) == HeavisideType::VARIANT)
                normal_term_ad =
                    normal_term_ad
                    + Math<T>::smooth_heaviside(
                        -X.row(0).dot(X.row(faces(a, 1))
                                          .cross(X.row(faces(a, 2)))
                                          .normalized()),
                        params.alpha_n, params.beta_n);
        }

        normal_term = normal_term_ad.getValue();
        hess = normal_term_ad.getHessian();
        grad = normal_term_ad.getGradient();
    }

    const double val = Math<double>::smooth_heaviside(normal_term - 1, 1., 0);
    const double grad_val =
        Math<double>::smooth_heaviside_grad(normal_term - 1, 1., 0);
    const double hess_val =
        Math<double>::smooth_heaviside_hess(normal_term - 1, 1., 0);

    return std::make_tuple(
        val, grad * grad_val,
        grad * hess_val * grad.transpose() + grad_val * hess);
}

bool Point3::smooth_point3_term_type(
    const Eigen::Matrix<double, -1, 3>& X,
    Eigen::ConstRef<RowVector3<double>> direc)
{
    otypes.set_size(edges.rows());

    const RowVector3<double> dn = direc.normalized();
    for (int a = 0; a < edges.rows(); a++) {
        const RowVector3<double> t = X.row(edges(a, 1)) - X.row(edges(a, 0));
        otypes.tangent_type(a) = otypes.compute_type(
            -dn.dot(t) / t.norm(), params.alpha_t, params.beta_t);
        if (otypes.tangent_type(a) == HeavisideType::ZERO)
            return false;
    }

    if (!orientable) {
        return true;
    }

    double normal_term = 0;
    for (int a = 0; a < faces.rows(); a++) {
        const RowVector3<double> t1 = X.row(faces(a, 1)) - X.row(faces(a, 0));
        const RowVector3<double> t2 = X.row(faces(a, 2)) - X.row(faces(a, 0));
        const double tmp = dn.dot(t1.cross(t2).normalized());
        otypes.normal_type(a) =
            otypes.compute_type(tmp, params.alpha_n, params.beta_n);
        normal_term +=
            Math<double>::smooth_heaviside(tmp, params.alpha_n, params.beta_n);
    }

    if (normal_term >= 1) {
        for (int a = 0; a < faces.rows(); a++)
            otypes.normal_type(a) = HeavisideType::ONE;
    }

    return normal_term > 0;
}

GradType<-1> Point3::smooth_point3_term_gradient(
    Eigen::ConstRef<RowVector3<double>> direc,
    Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> X,
    const SmoothContactParameters& params) const
{
    const int n_dofs = (X.rows() + 1) * 3;
    const int n_neighbor_dofs = n_neighbors * 3;

    const Eigen::Matrix<double, -1, 3, Eigen::RowMajor> tangents =
        X.bottomRows(n_neighbors).rowwise() - X.row(0);
    const Eigen::VectorXd tangents_vec =
        Eigen::Map<const Eigen::VectorXd>(tangents.data(), tangents.size());

    auto [dn, dn_grad] = normalize_vector_grad(direc);
    dn *= -1;
    dn_grad *= -1;

    const auto [tangent_term, tangent_grad] =
        smooth_point3_term_tangent_gradient(
            dn, tangents, params.alpha_t, params.beta_t);

    auto [normal_term, normal_grad] = smooth_point3_term_normal_gradient(
        dn, tangents, params.alpha_n, params.beta_n);

    double val = tangent_term * normal_term;

    // gradient wrt. [dn, tangents]
    Eigen::VectorXd grad_tmp =
        tangent_grad * normal_term + normal_grad * tangent_term;

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

HessianType<-1> Point3::smooth_point3_term_hessian(
    Eigen::ConstRef<RowVector3<double>> direc,
    Eigen::ConstRef<Eigen::Matrix<double, -1, 3>> X,
    const SmoothContactParameters& params) const
{
    const int n_dofs = (X.rows() + 1) * 3;
    const int n_neighbor_dofs = n_neighbors * 3;

    const Eigen::Matrix<double, -1, 3, Eigen::RowMajor> tangents =
        X.bottomRows(n_neighbors).rowwise() - X.row(0);
    const Eigen::VectorXd tangents_vec =
        Eigen::Map<const Eigen::VectorXd>(tangents.data(), tangents.size());

    auto [dn, dn_grad, dn_hess] = normalize_vector_hess(direc);
    dn *= -1;
    dn_grad *= -1;
    for (auto& mat : dn_hess)
        mat *= -1;

    const auto [tangent_term, tangent_grad, tangent_hess] =
        smooth_point3_term_tangent_hessian(
            dn, tangents, params.alpha_t, params.beta_t);

    auto [normal_term, normal_grad, normal_hess] =
        smooth_point3_term_normal_hessian(
            dn, tangents, params.alpha_n, params.beta_n);

    double val = tangent_term * normal_term;

    // gradient wrt. [dn, tangents]
    Eigen::VectorXd grad_tmp =
        tangent_grad * normal_term + normal_grad * tangent_term;

    // hessian wrt. [dn, tangents]
    Eigen::MatrixXd hess_tmp = tangent_hess * normal_term
        + normal_hess * tangent_term + tangent_grad * normal_grad.transpose()
        + normal_grad * tangent_grad.transpose();

    const double weight = tangents.squaredNorm() / 3.;
    hess_tmp *= weight;
    hess_tmp.bottomRightCorner(n_neighbor_dofs, n_neighbor_dofs)
        .diagonal()
        .array() += 2. / 3. * val;
    hess_tmp.bottomRows(n_neighbor_dofs) +=
        2. / 3. * tangents_vec * grad_tmp.transpose();
    hess_tmp.rightCols(n_neighbor_dofs) +=
        2. / 3. * grad_tmp * tangents_vec.transpose();

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

    hess.topLeftCorner<3, 3>() =
        dn_grad * hess_tmp.topLeftCorner<3, 3>() * dn_grad.transpose()
        + dn_hess[0] * grad_tmp(0) + dn_hess[1] * grad_tmp(1)
        + dn_hess[2] * grad_tmp(2);
    hess.bottomRightCorner(n_neighbor_dofs, n_neighbor_dofs) =
        hess_tmp.bottomRightCorner(n_neighbor_dofs, n_neighbor_dofs);

    hess.block(0, 6, 3, n_neighbor_dofs) =
        dn_grad * hess_tmp.topRightCorner(3, n_neighbor_dofs);
    hess.block(6, 0, n_neighbor_dofs, 3) =
        hess_tmp.bottomLeftCorner(n_neighbor_dofs, 3) * dn_grad;

    for (int k = 0, id = 3; k < n_neighbors; k++, id += 3) {
        hess.block<3, 3>(0, 3) -= dn_grad * hess_tmp.block<3, 3>(0, id);
        hess.block<3, 3>(3, 0) -= hess_tmp.block<3, 3>(id, 0) * dn_grad;
        for (int l = 0, jd = 3; l < n_neighbors; l++, jd += 3) {
            hess.block<3, 3>(3, 3) += hess_tmp.block<3, 3>(id, jd);
            hess.block<3, 3>(3, 3 + id) -= hess_tmp.block<3, 3>(jd, id);
            hess.block<3, 3>(3 + id, 3) -= hess_tmp.block<3, 3>(id, jd);
        }
    }

    return std::tuple(val, grad, hess);
}

template <typename scalar, int n_verts>
scalar Point3::smooth_point3_term(
    const Eigen::Matrix<scalar, n_verts, 3>& X,
    Eigen::ConstRef<RowVector3<scalar>> direc) const
{
    const RowVector3<scalar> dn = direc.normalized();
    scalar tangent_term(1.);
    scalar weight(0.);
    scalar normal_term(0.);
    for (int a = 0; a < edges.rows(); a++) {
        const RowVector3<scalar> t = X.row(edges(a, 1)) - X.row(edges(a, 0));
        if (otypes.tangent_type(a) == HeavisideType::VARIANT)
            tangent_term =
                tangent_term
                * Math<scalar>::smooth_heaviside(
                    -dn.dot(t) / t.norm(), params.alpha_t, params.beta_t);

        weight = weight + t.squaredNorm();
    }
    weight /= 3.;

    if (!orientable || otypes.normal_type(0) == HeavisideType::ONE)
        normal_term = scalar(1.);
    else {
        for (int a = 0; a < faces.rows(); a++) {
            const RowVector3<scalar> t1 =
                X.row(faces(a, 1)) - X.row(faces(a, 0));
            const RowVector3<scalar> t2 =
                X.row(faces(a, 2)) - X.row(faces(a, 0));
            normal_term = normal_term
                + Math<scalar>::smooth_heaviside(
                              dn.dot(t1.cross(t2).normalized()), params.alpha_n,
                              params.beta_n);
        }
        normal_term = Math<scalar>::smooth_heaviside(normal_term - 1, 1., 0);
    }

    return weight * normal_term * tangent_term;
}
} // namespace ipc
