#include "edge3.hpp"

#include <ipc/config.hpp>
#include <ipc/smooth_contact/primitives/autogen.hpp>
#include <ipc/utils/autodiff_types.hpp>

namespace ipc {
namespace {

    template <typename scalar>
    Eigen::Vector3<scalar> project(
        Eigen::ConstRef<Eigen::Vector3<scalar>> a,
        Eigen::ConstRef<Eigen::Vector3<scalar>> b)
    {
        return a - a.dot(b) * b;
    }

    bool smooth_edge3_term_type(
        Eigen::ConstRef<Eigen::Vector3d> dn,
        Eigen::ConstRef<Eigen::Vector3d> e0,
        Eigen::ConstRef<Eigen::Vector3d> e1,
        Eigen::ConstRef<Eigen::Vector3d> f0,
        Eigen::ConstRef<Eigen::Vector3d> f1,
        const SmoothContactParameters& params,
        OrientationTypes& otypes,
        const bool orientable)
    {
        otypes.set_size(2);

        const Eigen::Vector3d t0 =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f0, e0, e1)
                .normalized();
        const Eigen::Vector3d t1 =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f1, e0, e1)
                .normalized();

        {
            otypes.tangent_type(0) = OrientationTypes::compute_type(
                -dn.dot(t0), params.alpha_t, params.beta_t);
            otypes.tangent_type(1) = OrientationTypes::compute_type(
                -dn.dot(t1), params.alpha_t, params.beta_t);
            if (otypes.tangent_type(0) == HeavisideType::ZERO
                || otypes.tangent_type(1) == HeavisideType::ZERO) {
                return false;
            }
        }

        if (orientable) {
            const Eigen::Vector3d edge = (e0 - e1).normalized();
            const Eigen::Vector3d d = project<double>(dn, edge).normalized();
            otypes.normal_type(0) = OrientationTypes::compute_type(
                (d - t0).cross(d - t1).dot(edge), params.alpha_n,
                params.beta_n);
            otypes.normal_type(1) = otypes.normal_type(0);
        }

        return true;
    }

    double smooth_edge3_tangent_term(
        Eigen::ConstRef<Eigen::Vector3d> dn,
        Eigen::ConstRef<Eigen::Vector3d> e0,
        Eigen::ConstRef<Eigen::Vector3d> e1,
        Eigen::ConstRef<Eigen::Vector3d> f0,
        Eigen::ConstRef<Eigen::Vector3d> f1,
        const double alpha,
        const double beta,
        const OrientationTypes& otypes)
    {
        double tangent_term = 1.;
        if (otypes.tangent_type(0) != HeavisideType::ONE) {
            const Eigen::Vector3d t0 = PointEdgeDistance<
                double, 3>::point_line_closest_point_direction(f0, e0, e1);
            tangent_term *= opposite_direction_penalty(t0, -dn, alpha, beta);
        }
        if (otypes.tangent_type(1) != HeavisideType::ONE) {
            const Eigen::Vector3d t1 = PointEdgeDistance<
                double, 3>::point_line_closest_point_direction(f1, e0, e1);
            tangent_term *= opposite_direction_penalty(t1, -dn, alpha, beta);
        }

        return tangent_term;
    }

    GradType<15> smooth_edge3_tangent_term_gradient(
        Eigen::ConstRef<Eigen::Vector3d> dn,
        Eigen::ConstRef<Eigen::Vector3d> e0,
        Eigen::ConstRef<Eigen::Vector3d> e1,
        Eigen::ConstRef<Eigen::Vector3d> f0,
        Eigen::ConstRef<Eigen::Vector3d> f1,
        const double alpha,
        const double beta,
        const OrientationTypes& otypes)
    {
        Eigen::Vector2d vals;
        vals << 1., 1.;
        std::array<Vector<double, 15>, 2> grads;
        for (auto& g : grads) {
            g.setZero();
        }

        for (int d : { 0, 1 }) {
            const Eigen::Ref<const Eigen::Vector3d> f = d == 0 ? f0 : f1;
            if (otypes.tangent_type(d) != HeavisideType::ONE) {
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
                if (d == 1) {
                    indices.segment<3>(3).array() += 3;
                }

                grads[d](indices) = gradient_tmp;
            }
        }

        return std::make_tuple(
            vals.prod(), vals[0] * grads[1] + vals[1] * grads[0]);
    }

    HessianType<15> smooth_edge3_tangent_term_hessian(
        Eigen::ConstRef<Eigen::Vector3d> dn,
        Eigen::ConstRef<Eigen::Vector3d> e0,
        Eigen::ConstRef<Eigen::Vector3d> e1,
        Eigen::ConstRef<Eigen::Vector3d> f0,
        Eigen::ConstRef<Eigen::Vector3d> f1,
        const double alpha,
        const double beta,
        const OrientationTypes& otypes)
    {
        Eigen::Vector2d vals;
        vals << 1., 1.;
        std::array<Vector<double, 15>, 2> grads;
        std::array<Eigen::Matrix<double, 15, 15>, 2> hesses;
        for (auto& g : grads) {
            g.setZero();
        }
        for (auto& h : hesses) {
            h.setZero();
        }

        for (int d : { 0, 1 }) {
            const Eigen::Ref<const Eigen::Vector3d> f = d == 0 ? f0 : f1;
            if (otypes.tangent_type(d) != HeavisideType::ONE) {
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
                    + tmp_grad(0) * h[0] + tmp_grad(1) * h[1]
                    + tmp_grad(2) * h[2];

                // mixed 2nd derivatives
                hessian_tmp.block<3, 9>(0, 3) = tmp_hess.block<3, 3>(3, 0) * g;
                hessian_tmp.block<9, 3>(3, 0) =
                    g.transpose() * tmp_hess.block<3, 3>(0, 3);

                Vector<int, 12> indices;
                indices << 0, 1, 2, 9, 10, 11, 3, 4, 5, 6, 7, 8;
                if (d == 1) {
                    indices.segment<3>(3).array() += 3;
                }

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
        Eigen::ConstRef<Eigen::Vector3d> direc,
        Eigen::ConstRef<Eigen::Vector3d> e0,
        Eigen::ConstRef<Eigen::Vector3d> e1,
        Eigen::ConstRef<Eigen::Vector3d> f0,
        Eigen::ConstRef<Eigen::Vector3d> f1,
        const SmoothContactParameters& params,
        const OrientationTypes& otypes,
        const bool orientable)
    {
        const Eigen::Vector3d dn = direc.normalized();
        double tangent_term = smooth_edge3_tangent_term(
            dn, e0, e1, f0, f1, params.alpha_t, params.beta_t, otypes);
        double normal_term = !orientable
            ? 1
            : smooth_edge3_normal_term(
                  dn, e0, e1, f0, f1, params.alpha_n, params.beta_n, otypes);

        return (e1 - e0).squaredNorm() * tangent_term * normal_term;
    }

    GradType<15> smooth_edge3_term_gradient(
        Eigen::ConstRef<Eigen::Vector3d> direc,
        Eigen::ConstRef<Eigen::Vector3d> e0,
        Eigen::ConstRef<Eigen::Vector3d> e1,
        Eigen::ConstRef<Eigen::Vector3d> f0,
        Eigen::ConstRef<Eigen::Vector3d> f1,
        const SmoothContactParameters& params,
        const OrientationTypes& otypes,
        const bool orientable)
    {
        assert(otypes.size() == 2);

        const auto [dn, dn_grad] = normalize_vector_grad(direc);

        auto [t_term, t_grad] = smooth_edge3_tangent_term_gradient(
            dn, e0, e1, f0, f1, params.alpha_t, params.beta_t, otypes);

        double n_term = 1.;
        Vector15d n_grad = Vector15d::Zero();
        if (orientable) {
            std::tie(n_term, n_grad) = smooth_edge3_normal_term_gradient(
                dn, e0, e1, f0, f1, params.alpha_n, params.beta_n, otypes);
        }

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

    HessianType<15> smooth_edge3_term_hessian(
        Eigen::ConstRef<Eigen::Vector3d> direc,
        Eigen::ConstRef<Eigen::Vector3d> e0,
        Eigen::ConstRef<Eigen::Vector3d> e1,
        Eigen::ConstRef<Eigen::Vector3d> f0,
        Eigen::ConstRef<Eigen::Vector3d> f1,
        const SmoothContactParameters& params,
        const OrientationTypes& otypes,
        const bool orientable)
    {
        assert(otypes.size() == 2);

        const auto [dn, dn_grad, dn_hess] = normalize_vector_hess(direc);

        auto [t_term, t_grad, t_hess] = smooth_edge3_tangent_term_hessian(
            dn, e0, e1, f0, f1, params.alpha_t, params.beta_t, otypes);

        double n_term = 1.;
        Vector15d n_grad = Vector15d::Zero();
        Matrix15d n_hess = Matrix15d::Zero();
        if (orientable) {
            std::tie(n_term, n_grad, n_hess) = smooth_edge3_normal_term_hessian(
                dn, e0, e1, f0, f1, params.alpha_n, params.beta_n, otypes);
        }

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

        Matrix15d hessian = (weight * t_term) * n_hess
            + (weight * n_term) * t_hess
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

    /// @brief
    /// @tparam scalar
    /// @param dn from edge to point outside, normalized
    /// @param e0
    /// @param e1
    /// @param f0 face [f0, e0, e1]
    /// @param f1 face [f1, e1, e0]
    /// @param alpha
    /// @return
    template <typename scalar>
    scalar smooth_edge3_term_template(
        Eigen::ConstRef<Eigen::Vector3<scalar>> dn,
        Eigen::ConstRef<Eigen::Vector3<scalar>> e0,
        Eigen::ConstRef<Eigen::Vector3<scalar>> e1,
        Eigen::ConstRef<Eigen::Vector3<scalar>> f0,
        Eigen::ConstRef<Eigen::Vector3<scalar>> f1,
        const SmoothContactParameters& params,
        const OrientationTypes& otypes,
        const bool orientable)
    {
        scalar tangent_term = scalar(1.);
        const Eigen::Vector3<scalar> t0 =
            PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                f0, e0, e1)
                .normalized();
        if (otypes.tangent_type(0) != HeavisideType::ONE) {
            tangent_term = tangent_term
                * Math<scalar>::smooth_heaviside(
                               -dn.dot(t0), params.alpha_t, params.beta_t);
        }

        const Eigen::Vector3<scalar> t1 =
            PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                f1, e0, e1)
                .normalized();

        if (otypes.tangent_type(1) != HeavisideType::ONE) {
            tangent_term = tangent_term
                * Math<scalar>::smooth_heaviside(
                               -dn.dot(t1), params.alpha_t, params.beta_t);
        }

        scalar normal_term = scalar(1.);
        if (orientable && otypes.normal_type(0) != HeavisideType::ONE) {
            const Eigen::Vector3<scalar> edge = (e0 - e1).normalized();
            const Eigen::Vector3<scalar> d =
                project<scalar>(dn, edge).normalized();
            normal_term = Math<scalar>::smooth_heaviside(
                (d - t0).cross(d - t1).dot(edge), params.alpha_n,
                params.beta_n);
        }

        return (e1 - e0).squaredNorm() * tangent_term * normal_term;
    }

} // namespace
// d is a vector from any point on the edge to the point outside of the edge
Edge3::Edge3(
    const index_t id,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const VectorMax3d& d,
    const SmoothContactParameters& params)
    : Primitive(id, params)
{
    orientable =
        (mesh.is_orient_vertex(mesh.edges()(id, 0))
         && mesh.is_orient_vertex(mesh.edges()(id, 0)));

    std::array<index_t, 4> neighbors { { -1, -1, -1, -1 } };
    {
        std::array<index_t, 2> lf = { { mesh.edges_to_faces()(id, 0),
                                        mesh.edges_to_faces()(id, 1) } };
        has_neighbor_1 = lf[0] >= 0;
        has_neighbor_2 = lf[1] >= 0;

        if (has_neighbor_1) {
            const int j = 0;
            for (int i = 0; i < 3; i++) {
                const auto va = mesh.faces()(lf[j], i);

                if (va != mesh.edges()(id, 0) && va != mesh.edges()(id, 1)) {
                    neighbors[2 + j] = va;

                    // When there is a face adjacent to the edge, we orient the
                    // edge the same way as in the face
                    neighbors[0] = mesh.faces()(lf[j], (i + 1) % 3);
                    neighbors[1] = mesh.faces()(lf[j], (i + 2) % 3);
                    break;
                }
            }
            assert(neighbors[2] >= 0 && neighbors[0] >= 0 && neighbors[1] >= 0);
        }

        if (has_neighbor_2) {
            const int j = 1;
            for (int i = 0; i < 3; i++) {
                const auto va = mesh.faces()(lf[j], i);

                if (va != mesh.edges()(id, 0) && va != mesh.edges()(id, 1)) {
                    neighbors[2 + j] = va;

                    // When has_neighbor_1==false, we orient the edge the same
                    // way as in the second face.
                    if (!has_neighbor_1) {
                        neighbors[0] = mesh.faces()(lf[j], (i + 2) % 3);
                        neighbors[1] = mesh.faces()(lf[j], (i + 1) % 3);
                    } else {
                        assert(
                            neighbors[0] == mesh.faces()(lf[j], (i + 2) % 3));
                        assert(
                            neighbors[1] == mesh.faces()(lf[j], (i + 1) % 3));
                    }
                    break;
                }
            }
            assert(neighbors[3] >= 0);
        }
    }

    if (has_neighbor_1 && has_neighbor_2) {
        m_vertex_ids = std::vector<index_t>(
            neighbors.begin(), neighbors.begin() + neighbors.size());
        m_is_active = smooth_edge3_term_type(
            d.normalized(), vertices.row(m_vertex_ids[0]),
            vertices.row(m_vertex_ids[1]), vertices.row(m_vertex_ids[2]),
            vertices.row(m_vertex_ids[3]), params, otypes, orientable);
    } else if (has_neighbor_1 || has_neighbor_2) {
        m_vertex_ids = { { neighbors[0], neighbors[1],
                           has_neighbor_1 ? neighbors[2] : neighbors[3] } };

    } else {
        m_vertex_ids = { { neighbors[0], neighbors[1] } };
        m_is_active = true;
    }
}

int Edge3::n_vertices() const { return N_EDGE_NEIGHBORS_3D; }

double Edge3::potential(
    Eigen::ConstRef<Eigen::Vector3d> d, Eigen::ConstRef<Vector12d> x) const
{
#ifdef IPC_TOOLKIT_DEBUG_AUTODIFF
    return smooth_edge3_term_template<double>(
        d.normalized(), x.head<3>(), x.segment<3>(3), x.segment<3>(6),
        x.tail<3>(), m_params, otypes, orientable);
#else
    return smooth_edge3_term(
        d, x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), m_params,
        otypes, orientable);
#endif
}

Vector15d Edge3::grad(
    Eigen::ConstRef<Eigen::Vector3d> d, Eigen::ConstRef<Vector12d> x) const
{
#ifdef IPC_TOOLKIT_DEBUG_AUTODIFF
    Vector15d tmp;
    tmp << d, x;
    ScalarBase::setVariableCount(15);
    using T = ADGrad<15>;
    auto X = slice_positions<T, 5, 3>(tmp);
    return smooth_edge3_term_template<T>(
               X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3),
               X.row(4), m_params, otypes, orientable)
        .grad;
#else
    return std::get<1>(smooth_edge3_term_gradient(
        d, x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), m_params,
        otypes, orientable));
#endif
}

Matrix15d Edge3::hessian(
    Eigen::ConstRef<Eigen::Vector3d> d, Eigen::ConstRef<Vector12d> x) const
{
#ifdef IPC_TOOLKIT_DEBUG_AUTODIFF
    Vector15d tmp;
    tmp << d, x;
    ScalarBase::setVariableCount(15);
    using T = ADHessian<15>;
    auto X = slice_positions<T, 5, 3>(tmp);
    return smooth_edge3_term_template<T>(
               X.row(0) / X.row(0).norm(), X.row(1), X.row(2), X.row(3),
               X.row(4), m_params, otypes, orientable)
        .Hess;
#else
    return std::get<2>(smooth_edge3_term_hessian(
        d, x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(), m_params,
        otypes, orientable));
#endif
}

double smooth_edge3_normal_term(
    Eigen::ConstRef<Eigen::Vector3d> dn,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    Eigen::ConstRef<Eigen::Vector3d> f0,
    Eigen::ConstRef<Eigen::Vector3d> f1,
    const double alpha,
    const double beta,
    const OrientationTypes& otypes)
{
    if (otypes.normal_type(0) == HeavisideType::ONE
        || otypes.normal_type(1) == HeavisideType::ONE) {
        return 1.;
    }

    const Eigen::Vector3d t0 =
        PointEdgeDistance<double, 3>::point_line_closest_point_direction(
            f0, e0, e1)
            .normalized();
    const Eigen::Vector3d t1 =
        PointEdgeDistance<double, 3>::point_line_closest_point_direction(
            f1, e0, e1)
            .normalized();
    const Eigen::Vector3d edge = (e0 - e1).normalized();
    const Eigen::Vector3d d = project<double>(dn, edge).normalized();
    return Math<double>::smooth_heaviside(
        (d - t0).cross(d - t1).dot(edge), alpha, beta);
}

GradType<15> smooth_edge3_normal_term_gradient(
    Eigen::ConstRef<Eigen::Vector3d> dn,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    Eigen::ConstRef<Eigen::Vector3d> f0,
    Eigen::ConstRef<Eigen::Vector3d> f1,
    const double alpha,
    const double beta,
    const OrientationTypes& otypes)
{
    double val = 1.;
    Vector15d gradient = Vector15d::Zero();

    if (otypes.normal_type(0) == HeavisideType::ONE
        || otypes.normal_type(1) == HeavisideType::ONE) {
        return std::make_tuple(val, gradient);
    }

    {
        const auto t0n =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f0, e0, e1)
                .normalized();

        const auto t1n =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f1, e0, e1)
                .normalized();
        const auto edge = (e0 - e1).normalized();

        const Eigen::Vector3d d = project<double>(dn, edge).normalized();
        val = (d - t0n).cross(d - t1n).dot(edge);
    }

    autogen::edge_normal_term_gradient(
        dn[0], dn[1], dn[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2], f0[0],
        f0[1], f0[2], f1[0], f1[1], f1[2], gradient.data());

    gradient *= Math<double>::smooth_heaviside_grad(val, alpha, beta);
    val = Math<double>::smooth_heaviside(val, alpha, beta);
    return std::make_tuple(val, gradient);
}

HessianType<15> smooth_edge3_normal_term_hessian(
    Eigen::ConstRef<Eigen::Vector3d> dn,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    Eigen::ConstRef<Eigen::Vector3d> f0,
    Eigen::ConstRef<Eigen::Vector3d> f1,
    const double alpha,
    const double beta,
    const OrientationTypes& otypes)
{
    double val = 1.;
    Vector15d gradient = Vector15d::Zero();
    Matrix15d hessian = Matrix15d::Zero();

    if (otypes.normal_type(0) == HeavisideType::ONE
        || otypes.normal_type(1) == HeavisideType::ONE) {
        return std::make_tuple(val, gradient, hessian);
    }

    {
        const auto t0n =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f0, e0, e1)
                .normalized();

        const auto t1n =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                f1, e0, e1)
                .normalized();
        const auto edge = (e0 - e1).normalized();

        const Eigen::Vector3d d = project<double>(dn, edge).normalized();
        val = (d - t0n).cross(d - t1n).dot(edge);
    }

    autogen::edge_normal_term_gradient(
        dn[0], dn[1], dn[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2], f0[0],
        f0[1], f0[2], f1[0], f1[1], f1[2], gradient.data());

    autogen::edge_normal_term_hessian(
        dn[0], dn[1], dn[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2], f0[0],
        f0[1], f0[2], f1[0], f1[1], f1[2], hessian.data());

    const double tmp = Math<double>::smooth_heaviside_grad(val, alpha, beta);
    hessian = hessian * tmp
        + gradient * Math<double>::smooth_heaviside_hess(val, alpha, beta)
            * gradient.transpose();
    gradient *= tmp;
    val = Math<double>::smooth_heaviside(val, alpha, beta);

    // verify with autodiff
    // {
    //     ScalarBase::setVariableCount(15);
    //     using T = ADHessian<15>;
    //     Vector15d tmp;
    //     tmp << dn, e0, e1, f0, f1;
    //     auto X = slice_positions<T, 15, 1>(tmp);
    //     const Eigen::Vector3<T> n0 = (X.segment<3>(3) -
    //     X.segment<3>(9)).cross(X.segment<3>(6) - X.segment<3>(9)); const
    //     Eigen::Vector3<T> n1 = -(X.segment<3>(3) -
    //     X.tail<3>()).cross(X.segment<3>(6) - X.tail<3>()); T normal_term =
    //     Math<T>::smooth_heaviside(
    //         (Math<T>::smooth_heaviside(X.head<3>().dot(n0) / n0.norm(),
    //         alpha, beta)
    //         + Math<T>::smooth_heaviside(X.head<3>().dot(n1) / n1.norm(),
    //         alpha, beta)
    //         - 1), 1., 0);

    //     if ((normal_term.Hess - hessian).norm() > 1e-6 *
    //     hessian.norm()) {
    //         std::cout << "gradient: \n" <<
    //         normal_term.grad.transpose() << "\n vs \n" <<
    //         gradient.transpose() << std::endl; std::cout << "hessian
    //         mismatch: \n" << normal_term.Hess << "\n vs \n" <<
    //         hessian << std::endl;
    //     }
    // }

    return std::make_tuple(val, gradient, hessian);
}
} // namespace ipc
