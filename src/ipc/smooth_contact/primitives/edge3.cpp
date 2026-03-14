#include "edge3.hpp"

#include <ipc/config.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/utils/autodiff_types.hpp>

namespace ipc {

// ============================================================================
// Constructor
// ============================================================================

Edge3::Edge3(
    const index_t id,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<VectorMax3d> d,
    const SmoothContactParameters& params)
    : Primitive(id, params)
{
    orientable =
        (mesh.is_orient_vertex(mesh.edges()(id, 0))
         && mesh.is_orient_vertex(mesh.edges()(id, 1)));

    // Gather all face neighbors.
    // m_vertex_ids layout: [e0, e1, f0, f1, ..., f_{N-1}]
    // where f_i is the vertex opposite the edge in face i.
    const auto& adj_faces = mesh.edges_to_faces()[id];

    const index_t e0_id = mesh.edges()(id, 0);
    const index_t e1_id = mesh.edges()(id, 1);

    // Start with the two edge endpoints
    m_vertex_ids.clear();
    m_vertex_ids.push_back(e0_id);
    m_vertex_ids.push_back(e1_id);

    // For each adjacent face, find the opposite vertex and record local face
    // connectivity with per-face edge winding.
    // faces(i, :) = [local_ea, local_eb, local_fi]
    // where local_ea/eb are 0 or 1 (indices into m_vertex_ids for the two
    // edge endpoints), ordered to match the face winding.
    n_neighbors = 0;

    std::vector<Eigen::Vector3i> face_rows;

    for (const auto& fi : adj_faces) {
        // Find the opposite vertex in this face
        index_t opp_vertex = -1;
        int opp_local_in_face = -1;
        for (int k = 0; k < 3; k++) {
            const auto va = mesh.faces()(fi, k);
            if (va != e0_id && va != e1_id) {
                opp_vertex = va;
                opp_local_in_face = k;
                break;
            }
        }
        if (opp_vertex < 0) {
            continue; // degenerate face
        }

        // Add the opposite vertex to m_vertex_ids
        const int local_fi = static_cast<int>(m_vertex_ids.size());
        m_vertex_ids.push_back(opp_vertex);

        // Determine per-face local edge order from the face winding.
        // The face vertices after the opposite vertex are
        //   faces(fi, (opp+1)%3) and faces(fi, (opp+2)%3).
        const index_t face_e_a = mesh.faces()(fi, (opp_local_in_face + 1) % 3);
        // face_e_a is whichever of e0/e1 appears first after opp in the face
        const int local_ea = (face_e_a == e0_id) ? 0 : 1;
        const int local_eb = 1 - local_ea;

        face_rows.push_back(Eigen::Vector3i(local_ea, local_eb, local_fi));
        n_neighbors++;
    }

    // Build faces matrix
    faces.resize(n_neighbors, 3);
    for (int i = 0; i < n_neighbors; i++) {
        faces.row(i) = face_rows[i];
    }

    if (m_vertex_ids.size() > N_EDGE_NEIGHBORS_3D) {
        logger().error(
            "Too many face neighbors for edge3 primitive! {} > {}! Increase N_EDGE_NEIGHBORS_3D in common.hpp",
            m_vertex_ids.size(), N_EDGE_NEIGHBORS_3D);
    }

    m_is_active = smooth_edge3_term_type(vertices(m_vertex_ids, Eigen::all), d);
}

// ============================================================================
// Type determination (sets otypes)
// ============================================================================

bool Edge3::smooth_edge3_term_type(
    Eigen::ConstRef<Eigen::MatrixX3d> X,
    Eigen::ConstRef<Eigen::RowVector3d> direction)
{
    otypes.set_size(n_neighbors);

    const Eigen::RowVector3d dn = direction.normalized();
    const Eigen::RowVector3d e0 = X.row(0);
    const Eigen::RowVector3d e1 = X.row(1);

    // Tangent types: one per face neighbor
    for (int a = 0; a < n_neighbors; a++) {
        const int fi_local = faces(a, 2); // local index of opposite vertex
        const Eigen::RowVector3d t =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                X.row(fi_local), e0, e1);
        otypes.tangent_type(a) = OrientationTypes::compute_type(
            -dn.dot(t) / t.norm(), m_params.alpha_t, m_params.beta_t);
        if (otypes.tangent_type(a) == HeavisideType::ZERO) {
            return false;
        }
    }

    if (!orientable) {
        return true;
    }

    // Normal types: sum of per-face normal heavisides, then threshold.
    // For each face [e0, e1, fi], the outward face normal direction (from the
    // edge's perspective) is cross(e0 - fi, e1 - fi).
    // We check if dn aligns with this normal.
    double normal_sum = 0;
    for (int a = 0; a < n_neighbors; a++) {
        const int fi_local = faces(a, 2);
        const Eigen::RowVector3d t1 = X.row(faces(a, 0)) - X.row(fi_local);
        const Eigen::RowVector3d t2 = X.row(faces(a, 1)) - X.row(fi_local);
        const double tmp = dn.dot(t1.cross(t2).normalized());
        otypes.normal_type(a) = OrientationTypes::compute_type(
            tmp, m_params.alpha_n, m_params.beta_n);
        normal_sum += Math<double>::smooth_heaviside(
            tmp, m_params.alpha_n, m_params.beta_n);
    }

    if (normal_sum >= 1) {
        for (int a = 0; a < n_neighbors; a++) {
            otypes.normal_type(a) = HeavisideType::ONE;
        }
    }

    return normal_sum > 0;
}

// ============================================================================
// Templated term (for autodiff verification)
// ============================================================================

template <typename T, int n_verts>
T Edge3::smooth_edge3_term(
    Eigen::ConstRef<Eigen::Matrix<T, n_verts, 3>> X,
    Eigen::ConstRef<Eigen::RowVector3<T>> direction) const
{
    const Eigen::RowVector3<T> dn = direction.normalized();
    const Eigen::RowVector3<T> e0 = X.row(0);
    const Eigen::RowVector3<T> e1 = X.row(1);

    // Tangent term: product of heavisides
    T tangent_term(1.);
    for (int a = 0; a < n_neighbors; a++) {
        const int fi_local = faces(a, 2);
        const Eigen::RowVector3<T> t =
            PointEdgeDistance<T, 3>::point_line_closest_point_direction(
                X.row(fi_local), e0, e1);
        if (otypes.tangent_type(a) != HeavisideType::ONE) {
            tangent_term *= Math<T>::smooth_heaviside(
                -dn.dot(t) / t.norm(), m_params.alpha_t, m_params.beta_t);
        }
    }

    // Normal term
    T normal_term(1.);
    if (orientable && otypes.normal_type(0) != HeavisideType::ONE) {
        normal_term = T(0.);
        for (int a = 0; a < n_neighbors; a++) {
            const int fi_local = faces(a, 2);
            const Eigen::RowVector3<T> t1 =
                X.row(faces(a, 0)) - X.row(fi_local);
            const Eigen::RowVector3<T> t2 =
                X.row(faces(a, 1)) - X.row(fi_local);
            normal_term += Math<T>::smooth_heaviside(
                dn.dot(t1.cross(t2).normalized()), m_params.alpha_n,
                m_params.beta_n);
        }
        normal_term = Math<T>::smooth_heaviside(normal_term - 1, 1., 0);
    }

    // Weight: squared edge length
    const T weight = (e1 - e0).squaredNorm();

    return weight * tangent_term * normal_term;
}

// Explicit instantiations for autodiff
template double Edge3::smooth_edge3_term<double, Eigen::Dynamic>(
    Eigen::ConstRef<Eigen::MatrixX3d>,
    Eigen::ConstRef<Eigen::RowVector3d>) const;

// ============================================================================
// Tangent term gradient
// ============================================================================

GradientType<Eigen::Dynamic> Edge3::smooth_edge3_tangent_term_gradient(
    Eigen::ConstRef<Eigen::RowVector3d> dn,
    Eigen::ConstRef<Eigen::MatrixX3d> tangents,
    const double alpha,
    const double beta) const
{
    // tangents layout: rows are point_line_closest_point_direction for each
    // face neighbor. The gradient is w.r.t. [dn(3), tangents(nn*3)].
    const int nn = tangents.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Ones(nn);
    Eigen::VectorXd acc_val_1 = Eigen::VectorXd::Ones(nn);
    std::vector<Vector6d> tmp_grad(nn);
    for (int a = 0; a < nn; a++) {
        if (otypes.tangent_type(a) == HeavisideType::VARIANT) {
            std::tie(values(a), tmp_grad[a]) = opposite_direction_penalty_grad(
                tangents.row(a), dn, alpha, beta);
            for (int b = 0; b < nn; b++) {
                if (otypes.tangent_type(b) == HeavisideType::VARIANT
                    && b != a) {
                    acc_val_1(b) *= values(a);
                }
            }
        }
    }

    // Gradient w.r.t. [dn(3), tangent_0(3), tangent_1(3), ...]
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

// ============================================================================
// Tangent term hessian
// ============================================================================

HessianType<Eigen::Dynamic> Edge3::smooth_edge3_tangent_term_hessian(
    Eigen::ConstRef<Eigen::RowVector3d> dn,
    Eigen::ConstRef<Eigen::MatrixX3d> tangents,
    const double alpha,
    const double beta) const
{
    const int nn = tangents.rows();
    Eigen::VectorXd values = Eigen::VectorXd::Ones(nn);
    Eigen::VectorXd acc_val_1 = Eigen::VectorXd::Ones(nn);
    Eigen::MatrixXd acc_val_2 = Eigen::MatrixXd::Ones(nn, nn);
    std::vector<Vector6d> tmp_grad(nn);
    std::vector<Matrix6d> tmp_hess(nn);
    for (int a = 0; a < nn; a++) {
        if (otypes.tangent_type(a) == HeavisideType::VARIANT) {
            std::tie(values(a), tmp_grad[a], tmp_hess[a]) =
                opposite_direction_penalty_hess(
                    tangents.row(a), dn, alpha, beta);
            for (int b = 0; b < nn; b++) {
                if (b != a) {
                    acc_val_1(b) *= values(a);
                    for (int c = 0; c < nn; c++) {
                        if (a != c) {
                            acc_val_2(b, c) *= values(a);
                        }
                    }
                }
            }
        }
    }

    const int grad_size = (nn + 1) * 3;
    Eigen::VectorXd tangent_grad = Eigen::VectorXd::Zero(grad_size);
    Eigen::MatrixXd tangent_hess = Eigen::MatrixXd::Zero(grad_size, grad_size);
    Eigen::Vector3d tmp;
    for (int a = 0; a < nn; a++) {
        if (otypes.tangent_type(a) == HeavisideType::VARIANT) {
            const int id = (a + 1) * 3;
            tangent_grad.segment<3>(id) = tmp_grad[a].head<3>() * acc_val_1(a);
            tangent_grad.segment<3>(0) += tmp_grad[a].tail<3>() * acc_val_1(a);

            tmp.setZero();
            for (int b = 0; b < nn; b++) {
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

// ============================================================================
// Normal term gradient
// ============================================================================

GradientType<Eigen::Dynamic> Edge3::smooth_edge3_normal_term_gradient(
    Eigen::ConstRef<Eigen::RowVector3d> direction,
    Eigen::ConstRef<Eigen::MatrixX3d> X,
    const double alpha,
    const double beta) const
{
    // direction here is the (negated, normalized) direction: direction =
    // -d.normalized() The type check used dn_orig.dot(n_hat) where dn_orig =
    // d.normalized(). negative_orientation_penalty(t1, t2, d_arg, ...) =
    // H(d_arg . n_hat). We need d_arg = -direction = d.normalized() to match
    // the type check.
    //
    // X layout: row 0 = e0, row 1 = e1, rows 2..2+n_neighbors-1 = face verts
    const int n_dofs_local =
        static_cast<int>(X.rows()) * 3 + 3; // +3 for direction
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(n_dofs_local);
    if (!orientable || otypes.normal_type(0) == HeavisideType::ONE) {
        return std::make_tuple(1., grad);
    }

    double normal_sum = 0.;
    for (int a = 0; a < n_neighbors; a++) {
        if (otypes.normal_type(a) == HeavisideType::VARIANT) {
            const int fi_local = faces(a, 2);
            const Eigen::RowVector3d t1 = X.row(faces(a, 0)) - X.row(fi_local);
            const Eigen::RowVector3d t2 = X.row(faces(a, 1)) - X.row(fi_local);

            // negative_orientation_penalty(t1, t2, -direction, alpha, beta)
            // = H((-direction) . n_hat) = H(d.normalized() . n_hat)
            // which matches the type check.
            const int id_e0 =
                (faces(a, 0) + 1) * 3; // +1 because direction is at front
            const int id_e1 = (faces(a, 1) + 1) * 3;
            const int id_fi = (fi_local + 1) * 3;
            const auto [y, dy] = negative_orientation_penalty_grad(
                t1, t2, -direction, alpha, beta);

            normal_sum += y;

            // dy is w.r.t. [t1(3), t2(3), d_arg(3)] where d_arg = -direction
            // t1 = e_a0 - fi => d/d(e_a0) = +, d/d(fi) = -
            // t2 = e_a1 - fi => d/d(e_a1) = +, d/d(fi) = -
            // d_arg = -direction => d(d_arg)/d(direction) = -1
            grad.segment<3>(id_e0) += dy.head<3>();
            grad.segment<3>(id_fi) -= dy.head<3>();
            grad.segment<3>(id_e1) += dy.segment<3>(3);
            grad.segment<3>(id_fi) -= dy.segment<3>(3);
            grad.head<3>() -= dy.tail<3>(); // chain rule: d_arg = -direction
        }
    }

    const double val = Math<double>::smooth_heaviside(normal_sum - 1, 1., 0);
    const double grad_val =
        Math<double>::smooth_heaviside_grad(normal_sum - 1, 1., 0);

    return std::make_tuple(val, grad * grad_val);
}

// ============================================================================
// Normal term hessian
// ============================================================================

HessianType<Eigen::Dynamic> Edge3::smooth_edge3_normal_term_hessian(
    Eigen::ConstRef<Eigen::RowVector3d> direction,
    Eigen::ConstRef<Eigen::MatrixX3d> X,
    const double alpha,
    const double beta) const
{
    // direction here is -d.normalized() (already negated by caller).
    // We need H(d.normalized() . n_hat) = H((-direction) . n_hat).
    const int n_rows = static_cast<int>(X.rows());
    const int n_dofs_local = (n_rows + 1) * 3; // direction(3) + X(n_rows*3)
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(n_dofs_local);
    Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(n_dofs_local, n_dofs_local);
    if (!orientable || otypes.normal_type(0) == HeavisideType::ONE) {
        return std::make_tuple(1., grad, hess);
    }

    // Use autodiff for the inner sum (same approach as Point3)
    double normal_sum = 0;
    {
        using T = ADHessian<Eigen::Dynamic>;
        Eigen::VectorXd tmp(3 + n_rows * 3);
        tmp.head<3>() = direction;
        for (int i = 0; i < n_rows; i++) {
            tmp.segment<3>(3 + i * 3) = X.row(i);
        }

        ScalarBase::setVariableCount(tmp.size());
        const Eigen::Matrix<T, Eigen::Dynamic, DIM> XX =
            slice_positions<T, Eigen::Dynamic, DIM>(tmp);
        // XX row 0 = direction (= -d.normalized()), rows 1..n_rows = X rows

        T normal_sum_ad(0.);
        for (int a = 0; a < n_neighbors; a++) {
            if (otypes.normal_type(a) == HeavisideType::VARIANT) {
                const int fi_local = faces(a, 2);
                const Eigen::RowVector3<T> t1 =
                    XX.row(faces(a, 0) + 1) - XX.row(fi_local + 1);
                const Eigen::RowVector3<T> t2 =
                    XX.row(faces(a, 1) + 1) - XX.row(fi_local + 1);
                // -XX.row(0) = -direction = d.normalized()
                // So this computes H(d.normalized() . n_hat)
                normal_sum_ad = normal_sum_ad
                    + Math<T>::smooth_heaviside(
                                    -XX.row(0).dot(t1.cross(t2).normalized()),
                                    m_params.alpha_n, m_params.beta_n);
            }
        }

        normal_sum = normal_sum_ad.val;
        grad = normal_sum_ad.grad;
        hess = normal_sum_ad.Hess;
    }

    const double val = Math<double>::smooth_heaviside(normal_sum - 1, 1., 0);
    const double grad_val =
        Math<double>::smooth_heaviside_grad(normal_sum - 1, 1., 0);
    const double hess_val =
        Math<double>::smooth_heaviside_hess(normal_sum - 1, 1., 0);

    return std::make_tuple(
        val, grad * grad_val,
        grad * hess_val * grad.transpose() + grad_val * hess);
}

// ============================================================================
// Combined gradient
// ============================================================================

GradientType<Eigen::Dynamic> Edge3::smooth_edge3_term_gradient(
    Eigen::ConstRef<Eigen::RowVector3d> direction,
    Eigen::ConstRef<Eigen::MatrixX3d> X) const
{
    // X layout: row 0 = e0, row 1 = e1, rows 2..2+n_neighbors-1 = fi
    const int n_verts = static_cast<int>(X.rows());
    const int n_dofs = (n_verts + 1) * 3; // direction(3) + X(n_verts*3)

    // Compute tangent directions: t_a = point_line_closest_point_direction(fi,
    // e0, e1)
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> tangents(
        n_neighbors, 3);
    for (int a = 0; a < n_neighbors; a++) {
        tangents.row(a) =
            PointEdgeDistance<double, 3>::point_line_closest_point_direction(
                X.row(faces(a, 2)), X.row(0), X.row(1));
    }

    auto [dn, dn_grad] = normalization_and_jacobian(direction);
    // Negate dn so that:
    //   tangent term uses opposite_direction_penalty(t, dn, ...) = H(dn·t/|t|)
    //     = H(-original_dn·t/|t|), matching the type check.
    //   normal term gradient receives dn = -original_dn, and internally
    //     passes -dn = original_dn as d_arg to negative_orientation_penalty,
    //     giving H(original_dn · n_hat), matching the type check.
    dn *= -1;
    dn_grad *= -1;

    // Tangent term: gradient w.r.t. [dn(3), tangents(nn*3)]
    const auto [tangent_term, tangent_grad_raw] =
        smooth_edge3_tangent_term_gradient(
            dn, tangents, m_params.alpha_t, m_params.beta_t);

    // Normal term: gradient w.r.t. [dn(3), X(n_verts*3)]
    auto [normal_term, normal_grad_raw] = smooth_edge3_normal_term_gradient(
        dn, X, m_params.alpha_n, m_params.beta_n);

    // We need to assemble the gradient w.r.t. [direction(3), X(n_verts*3)].
    //
    // The tangent term depends on [dn, tangents], where tangents depend on
    // [e0, e1, f0, ..., f_{N-1}] via point_line_closest_point_direction.
    // We use autodiff to get the tangent derivatives w.r.t. the vertex
    // positions.
    //
    // The normal term already returns gradient w.r.t. [direction, X].
    //
    // Actually, to keep things simpler and correct, let's use the
    // tangent_grad_raw which is w.r.t. [dn, tangent_0, tangent_1, ...] and
    // chain-rule through the tangent directions.

    // Chain rule for tangent term: tangents depend on X through
    // point_line_closest_point_direction(fi, e0, e1)
    // tangent_grad_raw is size (nn+1)*3, laid out as [dn(3), t0(3), t1(3), ...]
    //
    // We need the Jacobian of each tangent w.r.t. [fi, e0, e1] (9 DOFs -> 3
    // DOFs). Use the precomputed gradient from PointEdgeDistanceDerivatives.

    // Expand tangent gradient to [direction(3), X(n_verts*3)] space
    Eigen::VectorXd tangent_grad = Eigen::VectorXd::Zero(n_dofs);
    {
        // dn part -> direction part via chain rule
        tangent_grad.head<3>() = dn_grad * tangent_grad_raw.head<3>();

        // tangent_a part -> [fi, e0, e1] parts via chain rule
        for (int a = 0; a < n_neighbors; a++) {
            const int fi_local = faces(a, 2);
            const Eigen::Vector3d tg_a =
                tangent_grad_raw.segment<3>((a + 1) * 3);

            // Get Jacobian of point_line_closest_point_direction w.r.t. [fi,
            // e0, e1]
            const auto [t_val, t_jac] = PointEdgeDistanceDerivatives<3>::
                point_line_closest_point_direction_grad(
                    X.row(fi_local), X.row(0), X.row(1));
            // t_jac is 3x9 matrix: d(tangent)/d([fi, e0, e1])
            const Eigen::Vector<double, 9> contrib = t_jac.transpose() * tg_a;

            // fi is at position fi_local in X => offset (fi_local+1)*3 in the
            // [direction, X] vector (since direction takes the first 3)
            const int id_fi = (fi_local + 1) * 3;
            const int id_e0 = (0 + 1) * 3; // = 3
            const int id_e1 = (1 + 1) * 3; // = 6

            tangent_grad.segment<3>(id_fi) += contrib.head<3>();
            tangent_grad.segment<3>(id_e0) += contrib.segment<3>(3);
            tangent_grad.segment<3>(id_e1) += contrib.tail<3>();
        }
    }

    // Normal gradient is already w.r.t. [direction(3), X(n_verts*3)], but the
    // direction part needs chain rule through normalization
    Eigen::VectorXd normal_grad = normal_grad_raw;
    normal_grad.head<3>() = dn_grad * normal_grad_raw.head<3>();

    // Combined: val = weight * tangent_term * normal_term
    double val = tangent_term * normal_term;

    Eigen::VectorXd grad_tmp =
        tangent_grad * normal_term + normal_grad * tangent_term;

    // Weight = (e1 - e0).squaredNorm()
    const Eigen::RowVector3d edge = X.row(1) - X.row(0);
    const double weight = edge.squaredNorm();
    grad_tmp *= weight;
    // Derivative of weight w.r.t. e0 and e1
    // d/d(e0) (e1-e0)^2 = 2*(e0-e1), d/d(e1) (e1-e0)^2 = 2*(e1-e0)
    const int id_e0 = 3; // offset in [direction, X]
    const int id_e1 = 6;
    grad_tmp.segment<3>(id_e0) += (2 * val) * (-edge).transpose();
    grad_tmp.segment<3>(id_e1) += (2 * val) * edge.transpose();

    val *= weight;

    return std::make_tuple(val, grad_tmp);
}

// ============================================================================
// Combined hessian
// ============================================================================

HessianType<Eigen::Dynamic> Edge3::smooth_edge3_term_hessian(
    Eigen::ConstRef<Eigen::RowVector3d> direction,
    Eigen::ConstRef<Eigen::MatrixX3d> X) const
{
    // Use full autodiff for the hessian. This is simpler and avoids
    // complex chain-rule assembly for the tangent term through
    // point_line_closest_point_direction.
    const int n_verts = static_cast<int>(X.rows());
    const int n_dofs = (n_verts + 1) * 3;

    using T = ADHessian<Eigen::Dynamic>;
    Eigen::VectorXd tmp(n_dofs);
    tmp.head<3>() = direction;
    for (int i = 0; i < n_verts; i++) {
        tmp.segment<3>(3 + i * 3) = X.row(i);
    }
    ScalarBase::setVariableCount(n_dofs);
    const Eigen::Matrix<T, Eigen::Dynamic, DIM> XX =
        slice_positions<T, Eigen::Dynamic, DIM>(tmp);

    // XX: row 0 = direction, row 1 = e0, row 2 = e1, rows 3.. = fi
    T result =
        smooth_edge3_term<T, Eigen::Dynamic>(XX.bottomRows(n_verts), XX.row(0));

    return std::make_tuple(result.val, result.grad, result.Hess);
}

// ============================================================================
// Public interface: potential, grad, hessian
// ============================================================================

double Edge3::potential(
    Eigen::ConstRef<Eigen::Vector3d> d,
    Eigen::ConstRef<VectorMax<double, MAX_SIZE>> x) const
{
    const Eigen::Matrix<double, Eigen::Dynamic, DIM> X =
        slice_positions<double, Eigen::Dynamic, DIM>(x);
    return smooth_edge3_term<double, Eigen::Dynamic>(X, d);
}

VectorMax<double, Edge3::MAX_SIZE + Edge3::DIM> Edge3::grad(
    Eigen::ConstRef<Eigen::Vector3d> d,
    Eigen::ConstRef<VectorMax<double, MAX_SIZE>> x) const
{
#ifdef IPC_TOOLKIT_DEBUG_AUTODIFF
    using T = ADGrad<Eigen::Dynamic>;
    Eigen::VectorXd tmp(x.size() + d.size());
    tmp << d, x;
    ScalarBase::setVariableCount(tmp.size());
    const Eigen::Matrix<T, Eigen::Dynamic, DIM> X =
        slice_positions<T, Eigen::Dynamic, DIM>(tmp);
    return smooth_edge3_term<T, Eigen::Dynamic>(
               X.bottomRows(X.rows() - 1), X.row(0))
        .grad;
#else
    const Eigen::Matrix<double, Eigen::Dynamic, DIM> X =
        slice_positions<double, Eigen::Dynamic, DIM>(x);
    const auto [val, grad] = smooth_edge3_term_gradient(d, X);
    return grad;
#endif
}

MatrixMax<double, Edge3::MAX_SIZE + Edge3::DIM, Edge3::MAX_SIZE + Edge3::DIM>
Edge3::hessian(
    Eigen::ConstRef<Eigen::Vector3d> d,
    Eigen::ConstRef<VectorMax<double, MAX_SIZE>> x) const
{
#ifdef IPC_TOOLKIT_DEBUG_AUTODIFF
    using T = ADHessian<Eigen::Dynamic>;
    Eigen::VectorXd tmp(x.size() + d.size());
    tmp << d, x;
    ScalarBase::setVariableCount(tmp.size());
    const Eigen::Matrix<T, Eigen::Dynamic, DIM> X =
        slice_positions<T, Eigen::Dynamic, DIM>(tmp);
    return smooth_edge3_term<T, Eigen::Dynamic>(
               X.bottomRows(X.rows() - 1), X.row(0))
        .Hess;
#else
    const auto X = slice_positions<double, Eigen::Dynamic, DIM>(x);
    const auto [val, grad, hess] = smooth_edge3_term_hessian(d, X);
    return hess;
#endif
}

} // namespace ipc