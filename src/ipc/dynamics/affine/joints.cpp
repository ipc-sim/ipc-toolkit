#include "joints.hpp"

#include <ipc/utils/logger.hpp>

#include <numeric>
#include <stdexcept>

namespace ipc::affine {

JointConstraints::JointConstraints(
    const std::shared_ptr<const rigid::RigidBodies>& bodies,
    const std::vector<rigid::Pose>& initial_poses)
    : m_bodies(bodies)
    , m_initial_poses(initial_poses)
{
    assert(m_bodies->dim() == 3); // TODO: Support 2D affine bodies
    assert(initial_poses.size() == m_bodies->num_bodies());
}

Eigen::Vector3d JointConstraints::world_to_material(
    const size_t body, Eigen::ConstRef<Eigen::Vector3d> world_point) const
{
    assert(body < m_bodies->num_bodies());
    const rigid::Pose& pose = m_initial_poses[body];
    // x̄ = A⁻¹ (p_world − p); A ∈ SO(3) initially, so A⁻¹ = Aᵀ
    return pose.rotation_matrix().transpose() * (world_point - pose.position);
}

void JointConstraints::add_point_coefficients(
    Eigen::Ref<Eigen::RowVectorXd> row,
    const size_t body,
    Eigen::ConstRef<Eigen::Vector3d> material_point,
    const int k,
    const double sign) const
{
    // p(x)_k = x[12b + k] + Σⱼ x[12b + 3 + k + 3j] x̄ⱼ
    row(12 * body + k) += sign;
    for (int j = 0; j < 3; j++) {
        row(12 * body + 3 + k + 3 * j) += sign * material_point(j);
    }
}

void JointConstraints::add_point_connection(
    const size_t body_a,
    const size_t body_b,
    Eigen::ConstRef<Eigen::Vector3d> world_anchor)
{
    assert(!m_finalized);
    assert(body_a != body_b);
    const Eigen::Vector3d xa = world_to_material(body_a, world_anchor);
    const Eigen::Vector3d xb = world_to_material(body_b, world_anchor);
    const int n = 12 * m_bodies->num_bodies();
    for (int k = 0; k < 3; k++) {
        Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(n);
        add_point_coefficients(row, body_a, xa, k, 1.0);
        add_point_coefficients(row, body_b, xb, k, -1.0);
        m_rows.push_back(row);
        m_rhs.push_back(0.0);
    }
}

void JointConstraints::add_fixed_point(
    const size_t body, Eigen::ConstRef<Eigen::Vector3d> world_anchor)
{
    assert(!m_finalized);
    const Eigen::Vector3d xb = world_to_material(body, world_anchor);
    const int n = 12 * m_bodies->num_bodies();
    for (int k = 0; k < 3; k++) {
        Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(n);
        add_point_coefficients(row, body, xb, k, 1.0);
        m_rows.push_back(row);
        m_rhs.push_back(world_anchor(k));
    }
}

void JointConstraints::add_hinge(
    const size_t body_a,
    const size_t body_b,
    Eigen::ConstRef<Eigen::Vector3d> world_axis_p0,
    Eigen::ConstRef<Eigen::Vector3d> world_axis_p1)
{
    add_point_connection(body_a, body_b, world_axis_p0);
    add_point_connection(body_a, body_b, world_axis_p1);
}

void JointConstraints::add_sliding_plane(
    const size_t body,
    Eigen::ConstRef<Eigen::Vector3d> world_anchor,
    Eigen::ConstRef<Eigen::Vector3d> normal)
{
    assert(!m_finalized);
    const Eigen::Vector3d xb = world_to_material(body, world_anchor);
    const Eigen::Vector3d n_hat = normal.normalized();
    const int n = 12 * m_bodies->num_bodies();
    Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(n);
    for (int k = 0; k < 3; k++) {
        add_point_coefficients(row, body, xb, k, n_hat(k));
    }
    m_rows.push_back(row);
    m_rhs.push_back(n_hat.dot(world_anchor));
}

void JointConstraints::add_fixed_body(const size_t body)
{
    assert(!m_finalized);
    const int n = 12 * m_bodies->num_bodies();
    const rigid::Pose& pose = m_initial_poses[body];
    const Eigen::Matrix3d A = pose.rotation_matrix();
    for (int d = 0; d < 12; d++) {
        Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(n);
        row(12 * body + d) = 1.0;
        m_rows.push_back(row);
        m_rhs.push_back(
            d < 3 ? pose.position(d) : A.reshaped()(d - 3)); // column-major
    }
}

void JointConstraints::finalize()
{
    if (m_finalized) {
        return;
    }

    const int n = 12 * int(m_bodies->num_bodies());
    const int m = int(m_rows.size());
    assert(m < n);

    // Assemble C (m×n) and s
    Eigen::MatrixXd C(m, n);
    Eigen::VectorXd s(m);
    for (int i = 0; i < m; i++) {
        C.row(i) = m_rows[i];
        s(i) = m_rhs[i];
    }

    // Gaussian elimination with column pivoting: make C upper triangular in a
    // permuted DOF order so that V = [C; 0 I] is a full-rank upper-triangular
    // basis extension [Chen et al. 2022, §4.1].
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);

    for (int i = 0; i < m; i++) {
        // Eliminate the previous pivots from row i
        for (int k = 0; k < i; k++) {
            const double f = C(i, k) / C(k, k);
            if (f != 0) {
                C.row(i) -= f * C.row(k);
                s(i) -= f * s(k);
            }
        }
        // Swap the largest remaining coefficient onto the diagonal
        int j_star;
        const double pivot =
            C.row(i).tail(n - i).cwiseAbs().maxCoeff(&j_star);
        j_star += i;
        if (pivot < 1e-12) {
            logger().error(
                "joint constraint row {} is linearly dependent; joints are "
                "redundant",
                i);
            throw std::runtime_error("redundant joint constraints");
        }
        if (j_star != i) {
            C.col(i).swap(C.col(j_star));
            std::swap(perm[i], perm[j_star]);
        }
    }

    // V_perm = [C̃; 0 I] (upper triangular) and U_perm = V_perm⁻¹
    Eigen::MatrixXd V_perm = Eigen::MatrixXd::Identity(n, n);
    V_perm.topRows(m) = C;

    const Eigen::MatrixXd U_perm = V_perm.triangularView<Eigen::Upper>().solve(
        Eigen::MatrixXd::Identity(n, n));

    // Fold the DOF permutation in: x_perm[i] = x[perm[i]]
    m_U.setZero(n, n);
    m_V.setZero(n, n);
    for (int i = 0; i < n; i++) {
        m_U.row(perm[i]) = U_perm.row(i);
        m_V.col(perm[i]) = V_perm.col(i);
    }

    m_s = s;
    m_finalized = true;
}

Eigen::VectorXd
JointConstraints::to_reduced(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    assert(m_finalized);
    return m_V * x;
}

Eigen::VectorXd
JointConstraints::to_full(Eigen::ConstRef<Eigen::VectorXd> z) const
{
    assert(m_finalized);
    return m_U * z;
}

Eigen::VectorXd
JointConstraints::residual(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    Eigen::VectorXd r(m_rows.size());
    for (size_t i = 0; i < m_rows.size(); i++) {
        r(i) = m_rows[i].dot(x) - m_rhs[i];
    }
    return r;
}

} // namespace ipc::affine
