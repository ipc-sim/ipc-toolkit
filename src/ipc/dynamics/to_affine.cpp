#include "to_affine.hpp"

#include <ipc/dynamics/rigid/pose.hpp>

namespace ipc::dynamics {

namespace {
    /// Extract the dense (n×n) diagonal block starting at @p offset from a
    /// (column-major) sparse matrix that is block diagonal per body.
    /// @note n ≤ affine_ndof ≤ 12.
    MatrixMax12d extract_diagonal_block(
        const Eigen::SparseMatrix<double>& H, const int offset, const int n)
    {
        MatrixMax12d block = MatrixMax12d::Zero(n, n);
        for (int c = 0; c < n; ++c) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(H, offset + c);
                 it; ++it) {
                const int r = int(it.row()) - offset;
                if (r >= 0 && r < n) {
                    block(r, c) = it.value();
                }
            }
        }
        return block;
    }

    /// Assemble a block-diagonal sparse matrix from per-body dense blocks.
    Eigen::SparseMatrix<double> assemble_block_diagonal(
        const std::vector<MatrixMax12d>& blocks, const int block_ndof)
    {
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(blocks.size() * block_ndof * block_ndof);
        for (size_t i = 0; i < blocks.size(); ++i) {
            const int offset = int(i) * block_ndof;
            for (int c = 0; c < block_ndof; ++c) {
                for (int r = 0; r < block_ndof; ++r) {
                    const double v = blocks[i](r, c);
                    if (v != 0.0) {
                        triplets.emplace_back(offset + r, offset + c, v);
                    }
                }
            }
        }
        Eigen::SparseMatrix<double> H(
            blocks.size() * block_ndof, blocks.size() * block_ndof);
        H.setFromTriplets(triplets.begin(), triplets.end());
        return H;
    }
} // namespace

// ---- Rigid to-affine map ---------------------------------------------

Eigen::VectorXd
RigidToAffine::to_affine(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const int d = m_dim, rndof = reduced_ndof(), andof = affine_ndof();
    assert(x.size() == int(m_num_bodies) * rndof);

    Eigen::VectorXd y(andof * m_num_bodies);
    for (size_t i = 0; i < m_num_bodies; ++i) {
        affine::Pose pose;
        const VectorMax3d theta = x.segment(i * rndof + d, rot_ndof());
        pose.set_rotation_vector(theta);
        y.segment(i * andof, d) = x.segment(i * rndof, d);
        y.segment(i * andof + d, d * d) = pose.rotation.reshaped();
    }
    return y;
}

Eigen::VectorXd
RigidToAffine::from_affine(Eigen::ConstRef<Eigen::VectorXd> y) const
{
    const int d = m_dim, rndof = reduced_ndof(), andof = affine_ndof();
    assert(y.size() == int(m_num_bodies) * andof);

    Eigen::VectorXd x(rndof * m_num_bodies);
    for (size_t i = 0; i < m_num_bodies; ++i) {
        affine::Pose pose;
        pose.rotation = y.segment(i * andof + d, d * d).reshaped(d, d);
        x.segment(i * rndof, d) = y.segment(i * andof, d);
        x.segment(i * rndof + d, rot_ndof()) = pose.rotation_vector();
    }
    return x;
}

std::vector<affine::Pose>
RigidToAffine::poses(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const int d = m_dim, rndof = reduced_ndof();
    std::vector<affine::Pose> out(m_num_bodies);
    for (size_t i = 0; i < m_num_bodies; ++i) {
        out[i].position = x.segment(i * rndof, d);
        const VectorMax3d theta = x.segment(i * rndof + d, rot_ndof());
        out[i].set_rotation_vector(theta);
    }
    return out;
}

Eigen::VectorXd RigidToAffine::apply_gradient(
    Eigen::ConstRef<Eigen::VectorXd> x,
    Eigen::ConstRef<Eigen::VectorXd> g_affine) const
{
    const int d = m_dim, rot = rot_ndof(), rndof = reduced_ndof(),
              andof = affine_ndof();
    assert(g_affine.size() == int(m_num_bodies) * andof);

    Eigen::VectorXd g(rndof * m_num_bodies);
    for (size_t i = 0; i < m_num_bodies; ++i) {
        const VectorMax3d theta = x.segment(i * rndof + d, rot);
        const MatrixMax<double, 9, 3> dQ =
            rigid::rotation_to_matrix_jacobian(theta);

        // Position block: dy/dx = I
        g.segment(i * rndof, d) = g_affine.segment(i * andof, d);
        // Rotation block: (dvec(Q)/dθ)ᵀ · g_A
        g.segment(i * rndof + d, rot) =
            dQ.transpose() * g_affine.segment(i * andof + d, d * d);
    }
    return g;
}

Eigen::SparseMatrix<double> RigidToAffine::apply_hessian(
    Eigen::ConstRef<Eigen::VectorXd> x,
    Eigen::ConstRef<Eigen::VectorXd> g_affine,
    const Eigen::SparseMatrix<double>& H_affine,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const int d = m_dim, rot = rot_ndof(), rndof = reduced_ndof(),
              andof = affine_ndof();
    assert(H_affine.rows() == int(m_num_bodies) * andof);

    std::vector<MatrixMax12d> blocks(m_num_bodies);
    for (size_t i = 0; i < m_num_bodies; ++i) {
        const MatrixMax12d Ha =
            extract_diagonal_block(H_affine, int(i) * andof, andof);

        const VectorMax3d theta = x.segment(i * rndof + d, rot);
        const MatrixMax<double, 9, 3> dQ =
            rigid::rotation_to_matrix_jacobian(theta);
        const MatrixMax<double, 9, 9> d2Q =
            rigid::rotation_to_matrix_hessian(theta);

        MatrixMax6d H = MatrixMax6d::Zero(rndof, rndof);
        // pos-pos
        H.topLeftCorner(d, d) = Ha.topLeftCorner(d, d);
        // pos-rot and rot-pos (zero for our separable terms, general here)
        H.block(0, d, d, rot) = Ha.block(0, d, d, d * d) * dQ;
        H.block(d, 0, rot, d) = dQ.transpose() * Ha.block(d, 0, d * d, d);
        // rot-rot: (dQ)ᵀ Haa (dQ) + Σₖ (g_A)ₖ d²Qₖ/dθ²
        MatrixMax3d H_rot = dQ.transpose() * Ha.block(d, d, d * d, d * d) * dQ;
        H_rot += (d2Q.transpose() * g_affine.segment(i * andof + d, d * d))
                     .reshaped(rot, rot);
        // Project only the rotation block (position block is PSD by
        // construction)
        H.block(d, d, rot, rot) = project_to_psd(H_rot, project_hessian_to_psd);

        blocks[i] = H;
    }
    return assemble_block_diagonal(blocks, rndof);
}

// ---- Affine (identity) to-affine map ---------------------------------

std::vector<affine::Pose>
AffineToAffine::poses(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const int d = m_dim, andof = affine_ndof();
    std::vector<affine::Pose> out(m_num_bodies);
    for (size_t i = 0; i < m_num_bodies; ++i) {
        out[i].position = x.segment(i * andof, d);
        out[i].rotation = x.segment(i * andof + d, d * d).reshaped(d, d);
    }
    return out;
}

Eigen::SparseMatrix<double> AffineToAffine::apply_hessian(
    Eigen::ConstRef<Eigen::VectorXd> x,
    Eigen::ConstRef<Eigen::VectorXd> g_affine,
    const Eigen::SparseMatrix<double>& H_affine,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    // x ≡ y, so the Hessian is unchanged apart from the (optional) PSD
    // projection of the affine-matrix block of each body.
    if (project_hessian_to_psd == PSDProjectionMethod::NONE) {
        return H_affine;
    }

    const int d = m_dim, andof = affine_ndof();
    std::vector<MatrixMax12d> blocks(m_num_bodies);
    for (size_t i = 0; i < m_num_bodies; ++i) {
        MatrixMax12d Ha =
            extract_diagonal_block(H_affine, int(i) * andof, andof);
        // Project the vec(A) block (position block is PSD by construction).
        const MatrixMax<double, 9, 9> A_block =
            Ha.bottomRightCorner(d * d, d * d);
        Ha.bottomRightCorner(d * d, d * d) =
            project_to_psd(A_block, project_hessian_to_psd);
        blocks[i] = Ha;
    }
    return assemble_block_diagonal(blocks, andof);
}

} // namespace ipc::dynamics
