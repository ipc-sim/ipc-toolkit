#pragma once

#include <ipc/dynamics/affine/pose.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <vector>

namespace ipc::dynamics {

/// @brief The map from the optimization DOFs x to the affine coordinates
/// y = [p; vec(A) column-major] per body, and the chain rule that pulls a
/// potential's affine-coordinate derivatives back to x.
///
/// Two parameterizations share the same affine-coordinate potentials and differ
/// only here:
///   - Rigid (RBD): x = [p; θ] per body, A = Q(θ) (a change of variables that
///     enforces A ∈ SO(dim) exactly).
///   - Affine (ABD): x = y (the identity), A optimized directly.
///
/// The affine and vertex representations are block diagonal per body; every
/// method operates block-by-body on stacked vectors/Hessians.
class ToAffine {
public:
    virtual ~ToAffine() = default;

    int dim() const { return m_dim; }
    size_t num_bodies() const { return m_num_bodies; }
    int pos_ndof() const { return m_dim; }
    /// @brief Reduced rotation DOFs per body (1|3 for rigid, dim² for affine).
    virtual int rot_ndof() const = 0;
    int reduced_ndof() const { return pos_ndof() + rot_ndof(); }
    int affine_ndof() const { return m_dim + m_dim * m_dim; }
    /// @brief Whether x ≡ y (affine parameterization).
    virtual bool is_identity() const { return false; }

    /// @brief Affine coordinates y(x): [p; vec(A)] per body.
    virtual Eigen::VectorXd
    to_affine(Eigen::ConstRef<Eigen::VectorXd> x) const = 0;

    /// @brief Reduced DOFs x(y). For rigid this is the log map (canonical θ).
    virtual Eigen::VectorXd
    from_affine(Eigen::ConstRef<Eigen::VectorXd> y) const = 0;

    /// @brief Per-body affine poses at x (exact for rigid: A = R(θ)).
    virtual std::vector<affine::Pose>
    poses(Eigen::ConstRef<Eigen::VectorXd> x) const = 0;

    /// @brief Pull an affine-coordinate gradient back to x: (dy/dx)ᵀ g.
    /// @param x The reduced DOFs.
    /// @param g_affine The gradient in affine coordinates (size affine_ndof·N).
    /// @return The gradient in reduced coordinates (size reduced_ndof·N).
    virtual Eigen::VectorXd apply_gradient(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> g_affine) const = 0;

    /// @brief Pull an affine-coordinate Hessian back to x:
    /// (dy/dx)ᵀ H (dy/dx) + Σₖ (g)ₖ d²yₖ/dx², projecting the rotation block to
    /// PSD once (per body) with @p project_hessian_to_psd.
    /// @note Assumes H_affine is block diagonal per body.
    virtual Eigen::SparseMatrix<double> apply_hessian(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> g_affine,
        const Eigen::SparseMatrix<double>& H_affine,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const = 0;

protected:
    ToAffine(const int dim, const size_t num_bodies)
        : m_dim(dim)
        , m_num_bodies(num_bodies)
    {
    }

    int m_dim;
    size_t m_num_bodies;
};

/// @brief The rigid to-affine map x = [p; θ] → y = [p; vec(Q(θ))].
class RigidToAffine : public ToAffine {
public:
    RigidToAffine(const int dim, const size_t num_bodies)
        : ToAffine(dim, num_bodies)
    {
    }

    int rot_ndof() const override { return m_dim == 2 ? 1 : 3; }

    Eigen::VectorXd
    to_affine(Eigen::ConstRef<Eigen::VectorXd> x) const override;
    Eigen::VectorXd
    from_affine(Eigen::ConstRef<Eigen::VectorXd> y) const override;
    std::vector<affine::Pose>
    poses(Eigen::ConstRef<Eigen::VectorXd> x) const override;

    Eigen::VectorXd apply_gradient(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> g_affine) const override;

    Eigen::SparseMatrix<double> apply_hessian(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> g_affine,
        const Eigen::SparseMatrix<double>& H_affine,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const override;
};

/// @brief The identity to-affine map x ≡ y (affine bodies).
class AffineToAffine : public ToAffine {
public:
    AffineToAffine(const int dim, const size_t num_bodies)
        : ToAffine(dim, num_bodies)
    {
    }

    int rot_ndof() const override { return m_dim * m_dim; }
    bool is_identity() const override { return true; }

    Eigen::VectorXd to_affine(Eigen::ConstRef<Eigen::VectorXd> x) const override
    {
        return x;
    }
    Eigen::VectorXd
    from_affine(Eigen::ConstRef<Eigen::VectorXd> y) const override
    {
        return y;
    }
    std::vector<affine::Pose>
    poses(Eigen::ConstRef<Eigen::VectorXd> x) const override;

    Eigen::VectorXd apply_gradient(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> g_affine) const override
    {
        return g_affine;
    }

    Eigen::SparseMatrix<double> apply_hessian(
        Eigen::ConstRef<Eigen::VectorXd> x,
        Eigen::ConstRef<Eigen::VectorXd> g_affine,
        const Eigen::SparseMatrix<double>& H_affine,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const override;
};

} // namespace ipc::dynamics
