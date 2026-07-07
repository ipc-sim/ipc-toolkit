#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::affine {

/// @brief Linear equality joint constraints for affine bodies.
///
/// A material point p(x) = A x̄ + p is *linear* in the affine DOFs
/// [p; vec(A)], so the common articulation joints (point connection, fixed
/// point, hinge, sliding plane) are linear equality constraints Cx = s
/// [Chen et al. 2022, §4.1]. They are enforced exactly via a change of
/// variables x = Uz: build an upper-triangular basis extension V of the
/// constraint rows (Gaussian elimination with column pivoting), set U = V⁻¹,
/// and pin the first m entries of z = Vx to s (Dirichlet-style) during the
/// Newton solve.
///
/// Anchor points are given in world space at the initial configuration and
/// converted to material coordinates using the initial poses.
class JointConstraints {
public:
    /// @param bodies The bodies (for DOF counts and dimensions).
    /// @param initial_poses The initial poses used to convert world-space
    ///     anchors to material coordinates.
    JointConstraints(
        const std::shared_ptr<const rigid::RigidBodies>& bodies,
        const std::vector<rigid::Pose>& initial_poses);

    /// @brief Glue a point of body_a to a point of body_b (3 rows).
    /// @param world_anchor The common point in world space (initially).
    void add_point_connection(
        const size_t body_a,
        const size_t body_b,
        Eigen::ConstRef<Eigen::Vector3d> world_anchor);

    /// @brief Fix a material point of a body at its initial world position (3 rows).
    void add_fixed_point(
        const size_t body, Eigen::ConstRef<Eigen::Vector3d> world_anchor);

    /// @brief Hinge between two bodies about the axis through two points
    /// (two point connections; 6 rows).
    void add_hinge(
        const size_t body_a,
        const size_t body_b,
        Eigen::ConstRef<Eigen::Vector3d> world_axis_p0,
        Eigen::ConstRef<Eigen::Vector3d> world_axis_p1);

    /// @brief Constrain a material point of a body to the plane through
    /// world_anchor with the given normal (1 row).
    void add_sliding_plane(
        const size_t body,
        Eigen::ConstRef<Eigen::Vector3d> world_anchor,
        Eigen::ConstRef<Eigen::Vector3d> normal);

    /// @brief Fix all 12 DOFs of a body at their initial values (12 rows).
    void add_fixed_body(const size_t body);

    /// @brief Number of constraint rows m.
    size_t num_constraints() const { return m_rows.size(); }

    bool empty() const { return m_rows.empty(); }

    /// @brief Build the change-of-variable matrices. Call after adding all
    /// joints (the Simulator calls this on construction).
    void finalize();

    // -- Post-finalize interface ----------------------------------------------

    /// @brief Map full DOFs to reduced coordinates: z = Vx.
    Eigen::VectorXd to_reduced(Eigen::ConstRef<Eigen::VectorXd> x) const;

    /// @brief Map reduced coordinates to full DOFs: x = Uz.
    Eigen::VectorXd to_full(Eigen::ConstRef<Eigen::VectorXd> z) const;

    /// @brief The change-of-variable matrix U (x = Uz).
    const Eigen::MatrixXd& U() const { return m_U; } // NOLINT

    /// @brief The pinned values s of the first m reduced coordinates.
    const Eigen::VectorXd& rhs() const { return m_s; }

    /// @brief Evaluate the constraint residual Cx − s.
    Eigen::VectorXd residual(Eigen::ConstRef<Eigen::VectorXd> x) const;

private:
    /// @brief Add the coefficients of world-point coordinate k of a material
    /// point x̄ on a body to a constraint row.
    void add_point_coefficients(
        Eigen::Ref<Eigen::RowVectorXd> row,
        const size_t body,
        Eigen::ConstRef<Eigen::Vector3d> material_point,
        const int k,
        const double sign) const;

    /// @brief Convert a world-space anchor to body material coordinates.
    Eigen::Vector3d world_to_material(
        const size_t body, Eigen::ConstRef<Eigen::Vector3d> world_point) const;

    std::shared_ptr<const rigid::RigidBodies> m_bodies;
    std::vector<rigid::Pose> m_initial_poses;

    /// Constraint rows of C (each of size ndof) and right-hand sides
    std::vector<Eigen::RowVectorXd> m_rows;
    std::vector<double> m_rhs;

    bool m_finalized = false;
    Eigen::MatrixXd m_U; ///< x = Uz (permutation folded in)
    Eigen::MatrixXd m_V; ///< z = Vx (permutation folded in)
    Eigen::VectorXd m_s; ///< pinned values of z.head(m)
};

} // namespace ipc::affine
