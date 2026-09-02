#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <string>

namespace ipc {

class HessianAssembler; // forward declaration (see utils/hessian_assembler.hpp)

/// @brief Base class for potentials.
/// @tparam TCollisions The type of the collisions.
template <class TCollisions> class Potential {
protected:
    using TCollision = typename TCollisions::value_type;
    /// @brief Maximum degrees of freedom per collision
    static constexpr int STENCIL_NDOF = 3 * TCollision::STENCIL_SIZE;
    using VectorMaxNd = VectorMax<double, STENCIL_NDOF>;
    using MatrixMaxNd = MatrixMax<double, STENCIL_NDOF, STENCIL_NDOF>;

public:
    Potential() = default;
    virtual ~Potential() = default;

    /// @brief The name of this potential (used for profiling).
    virtual std::string name() const = 0;

    // -- Cumulative methods ---------------------------------------------------

    /// @brief Compute the potential for a set of collisions.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @returns The potential for a set of collisions.
    double operator()(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the gradient of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param in_full_dof If true, return the gradient in full-mesh DOF (equivalent to `mesh.to_full_dof(gradient)`, but assembled directly in full DOF when possible, avoiding the extra map).
    /// @returns The gradient of the potential w.r.t. X. This will have a size of X.size() (or `mesh.full_ndof()` if in_full_dof).
    Eigen::VectorXd gradient(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const bool in_full_dof = false) const;

    /// @brief Compute the hessian of the potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @param in_full_dof If true, return the Hessian in full-mesh DOF (equivalent to `mesh.to_full_dof(hessian)`, but assembled directly in full DOF when possible, avoiding the two sparse-matrix products).
    /// @returns The Hessian of the potential w.r.t. X. This will have a size of X.size() by X.size() (or `mesh.full_ndof()` square if in_full_dof).
    Eigen::SparseMatrix<double> hessian(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE,
        const bool in_full_dof = false) const;

    /// @brief Assemble the Hessian of the potential using a custom assembler.
    ///
    /// Evaluates the local Hessian of every collision (in parallel) and feeds
    /// each to `assembler` (see HessianAssembler). This decouples the local
    /// derivative evaluation from the global matrix construction; hessian()
    /// is a thin wrapper around this using a TripletHessianAssembler.
    ///
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param X Degrees of freedom of the collision mesh (e.g., vertices or velocities).
    /// @param assembler The assembler that accumulates the local Hessians.
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @param in_full_dof If true, stencil vertex IDs are remapped to full-mesh vertex IDs (requires `mesh.is_selection_dof_map()`; throws otherwise).
    void assemble_hessian(
        const TCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        HessianAssembler& assembler,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE,
        const bool in_full_dof = false) const;

    // -- Single collision methods ---------------------------------------------

    /// @brief Compute the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @return The potential.
    virtual double operator()(
        const TCollision& collision, Eigen::ConstRef<VectorMaxNd> x) const = 0;

    /// @brief Compute the gradient of the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @return The gradient of the potential.
    virtual VectorMaxNd gradient(
        const TCollision& collision, Eigen::ConstRef<VectorMaxNd> x) const = 0;

    /// @brief Compute the hessian of the potential for a single collision.
    /// @param collision The collision.
    /// @param x The collision stencil's degrees of freedom.
    /// @param project_hessian_to_psd Whether to project the hessian to the positive semi-definite cone.
    /// @return The hessian of the potential.
    virtual MatrixMaxNd hessian(
        const TCollision& collision,
        Eigen::ConstRef<VectorMaxNd> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const = 0;
};

} // namespace ipc