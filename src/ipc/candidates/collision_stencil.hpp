#pragma once

#include <ipc/config.hpp>
#include <ipc/candidates/stencil_mixin.hpp>
#include <ipc/ccd/default_narrow_phase_ccd.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <array>

namespace ipc {

/// @brief A stencil representing a collision between at most four vertices.
///
/// The runtime-polymorphic stencil, used where the stencil type is not known
/// statically. The operations built on top of the primitives below live in
/// StencilMixin, which this shares with the non-polymorphic candidate types.
class CollisionStencil : public StencilMixin<CollisionStencil> {
    // So the mixin can reach the protected unnormalized-normal primitives.
    friend class StencilMixin<CollisionStencil>;

public:
    virtual ~CollisionStencil() = default;

    /// @brief Get this as a child collision type.
    /// @tparam CollisionStencilType The child collision type to get this as.
    /// @return This as the child collision type.
    /// @throws std::bad_cast if this is not of type CollisionStencilType.
    template <typename CollisionStencilType> CollisionStencilType& as()
    {
        static_assert(
            std::is_base_of<CollisionStencil, CollisionStencilType>::value,
            "Template argument CollisionStencilType must inherit from CollisionStencil");
        return dynamic_cast<CollisionStencilType&>(*this);
    }

    /// @brief Get this as a child collision type.
    /// @tparam CollisionStencilType The child collision type to get this as.
    /// @return This as the child collision type.
    /// @throws std::bad_cast if this is not of type CollisionStencilType.
    template <typename CollisionStencilType>
    const CollisionStencilType& as() const
    {
        static_assert(
            std::is_base_of<CollisionStencil, CollisionStencilType>::value,
            "Template argument CollisionStencilType must inherit from CollisionStencil");
        return dynamic_cast<const CollisionStencilType&>(*this);
    }

    /// @brief Get the number of vertices in the collision stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the vertex IDs of the collision stencil.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex IDs of the collision stencil. Elements i > num_vertices() are -1.
    virtual std::array<index_t, STENCIL_SIZE> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const = 0;

    // The mesh-based overloads live in the mixin; naming the positions-based
    // ones below would otherwise hide them.
    using StencilMixin<CollisionStencil>::compute_distance;
    using StencilMixin<CollisionStencil>::compute_distance_gradient;
    using StencilMixin<CollisionStencil>::compute_distance_hessian;
    using StencilMixin<CollisionStencil>::compute_coefficients;

    // ----------------------------------------------------------------------
    // NOTE: The following functions take stencil vertices as output by dof()
    // ----------------------------------------------------------------------

    /// @brief Compute the distance of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @note positions can be computed as stencil.dof(vertices, edges, faces)
    /// @return Distance of the stencil.
    virtual double
    compute_distance(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.
    /// @param positions Stencil's vertex positions.
    /// @note positions can be computed as stencil.dof(vertices, edges, faces)
    /// @return Distance gradient of the stencil w.r.t. the stencil's vertex positions.
    virtual VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    /// @param positions Stencil's vertex positions.
    /// @note positions can be computed as stencil.dof(vertices, edges, faces)
    /// @return Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    virtual MatrixMax12d
    compute_distance_hessian(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the coefficients of the stencil s.t. d(x) = ‖∑ cᵢ xᵢ‖².
    /// @param positions Stencil's vertex positions.
    /// @return Coefficients of the stencil.
    virtual VectorMax4d
    compute_coefficients(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Perform narrow-phase CCD on the candidate.
    /// @param[in] vertices_t0 Stencil vertices at the start of the time step.
    /// @param[in] vertices_t1 Stencil vertices at the end of the time step.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] min_distance Minimum separation distance between primitives.
    /// @param[in] tmax Maximum time (normalized) to look for collisions.
    /// @param[in] narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @return If the candidate had a collision over the time interval.
    virtual bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const = 0;

protected:
    /// @brief Compute the unnormalized normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @return Unnormalized normal of the stencil.
    virtual VectorMax3d compute_unnormalized_normal(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the Jacobian of the unnormalized normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @return Jacobian of the unnormalized normal of the stencil.
    virtual MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;
};

} // namespace ipc