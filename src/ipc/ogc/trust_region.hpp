#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc::ogc {

/// @brief A trust region for filtering optimization steps.
/// @see "Offset Geometric Contact" by Chen et al. [2025]
struct TrustRegion {
    /// @brief Center the trust region around the current position.
    Eigen::MatrixXd trust_region_centers;

    /// @brief Trust region radii for each vertex.
    Eigen::VectorXd trust_region_radii;

    /// @brief Scaling factor for relaxing the trust region radii.
    /// @note This should be in (0, 1).
    /// @note This is referred to as \f\(2\gamma_p\f\) in the paper.
    double relaxed_radius_scaling = 0.9;

    /// @brief Threshold for updating the trust region.
    /// If more than this fraction of vertices are restricted by the trust
    /// region, we update the trust region.
    /// @note This is referred to as \f\(\gamma_e\f\) in the paper.
    double update_threshold = 0.01;

    /// @brief Inflation radius for the trust region.
    /// This is computed each time step based on the predicted motion.
    double trust_region_inflation_radius;

    /// @brief If true, the trust region will be updated on the next call to `update_if_needed`.
    bool should_update_trust_region = true;

    TrustRegion() = delete;

    /// @brief Construct a new Trust Region object.
    /// @param dhat The offset distance for contact.
    explicit TrustRegion(double dhat);

    /// @brief Warm start the time step by moving towards the predicted positions.
    /// This also initializes the trust region.
    /// @param mesh The collision mesh.
    /// @param x Current vertex positions. (Will be modified)
    /// @param pred_x Predicted vertex positions.
    /// @param collisions Collisions to be initialized.
    /// @param dhat The offset distance for contact.
    /// @param min_distance Minimum distance between elements.
    /// @param broad_phase Broad phase collision detection.
    void warm_start_time_step(
        const CollisionMesh& mesh,
        Eigen::Ref<Eigen::MatrixXd> x,
        Eigen::ConstRef<Eigen::MatrixXd> pred_x,
        NormalCollisions& collisions,
        const double dhat,
        const double min_distance = 0.0,
        BroadPhase* broad_phase = nullptr);

    /// @brief Update the trust region based on the current positions.
    /// @param[in] mesh The collision mesh.
    /// @param[in] x Current vertex positions.
    /// @param[out] collisions Collisions to be updated.
    /// @param[in] min_distance Minimum distance between elements.
    /// @param[in] broad_phase Broad phase collision detection.
    void update(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> x,
        NormalCollisions& collisions,
        const double min_distance = 0.0,
        BroadPhase* broad_phase = nullptr);

    /// @brief Update the trust region if should_update_trust_region is true.
    /// @param[in] mesh The collision mesh.
    /// @param[in] x Current vertex positions.
    /// @param[out] collisions Collisions to be updated.
    /// @param[in] min_distance Minimum distance between elements.
    /// @param[in] broad_phase Broad phase collision detection.
    void update_if_needed(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> x,
        NormalCollisions& collisions,
        const double min_distance = 0.0,
        BroadPhase* broad_phase = nullptr)
    {
        if (should_update_trust_region) {
            update(mesh, x, collisions, min_distance, broad_phase);
            should_update_trust_region = false;
        }
    }

    /// @brief Filter the optimization step dx to stay within the trust region.
    /// @param mesh The collision mesh.
    /// @param x Current vertex positions.
    /// @param dx Proposed vertex displacements.
    /// @note Sets should_update_trust_region to true if the trust region should be updated on the next iteration.
    void filter_step(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> x,
        Eigen::Ref<Eigen::MatrixXd> dx);

    /// @brief Filter the optimization step dx using Planar-DAT (Divide and Truncate).
    ///
    /// For each collision pair, computes a direction-aware division plane and
    /// truncates only the component of displacement toward that plane. This
    /// eliminates the artificial damping and deadlock of the isotropic
    /// filter_step while retaining the penetration-free guarantee.
    ///
    /// @see "Divide and Truncate: A Penetration and Inversion Free Framework
    ///      for Coupled Multi-physics Systems" [ACM SIGGRAPH 2026]
    ///
    /// @param mesh The collision mesh.
    /// @param x Current vertex positions.
    /// @param dx Proposed vertex displacements (modified in-place).
    /// @param collisions Active collision pairs (e.g. from update()).
    void planar_filter_step(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> x,
        Eigen::Ref<Eigen::MatrixXd> dx,
        const NormalCollisions& collisions);

private:
    /// @brief Compute the truncation ratio for a vertex given a division plane.
    ///
    /// Finds the parameter t_i such that x_u + t_i * dx_u lies on the plane
    /// (n, p). If the crossing is in the future and within one step, returns
    /// relaxation_ratio * t_i; otherwise returns 1 (no truncation).
    ///
    /// @param x_u Current position of vertex u.
    /// @param dx_u Proposed displacement of vertex u.
    /// @param n Division plane normal (unit vector).
    /// @param p A point on the division plane.
    /// @return Truncation ratio in (0, 1].
    double planar_truncation_ratio(
        Eigen::ConstRef<VectorMax3d> x_u,
        Eigen::ConstRef<VectorMax3d> dx_u,
        Eigen::ConstRef<VectorMax3d> n,
        Eigen::ConstRef<VectorMax3d> p) const;
};

} // namespace ipc::ogc