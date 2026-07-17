#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/collision_mesh.hpp>
#include <ipc/barrier/barrier_type.hpp>
#include <ipc/collisions/normal/cuda/device_normal_collisions.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc::cuda {

/// @brief GPU implementation of the barrier collision potential.
///
/// Mirrors the public API of ipc::BarrierPotential, evaluating the potential
/// (and its derivatives) with CUDA kernels. Two call styles are provided:
/// - Drop-in overloads taking (NormalCollisions, CollisionMesh, vertices)
///   which upload to the device internally (convenient, but re-uploads on
///   every call).
/// - Fast-path overloads taking a prebuilt DeviceNormalCollisions and
///   DeviceVertices, for reuse across Newton iterations (build the collisions
///   once per time step, DeviceVertices::update() each iteration).
///
/// v1 limitations (throws a clear error if violated):
/// - 3D only.
/// - No plane-vertex collisions.
/// - Only BarrierType::CLAMPED_LOG and BarrierType::NORMALIZED_CLAMPED_LOG.
/// - No shape derivatives, force magnitude, or Gauss-Newton helpers.
class BarrierPotential {
public:
    /// @brief Construct a GPU barrier potential.
    /// @param dhat The activation distance of the barrier.
    /// @param stiffness The stiffness of the barrier.
    /// @param use_physical_barrier Whether to use the physical barrier.
    /// @param barrier_type The barrier function to use.
    BarrierPotential(
        const double dhat,
        const double stiffness,
        const bool use_physical_barrier = false,
        const BarrierType barrier_type = BarrierType::CLAMPED_LOG);

    /// @brief Get the activation distance of the barrier.
    double dhat() const { return m_dhat; }

    /// @brief Set the activation distance of the barrier.
    /// @param dhat The activation distance of the barrier.
    void set_dhat(const double dhat)
    {
        assert(dhat > 0);
        m_dhat = dhat;
    }

    /// @brief Get the stiffness of the barrier.
    double stiffness() const { return m_stiffness; }

    /// @brief Set the stiffness of the barrier.
    /// @param stiffness The stiffness of the barrier.
    void set_stiffness(const double stiffness)
    {
        assert(stiffness > 0);
        m_stiffness = stiffness;
    }

    /// @brief Get whether to use the physical barrier.
    bool use_physical_barrier() const { return m_use_physical_barrier; }

    /// @brief Set use physical barrier flag.
    /// @param use_physical_barrier Whether to use the physical barrier.
    void set_use_physical_barrier(const bool use_physical_barrier)
    {
        m_use_physical_barrier = use_physical_barrier;
    }

    /// @brief Get the barrier function type.
    BarrierType barrier_type() const { return m_barrier_type; }

    /// @brief Set the barrier function type.
    /// @param barrier_type The barrier function to use.
    /// @throws If the barrier type is not supported on the GPU.
    void set_barrier_type(const BarrierType barrier_type);

    // ------------------------------------------------------------------
    // Drop-in API (uploads to the device internally)

    /// @brief Compute the barrier potential for a set of collisions.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh (rowwise, |V| × 3).
    /// @return The sum of the potential over all collisions.
    double operator()(
        const NormalCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    /// @brief Compute the gradient of the barrier potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh (rowwise, |V| × 3).
    /// @return The gradient of the potential w.r.t. the vertex positions.
    Eigen::VectorXd gradient(
        const NormalCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;

    /// @brief Compute the hessian of the barrier potential.
    /// @param collisions The set of collisions.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh (rowwise, |V| × 3).
    /// @param project_hessian_to_psd Method of projecting the local hessians
    /// to the positive semi-definite cone.
    /// @return The hessian of the potential w.r.t. the vertex positions.
    Eigen::SparseMatrix<double> hessian(
        const NormalCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    // ------------------------------------------------------------------
    // Fast-path API (device-resident inputs)

    /// @brief Compute the barrier potential for a device-resident collision
    /// set.
    /// @param collisions The device-resident collisions.
    /// @param vertices The device-resident vertex positions.
    /// @return The sum of the potential over all collisions.
    double operator()(
        const DeviceNormalCollisions& collisions,
        const DeviceVertices& vertices) const;

    /// @brief Compute the gradient of the barrier potential for a
    /// device-resident collision set.
    /// @param collisions The device-resident collisions.
    /// @param vertices The device-resident vertex positions.
    /// @return The gradient of the potential w.r.t. the vertex positions.
    Eigen::VectorXd gradient(
        const DeviceNormalCollisions& collisions,
        const DeviceVertices& vertices) const;

    /// @brief Compute the hessian of the barrier potential for a
    /// device-resident collision set.
    /// @param collisions The device-resident collisions.
    /// @param vertices The device-resident vertex positions.
    /// @param project_hessian_to_psd Method of projecting the local hessians
    /// to the positive semi-definite cone (applied on the CPU during
    /// assembly).
    /// @return The hessian of the potential w.r.t. the vertex positions.
    Eigen::SparseMatrix<double> hessian(
        const DeviceNormalCollisions& collisions,
        const DeviceVertices& vertices,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

private:
    /// @brief The activation distance of the barrier.
    double m_dhat;

    /// @brief The stiffness of the barrier.
    double m_stiffness;

    /// @brief Whether to use the physical barrier.
    bool m_use_physical_barrier = false;

    /// @brief The barrier function type used to compute the potential.
    BarrierType m_barrier_type = BarrierType::CLAMPED_LOG;
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
