#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>

#include <Eigen/Core>

#include <array>
#include <memory>
#include <vector>

namespace ipc::cuda {

/// @brief Device-resident vertex positions (3D only).
///
/// Reusable across calls: build once, then call update() to re-upload new
/// positions (e.g. once per Newton iteration) without reallocating.
class DeviceVertices {
public:
    /// @brief Upload vertex positions to the device.
    /// @param vertices Vertex positions (rowwise, |V| × 3).
    explicit DeviceVertices(Eigen::ConstRef<Eigen::MatrixXd> vertices);

    ~DeviceVertices();

    DeviceVertices(DeviceVertices&&) noexcept;
    DeviceVertices& operator=(DeviceVertices&&) noexcept;
    DeviceVertices(const DeviceVertices&) = delete;
    DeviceVertices& operator=(const DeviceVertices&) = delete;

    /// @brief Re-upload vertex positions (no reallocation if the size is
    /// unchanged).
    /// @param vertices Vertex positions (rowwise, |V| × 3).
    void update(Eigen::ConstRef<Eigen::MatrixXd> vertices);

    /// @brief Get the number of vertices.
    index_t num_vertices() const { return m_num_vertices; }

    /// @brief Get the dimension of the vertices (always 3).
    int dim() const { return 3; }

    // Pimpl pattern to keep CUDA types out of this header. The Impl is
    // defined in device_normal_collisions_impl.cuh for use by the
    // ipc::cuda implementation files (.cu) only.
    struct Impl;
    const Impl& impl() const;

private:
    std::unique_ptr<Impl> m_impl;
    index_t m_num_vertices = 0;
};

/// @brief Device-resident normal collisions in a flat structure-of-arrays
/// layout with contiguous per-type ranges (vv | ev | ee | fv), matching the
/// flattened order of NormalCollisions::operator[].
///
/// Build once per collision set (e.g. once per time step) and reuse across
/// potential evaluations.
class DeviceNormalCollisions {
public:
    /// @brief Gather and upload a built collision set to the device.
    /// @param collisions The (CPU-built) normal collisions.
    /// @param mesh The collision mesh (for edge/face connectivity).
    /// @throws If the collisions contain plane-vertex collisions or the mesh
    /// is not 3D (unsupported on the GPU).
    DeviceNormalCollisions(
        const NormalCollisions& collisions, const CollisionMesh& mesh);

    ~DeviceNormalCollisions();

    DeviceNormalCollisions(DeviceNormalCollisions&&) noexcept;
    DeviceNormalCollisions& operator=(DeviceNormalCollisions&&) noexcept;
    DeviceNormalCollisions(const DeviceNormalCollisions&) = delete;
    DeviceNormalCollisions& operator=(const DeviceNormalCollisions&) = delete;

    /// @brief Get the total number of collisions.
    size_t size() const;

    /// @brief Get if the collision set is empty.
    bool empty() const { return size() == 0; }

    /// @brief Host-side mirror of the vertex ids of every collision (in the
    /// flattened order vv | ev | ee | fv; unused entries are -1). Used for
    /// CPU-side sparse matrix assembly.
    const std::vector<std::array<index_t, 4>>& host_vertex_ids() const;

    // Pimpl pattern to keep CUDA types out of this header. The Impl is
    // defined in device_normal_collisions_impl.cuh for use by the
    // ipc::cuda implementation files (.cu) only.
    struct Impl;
    const Impl& impl() const;

private:
    std::unique_ptr<Impl> m_impl;
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
