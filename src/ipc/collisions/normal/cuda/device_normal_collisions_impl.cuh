// Definitions of the pimpl structs of DeviceVertices and
// DeviceNormalCollisions. This header is CUDA-only and must be included from
// the ipc::cuda implementation files (.cu) exclusively.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/collisions/normal/cuda/device_normal_collisions.hpp>
#include <ipc/distance/distance_type.hpp>

#include <thrust/device_vector.h>

namespace ipc::cuda {

struct DeviceVertices::Impl {
    /// @brief Vertex positions (rowwise |V| × 3, row-major).
    thrust::device_vector<double> positions;
};

struct DeviceNormalCollisions::Impl {
    // Flat SoA arrays over all collisions in the flattened order
    // vv | ev | ee | fv (matching NormalCollisions::operator[]).
    // Unused vertex id entries are -1.

    /// @brief First vertex id of each collision stencil.
    thrust::device_vector<index_t> vertex_id_0;
    /// @brief Second vertex id of each collision stencil.
    thrust::device_vector<index_t> vertex_id_1;
    /// @brief Third vertex id of each collision stencil (-1 for vv).
    thrust::device_vector<index_t> vertex_id_2;
    /// @brief Fourth vertex id of each collision stencil (-1 for vv/ev).
    thrust::device_vector<index_t> vertex_id_3;

    /// @brief Quadrature weight of each collision (may be negative).
    thrust::device_vector<double> weight;
    /// @brief Minimum separation distance of each collision.
    thrust::device_vector<double> dmin;

    // Edge-edge-only arrays (size = n_ee; indexed by i - ee_begin()).

    /// @brief Mollifier activation threshold of each edge-edge collision.
    thrust::device_vector<double> ee_eps_x;
    /// @brief Runtime EdgeEdgeDistanceType of each edge-edge collision.
    thrust::device_vector<uint8_t> ee_dtype;

    // Per-type range sizes (flattened order vv | ev | ee | fv).

    size_t n_vv = 0; ///< Number of vertex-vertex collisions.
    size_t n_ev = 0; ///< Number of edge-vertex collisions.
    size_t n_ee = 0; ///< Number of edge-edge collisions.
    size_t n_fv = 0; ///< Number of face-vertex collisions.

    size_t vv_begin() const { return 0; }
    size_t ev_begin() const { return n_vv; }
    size_t ee_begin() const { return n_vv + n_ev; }
    size_t fv_begin() const { return n_vv + n_ev + n_ee; }
    size_t size() const { return n_vv + n_ev + n_ee + n_fv; }

    /// @brief Host-side mirror of the vertex ids (for CPU-side assembly).
    std::vector<std::array<index_t, 4>> host_vertex_ids;
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
