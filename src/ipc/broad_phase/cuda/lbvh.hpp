#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/broad_phase/lbvh.hpp> // ipc::LBVH::Node / Nodes / RightmostLeaves

#include <memory>

namespace ipc::cuda {

/// @brief GPU Linear Bounding Volume Hierarchy (LBVH) broad phase.
///
/// A first-class GPU counterpart to ipc::LBVH: it builds the vertex/edge/face
/// AABBs and their BVHs, and runs the traversal and mesh-connectivity
/// filtering, entirely on the device. Construction uses Morton codes + the
/// Apetrei 2014 single-pass bottom-up build and reuses the 32-byte
/// ipc::LBVH::Node layout, so the device tree can be copied back to the host
/// and validated against — or traversed by — the CPU code.
///
/// Detection runs the BVH descent (AABB overlap + triangular dedup) and the
/// connectivity (shared-vertex) exclusion on the device. The user vertex filter
/// (can_vertices_collide) is honored on the device only when it is the default
/// accept-all filter; a non-trivial filter is applied on the host while
/// materializing the device-emitted (connectivity-filtered) candidates. Either
/// way the output matches the CPU ipc::LBVH exactly for any filter.
class LBVH : public ipc::BroadPhase {
public:
    LBVH();
    ~LBVH();

    LBVH(LBVH&&) noexcept;
    LBVH& operator=(LBVH&&) noexcept;
    LBVH(const LBVH&) = delete;
    LBVH& operator=(const LBVH&) = delete;

    /// @brief Non-owning view of device-resident candidate pairs (SoA). The
    /// pointers address device memory owned by this LBVH and are valid until
    /// the next detect_*_device() call on the same type or clear().
    struct DeviceCandidateView {
        const int32_t* a = nullptr; ///< Device pointer to the first ids.
        const int32_t* b = nullptr; ///< Device pointer to the second ids.
        size_t size = 0;            ///< Number of candidate pairs.
    };

    /// @brief Get the name of the broad phase method.
    std::string name() const override { return "LBVH (CUDA)"; }

    using ipc::BroadPhase::build;

    /// @brief Build the broad phase for static collision detection.
    /// @param vertices Vertex positions (rowwise, |V| × 3).
    /// @param edges Collision mesh edges.
    /// @param faces Collision mesh faces.
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0) override;

    /// @brief Build the broad phase for continuous collision detection.
    /// @param vertices_t0 Starting vertex positions (rowwise, |V| × 3).
    /// @param vertices_t1 Ending vertex positions (rowwise, |V| × 3).
    /// @param edges Collision mesh edges.
    /// @param faces Collision mesh faces.
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
        Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const double inflation_radius = 0) override;

    /// @brief Build the broad phase from precomputed host vertex AABBs.
    /// The vertex boxes are uploaded; edge/face boxes and all BVHs are built on
    /// the device.
    /// @param vertex_boxes Precomputed vertex AABBs.
    /// @param edges Collision mesh edges.
    /// @param faces Collision mesh faces.
    /// @param dim Dimension of the simulation (2 or 3).
    void build(
        const AABBs& vertex_boxes,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const uint8_t dim) override;

    /// @brief Clear any built data.
    void clear() override;

    // ------------------------------------------------------------------
    // BroadPhase interface (host-materializing). The BVH descent (AABB overlap
    // + triangular dedup) and the mesh-connectivity (shared-vertex) exclusion
    // both run on the device. The user vertex filter is applied on the host
    // only when it is not accept-all (see can_*_collide); the output matches
    // the CPU ipc::LBVH exactly for any filter.

    void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const override;
    void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const override;
    void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const override;
    void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const override;
    void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const override;
    void detect_face_face_candidates(
        std::vector<FaceFaceCandidate>& candidates) const override;

    // ------------------------------------------------------------------
    // Device-resident candidate accessors (GPU-native pipeline). Each runs the
    // filtered traversal and returns a view of the connectivity-filtered pairs
    // left on the device. For the default (accept-all) vertex filter the view
    // is the exact candidate set; otherwise it is a connectivity-filtered
    // superset the caller must trim with the user vertex filter.

    DeviceCandidateView detect_vertex_vertex_candidates_device() const;
    DeviceCandidateView detect_edge_vertex_candidates_device() const;
    DeviceCandidateView detect_edge_edge_candidates_device() const;
    DeviceCandidateView detect_face_vertex_candidates_device() const;
    DeviceCandidateView detect_edge_face_candidates_device() const;
    DeviceCandidateView detect_face_face_candidates_device() const;

    // ------------------------------------------------------------------
    // Sizes (cheap; no device->host node copy).

    /// @brief Number of nodes in the vertex BVH (2 * n_leaves - 1, or 0).
    size_t num_vertex_nodes() const;
    /// @brief Number of nodes in the edge BVH (2 * n_leaves - 1, or 0).
    size_t num_edge_nodes() const;
    /// @brief Number of nodes in the face BVH (2 * n_leaves - 1, or 0).
    size_t num_face_nodes() const;

    // ------------------------------------------------------------------
    // Debug / validation: copy the device trees back to the host.

    /// @brief Copy the vertex BVH back to the host.
    void vertex_nodes_to_host(
        ipc::LBVH::Nodes& nodes,
        ipc::LBVH::RightmostLeaves& rightmost_leaves) const;
    /// @brief Copy the edge BVH back to the host.
    void edge_nodes_to_host(
        ipc::LBVH::Nodes& nodes,
        ipc::LBVH::RightmostLeaves& rightmost_leaves) const;
    /// @brief Copy the face BVH back to the host.
    void face_nodes_to_host(
        ipc::LBVH::Nodes& nodes,
        ipc::LBVH::RightmostLeaves& rightmost_leaves) const;

    // Pimpl pattern to keep CUDA types out of this header. The Impl is
    // defined in lbvh_impl.cuh for use by the ipc::cuda implementation files
    // (.cu) only.
    struct Impl;
    const Impl& impl() const;

protected:
    // Host-side collision filters, used to trim the device-emitted candidates
    // only when the user vertex filter is not accept-all (the device already
    // excludes shared-vertex pairs). Mirror ipc::LBVH.
    bool can_edge_vertex_collide(size_t ei, size_t vi) const override;
    bool can_edges_collide(size_t eai, size_t ebi) const override;
    bool can_face_vertex_collide(size_t fi, size_t vi) const override;
    bool can_edge_face_collide(size_t ei, size_t fi) const override;
    bool can_faces_collide(size_t fai, size_t fbi) const override;

private:
    std::unique_ptr<Impl> m_impl;
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
