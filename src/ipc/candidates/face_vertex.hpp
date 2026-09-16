#pragma once

#include <ipc/candidates/candidate_vector.hpp>
#include <ipc/candidates/stencil_adapter.hpp>
#include <ipc/candidates/stencil_mixin.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

/// @brief A candidate for face-vertex collision detection.
class FaceVertexCandidate : public StencilMixin<FaceVertexCandidate> {
    friend class StencilMixin<FaceVertexCandidate>;

public:
    /// @brief Construct a candidate with indeterminate edge IDs.
    /// @note Keeping this trivial is what makes the type trivially copyable.
    FaceVertexCandidate() = default;

    FaceVertexCandidate(index_t face_id, index_t vertex_id);

    // ------------------------------------------------------------------------
    // Stencil

    constexpr static int num_vertices() { return 4; }

    /// @brief Get the vertex IDs for the face-vertex pair
    /// @param edges The edge connectivity matrix
    /// @param faces The face connectivity matrix
    /// @return An array of vertex IDs in the order: [vi, f0i, f1i, f2i]
    std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return { { vertex_id, faces(face_id, 0), faces(face_id, 1),
                   faces(face_id, 2) } };
    }

    using StencilMixin<FaceVertexCandidate>::compute_coefficients;
    using StencilMixin<FaceVertexCandidate>::compute_distance;
    using StencilMixin<FaceVertexCandidate>::compute_distance_gradient;
    using StencilMixin<FaceVertexCandidate>::compute_distance_hessian;

    double compute_distance(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointTriangleDistanceType dtype =
            PointTriangleDistanceType::AUTO) const;

    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointTriangleDistanceType dtype =
            PointTriangleDistanceType::AUTO) const;

    MatrixMax12d compute_distance_hessian(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointTriangleDistanceType dtype =
            PointTriangleDistanceType::AUTO) const;

    VectorMax4d compute_coefficients(
        Eigen::ConstRef<VectorMax12d> positions,
        const PointTriangleDistanceType dtype =
            PointTriangleDistanceType::AUTO) const;

    // ------------------------------------------------------------------------

    bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const;

    // ------------------------------------------------------------------------

    constexpr static PointTriangleDistanceType known_dtype()
    {
        return PointTriangleDistanceType::AUTO;
    }

    bool operator==(const FaceVertexCandidate& other) const;
    bool operator!=(const FaceVertexCandidate& other) const;
    /// @brief Compare FaceVertexCandidate for sorting.
    bool operator<(const FaceVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexCandidate& fv)
    {
        return H::combine(std::move(h), fv.face_id, fv.vertex_id);
    }

    // ------------------------------------------------------------------------

    /// @brief ID of the face
    index_t face_id;
    /// @brief ID of the vertex
    index_t vertex_id;

protected:
    VectorMax3d
    compute_unnormalized_normal(Eigen::ConstRef<VectorMax12d> positions) const;

    MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const;
};

/// @brief FaceVertexCandidate with the runtime-polymorphic stencil interface.
///
/// The base of the corresponding collision types.
using FaceVertexStencil =
    DTypeStencilAdapter<FaceVertexCandidate, PointTriangleDistanceType>;

} // namespace ipc
