#pragma once

#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class FaceVertexCandidate : public ContinuousCollisionCandidate {
public:
    FaceVertexCandidate(long face_id, long vertex_id);

    // ------------------------------------------------------------------------
    // CollisionStencil

    int num_vertices() const override { return 4; };

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { vertex_id, faces(face_id, 0), faces(face_id, 1),
                   faces(face_id, 2) } };
    }

    using CollisionStencil::compute_coefficients;
    using CollisionStencil::compute_distance;
    using CollisionStencil::compute_distance_gradient;
    using CollisionStencil::compute_distance_hessian;

    double compute_distance(const VectorMax12d& positions) const override;

    VectorMax12d
    compute_distance_gradient(const VectorMax12d& positions) const override;

    MatrixMax12d
    compute_distance_hessian(const VectorMax12d& positions) const override;

    VectorMax4d
    compute_coefficients(const VectorMax12d& positions) const override;

    // ------------------------------------------------------------------------
    // ContinuousCollisionCandidate

    bool
    ccd(const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const override;

    // ------------------------------------------------------------------------

    virtual PointTriangleDistanceType known_dtype() const
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
    long face_id;
    /// @brief ID of the vertex
    long vertex_id;
};

} // namespace ipc
