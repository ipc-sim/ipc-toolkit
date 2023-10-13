#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class FaceVertexCandidate : public ContinuousCollisionCandidate {
public:
    FaceVertexCandidate(long face_id, long vertex_id);

    int num_vertices() const override { return 4; };

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { vertex_id, faces(face_id, 0), faces(face_id, 1),
                   faces(face_id, 2) } };
    }

    // ------------------------------------------------------------------------

    void print_ccd_query(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override;

    // ------------------------------------------------------------------------

    bool operator==(const FaceVertexCandidate& other) const;
    bool operator!=(const FaceVertexCandidate& other) const;
    /// @brief Compare FaceVertexCandidate for sorting.
    bool operator<(const FaceVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexCandidate& fv)
    {
        return H::combine(std::move(h), fv.face_id, fv.vertex_id);
    }

    long face_id;   ///< @brief ID of the face
    long vertex_id; ///< @brief ID of the vertex

    using CollisionStencil::compute_distance;
    using CollisionStencil::compute_distance_gradient;
    using CollisionStencil::compute_distance_hessian;
    using ContinuousCollisionCandidate::ccd;

protected:
    double compute_distance(const VectorMax12d& positions) const override;

    VectorMax12d
    compute_distance_gradient(const VectorMax12d& positions) const override;

    MatrixMax12d
    compute_distance_hessian(const VectorMax12d& positions) const override;

    bool
    ccd(const VectorMax12d& vertices_t0,
        const VectorMax12d& vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const override;

    virtual PointTriangleDistanceType known_dtype() const
    {
        return PointTriangleDistanceType::AUTO;
    }
};

} // namespace ipc
