#pragma once

#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class VertexVertexCandidate : public ContinuousCollisionCandidate {
public:
    VertexVertexCandidate(long vertex0_id, long vertex1_id);

    int num_vertices() const override { return 2; };

    /// @brief Get the indices of the vertices
    /// @param edges edge matrix of mesh
    /// @param faces face matrix of mesh
    /// @return List of vertex indices
    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { vertex0_id, vertex1_id, -1, -1 } };
    }

    // ------------------------------------------------------------------------

    std::ostream& write_ccd_query(
        std::ostream& out,
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override;

    // ------------------------------------------------------------------------

    bool operator==(const VertexVertexCandidate& other) const;
    bool operator!=(const VertexVertexCandidate& other) const;
    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const VertexVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexCandidate& vv)
    {
        long min_vi = std::min(vv.vertex0_id, vv.vertex1_id);
        long max_vi = std::max(vv.vertex0_id, vv.vertex1_id);
        return H::combine(std::move(h), min_vi, max_vi);
    }

    // ------------------------------------------------------------------------

    /// @brief ID of the first vertex
    long vertex0_id;
    /// @brief ID of the second vertex
    long vertex1_id;

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
};

} // namespace ipc
