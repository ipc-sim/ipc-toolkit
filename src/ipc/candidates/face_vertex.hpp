#pragma once

#include <ipc/candidates/continuous_collision_candidate.hpp>
#include <ipc/distance/distance_type.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

struct FaceVertexCandidate : ContinuousCollisionCandidate {
    FaceVertexCandidate(long face_index, long vertex_index);

    int num_vertices() const { return 4; };

    std::array<long, 4>
    vertex_indices(const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const
    {
        return { { vertex_index, //
                   F(face_index, 0), F(face_index, 1), F(face_index, 2) } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const PointTriangleDistanceType dtype =
            PointTriangleDistanceType::AUTO) const;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const PointTriangleDistanceType dtype =
            PointTriangleDistanceType::AUTO) const;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const PointTriangleDistanceType dtype =
            PointTriangleDistanceType::AUTO) const;

    // ------------------------------------------------------------------------

    bool
    ccd(const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const double tolerance = DEFAULT_CCD_TOLERANCE,
        const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS,
        const double conservative_rescaling =
            DEFAULT_CCD_CONSERVATIVE_RESCALING) const override;

    void print_ccd_query(
        const Eigen::MatrixXd& V0,
        const Eigen::MatrixXd& V1,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    // ------------------------------------------------------------------------

    bool operator==(const FaceVertexCandidate& other) const;
    bool operator!=(const FaceVertexCandidate& other) const;
    /// @brief Compare FaceVertexCandidate for sorting.
    bool operator<(const FaceVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexCandidate& fv)
    {
        return H::combine(std::move(h), fv.face_index, fv.vertex_index);
    }

    long face_index;
    long vertex_index;
};

} // namespace ipc
