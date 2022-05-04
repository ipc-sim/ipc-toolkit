#pragma once

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct FaceVertexConstraint : FaceVertexCandidate, CollisionConstraint {
    FaceVertexConstraint(long face_index, long vertex_index);
    FaceVertexConstraint(const FaceVertexCandidate& candidate);

    int num_vertices() const override { return 4; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, //
                   F(face_index, 0), F(face_index, 1), F(face_index, 2) } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;
};

} // namespace ipc

namespace std {
template <> struct hash<ipc::FaceVertexConstraint> {
    size_t operator()(ipc::FaceVertexConstraint const& fv) const noexcept;
};
} // namespace std
