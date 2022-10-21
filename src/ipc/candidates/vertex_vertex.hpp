#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

struct VertexVertexCandidate {
    VertexVertexCandidate(long vertex0_index, long vertex1_index);

    int num_vertices() const { return 2; };

    /// @brief Get the indices of the vertices
    /// @param E edge matrix of mesh
    /// @param F face matrix of mesh
    /// @return List of vertex indices
    std::array<long, 4>
    vertex_indices(const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const
    {
        return { { vertex0_index, vertex1_index, -1, -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const;

    VectorMax6d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const;

    MatrixMax6d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const;

    // ------------------------------------------------------------------------

    bool operator==(const VertexVertexCandidate& other) const;
    bool operator!=(const VertexVertexCandidate& other) const;
    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const VertexVertexCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexCandidate& vv)
    {
        long min_vi = std::min(vv.vertex0_index, vv.vertex1_index);
        long max_vi = std::max(vv.vertex0_index, vv.vertex1_index);
        return H::combine(std::move(h), min_vi, max_vi);
    }

    // ------------------------------------------------------------------------

    long vertex0_index;
    long vertex1_index;
};

} // namespace ipc
