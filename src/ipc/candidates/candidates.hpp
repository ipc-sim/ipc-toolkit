#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/face_vertex.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

struct Candidates {
    std::vector<EdgeVertexCandidate> ev_candidates;
    std::vector<EdgeEdgeCandidate> ee_candidates;
    std::vector<FaceVertexCandidate> fv_candidates;

    Candidates() { }

    size_t size() const;

    bool empty() const;

    void clear();

    ContinuousCollisionCandidate& operator[](size_t idx);
    const ContinuousCollisionCandidate& operator[](size_t idx) const;

    bool save_obj(
        const std::string& filename,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const;
};

} // namespace ipc
