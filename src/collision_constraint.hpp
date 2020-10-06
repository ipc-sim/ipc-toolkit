#pragma once

#include <ipc/spatial_hash/collision_candidate.hpp>

namespace ipc {

struct VertexVertexConstraint {
    VertexVertexConstraint(long vertex0_index, long vertex1_index);

    bool operator==(const VertexVertexConstraint& other) const;

    /// @brief Compare VertexVertexConstraints for sorting.
    bool operator<(const VertexVertexConstraint& other) const;

    long vertex0_index;
    long vertex1_index;
    long multiplicity = 1;
};

struct EdgeVertexConstraint : EdgeVertexCandidate {
    EdgeVertexConstraint(long edge_index, long vertex_index);
    EdgeVertexConstraint(const EdgeVertexCandidate& candidate);

    long multiplicity = 1;
};

struct EdgeEdgeConstraint : EdgeEdgeCandidate {
    EdgeEdgeConstraint(long edge0_index, long edge1_index, double eps_x);
    EdgeEdgeConstraint(const EdgeEdgeCandidate& candidate, double eps_x);

    double eps_x;
};

struct FaceVertexConstraint : FaceVertexCandidate {
    FaceVertexConstraint(long face_index, long vertex_index);
    FaceVertexConstraint(const FaceVertexCandidate& candidate);
};

struct Constraints {
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;

    size_t size() const;

    size_t num_constraints() const;

    void clear();
};

} // namespace ipc
