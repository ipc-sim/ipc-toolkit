#pragma once

#include <vector>

namespace ipc {

struct VertexVertexCandidate {
    VertexVertexCandidate(long vertex0_index, long vertex1_index);

    bool operator==(const VertexVertexCandidate& other) const;
    bool operator!=(const VertexVertexCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const VertexVertexCandidate& other) const;

    long vertex0_index;
    long vertex1_index;
};

struct EdgeVertexCandidate {
    EdgeVertexCandidate(long edge_index, long vertex_index);

    bool operator==(const EdgeVertexCandidate& other) const;
    bool operator!=(const EdgeVertexCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeVertexCandidates for sorting.
    bool operator<(const EdgeVertexCandidate& other) const;

    long edge_index;
    long vertex_index;
};

struct EdgeEdgeCandidate {
    EdgeEdgeCandidate(long edge0_index, long edge1_index);

    bool operator==(const EdgeEdgeCandidate& other) const;
    bool operator!=(const EdgeEdgeCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeEdgeCandidates for sorting.
    bool operator<(const EdgeEdgeCandidate& other) const;

    long edge0_index;
    long edge1_index;
};

/// @brief Candidate for <b>intersection</b> between edge and face.
///
/// Not included in Candidates because it is not a collision candidate.
struct EdgeFaceCandidate {
    EdgeFaceCandidate(long edge_index, long face_index);

    bool operator==(const EdgeFaceCandidate& other) const;
    bool operator!=(const EdgeFaceCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare EdgeFaceCandidate for sorting.
    bool operator<(const EdgeFaceCandidate& other) const;

    long edge_index;
    long face_index;
};

struct FaceVertexCandidate {
    FaceVertexCandidate(long face_index, long vertex_index);

    bool operator==(const FaceVertexCandidate& other) const;
    bool operator!=(const FaceVertexCandidate& other) const
    {
        return !(*this == other);
    }

    /// @brief Compare FaceVertexCandidate for sorting.
    bool operator<(const FaceVertexCandidate& other) const;

    long face_index;
    long vertex_index;
};

struct Candidates {
    std::vector<EdgeVertexCandidate> ev_candidates;
    std::vector<EdgeEdgeCandidate> ee_candidates;
    std::vector<FaceVertexCandidate> fv_candidates;

    size_t size() const;

    bool empty() const;

    void clear();
};

} // namespace ipc
