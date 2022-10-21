#pragma once

#include <ipc/ccd/ccd.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

/// @brief Candidate for <b>intersection</b> between edge and face.
///
/// Not included in Candidates because it is not a collision candidate.
struct EdgeFaceCandidate {
    EdgeFaceCandidate(long edge_index, long face_index);

    bool operator==(const EdgeFaceCandidate& other) const;
    bool operator!=(const EdgeFaceCandidate& other) const;
    /// @brief Compare EdgeFaceCandidate for sorting.
    bool operator<(const EdgeFaceCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeFaceCandidate& fv)
    {
        return H::combine(std::move(h), fv.edge_index, fv.face_index);
    }

    long edge_index;
    long face_index;
};

} // namespace ipc
