#pragma once

#include <ipc/ccd/ccd.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

/// @brief Candidate for <b>intersection</b> between edge and face.
///
/// Not included in Candidates because it is not a collision candidate.
class EdgeFaceCandidate {
public:
    EdgeFaceCandidate(long edge_id, long face_id);

    bool operator==(const EdgeFaceCandidate& other) const;
    bool operator!=(const EdgeFaceCandidate& other) const;
    /// @brief Compare EdgeFaceCandidate for sorting.
    bool operator<(const EdgeFaceCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const EdgeFaceCandidate& fv)
    {
        return H::combine(std::move(h), fv.edge_id, fv.face_id);
    }

    long edge_id; ///< @brief ID of the edge
    long face_id; ///< @brief ID of the face
};

} // namespace ipc
