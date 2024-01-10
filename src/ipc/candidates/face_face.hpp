#pragma once

#include <ipc/ccd/ccd.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

/// @brief Candidate for <b>intersection</b> between face and face.
///
/// Not included in Candidates because it is not a collision candidate.
class FaceFaceCandidate {
public:
    FaceFaceCandidate(long _face0_id, long _face1_id)
    : face0_id(_face0_id), face1_id(_face1_id)
    { }

    bool operator==(const FaceFaceCandidate& other) const
    {
        return face0_id == other.face0_id && face1_id == other.face1_id;
    }
    bool operator!=(const FaceFaceCandidate& other) const
    {
        return !(*this == other);
    }
    /// @brief Compare FaceFaceCandidate for sorting.
    bool operator<(const FaceFaceCandidate& other) const
    {
        if (face0_id == other.face0_id) {
            return face1_id < other.face1_id;
        }
        return face0_id < other.face0_id;
    }

    template <typename H>
    friend H AbslHashValue(H h, const FaceFaceCandidate& ff)
    {
        return H::combine(std::move(h), ff.face0_id, ff.face1_id);
    }

    /// @brief ID of the face
    long face0_id, face1_id;
};

} // namespace ipc
