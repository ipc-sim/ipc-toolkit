#pragma once

#include <ipc/ccd/ccd.hpp>
// #include <ipc/candidates/collision_stencil.hpp>
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
    virtual ~FaceFaceCandidate() = default;

    bool operator==(const FaceFaceCandidate& other) const
    {
        return (face0_id == other.face0_id && face1_id == other.face1_id) || 
                (face0_id == other.face1_id && face1_id == other.face0_id);
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

    /// @brief ID of the face
    long face0_id, face1_id;
};

} // namespace ipc
