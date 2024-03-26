#pragma once

#include <Eigen/Core>

namespace ipc {

/// @brief Candidate for collision between two faces.
/// @note This collision candidate is not needed for flat triangles because face-vertex and edge-edge collisions are sufficient.
/// @note This may be useful for nonlinear triangles in the future.
class FaceFaceCandidate {
public:
    FaceFaceCandidate(long face0_id, long face1_id);

    bool operator==(const FaceFaceCandidate& other) const;
    bool operator!=(const FaceFaceCandidate& other) const;
    /// @brief Compare FaceFaceCandidate for sorting.
    bool operator<(const FaceFaceCandidate& other) const;

    template <typename H>
    friend H AbslHashValue(H h, const FaceFaceCandidate& ff)
    {
        long min_fi = std::min(ff.face0_id, ff.face1_id);
        long max_fi = std::max(ff.face0_id, ff.face1_id);
        return H::combine(std::move(h), min_fi, max_fi);
    }

    /// @brief ID of the first face.
    long face0_id;
    /// @brief ID of the second face.
    long face1_id;
};

} // namespace ipc
