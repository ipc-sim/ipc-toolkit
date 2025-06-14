#include "face_face.hpp"

namespace ipc {

FaceFaceCandidate::FaceFaceCandidate(index_t _face0_id, index_t _face1_id)
    : face0_id(_face0_id)
    , face1_id(_face1_id)
{
}

bool FaceFaceCandidate::operator==(const FaceFaceCandidate& other) const
{
    // (i, j) == (i, j) || (i, j) == (j, i)
    return (this->face0_id == other.face0_id
            && this->face1_id == other.face1_id)
        || (this->face0_id == other.face1_id
            && this->face1_id == other.face0_id);
}

bool FaceFaceCandidate::operator!=(const FaceFaceCandidate& other) const
{
    return !(*this == other);
}

bool FaceFaceCandidate::operator<(const FaceFaceCandidate& other) const
{
    index_t this_min = std::min(this->face0_id, this->face1_id);
    index_t other_min = std::min(other.face0_id, other.face1_id);
    if (this_min == other_min) {
        return std::max(this->face0_id, this->face1_id)
            < std::max(other.face0_id, other.face1_id);
    }
    return this_min < other_min;
}

} // namespace ipc
