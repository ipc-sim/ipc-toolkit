#include "edge_face.hpp"

namespace ipc {

EdgeFaceCandidate::EdgeFaceCandidate(long edge_id, long face_id)
    : edge_id(edge_id)
    , face_id(face_id)
{
}

bool EdgeFaceCandidate::operator==(const EdgeFaceCandidate& other) const
{
    return edge_id == other.edge_id && face_id == other.face_id;
}

bool EdgeFaceCandidate::operator!=(const EdgeFaceCandidate& other) const
{
    return !(*this == other);
}

bool EdgeFaceCandidate::operator<(const EdgeFaceCandidate& other) const
{
    if (edge_id == other.edge_id) {
        return face_id < other.face_id;
    }
    return edge_id < other.edge_id;
}

} // namespace ipc
