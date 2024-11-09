#pragma once

#include <ipc/candidates/face_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>

namespace ipc {

class FaceVertexNormalCollision : public FaceVertexCandidate,
                                  public NormalCollision {
public:
    using FaceVertexCandidate::FaceVertexCandidate;

    FaceVertexNormalCollision(const FaceVertexCandidate& candidate)
        : FaceVertexCandidate(candidate)
    {
    }

    FaceVertexNormalCollision(
        const long _face_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : FaceVertexCandidate(_face_id, _vertex_id)
        , NormalCollision(_weight, _weight_gradient)
    {
    }

    PointTriangleDistanceType known_dtype() const override
    {
        // The distance type is known because of NormalCollisions::build()
        return PointTriangleDistanceType::P_T;
    }

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexNormalCollision& fv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const FaceVertexCandidate&>(fv));
    }
};

} // namespace ipc
