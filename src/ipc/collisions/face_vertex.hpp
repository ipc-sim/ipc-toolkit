#pragma once

#include <ipc/candidates/face_vertex.hpp>
#include <ipc/collisions/collision.hpp>

namespace ipc {

class FaceVertexCollision : public FaceVertexCandidate, public Collision {
public:
    using FaceVertexCandidate::FaceVertexCandidate;

    FaceVertexCollision(const FaceVertexCandidate& candidate)
        : FaceVertexCandidate(candidate)
    {
    }

    FaceVertexCollision(
        const long _face_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : FaceVertexCandidate(_face_id, _vertex_id)
        , Collision(_weight, _weight_gradient)
    {
    }

    PointTriangleDistanceType known_dtype() const override
    {
        // The distance type is known because of Collisions::build()
        return PointTriangleDistanceType::P_T;
    }

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexCollision& fv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const FaceVertexCandidate&>(fv));
    }
};

} // namespace ipc
