#pragma once

#include <ipc/candidates/face_vertex.hpp>
#include <ipc/collisions/collision_constraint.hpp>

namespace ipc {

class FaceVertexConstraint : public FaceVertexCandidate,
                             public CollisionConstraint {
public:
    using FaceVertexCandidate::FaceVertexCandidate;

    FaceVertexConstraint(const FaceVertexCandidate& candidate)
        : FaceVertexCandidate(candidate)
    {
    }

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexConstraint& fv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const FaceVertexCandidate&>(fv));
    }

protected:
    PointTriangleDistanceType known_dtype() const override
    {
        // The distance type is known because of Constraints::build()
        return PointTriangleDistanceType::P_T;
    }
};

} // namespace ipc
