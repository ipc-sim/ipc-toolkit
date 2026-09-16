#pragma once

#include <ipc/candidates/face_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>

namespace ipc {

class FaceVertexNormalCollision : public FaceVertexStencil,
                                  public NormalCollision {
public:
    using FaceVertexStencil::FaceVertexStencil;

    FaceVertexNormalCollision(const FaceVertexCandidate& candidate)
        : FaceVertexStencil(candidate)
    {
    }

    FaceVertexNormalCollision(
        const index_t _face_id,
        const index_t _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : FaceVertexStencil(_face_id, _vertex_id)
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
