#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>

namespace ipc {

class EdgeVertexNormalCollision : public EdgeVertexCandidate,
                                  public NormalCollision {
public:
    using EdgeVertexCandidate::EdgeVertexCandidate;

    EdgeVertexNormalCollision(const EdgeVertexCandidate& candidate)
        : EdgeVertexCandidate(candidate)
    {
    }

    EdgeVertexNormalCollision(
        const long _edge_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : EdgeVertexCandidate(_edge_id, _vertex_id)
        , NormalCollision(_weight, _weight_gradient)
    {
    }

    PointEdgeDistanceType known_dtype() const override
    {
        // The distance type is known because of NormalCollisions::build()
        return PointEdgeDistanceType::P_E;
    }

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexNormalCollision& ev)
    {
        return AbslHashValue(
            std::move(h), static_cast<const EdgeVertexCandidate&>(ev));
    }
};

} // namespace ipc
