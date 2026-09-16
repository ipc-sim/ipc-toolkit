#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>

namespace ipc {

class EdgeVertexNormalCollision : public EdgeVertexStencil,
                                  public NormalCollision {
public:
    using EdgeVertexStencil::EdgeVertexStencil;

    EdgeVertexNormalCollision(const EdgeVertexCandidate& candidate)
        : EdgeVertexStencil(candidate)
    {
    }

    EdgeVertexNormalCollision(
        const index_t _edge_id,
        const index_t _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : EdgeVertexStencil(_edge_id, _vertex_id)
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
