#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/collisions/collision_constraint.hpp>

namespace ipc {

class EdgeVertexConstraint : public EdgeVertexCandidate,
                             public CollisionConstraint {
public:
    using EdgeVertexCandidate::EdgeVertexCandidate;

    EdgeVertexConstraint(const EdgeVertexCandidate& candidate)
        : EdgeVertexCandidate(candidate)
    {
    }

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexConstraint& ev)
    {
        return AbslHashValue(
            std::move(h), static_cast<const EdgeVertexCandidate&>(ev));
    }

protected:
    PointEdgeDistanceType known_dtype() const override
    {
        // The distance type is known because of Constraints::build()
        return PointEdgeDistanceType::P_E;
    }
};

} // namespace ipc
