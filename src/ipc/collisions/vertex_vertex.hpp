#pragma once

#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class VertexVertexConstraint : public VertexVertexCandidate,
                               public CollisionConstraint {
public:
    using VertexVertexCandidate::VertexVertexCandidate;

    VertexVertexConstraint(const VertexVertexCandidate& candidate)
        : VertexVertexCandidate(candidate)
    {
    }

    VertexVertexConstraint(
        const long vertex0_id,
        const long vertex1_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
        : VertexVertexCandidate(vertex0_id, vertex1_id)
        , CollisionConstraint(weight, weight_gradient)
    {
    }

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexConstraint& vv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const VertexVertexCandidate&>(vv));
    }
};

} // namespace ipc
