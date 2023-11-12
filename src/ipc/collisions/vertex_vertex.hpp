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
        const long _vertex0_id,
        const long _vertex1_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : VertexVertexCandidate(_vertex0_id, _vertex1_id)
        , CollisionConstraint(_weight, _weight_gradient)
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
