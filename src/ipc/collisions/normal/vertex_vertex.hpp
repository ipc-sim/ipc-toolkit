#pragma once

#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

class VertexVertexNormalCollision : public VertexVertexStencil,
                                    public NormalCollision {
public:
    using VertexVertexStencil::VertexVertexStencil;

    VertexVertexNormalCollision(const VertexVertexCandidate& candidate)
        : VertexVertexStencil(candidate)
    {
    }

    VertexVertexNormalCollision(
        const index_t _vertex0_id,
        const index_t _vertex1_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : VertexVertexStencil(_vertex0_id, _vertex1_id)
        , NormalCollision(_weight, _weight_gradient)
    {
    }

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexNormalCollision& vv)
    {
        return AbslHashValue(
            std::move(h), static_cast<const VertexVertexCandidate&>(vv));
    }
};

} // namespace ipc
