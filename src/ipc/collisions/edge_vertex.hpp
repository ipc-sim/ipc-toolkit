#pragma once

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/collisions/collision.hpp>

namespace ipc {

class EdgeVertexCollision : public EdgeVertexCandidate, public Collision {
public:
    using EdgeVertexCandidate::EdgeVertexCandidate;

    EdgeVertexCollision(const EdgeVertexCandidate& candidate)
        : EdgeVertexCandidate(candidate)
    {
    }

    EdgeVertexCollision(
        const long _edge_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : EdgeVertexCandidate(_edge_id, _vertex_id)
        , Collision(_weight, _weight_gradient)
    {
    }

    PointEdgeDistanceType known_dtype() const override
    {
        // The distance type is known because of Collisions::build()
        return PointEdgeDistanceType::P_E;
    }

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexCollision& ev)
    {
        return AbslHashValue(
            std::move(h), static_cast<const EdgeVertexCandidate&>(ev));
    }


    double operator()(const VectorMax12d& positions, 
        const ParameterType &params) const override;

    VectorMax12d gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const override;

    MatrixMax12d hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

};

} // namespace ipc
