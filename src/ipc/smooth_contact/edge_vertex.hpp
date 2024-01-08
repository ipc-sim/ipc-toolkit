#pragma once
#include <ipc/collisions/edge_vertex.hpp>

namespace ipc {

class SmoothEdgeVertexCollision : public EdgeVertexCollision {
public:
    using EdgeVertexCollision::EdgeVertexCollision;

    SmoothEdgeVertexCollision(const EdgeVertexCandidate& candidate)
        : EdgeVertexCollision(candidate)
    {
    }

    SmoothEdgeVertexCollision(
        const long _edge_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : EdgeVertexCollision(_edge_id, _vertex_id, _weight, _weight_gradient)
    {
    }

    SmoothEdgeVertexCollision(
        const long _edge_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient,
        const double _eps)
        : EdgeVertexCollision(_edge_id, _vertex_id, _weight, _weight_gradient), local_eps(_eps)
    {
    }

    PointEdgeDistanceType known_dtype() const override
    {
        return PointEdgeDistanceType::AUTO;
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

    protected:
        double local_eps = -1;
};

} // namespace ipc
