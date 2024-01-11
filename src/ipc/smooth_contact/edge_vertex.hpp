#pragma once
#include <ipc/collisions/edge_vertex.hpp>

namespace ipc {

class SmoothEdgeVertexCollision : public EdgeVertexCollision {
public:
    SmoothEdgeVertexCollision(
        const long _edge_id,
        const long _vertex_id,
        const double _weight,
        const Eigen::SparseVector<double>& _weight_gradient)
        : EdgeVertexCollision(_edge_id, _vertex_id, _weight, _weight_gradient)
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
};

} // namespace ipc
