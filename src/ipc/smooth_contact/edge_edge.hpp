#pragma once

#include <ipc/collisions/edge_edge.hpp>

namespace ipc {

class SmoothEdgeEdgeCollision : public EdgeEdgeCollision {
public:
    using EdgeEdgeCollision::EdgeEdgeCollision;

    double operator()(const VectorMax12d& positions, 
        const ParameterType &params) const override;

    VectorMax12d gradient(
        const VectorMax12d& positions, 
        const ParameterType &params) const override;

    MatrixMax12d hessian(
        const VectorMax12d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

    bool is_mollified() const override { return false; }
};

}