#pragma once

#include "smooth_collision.hpp"

namespace ipc {

class SmoothEdgeEdgeCollision : public SmoothCollision<4> {
    constexpr static int dim = 2;
public:
    SmoothEdgeEdgeCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh)
    : SmoothCollision(primitive0_, primitive1_, mesh)
    { 
        vertices = vertex_ids(mesh.edges(), mesh.faces());
    }

    int num_vertices() const override
    {
        return 4;
    }

    std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { edges(primitive0, 0), edges(primitive0, 1),
                   edges(primitive1, 0), edges(primitive1, 1) } };
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

    void set_adaptive_dhat(const CollisionMesh &mesh, const double &dhat) override
    {
        dhat0 = std::min(dhat, mesh.min_distance_in_rest_config(primitive0));
        dhat1 = std::min(dhat, mesh.min_distance_in_rest_config(primitive1));
    }

private:
    Vector12d positions_to_3d(const Vector8d& positions) const;
    
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector8d& positions, ParameterType params) const;
};

}