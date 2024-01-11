#pragma once

#include <ipc/collisions/edge_edge.hpp>
#include <ipc/collision_mesh.hpp>

namespace ipc {

template <int dim_>
class SmoothEdgeEdgeCollision : public EdgeEdgeCollision {
public:
    // using EdgeEdgeCollision::EdgeEdgeCollision;
    SmoothEdgeEdgeCollision(
        const long _edge0_id,
        const long _edge1_id,
        const CollisionMesh &mesh)
    : EdgeEdgeCollision(_edge0_id, _edge1_id, 0.)
    { 
        vertices = vertex_ids(mesh.edges(), mesh.faces());
    }
    SmoothEdgeEdgeCollision(
        const long _edge0_id,
        const long _edge1_id,
        const double _eps_x,
        std::array<long, 4> _vertices)
    : EdgeEdgeCollision(_edge0_id, _edge1_id, _eps_x), vertices(_vertices)
    { }

    double compute_distance(const VectorMax12d& positions) const override;

    VectorMax12d
    compute_distance_gradient(const VectorMax12d& positions) const override;

    MatrixMax12d
    compute_distance_hessian(const VectorMax12d& positions) const override;

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

private:
    Vector12d positions_to_3d(const VectorMax12d& positions) const;
    
    template <typename scalar> 
    scalar evaluate_quadrature(const VectorMax12d& positions, const ParameterType &params) const;

    std::array<long, 4> vertices;
};

}