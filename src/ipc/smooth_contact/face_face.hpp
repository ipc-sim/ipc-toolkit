#pragma once

#include "smooth_collision.hpp"
#include <iostream>

namespace ipc {

class SmoothFaceFaceCollision : public SmoothCollision<6> {
public:
    SmoothFaceFaceCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh)
    : SmoothCollision<6>(primitive0_, primitive1_, mesh)
    { 
        vertices = vertex_ids(mesh.edges(), mesh.faces());
    }
    virtual ~SmoothFaceFaceCollision() { }

    std::array<long, 6> vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const override;

    double operator()(const Vector<double, -1, 18>& positions, 
        const ParameterType &params) const override;

    Vector<double, -1, 18> gradient(
        const Vector<double, -1, 18>& positions, 
        const ParameterType &params) const override;

    MatrixMax<double, 18, 18> hessian(
        const Vector<double, -1, 18>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

    void set_adaptive_dhat(const CollisionMesh &mesh, const double &dhat) override
    {
        // dhat0 = std::min(dhat, mesh.min_distance_in_rest_config(edge0_id));
        // dhat1 = std::min(dhat, mesh.min_distance_in_rest_config(edge1_id));
    }

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector<double, -1, 18>& positions, const ParameterType &params) const;
};

}