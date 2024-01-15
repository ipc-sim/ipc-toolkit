#pragma once

#include "smooth_collision.hpp"
#include <iostream>

namespace ipc {

class SmoothFaceFaceCollision : public SmoothCollision<8> {
public:
    SmoothFaceFaceCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh);
    virtual ~SmoothFaceFaceCollision() { }

    int num_vertices() const override
    {
        return 6;
    }

    std::array<long, 8> vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const override;

    double operator()(const Vector<double, -1, 24>& positions, 
        const ParameterType &params) const override;

    Vector<double, -1, 24> gradient(
        const Vector<double, -1, 24>& positions, 
        const ParameterType &params) const override;

    MatrixMax<double, 24, 24> hessian(
        const Vector<double, -1, 24>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

    void set_adaptive_dhat(const CollisionMesh &mesh, const double &dhat) override
    {
        // dhat0 = std::min(dhat, mesh.min_distance_in_rest_config(edge0_id));
        // dhat1 = std::min(dhat, mesh.min_distance_in_rest_config(edge1_id));
    }

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector<double, 18>& positions, const ParameterType &params) const;
};

}