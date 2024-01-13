#pragma once

#include "smooth_collision.hpp"

namespace ipc {

class SmoothEdgeEdge3Collision : public SmoothCollision<8> {
public:
    SmoothEdgeEdge3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh);
    virtual ~SmoothEdgeEdge3Collision() { }

    int num_vertices() const override
    {
        return 8;
    }

    std::array<long, 8> vertex_ids(
        const Eigen::MatrixXi& _edges, const Eigen::MatrixXi& _faces) const override;

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
private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector<double, 24>& positions, const ParameterType &params) const;

    Eigen::Matrix<int, 4, 3> face_to_vertex; // stores the local vertex ids for each vertex on each face
};

}