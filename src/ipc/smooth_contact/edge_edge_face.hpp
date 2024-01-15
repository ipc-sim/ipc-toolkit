#pragma once

#include "smooth_collision.hpp"
#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/math.hpp>

namespace ipc {

class SmoothEdgeEdge3Collision : public SmoothCollision<8> {
    SmoothEdgeEdge3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh);
public:
    SmoothEdgeEdge3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const Eigen::MatrixXd &V);
    virtual ~SmoothEdgeEdge3Collision() 
    {
    }

    int num_vertices() const override
    {
        return 8;
    }

    std::array<long, 8> vertex_ids(
        const Eigen::MatrixXi& _edges, const Eigen::MatrixXi& _faces) const override;
    
    double compute_distance(const Vector<double, -1, 24>& positions) const override;

    double operator()(const Vector<double, -1, 24>& positions, 
        const ParameterType &params) const override;

    Vector<double, -1, 24> gradient(
        const Vector<double, -1, 24>& positions, 
        const ParameterType &params) const override;

    MatrixMax<double, 24, 24> hessian(
        const Vector<double, -1, 24>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector<double, 24>& positions, const ParameterType &params) const;

    bool compute_types(const Vector<double, 24>& positions, const ParameterType &params); // return true if the potential is nonzero, return false if the potential is zero and can be skipped

    Eigen::Matrix<int, 4, 3> face_to_vertex; // stores the local vertex ids for each vertex on each face

    EdgeEdgeDistanceType dtype;
    std::array<PointEdgeDistanceType, 4> edge_dtypes;
    std::array<HEAVISIDE_TYPE, 4> tangent_types;
    std::array<HEAVISIDE_TYPE, 4> normal_types;
};

}