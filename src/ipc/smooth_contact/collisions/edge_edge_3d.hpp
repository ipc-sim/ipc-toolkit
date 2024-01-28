#pragma once

#include "smooth_collision.hpp"
#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/math.hpp>

namespace ipc {

class SmoothEdgeEdge3Collision : public SmoothCollision<max_vert_3d> {
public:
    SmoothEdgeEdge3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const std::array<double, 2> &dhats_,
        const Eigen::MatrixXd &V);
    virtual ~SmoothEdgeEdge3Collision() 
    {
    }

    int ndofs() const override
    {
        return num_vertices() * 3;
    }

    int num_vertices() const override
    {
        return 8;
    }

    double compute_distance(const Vector<double, -1, 3*max_vert_3d>& positions) const override;

    double operator()(const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const override;

    Vector<double, -1, 3*max_vert_3d> gradient(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params) const override;

    MatrixMax<double, 3*max_vert_3d, 3*max_vert_3d> hessian(
        const Vector<double, -1, 3*max_vert_3d>& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

    std::string name() const override { return "edge-edge"; }

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector<double, 24>& positions, ParameterType params) const;

    bool compute_types(const Vector<double, 24>& positions, ParameterType params); // return true if the potential is nonzero, return false if the potential is zero and can be skipped

    Eigen::Matrix<int, 4, 3> face_to_vertex; // stores the local vertex ids for each vertex on each face

    EdgeEdgeDistanceType dtype;
    std::array<ORIENTATION_TYPES, 2> otypes;
    std::array<HEAVISIDE_TYPE, 4> mtypes;
};

}