#pragma once

#include "smooth_collision.hpp"
#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/math.hpp>

namespace ipc {

class SmoothEdgeVertex3Collision : public SmoothCollision<max_vert_3d> {
public:
    SmoothEdgeVertex3Collision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const std::array<double, 2> &dhats_,
        const Eigen::MatrixXd &V);
    virtual ~SmoothEdgeVertex3Collision() 
    {
    }

    int ndofs() const override
    {
        return num_vertices() * 3;
    }

    int num_vertices() const override
    {
        return n_neighbors + 5;
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
    
    std::string name() const override { return "edge-vert"; }

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Eigen::VectorXd& positions, ParameterType params) const;

    bool compute_types(const Eigen::VectorXd& positions, ParameterType params); // return true if the potential is nonzero, return false if the potential is zero and can be skipped

    int n_neighbors;

    PointEdgeDistanceType dtype;
    ORIENTATION_TYPES otypes;
};

}