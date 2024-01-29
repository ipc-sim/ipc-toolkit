#pragma once

#include "smooth_collision.hpp"
#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/math.hpp>

namespace ipc {

class SmoothFaceVertexCollision : public SmoothCollision<max_vert_3d> {
public:
    constexpr static int max_size = SmoothCollision<max_vert_3d>::max_size;
    SmoothFaceVertexCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const double &dhat,
        const Eigen::MatrixXd &V);
    virtual ~SmoothFaceVertexCollision() 
    {
    }

    int ndofs() const override
    {
        return num_vertices() * 3;
    }

    int num_vertices() const override
    {
        return n_neighbors + 4;
    }

    double compute_distance(const Vector<double, -1, max_size>& positions) const override;

    double operator()(const Vector<double, -1, max_size>& positions, 
        const ParameterType &params) const override;

    Vector<double, -1, max_size> gradient(
        const Vector<double, -1, max_size>& positions, 
        const ParameterType &params) const override;

    MatrixMax<double, max_size, max_size> hessian(
        const Vector<double, -1, max_size>& positions, 
        const ParameterType &params) const override;
    
    std::string name() const override { return "face-vert"; }

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Eigen::VectorXd& positions, ParameterType params) const;

    bool compute_types(const Eigen::VectorXd& positions, ParameterType params); // return true if the potential is nonzero, return false if the potential is zero and can be skipped

    int n_neighbors;

    PointTriangleDistanceType dtype;
    ORIENTATION_TYPES otypes;
};

}