#pragma once

#include "smooth_collision.hpp"
#include <ipc/utils/math.hpp>
#include <ipc/distance/distance_type.hpp>
#include <iostream>

namespace ipc {

class SmoothFaceFaceCollision : public SmoothCollision<8> {
    SmoothFaceFaceCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh);
public:
    SmoothFaceFaceCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const Eigen::MatrixXd &V);
    virtual ~SmoothFaceFaceCollision() { }

    int num_vertices() const override
    {
        return 6;
    }

    std::array<long, 8> vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const override;

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
    scalar evaluate_quadrature(const Vector<double, 18>& positions, const ParameterType &params) const;

    bool compute_types(const Vector<double, 18>& positions, const ParameterType &params); // return true if the potential is nonzero, return false if the potential is zero and can be skipped

    std::array<PointTriangleDistanceType, 6> dtypes;
    std::array<HEAVISIDE_TYPE, 6> normal_types;
};

}