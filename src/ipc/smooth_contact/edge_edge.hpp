#pragma once

#include "smooth_collision.hpp"
#include <ipc/utils/math.hpp>

namespace ipc {

class SmoothEdgeEdgeCollision : public SmoothCollision<6> {
    constexpr static int dim = 2;
public:
    using Super = SmoothCollision<6>;
    
    SmoothEdgeEdgeCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const std::array<double, 2> &dhats_,
        const Eigen::MatrixXd &V);

    int num_vertices() const override
    {
        return 4;
    }

    std::array<long, 6> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return { { edges(primitive0, 0), edges(primitive0, 1),
                   edges(primitive1, 0), edges(primitive1, 1), -1, -1 } };
    }

    double compute_distance(const Vector<double, -1, 18>& positions) const override;

    double operator()(const VectorMax18d& positions, 
        const ParameterType &params) const override;

    VectorMax18d gradient(
        const VectorMax18d& positions, 
        const ParameterType &params) const override;

    MatrixMax18d hessian(
        const VectorMax18d& positions, 
        const ParameterType &params,
        const bool project_hessian_to_psd = false) const override;

private:
    template <typename scalar> 
    scalar evaluate_quadrature(const Vector8d& positions, ParameterType params) const;

    bool compute_types(const Vector8d& positions, const ParameterType &params); // return true if the potential is nonzero, return false if the potential is zero and can be skipped

    std::array<HEAVISIDE_TYPE, 4> normal_types;
};

}