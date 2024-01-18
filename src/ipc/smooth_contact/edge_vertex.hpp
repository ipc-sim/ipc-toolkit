#pragma once

#include "smooth_collision.hpp"
#include <ipc/utils/math.hpp>

namespace ipc {

class SmoothEdgeVertexCollision : public SmoothCollision<6> {
    constexpr static int dim = 2;
public:
    using Super = SmoothCollision<6>;
    
    SmoothEdgeVertexCollision(
        long primitive0_,
        long primitive1_,
        const CollisionMesh &mesh,
        const ParameterType &param,
        const std::array<double, 2> &dhats_,
        const Eigen::MatrixXd &V);

    int ndofs() const override
    {
        return 10;
    }

    int num_vertices() const override
    {
        return 5;
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
    scalar evaluate_quadrature(const Vector10d& positions, ParameterType params) const;

    bool compute_types(const Vector10d& positions, const ParameterType &params); // return true if the potential is nonzero, return false if the potential is zero and can be skipped
};

}