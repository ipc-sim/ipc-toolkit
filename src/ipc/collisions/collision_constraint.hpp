#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class CollisionConstraint : virtual public CollisionStencil {
public:
    virtual ~CollisionConstraint() { }

    virtual double compute_potential(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat) const;

    virtual VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat) const;

    virtual MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const bool project_hessian_to_psd) const;

    double minimum_distance = 0;
    double weight = 1;
    Eigen::SparseVector<double> weight_gradient;
};

} // namespace ipc
