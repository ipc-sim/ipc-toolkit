#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class CollisionConstraint : virtual public CollisionStencil {
public:
    CollisionConstraint() = default;

    CollisionConstraint(
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient);

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

    /// Compute the derivative of the potential gradient wrt the shape.
    virtual void compute_shape_derivative(
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        std::vector<Eigen::Triplet<double>>& triplets) const;

    double dmin = 0;
    double weight = 1;
    Eigen::SparseVector<double> weight_gradient;

protected:
    /// Compute (∇ₓw)(∇ᵤb)ᵀ
    virtual void compute_shape_derivative_first_term(
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        std::vector<Eigen::Triplet<double>>& triplets) const;

    /// Compute w ∇ₓ∇ᵤb
    virtual MatrixMax12d compute_shape_derivative_second_term(
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat) const;
};

} // namespace ipc
