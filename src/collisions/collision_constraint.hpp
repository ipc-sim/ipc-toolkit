#pragma once

#include <array>

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct CollisionConstraint {
public:
    virtual ~CollisionConstraint() { }

    virtual int num_vertices() const = 0;
    virtual std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const = 0;

    virtual double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const;

    virtual VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const;

    virtual MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const bool project_hessian_to_psd) const;

    double minimum_distance = 0;
    double weight = 1;
#ifdef IPC_TOOLKIT_COMPUTE_SHAPE_DERIVATIVE
    Eigen::SparseVector<double> weight_gradient;
#endif
};

} // namespace ipc
