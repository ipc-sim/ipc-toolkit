// Fucntions for computing the initial and updated barrier stiffnesses.

#include <ipc/barrier/adaptive_stiffness.hpp>

#include <ipc/barrier/barrier.hpp>
#include <ipc/ipc.hpp>

namespace ipc {

double world_bbox_diagonal(const Eigen::MatrixXd& V)
{
    Eigen::VectorXd min = V.colwise().minCoeff();
    Eigen::VectorXd max = V.colwise().maxCoeff();
    return (max - min).norm();
}

double intial_barrier_stiffness(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    double average_mass,
    const Eigen::VectorXd& grad_energy,
    double& max_barrier_stiffness,
    double min_barrier_stiffness_scale)
{
    double dhat_squared = dhat * dhat;

    double diag = world_bbox_diagonal(V);

    // Find a good initial value for κ
    double d0 = 1e-8 * diag;
    d0 *= d0;
    if (d0 >= dhat_squared) {
        d0 = 0.5 * dhat_squared; // TODO: this is untested
    }
    double min_barrier_stiffness = barrier_hessian(d0, dhat_squared) * 4 * d0;
    min_barrier_stiffness =
        min_barrier_stiffness_scale * average_mass / min_barrier_stiffness;
    assert(std::isfinite(min_barrier_stiffness));

    max_barrier_stiffness = 100 * min_barrier_stiffness;

    Constraints constraint_set;
    construct_constraint_set(V_rest, V, E, F, dhat, constraint_set);
    int num_active_barriers = constraint_set.num_constraints();
    Eigen::VectorXd grad_barrier =
        compute_barrier_potential_gradient(V, E, F, constraint_set, dhat);

    double kappa = 1.0;
    if (num_active_barriers > 0 && grad_barrier.squaredNorm() > 0) {
        // If this value is negative it will be clamped to κ_min anyways
        kappa = -grad_barrier.dot(grad_energy) / grad_barrier.squaredNorm();
        assert(std::isfinite(kappa));
    }

    return std::min(
        max_barrier_stiffness, std::max(min_barrier_stiffness, kappa));
}

void update_barrier_stiffness(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    double prev_min_distance,
    double& min_distance,
    double max_barrier_stiffness,
    double& barrier_stiffness,
    double dhat_epsilon_scale)
{
    double dhat_squared = dhat * dhat;

    double diag = world_bbox_diagonal(V);

    // Adaptive κ
    Constraints constraint_set;
    construct_constraint_set(
        /*V_rest=*/V, V, E, F, dhat_squared, constraint_set);
    min_distance = compute_minimum_distance(V, E, F, constraint_set);

    // Is the barrier having a difficulty pushing the bodies apart?
    double dhat_epsilon = dhat_epsilon_scale * diag;
    if (prev_min_distance < dhat_epsilon && min_distance < dhat_epsilon
        && min_distance < prev_min_distance) {
        // Then increase the barrier stiffness.
        barrier_stiffness =
            std::min(max_barrier_stiffness, 2 * barrier_stiffness);
    }
}

} // namespace ipc
