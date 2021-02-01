// Fucntions for computing the initial and updated barrier stiffnesses.

#include <ipc/barrier/adaptive_stiffness.hpp>

#include <ipc/barrier/barrier.hpp>
#include <ipc/ipc.hpp>

namespace ipc {

double world_bbox_diagonal(const Eigen::MatrixXd& V)
{
    Eigen::VectorX3d min = V.colwise().minCoeff();
    Eigen::VectorX3d max = V.colwise().maxCoeff();
    return (max - min).norm();
}

double initial_barrier_stiffness(
    double bbox_diagonal,
    double dhat,
    double average_mass,
    const Eigen::VectorXd& grad_energy,
    const Eigen::VectorXd& grad_barrier,
    double& max_barrier_stiffness,
    double min_barrier_stiffness_scale,
    double dmin)
{
    assert(average_mass > 0 && min_barrier_stiffness_scale > 0);
    assert(bbox_diagonal > 0);

    double dhat_squared = dhat * dhat;
    double dmin_squared = dmin * dmin;

    // Find a good initial value for κ
    double d0 = 1e-8 * bbox_diagonal + dmin;
    d0 *= d0;
    if (d0 - dmin_squared >= 2 * dmin * dhat + dhat_squared) {
        d0 = dmin * dhat + 0.5 * dhat_squared; // TODO: this is untested
    }
    double min_barrier_stiffness =
        barrier_hessian(d0 - dmin_squared, 2 * dmin * dhat + dhat_squared) * 4
        * d0;
    min_barrier_stiffness =
        min_barrier_stiffness_scale * average_mass / min_barrier_stiffness;
    assert(std::isfinite(min_barrier_stiffness));

    max_barrier_stiffness = 100 * min_barrier_stiffness;

    double kappa = 1.0;
    if (grad_barrier.squaredNorm() > 0) {
        // If this value is negative it will be clamped to κ_min anyways
        kappa = -grad_barrier.dot(grad_energy) / grad_barrier.squaredNorm();
        assert(std::isfinite(kappa));
    }

    return std::min(
        max_barrier_stiffness, std::max(min_barrier_stiffness, kappa));
}

// Use the mesh to compute the barrier gradient and world bounding box.
double initial_barrier_stiffness(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    double average_mass,
    const Eigen::VectorXd& grad_energy,
    double& max_barrier_stiffness,
    double min_barrier_stiffness_scale,
    double dmin)
{
    double diag = world_bbox_diagonal(V);

    Constraints constraint_set;
    construct_constraint_set(V_rest, V, E, F, dhat, constraint_set, dmin);
    Eigen::VectorXd grad_barrier =
        compute_barrier_potential_gradient(V, E, F, constraint_set, dhat);

    return initial_barrier_stiffness(
        diag, dhat, average_mass, grad_energy, grad_barrier,
        max_barrier_stiffness, min_barrier_stiffness_scale, dmin);
}

// Adaptive κ
void update_barrier_stiffness(
    double prev_min_distance,
    double min_distance,
    double max_barrier_stiffness,
    double& barrier_stiffness,
    double dhat_epsilon_scale,
    double bbox_diagonal,
    double dmin)
{
    // Is the barrier having a difficulty pushing the bodies apart?
    double dhat_epsilon = dhat_epsilon_scale * (bbox_diagonal + dmin);
    dhat_epsilon *= dhat_epsilon;
    if (prev_min_distance < dhat_epsilon && min_distance < dhat_epsilon
        && min_distance < prev_min_distance) {
        // Then increase the barrier stiffness.
        barrier_stiffness =
            std::min(max_barrier_stiffness, 2 * barrier_stiffness);
    }
}

// Adaptive κ
void update_barrier_stiffness(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    double prev_min_distance,
    double& min_distance,
    double max_barrier_stiffness,
    double& barrier_stiffness,
    double dhat_epsilon_scale,
    double dmin)
{
    Constraints constraint_set;
    construct_constraint_set(/*V_rest=*/V, V, E, F, dhat, constraint_set, dmin);
    // Use a temporay variable in case &prev_min_distance == &min_distance
    double current_min_distance =
        compute_minimum_distance(V, E, F, constraint_set);

    return update_barrier_stiffness(
        prev_min_distance, current_min_distance, max_barrier_stiffness,
        barrier_stiffness, dhat_epsilon_scale, world_bbox_diagonal(V), dmin);

    min_distance = current_min_distance;
}

} // namespace ipc
