#include <ipc/friction/normal_force_magnitude.hpp>

#include <ipc/barrier/barrier.hpp>

namespace ipc {

double compute_normal_force_magnitude(
    double distance_squared, double dhat, double barrier_stiffness, double dmin)
{
    double grad_b = barrier_gradient(
        distance_squared - dmin * dmin, 2 * dmin * dhat + dhat * dhat);
    grad_b *= barrier_stiffness;
    return -grad_b * 2 * sqrt(distance_squared); // / (h * h) eliminated here
}

VectorMax12d compute_normal_force_magnitude_gradient(
    double distance_squared,
    const Eigen::VectorXd& distance_squared_gradient,
    double dhat,
    double barrier_stiffness,
    double dmin)
{
    double arg_d = distance_squared - dmin * dmin;
    double arg_dhat = 2 * dmin * dhat + dhat * dhat;
    double distance = sqrt(distance_squared);

    return -barrier_stiffness
        * (barrier_hessian(arg_d, arg_dhat) * 2 * distance
           + barrier_gradient(arg_d, arg_dhat) / distance)
        * distance_squared_gradient;
}

} // namespace ipc
