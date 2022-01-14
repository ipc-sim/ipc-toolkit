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


} // namespace ipc
