#include "friction_collision.hpp"

#include <ipc/friction/normal_force_magnitude.hpp>

#include <ipc/config.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

void FrictionCollision::init(
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Barrier& barrier,
    const double dhat,
    const double barrier_stiffness,
    const double dmin)
{
    // do this to initialize dim()
    const int dim = positions.cols();
    tangent_basis.resize(dim, dim - 1);

    const VectorMax12d pos = dof(positions, edges, faces);
    closest_point = compute_closest_point(pos);
    tangent_basis = compute_tangent_basis(pos);
    normal_force_magnitude = compute_normal_force_magnitude(
        pos, barrier, dhat, barrier_stiffness, dmin);
}

double FrictionCollision::compute_normal_force_magnitude(
    const VectorMax12d& positions,
    const Barrier& barrier,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude(
        barrier, compute_distance(positions), dhat, barrier_stiffness, dmin);
}

VectorMax12d FrictionCollision::compute_normal_force_magnitude_gradient(
    const VectorMax12d& positions,
    const Barrier& barrier,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude_gradient(
        barrier, compute_distance(positions),
        compute_distance_gradient(positions), dhat, barrier_stiffness, dmin);
}

} // namespace ipc
