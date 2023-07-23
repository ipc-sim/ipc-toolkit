#include "continuous_collision_candidate.hpp"

namespace ipc {

bool ContinuousCollisionCandidate::ccd(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double& toi,
    const double min_distance,
    const double tmax,
    const double tolerance,
    const long max_iterations,
    const double conservative_rescaling) const
{
    return ccd(
        dof(vertices_t0, edges, faces), dof(vertices_t1, edges, faces), toi,
        min_distance, tmax, tolerance, max_iterations, conservative_rescaling);
}

} // namespace ipc
