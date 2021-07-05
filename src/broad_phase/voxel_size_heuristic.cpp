#include <ipc/broad_phase/voxel_size_heuristic.hpp>

namespace ipc {

double suggest_good_voxel_size(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, double inflation_radius)
{
    double edge_len = average_edge_length(V, V, E);
    return 2 * edge_len + inflation_radius;
}

double suggest_good_voxel_size(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    double inflation_radius)
{
    double edge_len = average_edge_length(V0, V1, E);
    // double disp_len = average_displacement_length(V1 - V0);
    // return 2 * std::max(edge_len, disp_len) + inflation_radius;
    return 2 * edge_len + inflation_radius;
}

/// @brief Compute the average edge length of a mesh.
double average_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E)
{
    double avg = 0;
    for (unsigned i = 0; i < E.rows(); i++) {
        avg += (V0.row(E(i, 0)) - V0.row(E(i, 1))).norm();
        avg += (V1.row(E(i, 0)) - V1.row(E(i, 1))).norm();
    }
    return avg / (2 * E.rows());
}

/// @brief Compute the average displacement length.
double average_displacement_length(const Eigen::MatrixXd& U)
{
    return U.rowwise().norm().sum() / U.rows();
}

} // namespace ipc
