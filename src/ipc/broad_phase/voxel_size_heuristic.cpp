#include "voxel_size_heuristic.hpp"

#include <ipc/utils/logger.hpp>

#include <igl/median.h>

namespace ipc {

double suggest_good_voxel_size(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, double inflation_radius)
{
    // double edge_len_std_deviation;
    // double edge_len = mean_edge_length(V, V, E, edge_len_std_deviation);
    double edge_len = median_edge_length(V, V, E);
    double voxel_size = 2 * edge_len + inflation_radius;
    // double voxel_size = edge_len + edge_len_std_deviation + inflation_radius;
    if (voxel_size <= 0) { // this case should not happen in real simulations
        voxel_size = std::numeric_limits<double>::max();
    }
    assert(std::isfinite(voxel_size));
    logger().trace(
        "suggesting voxel size of {} (avg_edge_len={})", voxel_size, edge_len);
    return voxel_size;
}

double suggest_good_voxel_size(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    double inflation_radius)
{
    // double edge_len_std_deviation;
    // double edge_len = mean_edge_length(V0, V1, E, edge_len_std_deviation);
    // double disp_len_std_deviation;
    // double disp_len = mean_displacement_length(V1 - V0,
    // disp_len_std_deviation);
    double edge_len = median_edge_length(V0, V1, E);
    double disp_len = median_displacement_length(V1 - V0);
    // double voxel_size = 2 * edge_len + inflation_radius;
    double voxel_size = 2 * std::max(edge_len, disp_len) + inflation_radius;
    // double voxel_size = std::max(
    //                         edge_len + edge_len_std_deviation,
    //                         disp_len + disp_len_std_deviation)
    //     + inflation_radius;
    if (voxel_size <= 0) { // this case should not happen in real simulations
        voxel_size = std::numeric_limits<double>::max();
    }
    assert(std::isfinite(voxel_size));
    logger().trace(
        "suggesting voxel size of {} (avg_edge_len={} avg_disp_len={})",
        voxel_size, edge_len, disp_len);
    return voxel_size;
}

/// @brief Compute the mean edge length of a mesh.
double mean_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    double& std_deviation)
{
    if (E.rows() == 0) {
        std_deviation = 0;
        return 0;
    }

    double sum = 0;
    for (int i = 0; i < E.rows(); i++) {
        sum += (V0.row(E(i, 0)) - V0.row(E(i, 1))).norm();
        sum += (V1.row(E(i, 0)) - V1.row(E(i, 1))).norm();
    }
    const double mean = sum / (2 * E.rows());

    for (int i = 0; i < E.rows(); i++) {
        std_deviation +=
            std::pow((V0.row(E(i, 0)) - V0.row(E(i, 1))).norm() - mean, 2);
        std_deviation +=
            std::pow((V1.row(E(i, 0)) - V1.row(E(i, 1))).norm() - mean, 2);
    }
    std_deviation = sqrt(std_deviation / (2 * E.rows()));

    return mean;
}

/// @brief Compute the mean displacement length.
double mean_displacement_length(const Eigen::MatrixXd& U, double& std_deviation)
{
    const double mean = U.rowwise().norm().mean();
    std_deviation =
        sqrt((U.rowwise().norm().array() - mean).pow(2).sum() / U.rows());
    return mean;
}

/// @brief Compute the median edge length of a mesh.
double median_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E)
{
    if (E.rows() == 0) {
        return 0;
    }

    Eigen::VectorXd lengths(2 * E.rows());
    for (int i = 0; i < E.rows(); i++) {
        lengths[2 * i + 0] = (V0.row(E(i, 0)) - V0.row(E(i, 1))).norm();
        lengths[2 * i + 1] = (V1.row(E(i, 0)) - V1.row(E(i, 1))).norm();
    }

    double median = -1;
    const bool success = igl::median(lengths, median);
    assert(success);
    return median;
}

/// @brief Compute the median displacement length.
double median_displacement_length(const Eigen::MatrixXd& U)
{
    double median = -1;
    const bool success = igl::median(U.rowwise().norm(), median);
    assert(success);
    return median;
}

} // namespace ipc
