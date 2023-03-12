#include "voxel_size_heuristic.hpp"

#include <ipc/utils/logger.hpp>

#include <igl/median.h>

namespace ipc {

double suggest_good_voxel_size(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& edges,
    double inflation_radius)
{
    // double edge_len_std_deviation;
    // double edge_len =
    //     mean_edge_length(V, V, edges,
    //     edge_len_std_deviation);
    // double voxel_size = edge_len + edge_len_std_deviation + inflation_radius;

    double edge_len = median_edge_length(V, V, edges);
    double voxel_size = 2 * edge_len + inflation_radius;

    // double voxel_size =
    //     max_edge_length(V, V, edges) + inflation_radius;

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
    const Eigen::MatrixXi& edges,
    double inflation_radius)
{
    // double edge_len_std_deviation;
    // double edge_len = mean_edge_length(
    //     V0, V1, edges, edge_len_std_deviation);
    // double disp_len_std_deviation;
    // double disp_len = mean_displacement_length(
    //     V1 - V0, disp_len_std_deviation);
    // double voxel_size = std::max(
    //                         edge_len + edge_len_std_deviation,
    //                         disp_len + disp_len_std_deviation)
    //     + inflation_radius;

    double edge_len = median_edge_length(V0, V1, edges);
    double disp_len = median_displacement_length(V1 - V0);
    double voxel_size = 2 * std::max(edge_len, disp_len) + inflation_radius;

    // double voxel_size =
    //     std::max(
    //         max_edge_length(V0, V1, edges),
    //         max_displacement_length(V1 - V0))
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
    const Eigen::MatrixXi& edges,
    double& std_deviation)
{
    if (edges.rows() == 0) {
        std_deviation = 0;
        return 0;
    }

    double sum = 0;
    for (int i = 0; i < edges.rows(); i++) {
        const size_t e0i = edges(i, 0), e1i = edges(i, 1);
        sum += (V0.row(e0i) - V0.row(e1i)).norm();
        sum += (V1.row(e0i) - V1.row(e1i)).norm();
    }
    const double mean = sum / (2 * edges.rows());

    for (int i = 0; i < edges.rows(); i++) {
        const size_t e0i = edges(i, 0), e1i = edges(i, 1);
        std_deviation += std::pow((V0.row(e0i) - V0.row(e1i)).norm() - mean, 2);
        std_deviation += std::pow((V1.row(e0i) - V1.row(e1i)).norm() - mean, 2);
    }
    std_deviation = sqrt(std_deviation / (2 * edges.rows()));

    return mean;
}

double mean_displacement_length(
    const Eigen::MatrixXd& displacements, double& std_deviation)
{
    const double mean = displacements.rowwise().norm().mean();
    std_deviation = sqrt(
        (displacements.rowwise().norm().array() - mean).pow(2).sum()
        / displacements.rows());
    return mean;
}

double median_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges)
{
    if (edges.rows() == 0) {
        return 0;
    }

    Eigen::VectorXd lengths(2 * edges.rows());
    for (int i = 0; i < edges.rows(); i++) {
        const size_t e0i = edges(i, 0), e1i = edges(i, 1);
        lengths[2 * i + 0] = (V0.row(e0i) - V0.row(e1i)).norm();
        lengths[2 * i + 1] = (V1.row(e0i) - V1.row(e1i)).norm();
    }

    double median = -1;
    const bool success = igl::median(lengths, median);
    assert(success);
    return median;
}

double median_displacement_length(const Eigen::MatrixXd& displacements)
{
    double median = -1;
    const bool success = igl::median(displacements.rowwise().norm(), median);
    assert(success);
    return median;
}

double max_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& edges)
{
    double max_edge = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < edges.rows(); i++) {
        const size_t e0i = edges(i, 0), e1i = edges(i, 1);
        max_edge = std::max({
            max_edge,
            (V0.row(e0i) - V0.row(e1i)).norm(),
            (V1.row(e0i) - V1.row(e1i)).norm(),
        });
    }
    return max_edge;
}

double max_displacement_length(const Eigen::MatrixXd& displacements)
{
    return displacements.rowwise().norm().maxCoeff();
}

} // namespace ipc
