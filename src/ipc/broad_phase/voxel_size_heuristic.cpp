#include "voxel_size_heuristic.hpp"

#include <ipc/utils/logger.hpp>

#include <igl/median.h>

namespace ipc {

namespace {
    // Avoid unused variable warnings
    inline void check_success(bool success) { assert(success); }

    // Faster than std::pow(x, 2)
    inline double sqr(double x) { return x * x; }
} // namespace

double suggest_good_voxel_size(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    // double edge_len_std_deviation;
    // double edge_len =
    //     mean_edge_length(vertices, vertices, edges, edge_len_std_deviation);
    // double voxel_size = edge_len + edge_len_std_deviation + inflation_radius;

    double edge_len = median_edge_length(vertices, vertices, edges);
    double voxel_size = 2 * edge_len + inflation_radius;

    // double voxel_size =
    //     max_edge_length(vertices, vertices, edges) + inflation_radius;

    if (voxel_size <= 0) { // this case should not happen in real simulations
        voxel_size = std::numeric_limits<double>::max();
    }
    assert(std::isfinite(voxel_size));
    logger().trace(
        "suggesting voxel size of {:g} (avg_edge_len={:g})", voxel_size,
        edge_len);
    return voxel_size;
}

double suggest_good_voxel_size(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    // double edge_len_std_deviation;
    // double edge_len = mean_edge_length(
    //     vertices_t0, vertices_t1, edges, edge_len_std_deviation);
    // double disp_len_std_deviation;
    // double disp_len = mean_displacement_length(
    //     vertices_t1 - vertices_t0, disp_len_std_deviation);
    // double voxel_size = std::max(
    //                         edge_len + edge_len_std_deviation,
    //                         disp_len + disp_len_std_deviation)
    //     + inflation_radius;

    double edge_len = median_edge_length(vertices_t0, vertices_t1, edges);
    double disp_len = median_displacement_length(vertices_t1 - vertices_t0);
    double voxel_size = 2 * std::max(edge_len, disp_len) + inflation_radius;

    // double voxel_size = std::max(
    //                         max_edge_length(vertices_t0, vertices_t1, edges),
    //                         max_displacement_length(vertices_t1 -
    //                         vertices_t0))
    //     + inflation_radius;

    if (voxel_size <= 0) { // this case should not happen in real simulations
        voxel_size = std::numeric_limits<double>::max();
    }
    assert(std::isfinite(voxel_size));
    logger().trace(
        "suggesting voxel size of {:g} (avg_edge_len={:g} avg_disp_len={:g})",
        voxel_size, edge_len, disp_len);
    return voxel_size;
}

double mean_edge_length(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    double& std_deviation)
{
    if (edges.rows() == 0) {
        std_deviation = 0;
        return 0;
    }

    double sum = 0;
    for (int i = 0; i < edges.rows(); i++) {
        const int e0i = edges(i, 0), e1i = edges(i, 1);
        sum += (vertices_t0.row(e0i) - vertices_t0.row(e1i)).norm();
        sum += (vertices_t1.row(e0i) - vertices_t1.row(e1i)).norm();
    }
    const double mean = sum / (2 * edges.rows());

    std_deviation = 0;
    for (int i = 0; i < edges.rows(); i++) {
        const int e0i = edges(i, 0), e1i = edges(i, 1);
        std_deviation +=
            sqr((vertices_t0.row(e0i) - vertices_t0.row(e1i)).norm() - mean);
        std_deviation +=
            sqr((vertices_t1.row(e0i) - vertices_t1.row(e1i)).norm() - mean);
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
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges)
{
    if (edges.rows() == 0) {
        return 0;
    }

    Eigen::VectorXd lengths(2 * edges.rows());
    for (int i = 0; i < edges.rows(); i++) {
        const size_t e0i = edges(i, 0), e1i = edges(i, 1);
        lengths[2 * i + 0] =
            (vertices_t0.row(e0i) - vertices_t0.row(e1i)).norm();
        lengths[2 * i + 1] =
            (vertices_t1.row(e0i) - vertices_t1.row(e1i)).norm();
    }

    double median = -1;
    check_success(igl::median(lengths, median));
    return median;
}

double median_displacement_length(const Eigen::MatrixXd& displacements)
{
    double median = -1;
    check_success(igl::median(displacements.rowwise().norm(), median));
    return median;
}

double max_edge_length(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges)
{
    double max_edge = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < edges.rows(); i++) {
        const size_t e0i = edges(i, 0), e1i = edges(i, 1);
        max_edge = std::max({
            max_edge,
            (vertices_t0.row(e0i) - vertices_t0.row(e1i)).norm(),
            (vertices_t1.row(e0i) - vertices_t1.row(e1i)).norm(),
        });
    }
    return max_edge;
}

double max_displacement_length(const Eigen::MatrixXd& displacements)
{
    return displacements.rowwise().norm().maxCoeff();
}

} // namespace ipc
