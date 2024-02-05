#pragma once

#include <ipc/config.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Generate composite trapezoidal quadrature points and weights for unit line [0, 1]
/// @param N Number of quadrature points
/// @param pts Barycentric positions of quadrature points Nx1
/// @param weights Quadrature weights Nx1
void line_quadrature(
    const int N, Eigen::VectorXd& pts, Eigen::VectorXd& weights);

/// @brief Generate composite trapezoidal quadrature points and weights for unit triangle x+y<=1, x>=0, y>=0
/// @param N Control number of samples, actual number of samples is N*(N+1)/2
/// @param pts UV of quadrature points
/// @param weights Quadrature weights
/// @return
void triangle_quadrature(
    const int N, Eigen::Matrix<double, -1, 2>& pts, Eigen::VectorXd& weights);
} // namespace ipc