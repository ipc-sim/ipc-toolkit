#pragma once

#include <ipc/config.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

    /// @brief Generate quadrature points and weights for unit line [0, 1]
    /// @param order Order of quadrature rule, only linear for now
    /// @param N Number of quadrature points
    /// @param pts Barycentric positions of quadrature points Nx1
    /// @param weights Quadrature weights Nx1
    void line_quadrature(
        const int order, 
        const int N, 
        Eigen::VectorXd &pts, 
        Eigen::VectorXd &weights);
}