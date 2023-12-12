#include "quadrature.hpp"

namespace ipc {
    void line_quadrature(
        const int order, 
        const int N, 
        Eigen::VectorXd &pts, 
        Eigen::VectorXd &weights)
    {
        pts = Eigen::VectorXd::LinSpaced(N, 0, 1);
        weights = Eigen::VectorXd::Ones(N) / N;
    }
}