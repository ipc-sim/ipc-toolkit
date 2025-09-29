#include "quadrature.hpp"

namespace ipc {
void line_quadrature(
    const int N, Eigen::VectorXd& pts, Eigen::VectorXd& weights)
{
    pts = Eigen::VectorXd::LinSpaced(N, 0, 1);
    weights = Eigen::VectorXd::Ones(N) / (N - 1);
    weights(0) /= 2;
    weights(N - 1) /= 2;
}

void triangle_quadrature(
    const int N, Eigen::Matrix<double, -1, 2>& pts, Eigen::VectorXd& weights)
{
    const double h = 1.0 / (N - 1);
    const int N_pts = N * (N + 1) / 2;
    pts.setZero(N_pts, 2);
    weights.setZero(N_pts);
    for (int i = 0, id = 0; i < N; i++) {
        for (int j = 0; i + j < N; j++, id++) {
            pts.row(id) << i, j;
            int flag = static_cast<int>(i == 0) + static_cast<int>(j == 0)
                + static_cast<int>(i + j == N - 1);
            weights(id) = (flag == 2) ? 1.0 / 6.0
                : (flag == 1)         ? 1.0 / 2.0
                                      : 1.0;
        }
    }
    pts *= h;
    weights *= h * h;
}
} // namespace ipc