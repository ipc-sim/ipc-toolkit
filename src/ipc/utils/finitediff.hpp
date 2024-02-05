#pragma once

#include <ipc/config.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {
enum class FD_RULE { CENTRAL, LEFT, RIGHT };

void finite_gradient(
    const Eigen::VectorXd& x,
    const std::function<double(const Eigen::VectorXd&)>& f,
    Eigen::VectorXd& grad,
    FD_RULE rule = FD_RULE::CENTRAL,
    const double eps = 1e-7)
{
    grad.setZero(x.size());
    switch (rule) {
    case FD_RULE::CENTRAL:
        for (int i = 0; i < x.size(); i++)
            for (int d : { -1, 1 }) {
                auto y = x;
                y(i) += d * eps;
                grad(i) += d * f(y) / (2 * eps);
            }
        break;
    case FD_RULE::LEFT:
        for (int i = 0; i < x.size(); i++) {
            auto y = x;
            grad(i) += f(y) / eps;
            y(i) -= eps;
            grad(i) -= f(y) / eps;
        }
        break;
    case FD_RULE::RIGHT:
        for (int i = 0; i < x.size(); i++) {
            auto y = x;
            grad(i) -= f(y) / eps;
            y(i) += eps;
            grad(i) += f(y) / eps;
        }
        break;
    default:
        assert(false);
    }
}
} // namespace ipc