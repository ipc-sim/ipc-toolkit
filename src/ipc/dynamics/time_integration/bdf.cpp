#include "bdf.hpp"

#include <ipc/utils/logger.hpp>

#include <array>
#include <cassert>

namespace ipc::dynamics {

void BDF::set_target_order(const int target_order)
{
    if (target_order < 1 || target_order > 6) {
        log_and_throw_error("BDF order must be 1 ≤ n ≤ 6");
    }
    m_target_order = target_order;
}

const std::vector<double>& BDF::alphas(const int i)
{
    static const std::array<std::vector<double>, 6> alphas = { {
        { 1 },
        { 4.0 / 3.0, -1.0 / 3.0 },
        { 18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0 },
        { 48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0 },
        { 300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0,
          12.0 / 137.0 },
        { 360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0,
          72.0 / 147.0, -10.0 / 147.0 },
    } };
    assert(i >= 0 && i < alphas.size());
    return alphas[i];
}

double BDF::betas(const int i)
{
    static const std::array<double, 6> betas = { {
        1.0,
        2.0 / 3.0,
        6.0 / 11.0,
        12.0 / 25.0,
        60.0 / 137.0,
        60.0 / 147.0,
    } };
    assert(i >= 0 && i < betas.size());
    return betas[i];
}

Eigen::VectorXd BDF::weighted_sum_x_prevs() const
{
    const std::vector<double>& alpha = alphas(order() - 1);

    Eigen::VectorXd sum = alpha[0] * x_prev(0);
    for (int i = 1; i < order(); i++) {
        sum += alpha[i] * x_prev(i);
    }
    return sum;
}

Eigen::VectorXd BDF::weighted_sum_v_prevs() const
{
    const std::vector<double>& alpha = alphas(order() - 1);

    Eigen::VectorXd sum = alpha[0] * v_prev(0);
    for (int i = 1; i < order(); i++) {
        sum += alpha[i] * v_prev(i);
    }
    return sum;
}

Eigen::VectorXd BDF::predicted_positions() const
{
    return weighted_sum_x_prevs()
        + betas(order() - 1) * dt() * weighted_sum_v_prevs();
}

Eigen::VectorXd
BDF::compute_position(Eigen::ConstRef<Eigen::VectorXd> v) const
{
    return weighted_sum_x_prevs() + beta_dt() * v;
}

Eigen::VectorXd
BDF::compute_velocity(Eigen::ConstRef<Eigen::VectorXd> x) const
{
    return (x - weighted_sum_x_prevs()) / beta_dt();
}

Eigen::VectorXd
BDF::compute_acceleration(Eigen::ConstRef<Eigen::VectorXd> v) const
{
    return (v - weighted_sum_v_prevs()) / beta_dt();
}

double BDF::acceleration_scaling() const
{
    const double beta_dt = this->beta_dt();
    return beta_dt * beta_dt;
}

double BDF::dvdx(const unsigned prev_ti) const
{
    if (prev_ti == 0) {
        return 1 / beta_dt();
    }
    if (int(prev_ti) >= order()) {
        return 0;
    }
    return -alphas(order() - 1)[prev_ti] / beta_dt();
}

double BDF::beta_dt() const { return betas(order() - 1) * dt(); }

} // namespace ipc::dynamics
