#pragma once

#include <ipc/dynamics/time_integration/time_integrator.hpp>

namespace ipc::dynamics {

/// @brief Backward Differentiation Formulas (BDF) of order 1 ≤ n ≤ 6.
/// \f[
///     x^{t+1} = \left(\sum_{i=0}^{n-1} \alpha_i x^{t-i}\right) + \beta\Delta
///     t\, v^{t+1}, \quad v^{t+1} = \left(\sum_{i=0}^{n-1} \alpha_i
///     v^{t-i}\right) + \beta\Delta t\, a^{t+1}
/// \f]
/// BDF-1 is the standard implicit (backward) Euler method.
/// @see https://en.wikipedia.org/wiki/Backward_differentiation_formula
class BDF : public ImplicitTimeIntegrator {
public:
    BDF() = default;

    /// @brief Construct a BDF integrator of the given target order.
    explicit BDF(const int target_order) { set_target_order(target_order); }

    /// @brief Predicted positions
    /// \f$\hat{x} = \sum_i \alpha_i x^{t-i} + \beta\Delta t \sum_i \alpha_i
    /// v^{t-i}\f$.
    Eigen::VectorXd predicted_positions() const override;

    /// @brief \f$x = \sum_i \alpha_i x^{t-i} + \beta\Delta t\, v\f$.
    Eigen::VectorXd
    compute_position(Eigen::ConstRef<Eigen::VectorXd> v) const override;

    /// @brief \f$v = (x - \sum_i \alpha_i x^{t-i}) / (\beta\Delta t)\f$.
    Eigen::VectorXd
    compute_velocity(Eigen::ConstRef<Eigen::VectorXd> x) const override;

    /// @brief \f$a = (v - \sum_i \alpha_i v^{t-i}) / (\beta\Delta t)\f$.
    Eigen::VectorXd
    compute_acceleration(Eigen::ConstRef<Eigen::VectorXd> v) const override;

    /// @brief \f$(\beta\Delta t)^2\f$.
    double acceleration_scaling() const override;

    /// @brief \f$\beta\Delta t\f$.
    double velocity_scaling() const override { return beta_dt(); }

    /// @brief \f$\partial v/\partial x = 1/(\beta\Delta t)\f$,
    /// \f$\partial v/\partial x^{t-i} = -\alpha_i/(\beta\Delta t)\f$.
    double dvdx(const unsigned prev_ti = 0) const override;

    /// @brief The requested (maximum) BDF order.
    int target_order() const { return m_target_order; }

    /// @brief Set the requested BDF order (1 ≤ n ≤ 6).
    /// @note Must be set before init(); it determines the history buffer size.
    void set_target_order(const int target_order);

    /// @brief The current order (min of the target and the available history).
    int order() const
    {
        return std::min<int>(m_target_order, available_steps());
    }

    /// @brief \f$\sum_{i=0}^{n-1} \alpha_i x^{t-i}\f$.
    Eigen::VectorXd weighted_sum_x_prevs() const;

    /// @brief \f$\sum_{i=0}^{n-1} \alpha_i v^{t-i}\f$.
    Eigen::VectorXd weighted_sum_v_prevs() const;

    /// @brief \f$\beta\Delta t\f$.
    double beta_dt() const;

    /// @brief The \f$\alpha\f$ coefficients for order-`i+1` BDF.
    static const std::vector<double>& alphas(const int i);

    /// @brief The \f$\beta\f$ coefficient for order-`i+1` BDF.
    static double betas(const int i);

protected:
    unsigned max_steps() const override { return unsigned(m_target_order); }

private:
    /// @brief The requested (maximum) BDF order.
    int m_target_order = 1;
};

} // namespace ipc::dynamics
