#pragma once

#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc::dynamics {

/// @brief Implicit time integrator of the second-order ODE M ẍ = f(x).
///
/// The step is posed as the minimization of an incremental potential
/// \f[
///     E(x) = \tfrac12 (x - \hat{x})^\top M (x - \hat{x}) + s\, V(x),
/// \f]
/// where \f$\hat{x} =\f$ predicted_positions() and \f$s =\f$
/// acceleration_scaling(); the minimizer satisfies
/// \f$M(x-\hat{x}) = s\, f(x)\f$ with \f$f = -\nabla V\f$.
///
/// A ring buffer of the previous positions, velocities, and accelerations is
/// stored so multi-step schemes (e.g., BDF-n) can be implemented. The state is
/// affine-shaped ([p; vec(Q) column-major] per body); the pose() helpers decode
/// it into rigid::AffinePose.
class ImplicitTimeIntegrator {
public:
    virtual ~ImplicitTimeIntegrator() = default;

    /// @brief Initialize the integrator with the previous x, v, and a.
    /// @param x_prev Previous positions.
    /// @param v_prev Previous velocities.
    /// @param a_prev Previous accelerations.
    /// @param num_bodies Number of bodies (used to decode the affine-shaped
    ///     state into poses).
    /// @note For multi-step schemes, the scheme's parameters (e.g., the BDF
    ///     order) must be set before calling init() since they determine the
    ///     history buffer size.
    void init(
        Eigen::ConstRef<Eigen::VectorXd> x_prev,
        Eigen::ConstRef<Eigen::VectorXd> v_prev,
        Eigen::ConstRef<Eigen::VectorXd> a_prev,
        const size_t num_bodies);

    /// @brief Advance the stored history with a newly-accepted state x.
    /// @param x The new positions.
    void update(Eigen::ConstRef<Eigen::VectorXd> x);

    // -----------------------------------------------------------------------
    // Scheme-specific integration formulas
    // -----------------------------------------------------------------------

    /// @brief Predicted positions \f$\hat{x}\f$ used in the inertia term.
    virtual Eigen::VectorXd predicted_positions() const = 0;

    /// @brief Compute the current position given the current velocity.
    virtual Eigen::VectorXd
    compute_position(Eigen::ConstRef<Eigen::VectorXd> v) const = 0;

    /// @brief Compute the current velocity given the current position.
    virtual Eigen::VectorXd
    compute_velocity(Eigen::ConstRef<Eigen::VectorXd> x) const = 0;

    /// @brief Compute the current acceleration given the current velocity.
    virtual Eigen::VectorXd
    compute_acceleration(Eigen::ConstRef<Eigen::VectorXd> v) const = 0;

    /// @brief Scaling of the potential forces in the incremental potential
    /// (\f$(\beta\Delta t)^2\f$ for BDF; \f$\Delta t^2\f$ for implicit Euler).
    virtual double acceleration_scaling() const = 0;

    /// @brief First-order (velocity) force scaling (\f$\beta\Delta t\f$ for BDF).
    virtual double velocity_scaling() const = 0;

    /// @brief Derivative of the velocity with respect to the i-th previous
    /// positions (0 → current).
    virtual double dvdx(const unsigned prev_ti = 0) const = 0;

    // -----------------------------------------------------------------------
    // Time step
    // -----------------------------------------------------------------------

    /// @brief The time step size.
    double dt() const { return m_dt; }

    /// @brief Set the time step size (resets the accumulated history).
    void set_dt(const double dt);

    // -----------------------------------------------------------------------
    // History access (0 → most recent, 1 → previous, ...)
    // -----------------------------------------------------------------------

    /// @brief The i-th previous positions.
    Eigen::Map<const Eigen::VectorXd> x_prev(const int i = 0) const;
    /// @brief The i-th previous velocities.
    Eigen::Map<const Eigen::VectorXd> v_prev(const int i = 0) const;
    /// @brief The i-th previous accelerations.
    Eigen::Map<const Eigen::VectorXd> a_prev(const int i = 0) const;

    // -----------------------------------------------------------------------
    // Pose decoding of the affine-shaped state
    // -----------------------------------------------------------------------

    /// @brief Decode body i's affine pose from a state vector x.
    rigid::AffinePose
    pose(Eigen::ConstRef<Eigen::VectorXd> x, const size_t i) const;

    /// @brief Decode body i's affine pose from the most recent state.
    rigid::AffinePose pose(const size_t i) const { return pose(x_prev(0), i); }

    /// @brief Decode all bodies' poses from predicted_positions().
    std::vector<rigid::AffinePose> predicted_pose() const;

    /// @brief Number of bodies.
    size_t num_bodies() const { return m_num_bodies; }
    /// @brief Position DOFs per body (2 or 3).
    size_t pos_ndof() const { return m_pos_ndof; }
    /// @brief Rotation DOFs per body (1 or 9).
    size_t rot_ndof() const { return m_rot_ndof; }

protected:
    /// @brief The current number of stored history steps.
    int available_steps() const { return m_available_steps; }

    /// @brief The maximum number of history steps the scheme requires.
    virtual unsigned max_steps() const = 0;

private:
    /// @brief Ring buffers of the previous positions/velocities/accelerations
    /// (one column per stored step).
    Eigen::MatrixXd m_x_prevs;
    Eigen::MatrixXd m_v_prevs;
    Eigen::MatrixXd m_a_prevs;

    /// @brief Time step size.
    double m_dt = 1.0 / 60.0;

    /// @brief Column index of the most recent stored step.
    unsigned m_current_ptr = 0;

    /// @brief Number of stored history steps (≤ max_steps()).
    unsigned m_available_steps = 1;

    /// @brief Number of bodies and the per-body state layout.
    size_t m_num_bodies = 0;
    size_t m_pos_ndof = 3;
    size_t m_rot_ndof = 9;
};

} // namespace ipc::dynamics
