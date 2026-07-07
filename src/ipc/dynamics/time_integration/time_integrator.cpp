#include "time_integrator.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cassert>

namespace ipc::dynamics {

namespace {
    /// @brief Modulus that is always non-negative (unlike %).
    inline int true_mod(const int a, const int b) { return (a % b + b) % b; }
} // namespace

void ImplicitTimeIntegrator::init(
    Eigen::ConstRef<Eigen::VectorXd> x_prev,
    Eigen::ConstRef<Eigen::VectorXd> v_prev,
    Eigen::ConstRef<Eigen::VectorXd> a_prev,
    const size_t num_bodies)
{
    assert(x_prev.size() == v_prev.size());
    assert(v_prev.size() == a_prev.size());

    m_num_bodies = num_bodies;

    m_pos_ndof = 3, m_rot_ndof = 9; // Default to the 3D case
    // TODO: This can be broken if there are a multiple of 4 bodies in 2D. We
    //       should probably pass in the number of DOFs instead of inferring it.
    if (x_prev.size() != num_bodies * (m_pos_ndof + m_rot_ndof)) {
        m_pos_ndof = 2, m_rot_ndof = 1; // 2D case
    }
    // If neither body layout matches exactly, x_prev is not affine-shaped
    // (e.g., a generic state used only for the raw multistep formulas below);
    // pose()/predicted_pose() are only meaningful when the layout does fit,
    // so we don't require it here.

    m_available_steps = 1;
    m_current_ptr = 0;

    m_x_prevs.resize(x_prev.size(), max_steps());
    m_v_prevs.resize(v_prev.size(), max_steps());
    m_a_prevs.resize(a_prev.size(), max_steps());

    m_x_prevs.col(m_current_ptr) = x_prev;
    m_v_prevs.col(m_current_ptr) = v_prev;
    m_a_prevs.col(m_current_ptr) = a_prev;
}

void ImplicitTimeIntegrator::update(Eigen::ConstRef<Eigen::VectorXd> x)
{
    assert(x.size() == m_x_prevs.rows());

    const Eigen::VectorXd v = compute_velocity(x);
    const Eigen::VectorXd a = compute_acceleration(v);

    m_current_ptr = (m_current_ptr + 1) % max_steps();

    m_x_prevs.col(m_current_ptr) = x;
    m_v_prevs.col(m_current_ptr) = v;
    m_a_prevs.col(m_current_ptr) = a;

    m_available_steps = std::min<unsigned>(m_available_steps + 1, max_steps());
}

void ImplicitTimeIntegrator::set_dt(const double dt)
{
    m_dt = dt;
    m_available_steps = 1; // Reset the history: a new dt invalidates it.
}

Eigen::Map<const Eigen::VectorXd>
ImplicitTimeIntegrator::x_prev(const int i) const
{
    assert(i < available_steps());
    const int ptr = true_mod(int(m_current_ptr) - i, max_steps());
    return Eigen::Map<const Eigen::VectorXd>(
        m_x_prevs.data() + ptr * m_x_prevs.rows(), m_x_prevs.rows());
}

Eigen::Map<const Eigen::VectorXd>
ImplicitTimeIntegrator::v_prev(const int i) const
{
    assert(i < available_steps());
    const int ptr = true_mod(int(m_current_ptr) - i, max_steps());
    return Eigen::Map<const Eigen::VectorXd>(
        m_v_prevs.data() + ptr * m_v_prevs.rows(), m_v_prevs.rows());
}

Eigen::Map<const Eigen::VectorXd>
ImplicitTimeIntegrator::a_prev(const int i) const
{
    assert(i < available_steps());
    const int ptr = true_mod(int(m_current_ptr) - i, max_steps());
    return Eigen::Map<const Eigen::VectorXd>(
        m_a_prevs.data() + ptr * m_a_prevs.rows(), m_a_prevs.rows());
}

rigid::AffinePose ImplicitTimeIntegrator::pose(
    Eigen::ConstRef<Eigen::VectorXd> x, const size_t i) const
{
    rigid::AffinePose pose;

    pose.position = x.segment(i * (m_pos_ndof + m_rot_ndof), m_pos_ndof);
    if (m_rot_ndof == 1) {
        pose.rotation.resize(1, 1);
        pose.rotation(0, 0) = x(i * (m_pos_ndof + m_rot_ndof) + m_pos_ndof);
    } else {
        pose.rotation =
            x.segment(i * (m_pos_ndof + m_rot_ndof) + m_pos_ndof, m_rot_ndof)
                .reshaped(3, 3);
    }

    return pose;
}

std::vector<rigid::AffinePose> ImplicitTimeIntegrator::predicted_pose() const
{
    const Eigen::VectorXd x_hat = predicted_positions();

    std::vector<rigid::AffinePose> predicted(m_num_bodies);
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, m_num_bodies),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                predicted[i] = pose(x_hat, i);
            }
        });

    return predicted;
}

} // namespace ipc::dynamics
