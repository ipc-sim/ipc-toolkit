#include "kinematic_driver.hpp"

namespace ipc::demo {

rigid::Pose KinematicDriver::target(
    const rigid::Pose& current,
    const rigid::Pose& velocity,
    const double dt) const
{
    if (is_scripted()) {
        // Scripted absolute target pose (advanced each step).
        return m_scripted_poses.front();
    }

    // Integrate the body's prescribed velocity: p̂ = p + dt·v and, for the
    // rotation, θ̂ = θ + dt·ω in 2D (a linear space) or the exact exponential
    // Q̂ = exp(dt·[ω]×)·Q in 3D.
    rigid::Pose result = current;
    result.position += dt * velocity.position;

    if (current.rotation.size() == 1) {
        result.rotation(0) += dt * velocity.rotation(0);
    } else {
        result.rotation = rigid::rotation_matrix_to_vector(
            rigid::rotation_vector_to_matrix(dt * velocity.rotation)
            * rigid::rotation_vector_to_matrix(current.rotation));
    }

    return result;
}

} // namespace ipc::demo
