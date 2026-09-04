#include "rigid_trajectory.hpp"

namespace ipc::rigid {

RigidTrajectory::RigidTrajectory(
    Eigen::ConstRef<VectorMax3d> rest_position,
    const Pose& pose_t0,
    const Pose& pose_t1)
    : m_rest_position(rest_position)
    , m_pose_t0(pose_t0)
    , m_pose_t1(pose_t1)
    , m_radius(rest_position.norm())
    , m_dtheta_norm((pose_t1.rotation - pose_t0.rotation).norm())
{
    assert(pose_t0.position.size() == rest_position.size());
    assert(pose_t1.position.size() == rest_position.size());
    assert(pose_t0.rotation.size() == pose_t1.rotation.size());
}

VectorMax3d RigidTrajectory::operator()(const double t) const
{
    const Pose pose(
        (1 - t) * m_pose_t0.position + t * m_pose_t1.position,
        (1 - t) * m_pose_t0.rotation + t * m_pose_t1.rotation);
    return pose.rotation_matrix() * m_rest_position + pose.position;
}

double RigidTrajectory::max_distance_from_linear(
    const double t0, const double t1) const
{
    return rigid_chord_deviation(m_radius, (t1 - t0) * m_dtheta_norm);
}

} // namespace ipc::rigid
