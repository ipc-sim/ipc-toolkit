#include "body_potentials.hpp"

#include <ipc/dynamics/time_integration/time_integrator.hpp>

namespace ipc::dynamics {

BodyPotentials::BodyPotentials(
    const rigid::RigidBodies& bodies,
    std::shared_ptr<const ImplicitTimeIntegrator> time_integrator,
    std::shared_ptr<const ToAffine> to_affine,
    const double orthogonality_stiffness,
    const affine::AugmentedLagrangian::Params& al_params)
    : m_time_integrator(std::move(time_integrator))
    , m_to_affine(std::move(to_affine))
    , m_inertial(bodies, m_time_integrator)
    , m_al(al_params)
{
    // The orthogonality penalty keeps A near SO(dim); it is only meaningful for
    // the affine (identity) parameterization — rigid rotations are exact.
    if (m_to_affine->is_identity()) {
        m_orthogonality.emplace(orthogonality_stiffness);
    }
}

void BodyPotentials::update(const rigid::RigidBodies& bodies)
{
    m_inertial.update(bodies);

    // Start-of-step poses used to transform world torques into body space.
    std::vector<affine::Pose> poses(bodies.num_bodies());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        poses[i] = m_time_integrator->pose(i);
    }
    m_body_forces.update(bodies, poses);
}

Eigen::VectorXd BodyPotentials::affine_gradient(
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> y,
    const double scaling,
    const bool include_al) const
{
    Eigen::VectorXd g = m_inertial.gradient(bodies, y)
        + scaling * m_body_forces.gradient(bodies, y);
    if (m_orthogonality) {
        g += scaling * m_orthogonality->gradient(bodies, y);
    }
    if (include_al && m_al.active()) {
        g += scaling * m_al.gradient(bodies, y);
    }
    return g;
}

double BodyPotentials::energy(
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const double scaling,
    const bool include_al) const
{
    const Eigen::VectorXd y = m_to_affine->to_affine(x);

    double E = m_inertial(bodies, y) + scaling * m_body_forces(bodies, y);
    if (m_orthogonality) {
        E += scaling * (*m_orthogonality)(bodies, y);
    }
    if (include_al && m_al.active()) {
        E += scaling * m_al(bodies, y);
    }
    return E;
}

Eigen::VectorXd BodyPotentials::gradient(
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const double scaling,
    const bool include_al) const
{
    const Eigen::VectorXd y = m_to_affine->to_affine(x);
    return m_to_affine->apply_gradient(
        x, affine_gradient(bodies, y, scaling, include_al));
}

Eigen::SparseMatrix<double> BodyPotentials::hessian(
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const double scaling,
    const PSDProjectionMethod project_hessian_to_psd,
    const bool include_al) const
{
    const Eigen::VectorXd y = m_to_affine->to_affine(x);

    // Affine-coordinate Hessian (body forces are linear → zero Hessian). The
    // term Hessians are left unprojected; apply_hessian projects once.
    Eigen::SparseMatrix<double> H = m_inertial.mass_matrix();
    if (m_orthogonality) {
        H += scaling
            * m_orthogonality->hessian(bodies, y, PSDProjectionMethod::NONE);
    }
    if (include_al && m_al.active()) {
        H += scaling * m_al.hessian(bodies, y);
    }

    const Eigen::VectorXd g = affine_gradient(bodies, y, scaling, include_al);
    return m_to_affine->apply_hessian(x, g, H, project_hessian_to_psd);
}

void BodyPotentials::init_augmented_lagrangian(
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x0,
    const std::vector<rigid::Pose>& targets)
{
    m_al.init(bodies, m_to_affine->to_affine(x0), targets);
}

void BodyPotentials::update_augmented_lagrangian(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    m_al.update(bodies, m_to_affine->to_affine(x));
}

} // namespace ipc::dynamics
