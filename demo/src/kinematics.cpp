#include "kinematics.hpp"

#include <ipc/candidates/candidates.hpp>
#include <ipc/ccd/additive_ccd.hpp>
#include <ipc/ccd/nonlinear_ccd.hpp>
#include <ipc/dynamics/affine/affine_dof.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/rigid_candidates.hpp>

#include <cassert>

namespace ipc::demo {

Kinematics::Kinematics(const std::shared_ptr<const rigid::RigidBodies>& bodies)
    : m_bodies(bodies)
{
    assert(m_bodies != nullptr);
}

Kinematics::~Kinematics() = default;

namespace {

    /// @brief Rigid parameterization: x_i = [p; θ] (6-DOF in 3D with θ the
    /// rotation vector; 3-DOF in 2D with θ the scalar angle).
    class RigidKinematics : public Kinematics {
    public:
        explicit RigidKinematics(
            const std::shared_ptr<const rigid::RigidBodies>& bodies)
            : Kinematics(bodies)
        {
            m_to_affine = std::make_shared<dynamics::RigidToAffine>(
                bodies->dim(), bodies->num_bodies());
        }

        int pose_ndof() const { return m_bodies->dim() == 2 ? 3 : 6; }

        int ndof() const override
        {
            return pose_ndof() * int(m_bodies->num_bodies());
        }

        Eigen::MatrixXd
        world_vertices(Eigen::ConstRef<Eigen::VectorXd> x) const override
        {
            return m_bodies->vertices(x);
        }

        Eigen::VectorXd map_vertex_gradient(
            Eigen::ConstRef<Eigen::VectorXd> x,
            Eigen::ConstRef<Eigen::VectorXd> grad) const override
        {
            return m_bodies->to_rigid_dof(
                rigid::Pose::to_poses(x, m_bodies->dim()), grad);
        }

        Eigen::SparseMatrix<double> map_vertex_hessian(
            Eigen::ConstRef<Eigen::VectorXd> x,
            Eigen::ConstRef<Eigen::VectorXd> grad,
            const Eigen::SparseMatrix<double>& hess) const override
        {
            return m_bodies->to_rigid_dof(
                rigid::Pose::to_poses(x, m_bodies->dim()), grad, hess);
        }

        std::vector<affine::Pose>
        poses(Eigen::ConstRef<Eigen::VectorXd> x) const override
        {
            return m_to_affine->poses(x);
        }

        Eigen::VectorXd
        dof(const std::vector<rigid::Pose>& poses) const override
        {
            return rigid::Pose::from_poses(poses);
        }

        Eigen::VectorXd
        to_integrator_state(Eigen::ConstRef<Eigen::VectorXd> x) const override
        {
            // The integrator state is affine-shaped [p; vec(Q)] per body; the
            // log map keeps the optimization variable θ bounded on the way
            // back.
            return m_to_affine->to_affine(x);
        }

        Eigen::VectorXd
        from_integrator_state(Eigen::ConstRef<Eigen::VectorXd> X) const override
        {
            return m_to_affine->from_affine(X);
        }

        void update_candidates(
            Eigen::ConstRef<Eigen::VectorXd> x,
            const double inflation_radius,
            BroadPhase* broad_phase) override
        {
            // Discrete (distance) candidates through the rigid body-pair
            // broad phase: reuses each body's static rest-frame BVH instead
            // of rebuilding a world-space one [Ferguson et al. 2021].
            m_candidates.build(
                *m_bodies, rigid::Pose::to_poses(x, m_bodies->dim()),
                inflation_radius, broad_phase);
        }

        void update_candidates(
            Eigen::ConstRef<Eigen::VectorXd> x0,
            Eigen::ConstRef<Eigen::VectorXd> x1,
            const double inflation_radius,
            BroadPhase* broad_phase) override
        {
            m_candidates.build(
                *m_bodies, rigid::Pose::to_poses(x0, m_bodies->dim()),
                rigid::Pose::to_poses(x1, m_bodies->dim()), inflation_radius,
                broad_phase);
        }

        double max_step_size(
            Eigen::ConstRef<Eigen::VectorXd> x0,
            Eigen::ConstRef<Eigen::VectorXd> x1,
            const double inflation_radius,
            BroadPhase* broad_phase,
            const double min_distance) override
        {
            update_candidates(x0, x1, inflation_radius, broad_phase);
            return m_candidates.compute_collision_free_stepsize(
                *m_bodies, rigid::Pose::to_poses(x0, m_bodies->dim()),
                rigid::Pose::to_poses(x1, m_bodies->dim()), min_distance,
                m_ccd);
        }

        const Candidates& candidates() const override { return m_candidates; }

        void clear_candidates() override { m_candidates.clear(); }

    private:
        /// @brief Rigid collision candidates (rotational AABB inflation).
        rigid::RigidCandidates m_candidates;
        /// @brief Nonlinear CCD along the curved rigid trajectories.
        NonlinearCCD m_ccd;
    };

    /// @brief Affine parameterization: x_i = [p; vec(A) column-major]
    /// (12-DOF per body in 3D; 6-DOF in 2D).
    class AffineKinematics : public Kinematics {
    public:
        explicit AffineKinematics(
            const std::shared_ptr<const rigid::RigidBodies>& bodies)
            : Kinematics(bodies)
            , m_J_all(affine::affine_jacobian(*bodies))
        {
            m_to_affine = std::make_shared<dynamics::AffineToAffine>(
                bodies->dim(), bodies->num_bodies());
        }

        int body_ndof() const
        {
            return m_bodies->dim() + m_bodies->dim() * m_bodies->dim();
        }

        int ndof() const override
        {
            return body_ndof() * int(m_bodies->num_bodies());
        }

        Eigen::MatrixXd
        world_vertices(Eigen::ConstRef<Eigen::VectorXd> x) const override
        {
            return affine::vertices(*m_bodies, x);
        }

        Eigen::VectorXd map_vertex_gradient(
            Eigen::ConstRef<Eigen::VectorXd> x,
            Eigen::ConstRef<Eigen::VectorXd> grad) const override
        {
            return m_J_all.transpose() * grad;
        }

        Eigen::SparseMatrix<double> map_vertex_hessian(
            Eigen::ConstRef<Eigen::VectorXd> x,
            Eigen::ConstRef<Eigen::VectorXd> grad,
            const Eigen::SparseMatrix<double>& hess) const override
        {
            // V(x) is linear in the affine DOFs, so there is no second-order
            // term (grad is unused).
            return m_J_all.transpose() * hess * m_J_all;
        }

        std::vector<affine::Pose>
        poses(Eigen::ConstRef<Eigen::VectorXd> x) const override
        {
            return m_to_affine->poses(x);
        }

        Eigen::VectorXd
        dof(const std::vector<rigid::Pose>& poses) const override
        {
            const int dim = m_bodies->dim();
            Eigen::VectorXd x(body_ndof() * poses.size());
            for (size_t i = 0; i < poses.size(); ++i) {
                x.segment(body_ndof() * i, dim) = poses[i].position;
                x.segment(body_ndof() * i + dim, dim * dim) =
                    poses[i].rotation_matrix().reshaped();
            }
            return x;
        }

        Eigen::VectorXd
        to_integrator_state(Eigen::ConstRef<Eigen::VectorXd> x) const override
        {
            return m_to_affine->to_affine(x); // identity: DOFs are the state
        }

        Eigen::VectorXd
        from_integrator_state(Eigen::ConstRef<Eigen::VectorXd> X) const override
        {
            return m_to_affine->from_affine(X); // identity: DOFs are the state
        }

        void update_candidates(
            Eigen::ConstRef<Eigen::VectorXd> x,
            const double inflation_radius,
            BroadPhase* broad_phase) override
        {
            m_candidates.build(
                *m_bodies, world_vertices(x), inflation_radius, broad_phase);
        }

        void update_candidates(
            Eigen::ConstRef<Eigen::VectorXd> x0,
            Eigen::ConstRef<Eigen::VectorXd> x1,
            const double inflation_radius,
            BroadPhase* broad_phase) override
        {
            m_candidates.build(
                *m_bodies, world_vertices(x0), world_vertices(x1),
                inflation_radius, broad_phase);
        }

        double max_step_size(
            Eigen::ConstRef<Eigen::VectorXd> x0,
            Eigen::ConstRef<Eigen::VectorXd> x1,
            const double inflation_radius,
            BroadPhase* broad_phase,
            const double min_distance) override
        {
            const Eigen::MatrixXd V0 = world_vertices(x0);
            const Eigen::MatrixXd V1 = world_vertices(x1);

            m_candidates.build(
                *m_bodies, V0, V1, inflation_radius, broad_phase);

            // The per-iterate vertex trajectories are linear in the step
            // parameter, so linear CCD is exact [Lan et al. 2022, §4.1].
            return m_candidates.compute_collision_free_stepsize(
                *m_bodies, V0, V1, min_distance, m_ccd);
        }

        const Candidates& candidates() const override { return m_candidates; }

        void clear_candidates() override { m_candidates.clear(); }

    private:
        /// @brief Constant Jacobian dV/dx (V is linear in the affine DOFs).
        Eigen::SparseMatrix<double> m_J_all;
        /// @brief Collision candidates on the linear vertex trajectories.
        Candidates m_candidates;
        /// @brief Linear CCD (exact for affine per-iterate trajectories).
        AdditiveCCD m_ccd;
    };

} // namespace

std::unique_ptr<Kinematics> Kinematics::create_rigid(
    const std::shared_ptr<const rigid::RigidBodies>& bodies)
{
    return std::make_unique<RigidKinematics>(bodies);
}

std::unique_ptr<Kinematics> Kinematics::create_affine(
    const std::shared_ptr<const rigid::RigidBodies>& bodies)
{
    return std::make_unique<AffineKinematics>(bodies);
}

} // namespace ipc::demo
