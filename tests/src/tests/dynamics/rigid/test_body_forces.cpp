#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <finitediff.hpp>

#include <ipc/dynamics/rigid/body_forces.hpp>

using namespace ipc;
using namespace ipc::rigid;

namespace {

// Helper function to generate a RigidBody
RigidBodies rigid_bodies()
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));

    const double L = 0.5, W = 1.0, H = 2.0;
    V.col(0).array() *= L;
    V.col(1).array() *= W;
    V.col(2).array() *= H;

    const double density = GENERATE(1.0, 2.0, 3.0);

    std::vector<Pose> initial_poses;
    initial_poses.push_back(Pose::Zero(3)); // Initial pose at the origin

    RigidBodies bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { density }, initial_poses);

    bodies[0].set_external_force(
        Pose(Eigen::Vector3d::Random(), Eigen::Vector3d::Random()));

    return bodies;
}

std::shared_ptr<const ImplicitEuler> time_integrator(const RigidBodies& bodies)
{
    Eigen::VectorXd x = Eigen::VectorXd::Zero(12 * bodies.num_bodies());
    for (int i = 0; i < bodies.num_bodies(); ++i) {
        x.segment<9>(i * 12 + 3) = Eigen::Matrix3d::Identity().reshaped();
    }
    const Eigen::VectorXd v = Eigen::VectorXd::Random(12 * bodies.num_bodies());
    const Eigen::VectorXd a = Eigen::VectorXd::Random(12 * bodies.num_bodies());
    const double dt = 0.01; // Arbitrary time step

    return std::make_shared<const ImplicitEuler>(
        x, v, a, dt, bodies.num_bodies());
}

} // namespace

TEST_CASE(
    "BodyForces energy and derivative",
    "[body_forces][energy][gradient][hessian]")
{
    RigidBodies bodies = rigid_bodies();
    BodyForces body_forces(time_integrator(bodies));
    body_forces.update(bodies);

    VectorMax6d x = VectorMax6d::Random(6);
    VectorMax3d force = body_forces.forces()[0];
    MatrixMax3d torque = body_forces.torques()[0];

    double energy = body_forces(bodies[0], x, force, torque);

    // Since we don't have a ground truth, we can only check if the energy is a
    // valid number
    CHECK(std::isfinite(energy));

    VectorMax6d analytical_gradient;
    MatrixMax6d analytical_hessian;

    {
        analytical_gradient = body_forces.gradient(bodies[0], x, force, torque);

        // Numerical gradient calculation
        auto f = [&](const Eigen::VectorXd& x_arg) {
            return body_forces(
                bodies[0], x_arg, force,
                torque); // Capture other variables by value
        };
        Eigen::VectorXd numerical_gradient;
        fd::finite_gradient(x, f, numerical_gradient);

        // std::cout << "Analytical Gradient:\n" << analytical_gradient <<
        // "\n\n"; std::cout << "Numerical Gradient:\n" << numerical_gradient <<
        // "\n\n";

        CHECK(fd::compare_gradient(analytical_gradient, numerical_gradient));
    }

    {
        analytical_hessian = body_forces.hessian(bodies[0], x, force, torque);

        // Numerical hessian calculation
        auto f = [&](const Eigen::VectorXd& x_arg) {
            return body_forces.gradient(bodies[0], x_arg, force, torque);
        };
        Eigen::MatrixXd numerical_hessian;
        fd::finite_jacobian(x, f, numerical_hessian);

        // std::cout << "Analytical Hessian:\n" << analytical_hessian << "\n\n";
        // std::cout << "Numerical Hessian:\n" << numerical_hessian << "\n\n";

        // Compare analytical and numerical Hessians
        CHECK(fd::compare_jacobian(analytical_hessian, numerical_hessian));
    }

    // Newton direction
    {
        Eigen::VectorXd newton_direction =
            -analytical_hessian.ldlt().solve(analytical_gradient);
        // std::cout << "Newton Direction:\n" << newton_direction << "\n\n";

        if (newton_direction.dot(analytical_gradient) < 0.0) {
            CHECK(true);
        } else {
            analytical_hessian = body_forces.hessian(
                bodies[0], x, force, torque, PSDProjectionMethod::ABS);

            newton_direction =
                -analytical_hessian.ldlt().solve(analytical_gradient);

            // std::cout << "Newton Direction:\n" << newton_direction << "\n\n";

            CHECK(newton_direction.dot(analytical_gradient) < 0.0);
        }
    }
}

TEST_CASE(
    "BodyForces total energy and derivatives",
    "[rigid][body_forces][total_energy][total_gradient][total_hessian]")
{
    RigidBodies bodies = rigid_bodies();
    BodyForces body_forces(time_integrator(bodies));
    body_forces.set_gravity(Eigen::Vector3d(0, -9.81, 0));
    body_forces.update(bodies);

    VectorMax6d x = VectorMax6d::Random(6);

    REQUIRE(body_forces.forces().size() == 1);
    REQUIRE(body_forces.torques().size() == 1);

    const auto& force = body_forces.forces()[0];
    const auto& torque = body_forces.torques()[0];

    CHECK(body_forces(bodies, x) == body_forces(bodies[0], x, force, torque));

    CHECK(
        fd::compare_gradient(
            body_forces.gradient(bodies, x),
            body_forces.gradient(bodies[0], x, force, torque)));

    CHECK(
        body_forces.hessian(bodies, x, PSDProjectionMethod::ABS)
        == body_forces.hessian(
            bodies[0], x, force, torque, PSDProjectionMethod::ABS));
}