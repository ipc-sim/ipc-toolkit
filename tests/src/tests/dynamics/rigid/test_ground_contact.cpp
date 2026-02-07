#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/ground_contact.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <iostream>

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

    return bodies;
}

} // namespace

TEST_CASE(
    "GroundContact energy and derivative",
    "[ground_contact][energy][gradient][hessian]")
{
    RigidBodies bodies = rigid_bodies();
    GroundContact ground_contact(0.0);
    ground_contact.set_dhat(0.1);

    VectorMax6d x = VectorMax6d::Zero(6);
    SECTION("Translation only") { x(1) = 0.55; }
    SECTION("Rotation")
    {
        x(1) = 0.55;
        x(3) = 15.0 * igl::PI / 180.0; // 15 degrees roll
    }

    double energy = ground_contact(bodies, 0, Pose(x.segment(0, 6)));

    // Since we don't have a ground truth, we can only check if the energy is a
    // valid number
    CHECK(std::isfinite(energy));

    Eigen::VectorXd analytical_gradient;
    Eigen::MatrixXd analytical_hessian;

    {
        analytical_gradient =
            ground_contact.gradient(bodies, 0, Pose(x.segment(0, 6)));

        // Numerical gradient calculation
        auto f = [&](const Eigen::VectorXd& x_arg) {
            return ground_contact(bodies, 0, Pose(x_arg.segment(0, 6)));
        };
        Eigen::VectorXd numerical_gradient;
        fd::finite_gradient(x, f, numerical_gradient);

        CHECK(fd::compare_gradient(analytical_gradient, numerical_gradient));
        if (!fd::compare_gradient(analytical_gradient, numerical_gradient)) {
            std::cout << "Analytical Gradient:\n"
                      << analytical_gradient << "\n\n";
            std::cout << "Numerical Gradient:\n"
                      << numerical_gradient << "\n\n";
        }
    }

    {
        analytical_hessian =
            ground_contact.hessian(bodies, 0, Pose(x.segment(0, 6)));

        // Numerical hessian calculation
        auto f = [&](const Eigen::VectorXd& x_arg) {
            return ground_contact.gradient(
                bodies, 0, Pose(x_arg.segment(0, 6)));
        };
        Eigen::MatrixXd numerical_hessian;
        fd::finite_jacobian(x, f, numerical_hessian);

        // Compare analytical and numerical Hessians
        CHECK(fd::compare_jacobian(analytical_hessian, numerical_hessian));
        if (!fd::compare_jacobian(analytical_hessian, numerical_hessian)) {
            std::cout << "Analytical Hessian:\n"
                      << analytical_hessian << "\n\n";
            std::cout << "Numerical Hessian:\n" << numerical_hessian << "\n\n";
        }
    }

    // Newton direction
    {
        // Make the Hessian full rank
        analytical_hessian += 1e-8 * Eigen::MatrixXd::Identity(6, 6);

        Eigen::VectorXd newton_direction =
            -analytical_hessian.ldlt().solve(analytical_gradient);

        if (newton_direction.dot(analytical_gradient) < 0.0) {
            CHECK(true);
        } else {
            analytical_hessian = ground_contact.hessian(
                bodies, 0, Pose(x.segment(0, 6)), PSDProjectionMethod::ABS);

            newton_direction =
                -analytical_hessian.ldlt().solve(analytical_gradient);

            CHECK(newton_direction.dot(analytical_gradient) < 0.0);
        }
    }
}

TEST_CASE(
    "GroundContact total energy and derivatives",
    "[rigid][ground_contact][total_energy][total_gradient][total_hessian]")
{
    RigidBodies bodies = rigid_bodies();
    GroundContact ground_contact(0.0);
    ground_contact.set_dhat(0.1);

    VectorMax6d x = VectorMax6d::Random(6);

    CHECK(ground_contact(bodies, x) == ground_contact(bodies, 0, Pose(x)));

    CHECK(
        fd::compare_gradient(
            ground_contact.gradient(bodies, x),
            ground_contact.gradient(bodies, 0, Pose(x))));

    CHECK(
        ground_contact.hessian(bodies, x, PSDProjectionMethod::ABS)
        == ground_contact.hessian(
            bodies, 0, Pose(x), PSDProjectionMethod::ABS));
}