#include <catch2/catch_test_macros.hpp>

#include <tests/utils.hpp>

#include <finitediff.hpp>

#include <ipc/dynamics/affine/orthogonality_potential.hpp>

#include <iostream>

using namespace ipc;
using namespace ipc::rigid;
using namespace ipc::affine;

// Helper function to generate a RigidBody
auto rigid_bodies()
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
    initial_poses.push_back(Pose::Identity(3)); // Initial pose at the origin

    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { density }, initial_poses);

    (*bodies)[0].set_external_force(
        Pose(Eigen::Vector3d::Random(), Eigen::Vector3d::Random()));

    return bodies;
}

TEST_CASE("Orthogonality potential", "[affine]")
{
    auto bodies = rigid_bodies();

    OrthogonalityPotential V_perp(1);

    VectorMax12d x = VectorMax12d::Random(12);

    double energy = V_perp(*bodies, x);

    // Since we don't have a ground truth, we can only check if the energy is a
    // valid number
    CHECK(std::isfinite(energy));

    VectorMax12d analytical_gradient;
    MatrixMax12d analytical_hessian;

    {
        analytical_gradient = V_perp.gradient(*bodies, x);
        REQUIRE(analytical_gradient.squaredNorm() > 1e-8);

        // Compute the gradient using finite differences
        auto f = [&](const Eigen::VectorXd& x_arg) {
            return V_perp(*bodies, x_arg);
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
        analytical_hessian = V_perp.hessian(*bodies, x);

        // Numerical hessian calculation
        auto f = [&](const Eigen::VectorXd& x_arg) {
            return V_perp.gradient(*bodies, x_arg);
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
        Eigen::VectorXd newton_direction =
            -analytical_hessian.ldlt().solve(analytical_gradient);
        // std::cout << "Newton Direction:\n" << newton_direction << "\n\n";

        if (newton_direction.dot(analytical_gradient) < 0.0) {
            CHECK(true);
        } else {
            analytical_hessian =
                V_perp.hessian(*bodies, x, PSDProjectionMethod::ABS);

            newton_direction =
                -analytical_hessian.ldlt().solve(analytical_gradient);

            // std::cout << "Newton Direction:\n" << newton_direction << "\n\n";

            CHECK(newton_direction.dot(analytical_gradient) < 0.0);
        }
    }
}