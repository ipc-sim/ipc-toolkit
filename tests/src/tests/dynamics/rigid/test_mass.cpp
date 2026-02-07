// Test the mass utilities.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/mass.hpp>

#include <igl/moments.h>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Mass properties", "[rigid][mass]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;

    const double density = GENERATE(1.0, 2.0, 3.0);

    double expected_mass = -1;
    Eigen::Vector3d expected_center;
    Eigen::Matrix3d expected_inertia;
    SECTION("Cube")
    {
        REQUIRE(tests::load_mesh("cube.ply", V, E, F));

        const double L = 0.5, W = 1.0, H = 2.0;

        V.col(0).array() *= L;
        V.col(1).array() *= W;
        V.col(2).array() *= H;

        expected_mass = density * L * W * H;

        expected_center.setZero();

        const double Ixx = expected_mass * (W * W + H * H) / 12.0;
        const double Iyy = expected_mass * (L * L + H * H) / 12.0;
        const double Izz = expected_mass * (L * L + W * W) / 12.0;
        expected_inertia = Eigen::DiagonalMatrix<double, 3>(Ixx, Iyy, Izz);
    }
    SECTION("Bunny (Low-Poly)")
    {
        REQUIRE(tests::load_mesh("bunny (lowpoly).ply", V, E, F));
        double volume;
        igl::moments(V, F, volume, expected_center, expected_inertia);
        expected_center /= volume;
        expected_inertia.array() *= density;
        expected_mass = density * volume;
    }
    SECTION("Bunny")
    {
        REQUIRE(tests::load_mesh("bunny.ply", V, E, F));
        double volume;
        igl::moments(V, F, volume, expected_center, expected_inertia);
        expected_center /= volume;
        expected_inertia.array() *= density;
        expected_mass = density * volume;
    }
    SECTION("Bowl")
    {
        REQUIRE(tests::load_mesh("bowl.ply", V, E, F));
        double volume;
        igl::moments(V, F, volume, expected_center, expected_inertia);
        expected_center /= volume;
        expected_inertia.array() *= density;
        expected_mass = density * volume;
    }

    double total_mass;
    VectorMax3d center;
    MatrixMax3d inertia;
    compute_mass_properties(V, F, density, total_mass, center, inertia);

    REQUIRE(total_mass == Catch::Approx(expected_mass).margin(1e-6));
    {
        CAPTURE(density, center, expected_center);
        REQUIRE(center.isApprox(expected_center, 1e-6));
    }
    {
        CAPTURE(density, inertia, expected_inertia);
        REQUIRE(inertia.isApprox(expected_inertia, 1e-6));
    }
}