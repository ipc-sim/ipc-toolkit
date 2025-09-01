// Test the mass utilities.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/mass.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Mass properties", "[rigid][mass]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));

    const double L = 0.5, W = 1.0, H = 2.0;

    V.col(0).array() *= L;
    V.col(1).array() *= W;
    V.col(2).array() *= H;

    const double density = GENERATE(1.0, 2.0, 3.0);

    double total_mass;
    VectorMax3d center;
    MatrixMax3d inertia;
    compute_mass_properties(V, F, density, total_mass, center, inertia);

    const double m = density * L * W * H;
    const double Ixx = m * (W * W + H * H) / 12.0;
    const double Iyy = m * (L * L + H * H) / 12.0;
    const double Izz = m * (L * L + W * W) / 12.0;
    Eigen::Matrix3d expected_inertia =
        Eigen::DiagonalMatrix<double, 3>(Ixx, Iyy, Izz);

    CAPTURE(total_mass, center, inertia);
    REQUIRE(total_mass == Catch::Approx(m).margin(1e-6));
    REQUIRE(center.isApprox(Eigen::Vector3d(0.0, 0.0, 0.0), 1e-6));
    REQUIRE(inertia.isApprox(expected_inertia, 1e-6));
}