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

TEST_CASE("2D mass properties", "[2d][mass]")
{
    // A W×H rectangle outline with edge-Voronoi lumped masses.
    const double W = 2.0, H = 1.0, density = GENERATE(1.0, 100.0);

    Eigen::MatrixXd V(4, 2);
    V << -W / 2, -H / 2, W / 2, -H / 2, W / 2, H / 2, -W / 2, H / 2;
    Eigen::MatrixXi E(4, 2);
    E << 0, 1, 1, 2, 2, 3, 3, 0;

    double mass;
    ipc::VectorMax3d center;
    ipc::MatrixMax3d inertia;
    ipc::rigid::compute_mass_properties(V, E, density, mass, center, inertia);

    // Mass = density × perimeter; centered at the origin.
    CHECK(mass == Catch::Approx(density * 2 * (W + H)));
    CHECK(center.norm() == Catch::Approx(0).margin(1e-12));

    // Second moment ∫ρ x̄x̄ᵀ of 4 lumped corner masses m_i = ρ(W + H)/2 at
    // (±W/2, ±H/2).
    const double corner_mass = density * (W + H) / 2;
    REQUIRE(inertia.rows() == 2);
    CHECK(inertia(0, 0) == Catch::Approx(4 * corner_mass * W * W / 4));
    CHECK(inertia(1, 1) == Catch::Approx(4 * corner_mass * H * H / 4));
    CHECK(inertia(0, 1) == Catch::Approx(0).margin(1e-10));
}

TEST_CASE("2D point cloud mass properties", "[2d][mass]")
{
    Eigen::MatrixXd V(3, 2);
    V << 0, 0, 2, 0, 0, 2;
    const Eigen::MatrixXi E(0, 2);

    double mass;
    ipc::VectorMax3d center;
    ipc::MatrixMax3d inertia;
    ipc::rigid::compute_mass_properties(V, E, 1.0, mass, center, inertia);

    CHECK(mass == Catch::Approx(3.0)); // unit mass per point
    CHECK(center.x() == Catch::Approx(2.0 / 3.0));
    CHECK(center.y() == Catch::Approx(2.0 / 3.0));
    REQUIRE(inertia.rows() == 2);
    // Second moment about the COM: Σ (vᵢ − c)(vᵢ − c)ᵀ
    double sxx = 0, syy = 0, sxy = 0;
    for (int i = 0; i < 3; ++i) {
        const Eigen::Vector2d r = V.row(i).transpose() - center;
        sxx += r.x() * r.x();
        syy += r.y() * r.y();
        sxy += r.x() * r.y();
    }
    CHECK(inertia(0, 0) == Catch::Approx(sxx));
    CHECK(inertia(1, 1) == Catch::Approx(syy));
    CHECK(inertia(0, 1) == Catch::Approx(sxy));
}
