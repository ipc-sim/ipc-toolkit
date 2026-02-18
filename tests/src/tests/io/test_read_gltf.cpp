// Test the mass utilities.

#include <catch2/catch_test_macros.hpp>

#include <tests/config.hpp>

#include <ipc/io/read_gltf.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Read GLTF (5-cubes)", "[io][rigid]")
{
    std::shared_ptr<RigidBodies> bodies;
    std::vector<Pose> initial_poses;

    std::tie(bodies, initial_poses) = read_gltf(
        tests::DATA_DIR / "rigid-ipc/16-unit-tests/5-cubes.glb", true);

    REQUIRE(bodies != nullptr);

    CHECK(bodies->num_bodies() == 5);
    CHECK(initial_poses.size() == 5);
    CHECK(bodies->planes.size() == 1);
    CAPTURE(bodies->planes[0].normal, bodies->planes[0].origin);
    CHECK(bodies->planes[0].origin.isApprox(Eigen::Vector3d(0, -1, 0)));
    CHECK(bodies->planes[0].normal.isApprox(Eigen::Vector3d(0, 1, 0)));
}

TEST_CASE("Read GLTF (incline plane)", "[io][rigid]")
{
    std::shared_ptr<RigidBodies> bodies;
    std::vector<Pose> initial_poses;

    std::tie(bodies, initial_poses) = read_gltf(
        tests::DATA_DIR
            / "rigid-ipc/18-high-school-physics-friction-test-mu=0.5.glb",
        true);

    REQUIRE(bodies != nullptr);

    CHECK(bodies->num_bodies() == 1);
    CHECK(initial_poses.size() == 1);
    CHECK(bodies->planes.size() == 1);
    CAPTURE(bodies->planes[0].normal, bodies->planes[0].origin);
    CHECK(bodies->planes[0].origin.isApprox(Eigen::Vector3d(2, -1.565, 0)));

    const double theta = std::atan(0.5);
    Eigen::Vector3d n(std::sin(theta), std::cos(theta), 0);
    CAPTURE(theta, n);
    CHECK(bodies->planes[0].normal.isApprox(n, 1e-6));
}

TEST_CASE("Read GLTF (card house)", "[io][rigid]")
{
    std::shared_ptr<RigidBodies> bodies;
    std::vector<Pose> initial_poses;

    std::tie(bodies, initial_poses) = read_gltf(
        tests::DATA_DIR / "rigid-ipc/10-codimensional-card-house.glb", true);

    REQUIRE(bodies != nullptr);

    CHECK(bodies->num_bodies() == 17);
    CHECK(initial_poses.size() == bodies->num_bodies());
    CHECK(bodies->planes.size() == 1);
    CAPTURE(bodies->planes[0].normal, bodies->planes[0].origin);
    CHECK(bodies->planes[0].origin.isApprox(Eigen::Vector3d(0, 0, 0)));
    CHECK(bodies->planes[0].normal.isApprox(Eigen::Vector3d(0, 1, 0)));
}