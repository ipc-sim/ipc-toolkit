// Test the mass utilities.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <ipc/io/read_gltf.hpp>
#include <ipc/io/write_gltf.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Read GLTF (5-cubes)", "[io][rigid]")
{
    std::shared_ptr<RigidBodies> bodies;
    std::vector<Pose> initial_poses;

    std::tie(bodies, initial_poses) = read_gltf(
        (tests::DATA_DIR / "rigid-ipc/16-unit-tests/5-cubes.glb").string(),
        true);

    REQUIRE(bodies != nullptr);

    CHECK(bodies->num_bodies() == 5);
    CHECK(initial_poses.size() == 5);
    CHECK(bodies->planes.size() == 1);
    CAPTURE(bodies->planes[0].normal(), bodies->planes[0].offset());
    CHECK(bodies->planes[0].offset() == Catch::Approx(1));
    CHECK(bodies->planes[0].normal().isApprox(Eigen::Vector3d(0, 1, 0)));
}

TEST_CASE("Read GLTF (incline plane)", "[io][rigid]")
{
    std::shared_ptr<RigidBodies> bodies;
    std::vector<Pose> initial_poses;

    std::tie(bodies, initial_poses) = read_gltf(
        (tests::DATA_DIR
         / "rigid-ipc/18-high-school-physics-friction-test-mu=0.5.glb")
            .string(),
        true);

    REQUIRE(bodies != nullptr);

    CHECK(bodies->num_bodies() == 1);
    CHECK(initial_poses.size() == 1);
    CHECK(bodies->planes.size() == 1);
    CAPTURE(bodies->planes[0].normal(), bodies->planes[0].offset());
    CHECK(bodies->planes[0].offset() == Catch::Approx(0.50535359193774421));

    const double theta = std::atan(0.5);
    Eigen::Vector3d n(std::sin(theta), std::cos(theta), 0);
    CAPTURE(theta, n);
    CHECK(bodies->planes[0].normal().isApprox(n, 1e-6));
}

TEST_CASE("Read GLTF (card house)", "[io][rigid]")
{
    std::shared_ptr<RigidBodies> bodies;
    std::vector<Pose> initial_poses;

    std::tie(bodies, initial_poses) = read_gltf(
        (tests::DATA_DIR / "rigid-ipc/10-codimensional-card-house.glb")
            .string(),
        true);

    REQUIRE(bodies != nullptr);

    CHECK(bodies->num_bodies() == 17);
    CHECK(initial_poses.size() == bodies->num_bodies());
    CHECK(bodies->planes.size() == 1);
    CAPTURE(bodies->planes[0].normal(), bodies->planes[0].offset());
    CHECK(bodies->planes[0].offset() == 0);
    CHECK(bodies->planes[0].normal().isApprox(Eigen::Vector3d(0, 1, 0)));
}

TEST_CASE("Read GLTF names, zero-scale nodes, and DOF masks", "[io][rigid]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    // Three named cubes in a row; the last two share a name.
    std::vector<Pose> initial_poses(3, Pose::Identity(3));
    for (size_t i = 0; i < initial_poses.size(); i++) {
        initial_poses[i].position = Eigen::Vector3d(3.0 * i, 0.0, 0.0);
    }
    auto bodies = RigidBodies::build_from_meshes(
        { V, V, V }, { E, E, E }, { F, F, F }, { 1000.0, 1000.0, 1000.0 },
        initial_poses, /*convert_planes=*/false, /*is_dof_fixed=*/ {},
        { "Hidden", "Cube", "Cube" });

    // Write an ASCII glTF, then scale the first node to zero.
    const auto tmp = std::filesystem::temp_directory_path();
    const std::string base = (tmp / "ipctk_read_base.gltf").string();
    REQUIRE(write_gltf(
        base, *bodies, std::list<std::vector<Pose>> { initial_poses }, 0.01,
        /*embed_buffers=*/true, /*write_binary=*/false));
    nlohmann::json j;
    {
        std::ifstream f(base);
        f >> j;
    }
    j["nodes"][0]["scale"] = { 0.0, 0.0, 0.0 };
    // write_gltf stores the unique body names, so we restore the duplicate.
    j["nodes"][2]["name"] = "Cube";
    const std::string patched = (tmp / "ipctk_read_patched.gltf").string();
    {
        std::ofstream f(patched);
        f << j.dump();
    }

    VectorMax6b slider(6);
    slider << false, true, true, true, true, true;
    auto [bodies_in, poses_in] = read_gltf(
        patched, /*convert_planes=*/false,
        { { "Cube.1", slider }, { "NoSuchBody", slider } });
    REQUIRE(bodies_in != nullptr);

    // The zero-scale node is skipped. The others keep their node names, made
    // unique, and masks are matched against those unique names.
    REQUIRE(bodies_in->num_bodies() == 2);
    CHECK(poses_in.size() == 2);
    CHECK(bodies_in->body_name(0) == "Cube");
    CHECK(bodies_in->body_name(1) == "Cube.1");

    // DOF masks are applied by name; unlisted bodies are free.
    CHECK(!(*bodies_in)[0].is_dof_fixed().any());
    CHECK((*bodies_in)[1].is_dof_fixed() == slider);

    // Geometry only: every body gets the default density.
    CHECK((*bodies_in)[0].density() == Catch::Approx(1000.0));

    std::filesystem::remove(base);
    std::filesystem::remove(patched);
}
