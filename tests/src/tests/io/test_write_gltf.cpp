#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/io/read_gltf.hpp>
#include <ipc/io/write_gltf.hpp>

#include <igl/PI.h>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Write and read glTF round-trip", "[io][gltf]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    Eigen::MatrixXd V_box = V;
    V_box.col(0) *= 2.0; // Anisotropic body (nontrivial R₀)

    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[0].position = Eigen::Vector3d(0, 1, 0);
    initial_poses[1].position = Eigen::Vector3d(2, 0.5, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V_box, V }, { E, E }, { F, F }, { 1000.0, 500.0 }, initial_poses);

    // Synthesize a short pose history (falling and rotating)
    const double dt = 0.1;
    std::list<std::vector<Pose>> pose_history;
    for (int step = 0; step < 5; step++) {
        std::vector<Pose> poses = initial_poses;
        for (size_t i = 0; i < poses.size(); i++) {
            poses[i].position.y() -= 0.1 * step;
            poses[i].rotation = rotation_matrix_to_vector(
                Eigen::AngleAxisd(0.1 * step, Eigen::Vector3d::UnitZ())
                    .toRotationMatrix()
                * poses[i].rotation_matrix());
        }
        pose_history.push_back(poses);
    }

    const bool write_binary = GENERATE(true, false);
    const std::string filename =
        (std::filesystem::temp_directory_path()
         / (write_binary ? "ipctk_round_trip.glb" : "ipctk_round_trip.gltf"))
            .string();

    REQUIRE(write_gltf(
        filename, *bodies, pose_history, dt, /*embed_buffers=*/true,
        write_binary));

    // Round trip
    auto [bodies_in, poses_in] = read_gltf(filename);
    REQUIRE(bodies_in != nullptr);

    REQUIRE(bodies_in->num_bodies() == bodies->num_bodies());
    REQUIRE(poses_in.size() == bodies->num_bodies());

    for (size_t i = 0; i < bodies->num_bodies(); i++) {
        CHECK(bodies_in->body_num_vertices(i) == bodies->body_num_vertices(i));
        CHECK(bodies_in->body_num_faces(i) == bodies->body_num_faces(i));

        // The world-space vertices at the initial pose must round-trip
        // (through the R₀ baking on write and re-diagonalization on read).
        // Vertex order is preserved by the writer.
        const Eigen::MatrixXd V_out =
            bodies->body_vertices(i, pose_history.front()[i]);
        const Eigen::MatrixXd V_in = bodies_in->body_vertices(i, poses_in[i]);
        REQUIRE(V_in.rows() == V_out.rows());
        CHECK((V_in - V_out).cwiseAbs().maxCoeff() < 1e-5); // float buffers
    }

    std::filesystem::remove(filename);
}
