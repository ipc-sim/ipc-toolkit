#include <tests/config.hpp>

#ifdef IPC_TOOLKIT_WITH_PROFILER

#include <ipc/ipc.hpp>
#include <ipc/ccd/additive_ccd.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/utils/profiler.hpp>

#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

#include <catch2/catch_all.hpp>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

using namespace ipc;

TEST_CASE("Profiler", "[profiler]")
{
    profiler().reset();

    constexpr int sleep_time_ms = 100;
    constexpr int num_threads = 10;

    {
        IPC_TOOLKIT_PROFILE_BLOCK("Block 1");
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time_ms));
    }

    nlohmann::json data = profiler().data();

    CHECK(data.size() == 1);

    REQUIRE(data.contains("Block 1"));
    nlohmann::json block1 = data.at("Block 1");

    CHECK(block1.size() == 2); // count, time_ms
    CHECK(
        block1["time_ms"].get<double>()
        == Catch::Approx(sleep_time_ms).margin(10));
    CHECK(block1["count"].get<int>() == 1);

    // ---------------------------------------------------------------------

    tbb::parallel_for(0, num_threads, [&](int) {
        // for (int i = 0; i < num_threads; ++i) {
        IPC_TOOLKIT_PROFILE_BLOCK("Block 2");
        std::this_thread::sleep_for(
            std::chrono::milliseconds(sleep_time_ms / num_threads));
        // }
    });

    data = profiler().data();
    profiler().print();

    CHECK(data.size() == 2);

    REQUIRE(data.contains("Block 2"));
    nlohmann::json block2 = data.at("Block 2");
    CHECK(block2.size() == 2); // count, time_ms
    // Only the coordinator thread's iterations are recorded
    CHECK(block2["count"].get<int>() >= 1);

    // ---------------------------------------------------------------------

    auto foo = []() {
        IPC_TOOLKIT_PROFILE_BLOCK("Block 3");
        {
            IPC_TOOLKIT_PROFILE_BLOCK("Block 4");
            tbb::parallel_for(0, num_threads, [&](int) {
                // for (int i = 0; i < num_threads; ++i) {
                {
                    IPC_TOOLKIT_PROFILE_BLOCK("Block 5");
                    std::this_thread::sleep_for(
                        std::chrono::milliseconds(sleep_time_ms / num_threads));
                }

                {
                    IPC_TOOLKIT_PROFILE_BLOCK("Block 6");
                    std::this_thread::sleep_for(
                        std::chrono::milliseconds(sleep_time_ms / num_threads));
                }
                // }
            });
        }
    };

    foo();
    foo(); // Run it twice to check that counts are aggregated correctly

    data = profiler().data();
    profiler().print();

    CHECK(data.size() == 3);

    REQUIRE(data.contains("Block 3"));
    CHECK_FALSE(data.contains("Block 4"));
    nlohmann::json block3 = data.at("Block 3");
    CHECK(block3.size() == 3); // count, time_ms, Block 4
    CHECK(block3["count"].get<int>() == 2);

    REQUIRE(block3.contains("Block 4"));
    nlohmann::json block4 = block3.at("Block 4");

    REQUIRE(block4.contains("Block 5"));
    nlohmann::json block5 = block4.at("Block 5");
    CHECK(block5.size() == 2); // count, time_ms
    // foo() runs twice; coordinator records at least once per call
    CHECK(block5["count"].get<int>() >= 2);
    CHECK(block5["time_ms"].get<double>() < block4["time_ms"].get<double>());

    REQUIRE(block4.contains("Block 6"));
    nlohmann::json block6 = block4.at("Block 6");
    CHECK(block6.size() == 2); // count, time_ms
    CHECK(block6["count"].get<int>() >= 2);
    CHECK(block6["time_ms"].get<double>() < block4["time_ms"].get<double>());

    CHECK(
        block5["time_ms"].get<double>() + block6["time_ms"].get<double>()
        < block4["time_ms"].get<double>());
}

TEST_CASE("Profile full pipeline", "[!benchmark][profiler]")
{
    // const std::string mesh_t0 = "cloth_ball92.ply";
    // const std::string mesh_t1 = "cloth_ball93.ply";
    const std::string mesh_t0 = "rod-twist/3036.ply";
    const std::string mesh_t1 = "rod-twist/3037.ply";

    Eigen::MatrixXd V0_full, V1_full;
    Eigen::MatrixXi F0, F1;

    const bool loaded_t0 = igl::read_triangle_mesh(
        (ipc::tests::DATA_DIR / mesh_t0).string(), V0_full, F0);
    const bool loaded_t1 = igl::read_triangle_mesh(
        (ipc::tests::DATA_DIR / mesh_t1).string(), V1_full, F1);

    if (!loaded_t0 || !loaded_t1) {
        WARN("Skipping profiler test: rod-twist data not found");
        return;
    }

    REQUIRE(F0.rows() == F1.rows());
    REQUIRE(V0_full.rows() == V1_full.rows());

    Eigen::MatrixXi E;
    igl::edges(F0, E);

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0_full, E, F0);
    Eigen::MatrixXd V0 = mesh.vertices(V0_full);
    Eigen::MatrixXd V1 = mesh.vertices(V1_full);

    profiler().reset();

    const double dhat = 1e-3;
    const double mu = 0.3;
    const double eps_v = 1e-3;

    // CCD step size
    const double step = compute_collision_free_stepsize(
        mesh, V0, V1, 0.0, nullptr, AdditiveCCD());

    // Normal collisions and barrier potential
    NormalCollisions collisions;
    collisions.build(mesh, V0, dhat);

    BarrierPotential bp(dhat, /*stiffness=*/1e4);
    bp(collisions, mesh, V0);
    bp.gradient(collisions, mesh, V0);
    bp.hessian(collisions, mesh, V0);

    // Tangential (friction) collisions and friction potential
    TangentialCollisions tangential;
    tangential.build(mesh, V0, collisions, bp, mu);

    FrictionPotential fp(eps_v);
    fp(tangential, mesh, V0);
    fp.gradient(tangential, mesh, V0);
    fp.hessian(tangential, mesh, V0);

    const std::string output = "profiler_output.csv";
    profiler().write_csv(output);

    logger().info(
        "Profiler output written to: {}\n"
        "  vertices: {}, faces: {}, collisions: {}, step_size: {:.6g}\n",
        output, V0.rows(), F0.rows(), collisions.size(), step);
}

#endif