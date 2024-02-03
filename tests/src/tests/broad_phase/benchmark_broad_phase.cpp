#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/broad_phase/brute_force.hpp>

using namespace ipc;

TEST_CASE("Benchmark broad phase", "[!benchmark][broad_phase]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

#ifdef NDEBUG
    std::string filename =
        GENERATE(std::string("cube.obj"), std::string("bunny.obj"));
#else
    std::string filename = "cube.obj";
#endif
    REQUIRE(tests::load_mesh(filename, V0, E, F));

    VectorMax3d min_V = V0.colwise().minCoeff();
    VectorMax3d max_V = V0.colwise().maxCoeff();
    VectorMax3d center = (max_V + min_V) / 2;
    V0.rowwise() -= center.transpose();

    std::string testcase_name;
    SECTION("Squish")
    {
        double s = GENERATE(0.1, 0.0, -0.1);
        V1 = V0;
        V1.col(0) *= s;

        testcase_name = fmt::format("Squish s={}", s);
    }
    SECTION("Interobject")
    {
        VectorMax3d scale = max_V - min_V;
        for (int j = 0; j < V0.cols(); j++) {
            V0.col(j).array() /= scale(j);
        }

        int n = V0.rows();

        Eigen::MatrixXd V(2 * n, V0.cols());
        V.topRows(n) = V0;
        V.block(0, 0, n, 1).array() += 0.5;
        V.bottomRows(n) = V0;
        V.block(n, 0, n, 1).array() -= 0.5;
        V0 = V;

        double dx = GENERATE(0.0, 0.25, 0.75, 1.0);
        V1 = V0;
        V1.block(0, 0, n, 1).array() -= dx;
        V1.block(n, 0, n, 1).array() += dx;

        Eigen::MatrixXi F_stack(2 * F.rows(), F.cols());
        F_stack.topRows(F.rows()) = F;
        F_stack.bottomRows(F.rows()) = F.array() + n;
        F = F_stack;

        testcase_name = fmt::format("Interobject dx={}", dx);
    }

    double inflation_radius = 1e-2; // GENERATE(take(5, random(0.0, 0.1)));

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    const static std::vector<std::string> BP_names = {
        "BF", "HG", "SH", "BVH", "STQ", "GPU_STQ",
    };
    for (int i = 0; i < NUM_BROAD_PHASE_METHODS; i++) {
        BroadPhaseMethod method = static_cast<BroadPhaseMethod>(i);
        BENCHMARK(fmt::format("BP {} ({})", testcase_name, BP_names[i]))
        {
            Candidates candidates;
            candidates.build(mesh, V0, V1, inflation_radius, method);
        };
    }
}

TEST_CASE(
    "Benchmark broad phase on real data",
    "[!benchmark][broad_phase][real_data]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    std::string mesh_name_t0, mesh_name_t1;
    SECTION("Data 0")
    {
        mesh_name_t0 = "private/slow-broadphase-ccd/0.obj";
        mesh_name_t1 = "private/slow-broadphase-ccd/1.obj";
    }
    SECTION("Data 1")
    {
        mesh_name_t0 = "private/slow-broadphase-ccd/s0.obj";
        mesh_name_t1 = "private/slow-broadphase-ccd/s1.obj";
    }
    SECTION("Cloth-Ball")
    {
        mesh_name_t0 = "cloth_ball92.ply";
        mesh_name_t1 = "cloth_ball93.ply";
    }
    // SECTION("Squishy-Ball")
    // {
    //     mesh_name_t0 = "private/puffer-ball/20.ply";
    //     mesh_name_t1 = "private/puffer-ball/21.ply";
    // }

    if (!tests::load_mesh(mesh_name_t0, V0, E, F)
        || !tests::load_mesh(mesh_name_t1, V1, E, F)) {
        return; // Data is private
    }

    double inflation_radius = 1e-2; // GENERATE(take(5, random(0.0, 0.1)));

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    const static std::vector<std::string> BP_names = {
        "BF", "HG", "SH", "BVH", "STQ", "GPU_STQ",
    };
    for (int i = 1; i < NUM_BROAD_PHASE_METHODS; i++) {
        BroadPhaseMethod method = static_cast<BroadPhaseMethod>(i);
        BENCHMARK(fmt::format("BP Real Data ({})", BP_names[i]))
        {
            Candidates candidates;
            candidates.build(mesh, V0, V1, inflation_radius, method);
        };
    }
}
