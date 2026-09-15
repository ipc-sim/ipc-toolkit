// Per-stage breakdown of ipc::LBVH and ipc::cuda::LBVH, for the stacked-bar
// benchmark figure. Both broad phases are instrumented with
// IPC_TOOLKIT_PROFILE_BLOCK, so this only has to drive them and hand the
// profiler tree out; the stage names come from the library, not from here.
//
// Device stages are timed with CUDA events (see
// ipc/utils/cuda/cuda_profiler.cuh), which waits on each stage's stop event.
// That serializes stages the uninstrumented build overlaps, so the totals
// here run larger than the fused build and detect times. Use this for stage
// *shares*; take absolute times from a build without the profiler.
//
// Run with:
//   IPC_TOOLKIT_BENCH_OUTPUT=stages.json ./ipc_toolkit_tests "[lbvh_stages]"
//
// Environment:
//   IPC_TOOLKIT_BENCH_SAMPLES  timed calls per scene and phase (default 10)
//   IPC_TOOLKIT_BENCH_OUTPUT   write the results as JSON to this path

#include <ipc/config.hpp>

#if defined(IPC_TOOLKIT_WITH_CUDA) && defined(IPC_TOOLKIT_WITH_PROFILER)

#include <tests/gpu_utils.hpp>
#include <tests/utils.hpp>

#include <ipc/broad_phase/cuda/lbvh.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/utils/profiler.hpp>

#include <catch2/catch_test_macros.hpp>

#include <nlohmann/json.hpp>

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace ipc;

namespace {

struct Scene {
    std::string name, mesh_t0, mesh_t1;
};

const std::vector<Scene> SCENES = {
    { "Cloth-Funnel", "cloth-funnel/227.ply", "cloth-funnel/228.ply" },
    { "Armadillo-Rollers", "armadillo-rollers/326.ply",
      "armadillo-rollers/327.ply" },
    { "Rod-Twist", "rod-twist/3036.ply", "rod-twist/3037.ply" },
    { "Cloth-Ball", "cloth_ball92.ply", "cloth_ball93.ply" },
    { "N-Body-Simulation", "n-body-simulation/balls16_18.ply",
      "n-body-simulation/balls16_19.ply" },
    { "Puffer-Ball", "puffer-ball/20.ply", "puffer-ball/21.ply" },
};

int env_int(const char* name, const int fallback)
{
    const char* v = std::getenv(name);
    return (v != nullptr && *v != '\0') ? std::atoi(v) : fallback;
}

/// @brief Run `f` `samples` times with the profiler cleared first, and return
/// the resulting scope tree. Each scope carries its own accumulated time_ms
/// and count, so the consumer divides rather than this dividing for it.
template <typename F> nlohmann::json profile_calls(const int samples, F&& f)
{
    f(); // warm up: first-touch allocation, learned buffer capacities
    profiler().clear();
    for (int i = 0; i < samples; ++i) {
        f();
    }
    return profiler().data();
}

} // namespace

TEST_CASE("Benchmark LBVH stages", "[!benchmark][broad_phase][lbvh_stages]")
{
    tests::skip_if_no_cuda_device();

    constexpr double inflation_radius = 0;
    const int samples = env_int("IPC_TOOLKIT_BENCH_SAMPLES", 10);

    nlohmann::json report;
    report["samples"] = samples;
    report["scenes"] = nlohmann::json::array();

    for (const auto& [name, mesh_t0, mesh_t1] : SCENES) {
        Eigen::MatrixXd vertices_t0, vertices_t1;
        Eigen::MatrixXi edges, faces;
        REQUIRE(tests::load_mesh(mesh_t0, vertices_t0, edges, faces));
        REQUIRE(tests::load_mesh(mesh_t1, vertices_t1, edges, faces));

        nlohmann::json scene;
        scene["scene"] = name;
        scene["num_faces"] = faces.rows();
        scene["num_edges"] = edges.rows();
        scene["num_vertices"] = vertices_t0.rows();

        ipc::LBVH cpu;
        ipc::cuda::LBVH gpu;

        scene["cpu_build"] = profile_calls(samples, [&]() {
            cpu.build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
        });
        scene["gpu_build"] = profile_calls(samples, [&]() {
            gpu.build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
        });

        // Both trees are now built; time detection against them. The output
        // vector is fresh per call, as in the Catch2 detect benchmarks, so
        // the candidate storage is allocated here and not reused across
        // calls -- on Cloth-Ball that is 42 MB per call.
        size_t n_cpu = 0, n_gpu = 0;
        scene["cpu_detect"] = profile_calls(samples, [&]() {
            std::vector<EdgeEdgeCandidate> candidates;
            cpu.detect_edge_edge_candidates(candidates);
            n_cpu = candidates.size();
        });
        scene["gpu_detect"] = profile_calls(samples, [&]() {
            std::vector<EdgeEdgeCandidate> candidates;
            gpu.detect_edge_edge_candidates(candidates);
            n_gpu = candidates.size();
        });
        // The two broad phases must agree, or the stage split is comparing
        // different amounts of work.
        REQUIRE(n_gpu == n_cpu);
        scene["num_ee_candidates"] = n_cpu;

        report["scenes"].push_back(scene);
    }

    profiler().clear();

    const char* out = std::getenv("IPC_TOOLKIT_BENCH_OUTPUT");
    if (out != nullptr && *out != '\0') {
        std::ofstream file(out);
        REQUIRE(file.is_open());
        file << report.dump(2) << std::endl;
    } else {
        WARN("IPC_TOOLKIT_BENCH_OUTPUT unset; stage results not written");
    }
}

#endif
