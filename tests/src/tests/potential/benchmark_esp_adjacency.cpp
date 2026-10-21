// End-to-end ESP timings that track the cost of the candidate adjacency lists
// (Candidates::*_set), which ESP queries for every primitive it integrates
// over. The first test goes through the public ESP API only, so the very same
// benchmark runs against any implementation of those lists.
//
// Each scene benchmarks
//   - ESPCollisions::build  (convert_candidates_to_sets + the builder queries),
//   - ESPPotential value, gradient, and Hessian,
// and prints the results themselves, so two implementations can be checked
// against each other. They agree to roundoff, not bitwise: the sums run over
// hash containers whose iteration order changes from one process to the next.
//
// Run with (single-threaded is the primary measurement):
//   ./ipc_toolkit_tests "[esp_adjacency_benchmark]" --nthreads 1 \
//       --benchmark-samples 7
// A single build takes hundreds of milliseconds on the large scenes, so lower
// the sample count from Catch2's default of 100.
//
// Environment:
//   IPC_TOOLKIT_BENCH_SCENES   comma separated scene labels to run

#include <tests/config.hpp>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/esp/esp_potential.hpp>

#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

using namespace ipc;

namespace {

struct SceneSpec {
    std::string label;
    std::string mesh_name;
    double dhat;
    int quad_rule; // 0 = vertices, 1 = vertices + centroid
};

const std::vector<SceneSpec>& scene_specs()
{
    static const std::vector<SceneSpec> specs = {
        { "sphere", "wrapped_sphere.obj", 0.15, 0 },
        { "sphere-face", "wrapped_sphere.obj", 0.15, 1 },
        { "two-cubes", "two-cubes-close.ply", 0.1, 0 },
        { "bunny", "bunny.ply", 3e-2, 0 },
        { "cloth-ball", "cloth_ball92.ply", 1e-3, 0 },
        { "cloth-funnel", "cloth-funnel/227.ply", 4.0e-3, 0 },
    };
    return specs;
}

bool scene_selected(const std::string& label)
{
    const char* v = std::getenv("IPC_TOOLKIT_BENCH_SCENES");
    if (!v) {
        return true;
    }
    std::stringstream ss(v);
    for (std::string tok; std::getline(ss, tok, ',');) {
        if (tok == label) {
            return true;
        }
    }
    return false;
}

std::vector<SceneSpec> selected_scenes(const bool vertices_only = false)
{
    std::vector<SceneSpec> out;
    for (const SceneSpec& spec : scene_specs()) {
        if (scene_selected(spec.label)
            && (!vertices_only || spec.quad_rule == 0)) {
            out.push_back(spec);
        }
    }
    return out;
}

} // namespace

TEST_CASE(
    "ESP adjacency end-to-end timings", "[esp_adjacency_benchmark][!benchmark]")
{
    const auto scenes = selected_scenes();
    if (scenes.empty()) {
        SKIP("no scene selected");
    }
    const SceneSpec spec = GENERATE_COPY(from_range(scenes));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    const std::string path = (tests::DATA_DIR / spec.mesh_name).string();
    if (!igl::read_triangle_mesh(path, V, F)) {
        SKIP("missing " << path);
    }
    igl::edges(F, E);
    const CollisionMesh mesh(V, E, F);

    ESPParameters params(spec.dhat, 1., 0);
    if (spec.quad_rule == 1) {
        params.set_quad_rule(
            {
                { { { 1., 0., 0. } }, 0.25 },
                { { { 0., 1., 0. } }, 0.25 },
                { { { 0., 0., 1. } }, 0.25 },
                { { { 1. / 3, 1. / 3, 1. / 3 } }, 0.25 },
            });
    }
    const ESPPotential potential(params);

    ESPCollisions collisions;
    double value = 0;
    Eigen::VectorXd grad;
    Eigen::SparseMatrix<double> hess;

    BENCHMARK(fmt::format("{}: build", spec.label))
    {
        collisions.build(mesh, V, params);
        return collisions.size();
    };

    BENCHMARK(fmt::format("{}: value", spec.label))
    {
        value = potential(collisions, mesh, V);
        return value;
    };

    BENCHMARK(fmt::format("{}: gradient", spec.label))
    {
        grad = potential.gradient(collisions, mesh, V);
        return grad.size();
    };

    // The Hessian of a scene with 100k collisions takes minutes and says
    // nothing the gradient does not, so we skip it past this size.
    constexpr size_t MAX_HESSIAN_COLLISIONS = 20000;
    if (collisions.size() <= MAX_HESSIAN_COLLISIONS) {
        BENCHMARK(fmt::format("{}: hessian", spec.label))
        {
            hess = potential.hessian(collisions, mesh, V);
            return hess.nonZeros();
        };
    }

    std::printf(
        "%-12s verts=%zu faces=%zu colls=%zu P=%.17g |g|=%.17g sum(g)=%.17g "
        "|H|=%.17g\n",
        spec.label.c_str(), size_t(V.rows()), size_t(F.rows()),
        collisions.size(), value, grad.norm(), grad.sum(), hess.norm());
    std::fflush(stdout);
}

// The accessors in isolation: build the candidates once, then sweep every
// primitive id through each accessor and walk the result, as the ESP code does.
// This is an upper bound on what a faster container can buy -- it measures the
// queries with nothing else around them -- so the end-to-end numbers above are
// the ones to trust for the speedup.
TEST_CASE(
    "ESP adjacency accessors in isolation",
    "[esp_adjacency_benchmark][!benchmark]")
{
    const auto scenes = selected_scenes(/*vertices_only=*/true);
    if (scenes.empty()) {
        SKIP("no scene selected");
    }
    const SceneSpec spec = GENERATE_COPY(from_range(scenes));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    const std::string path = (tests::DATA_DIR / spec.mesh_name).string();
    if (!igl::read_triangle_mesh(path, V, F)) {
        SKIP("missing " << path);
    }
    igl::edges(F, E);
    const CollisionMesh mesh(V, E, F);

    Candidates candidates;
    candidates.build(mesh, V, spec.dhat, nullptr, /*all_types=*/true);

    const index_t nv = index_t(mesh.num_vertices());
    const index_t ne = index_t(mesh.num_edges());
    const index_t nf = index_t(mesh.num_faces());

    BENCHMARK(fmt::format("{}: convert", spec.label))
    {
        candidates.convert_candidates_to_sets();
        return candidates.size();
    };

    // The accessors read the lists that convert_candidates_to_sets built.
    candidates.convert_candidates_to_sets();

    const auto sweep = [&](const std::string& name, const index_t n,
                           auto&& query) {
        BENCHMARK(fmt::format("{}: {}", spec.label, name))
        {
            size_t sum = 0;
            for (index_t id = 0; id < n; ++id) {
                for (const index_t x : query(id)) {
                    sum += size_t(x) + 1;
                }
            }
            return sum;
        };
    };

    sweep("vv", nv, [&](index_t i) { return candidates.vv_set(mesh, i); });
    sweep("ve", nv, [&](index_t i) { return candidates.ve_set(i); });
    sweep("vf", nv, [&](index_t i) { return candidates.vf_set(i); });
    sweep("ev", ne, [&](index_t i) { return candidates.ev_set(mesh, i); });
    sweep("ee", ne, [&](index_t i) { return candidates.ee_set(mesh, i); });
    sweep("ef", ne, [&](index_t i) { return candidates.ef_set(mesh, i); });
    sweep("fv", nf, [&](index_t i) { return candidates.fv_set(mesh, i); });
    sweep("fe", nf, [&](index_t i) { return candidates.fe_set(mesh, i); });
    sweep("ff", nf, [&](index_t i) { return candidates.ff_set(mesh, i); });
}
