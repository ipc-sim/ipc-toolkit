// End-to-end ESP timings that track the cost of the candidate adjacency lists
// (Candidates::*_set), which ESP queries for every primitive it integrates
// over. The first test goes through the public ESP API only, so the very same
// benchmark runs against any implementation of those lists.
//
// Each scene reports the median wall time of
//   - ESPCollisions::build  (convert_candidates_to_sets + the builder queries),
//   - ESPPotential value, gradient, and Hessian,
// plus the results themselves, so two implementations can be checked against
// each other. They agree to roundoff, not bitwise: the sums run over hash
// containers whose iteration order changes from one process to the next.
//
// Run with (single-threaded is the primary measurement):
//   ./ipc_toolkit_tests "[esp_adjacency_benchmark]"
//
// Environment:
//   IPC_TOOLKIT_BENCH_SAMPLES  number of timed runs per cell (default 7)
//   IPC_TOOLKIT_BENCH_THREADS  TBB threads (default 1)
//   IPC_TOOLKIT_BENCH_SCENES   comma separated scene labels to run

#include <tests/config.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/esp/esp_potential.hpp>

#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

#include <tbb/global_control.h>

#include <algorithm>
#include <chrono>
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

double now_ms()
{
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

int env_int(const char* name, const int fallback)
{
    const char* v = std::getenv(name);
    return v ? std::max(1, std::atoi(v)) : fallback;
}

template <typename F> double median_ms(const int samples, F&& f)
{
    std::vector<double> t;
    for (int i = 0; i < samples; ++i) {
        const double t0 = now_ms();
        f();
        t.push_back(now_ms() - t0);
    }
    std::sort(t.begin(), t.end());
    return t[t.size() / 2];
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

} // namespace

TEST_CASE(
    "ESP adjacency end-to-end timings", "[esp_adjacency_benchmark][!benchmark]")
{
    const int samples = env_int("IPC_TOOLKIT_BENCH_SAMPLES", 7);
    const int threads = env_int("IPC_TOOLKIT_BENCH_THREADS", 1);
    tbb::global_control limit(
        tbb::global_control::max_allowed_parallelism, size_t(threads));

    std::printf("threads=%d samples=%d\n", threads, samples);
    std::printf(
        "%-12s %8s %8s %9s | %9s %9s %9s %9s | checksum\n", "scene", "verts",
        "faces", "colls", "build", "value", "grad", "hess");

    for (const SceneSpec& spec : scene_specs()) {
        if (!scene_selected(spec.label)) {
            continue;
        }

        Eigen::MatrixXd V;
        Eigen::MatrixXi E, F;
        const std::string path = (tests::DATA_DIR / spec.mesh_name).string();
        if (!igl::read_triangle_mesh(path, V, F)) {
            std::printf("%-12s missing %s\n", spec.label.c_str(), path.c_str());
            continue;
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
        const double t_build =
            median_ms(samples, [&] { collisions.build(mesh, V, params); });

        double value = 0;
        const double t_value =
            median_ms(samples, [&] { value = potential(collisions, mesh, V); });

        Eigen::VectorXd grad;
        const double t_grad = median_ms(
            samples, [&] { grad = potential.gradient(collisions, mesh, V); });

        // The Hessian of a scene with 100k collisions takes minutes and says
        // nothing the gradient does not, so we skip it past this size.
        constexpr size_t MAX_HESSIAN_COLLISIONS = 20000;
        Eigen::SparseMatrix<double> hess;
        double t_hess = -1;
        if (collisions.size() <= MAX_HESSIAN_COLLISIONS) {
            t_hess = median_ms(samples, [&] {
                hess = potential.hessian(collisions, mesh, V);
            });
        }

        std::printf(
            "%-12s %8zu %8zu %9zu | %9.2f %9.2f %9.2f %9.2f | "
            "P=%.17g |g|=%.17g sum(g)=%.17g |H|=%.17g\n",
            spec.label.c_str(), size_t(V.rows()), size_t(F.rows()),
            collisions.size(), t_build, t_value, t_grad, t_hess, value,
            grad.norm(), grad.sum(), hess.norm());
        std::fflush(stdout);
    }
    SUCCEED();
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
    const int samples = env_int("IPC_TOOLKIT_BENCH_SAMPLES", 7);
    const int threads = env_int("IPC_TOOLKIT_BENCH_THREADS", 1);
    tbb::global_control limit(
        tbb::global_control::max_allowed_parallelism, size_t(threads));

    std::printf(
        "%-12s %9s | %8s | %8s %8s %8s | %8s %8s %8s | %8s %8s %8s | "
        "checksum\n",
        "scene", "cands", "convert", "vv", "ve", "vf", "ev", "ee", "ef", "fv",
        "fe", "ff");

    for (const SceneSpec& spec : scene_specs()) {
        if (!scene_selected(spec.label) || spec.quad_rule != 0) {
            continue;
        }

        Eigen::MatrixXd V;
        Eigen::MatrixXi E, F;
        const std::string path = (tests::DATA_DIR / spec.mesh_name).string();
        if (!igl::read_triangle_mesh(path, V, F)) {
            continue;
        }
        igl::edges(F, E);
        const CollisionMesh mesh(V, E, F);

        Candidates candidates;
        candidates.build(mesh, V, spec.dhat, nullptr, /*all_types=*/true);
        const double t_convert = median_ms(
            samples, [&] { candidates.convert_candidates_to_sets(); });

        const index_t nv = index_t(mesh.num_vertices());
        const index_t ne = index_t(mesh.num_edges());
        const index_t nf = index_t(mesh.num_faces());

        size_t checksum = 0;
        const auto sweep = [&](const index_t n, auto&& query) {
            size_t sum = 0;
            const double t = median_ms(samples, [&] {
                sum = 0;
                for (index_t id = 0; id < n; ++id) {
                    for (const index_t x : query(id)) {
                        sum += size_t(x) + 1;
                    }
                }
            });
            checksum = checksum * 1000003 + sum;
            return t;
        };

        const double vv =
            sweep(nv, [&](index_t i) { return candidates.vv_set(mesh, i); });
        const double ve =
            sweep(nv, [&](index_t i) { return candidates.ve_set(i); });
        const double vf =
            sweep(nv, [&](index_t i) { return candidates.vf_set(i); });
        const double ev =
            sweep(ne, [&](index_t i) { return candidates.ev_set(mesh, i); });
        const double ee =
            sweep(ne, [&](index_t i) { return candidates.ee_set(mesh, i); });
        const double ef =
            sweep(ne, [&](index_t i) { return candidates.ef_set(mesh, i); });
        const double fv =
            sweep(nf, [&](index_t i) { return candidates.fv_set(mesh, i); });
        const double fe =
            sweep(nf, [&](index_t i) { return candidates.fe_set(mesh, i); });
        const double ff =
            sweep(nf, [&](index_t i) { return candidates.ff_set(mesh, i); });

        std::printf(
            "%-12s %9zu | %8.3f | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f | "
            "%8.3f %8.3f %8.3f | %zu\n",
            spec.label.c_str(), candidates.size(), t_convert, vv, ve, vf, ev,
            ee, ef, fv, fe, ff, checksum);
        std::fflush(stdout);
    }
    SUCCEED();
}
