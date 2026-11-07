// Sweep benchmark for Potential<T>::gradient's assembly, measuring cost
// against collision count, DOF count, and thread count. This is the sweep the
// crossover between the two assembly strategies was calibrated on; rerun it
// before retuning that condition.
//
// The sweeps use synthetic scenes (synthetic_contact_scene.hpp) because real
// meshes tie the collision count to the DOF count, which is exactly the axis
// that has to be varied independently to see a crossover. A Catch2 BENCHMARK
// case on the real assembly scenes spot-checks the synthetic conclusions.
//
// Run with (RELEASE build only — debug timings are meaningless):
//   ./ipc_toolkit_tests "[grad-assembly-sweep]"          # full sweep + CSV
//   ./ipc_toolkit_tests "[grad-assembly]" --benchmark-samples 10
//
// Environment:
//   IPC_GRAD_SWEEP=main|types|dofmap|all   subset selection (default: all)
//   IPC_GRAD_SWEEP_CSV=<path>              CSV output (default:
//                                          gradient_assembly_sweep.csv)

#include "assembly_scene.hpp"
#include "synthetic_contact_scene.hpp"

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/utils/eigen_ext.hpp>

#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_reduce.h>

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

using namespace ipc;
using namespace ipc::tests;

namespace {

/// @brief How the output DOF relate to the assembled DOF.
enum class DofMapMode : std::uint8_t {
    NONE, ///< in_full_dof = false: assemble and return in collision DOF.
    FOLD, ///< in_full_dof = true, selection map: assemble in full DOF.
    MAP,  ///< in_full_dof = true, displacement map: assemble then to_full_dof.
};

const char* to_string(const DofMapMode mode)
{
    switch (mode) {
    case DofMapMode::NONE:
        return "none";
    case DofMapMode::FOLD:
        return "fold";
    default:
        return "map";
    }
}

double now_ms()
{
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

/// @brief Evaluate every local gradient (DOF gather + gradient + vertex ids +
/// optional fold remap) without assembling anything: the "compute" work that
/// both assembly strategies share. The assemble share of a measurement is
/// estimated as its total minus this.
double local_compute_baseline(const AssemblyScene& scene, const bool fold)
{
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const CollisionMesh& mesh = scene.mesh();
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();
    const Eigen::MatrixXd& X = scene.vertices();

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double partial) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const VectorMax12d local_grad = potential.gradient(
                    collisions[i], collisions[i].dof(X, edges, faces));
                auto ids = collisions[i].vertex_ids(edges, faces);
                if (fold) {
                    for (auto& id : ids) {
                        if (id >= 0) {
                            id = mesh.to_full_vertex_id(id);
                        }
                    }
                }
                partial += local_grad(0) + double(ids[0]);
            }
            return partial;
        },
        std::plus<double>());
}

struct Timing {
    double median_ms = 0;
    double min_ms = 0;
    int samples = 0;
};

/// @brief Median/min wall time of f() with an adaptive sample count: at least
/// MIN_SAMPLES, up to MAX_SAMPLES or until BUDGET_MS is spent.
template <typename F> Timing time_call(F&& f)
{
    constexpr int MIN_SAMPLES = 3;
    constexpr int MAX_SAMPLES = 7;
    constexpr double BUDGET_MS = 1000.0;

    std::vector<double> samples;
    const double t_begin = now_ms();
    while (int(samples.size()) < MAX_SAMPLES
           && (int(samples.size()) < MIN_SAMPLES
               || now_ms() - t_begin < BUDGET_MS)) {
        const double t0 = now_ms();
        f();
        samples.push_back(now_ms() - t0);
    }
    std::sort(samples.begin(), samples.end());

    Timing t;
    t.samples = int(samples.size());
    t.min_ms = samples.front();
    t.median_ms = samples[samples.size() / 2];
    return t;
}

/// @brief Analytic peak extra memory of a branch beyond the output vector.
std::uint64_t peak_extra_bytes(
    const bool two_pass,
    const size_t num_collisions,
    const size_t out_ndof,
    const size_t num_threads)
{
    constexpr size_t STENCIL = 4;
    if (two_pass) {
        return num_collisions
            * (sizeof(VectorMax12d) + STENCIL * sizeof(index_t));
    }
    return std::uint64_t(num_threads) * out_ndof * 8;
}

/// @brief One measured cell of the sweep.
struct SweepRow {
    std::string sweep;
    std::string scene_label;
    bool two_pass;
    SyntheticContactType type;
    int dim;
    size_t num_collisions;
    size_t ndof;      ///< Collision-mesh DOF.
    size_t full_ndof; ///< Full-mesh DOF (= output size when folding/mapping).
    DofMapMode dof_map;
    size_t threads;
    Timing total;
    double compute_ms;  ///< Shared per-thread-count local-compute baseline.
    double assemble_ms; ///< max(total - compute, 0).
    std::uint64_t extra_bytes;
    double checksum;
};

void write_csv(const std::vector<SweepRow>& rows)
{
    const char* path_env = std::getenv("IPC_GRAD_SWEEP_CSV");
    const std::string path =
        (path_env != nullptr) ? path_env : "gradient_assembly_sweep.csv";

    std::ofstream csv(path);
    csv << "sweep,scene,branch,type,dim,num_collisions,ndof,full_ndof,"
           "dof_map,threads,samples,total_ms,total_min_ms,compute_ms,"
           "assemble_ms,extra_bytes,checksum\n";
    for (const SweepRow& r : rows) {
        csv << fmt::format(
            "{},{},{},{},{},{},{},{},{},{},{},{:.6f},{:.6f},{:.6f},{:.6f},"
            "{},{:.17g}\n",
            r.sweep, r.scene_label, r.two_pass ? "two-pass" : "thread-locals",
            ipc::tests::to_string(r.type), r.dim, r.num_collisions, r.ndof,
            r.full_ndof, to_string(r.dof_map), r.threads, r.total.samples,
            r.total.median_ms, r.total.min_ms, r.compute_ms, r.assemble_ms,
            r.extra_bytes, r.checksum);
    }
    fmt::print("\nWrote {} rows to {}\n", rows.size(), path);
}

/// @brief Measure one scene at one thread count, appending a row to `rows`.
void measure_scene(
    const std::string& sweep,
    const AssemblyScene& scene,
    const SyntheticContactSpec& spec,
    const DofMapMode dof_map,
    const size_t num_threads,
    std::vector<SweepRow>& rows)
{
    const bool in_full_dof = dof_map != DofMapMode::NONE;
    const size_t out_ndof =
        dof_map == DofMapMode::FOLD ? scene.full_ndof() : scene.ndof();

    tbb::global_control control(
        tbb::global_control::max_allowed_parallelism, num_threads);

    // Strategy-invariant local-compute baseline (median); the assemble share
    // is estimated as total minus this.
    const Timing compute = time_call([&] {
        (void)local_compute_baseline(scene, dof_map == DofMapMode::FOLD);
    });

    const auto run = [&]() -> Eigen::VectorXd {
        return scene.potential().gradient(
            scene.collisions(), scene.mesh(), scene.vertices(), in_full_dof);
    };

    // Warm up: the first call touches fresh pages in the per-thread buffers,
    // which would otherwise over-report.
    (void)run();
    const Eigen::VectorXd grad = run();

    SweepRow row;
    row.sweep = sweep;
    row.scene_label = scene.label();
    row.two_pass = assembly_is_two_pass(scene.num_collisions(), out_ndof);
    row.type = spec.type;
    row.dim = spec.dim;
    row.num_collisions = scene.num_collisions();
    row.ndof = scene.ndof();
    row.full_ndof = scene.full_ndof();
    row.dof_map = dof_map;
    row.threads = num_threads;
    row.total = time_call([&] { (void)run(); });
    row.compute_ms = compute.median_ms;
    row.assemble_ms = std::max(row.total.median_ms - compute.median_ms, 0.0);
    row.extra_bytes = peak_extra_bytes(
        row.two_pass, scene.num_collisions(), out_ndof, num_threads);
    row.checksum = grad.sum();
    rows.push_back(row);

    fmt::print(
        "{:<8} {:<28} {:>14} thr={:<3} total={:>9.4f}ms assemble={:>9.4f}ms "
        "sum={:.6g}\n",
        sweep, scene.label(), row.two_pass ? "two-pass" : "thread-locals",
        num_threads, row.total.median_ms, row.assemble_ms, row.checksum);
    std::fflush(stdout);
}

std::vector<size_t> sweep_thread_counts()
{
    const size_t hw = std::max(1U, std::thread::hardware_concurrency());
    std::vector<size_t> threads = { 1, 2, 4, 8, hw };
    std::sort(threads.begin(), threads.end());
    threads.erase(std::unique(threads.begin(), threads.end()), threads.end());
    threads.erase(
        std::remove_if(
            threads.begin(), threads.end(),
            [&](const size_t t) { return t > hw; }),
        threads.end());
    return threads;
}

bool sweep_enabled(const std::string& name)
{
    const char* env = std::getenv("IPC_GRAD_SWEEP");
    const std::string selection = (env != nullptr) ? env : "all";
    return selection == "all" || selection == name;
}

} // namespace

TEST_CASE("Gradient assembly sweep", "[.][grad-assembly-sweep]")
{
    const std::vector<size_t> threads = sweep_thread_counts();
    std::vector<SweepRow> rows;

    // --- Main grid: 3D face-vertex, collision DOF ---------------------------
    if (sweep_enabled("main")) {
        const std::vector<size_t> num_collisions = { 10,     100,    1'000,
                                                     5'000,  20'000, 100'000,
                                                     500'000 };
        const std::vector<size_t> ndofs = { 1'000, 10'000, 100'000, 300'000,
                                            1'000'000 };

        for (const size_t nc : num_collisions) {
            for (const size_t ndof : ndofs) {
                SyntheticContactSpec spec;
                spec.num_collisions = nc;
                spec.target_ndof = ndof;
                spec.dim = 3;
                spec.type = SyntheticContactType::FACE_VERTEX;
                spec.seed = 20260807;
                const AssemblyScene scene = build_synthetic_contact_scene(spec);
                for (const size_t t : threads) {
                    measure_scene(
                        "main", scene, spec, DofMapMode::NONE, t, rows);
                }
            }
        }
    }

    // --- Collision-type side sweep: 3D edge-edge and 2D edge-vertex ---------
    if (sweep_enabled("types")) {
        const std::vector<size_t> num_collisions = { 1'000, 20'000, 100'000 };
        const std::vector<size_t> ndofs = { 10'000, 300'000 };
        const std::vector<size_t> type_threads = { 1, threads.back() };

        for (const bool two_d : { false, true }) {
            for (const size_t nc : num_collisions) {
                for (const size_t ndof : ndofs) {
                    SyntheticContactSpec spec;
                    spec.num_collisions = nc;
                    spec.target_ndof = ndof;
                    spec.dim = two_d ? 2 : 3;
                    spec.type = two_d ? SyntheticContactType::EDGE_VERTEX
                                      : SyntheticContactType::EDGE_EDGE;
                    spec.seed = 20260807;
                    const AssemblyScene scene =
                        build_synthetic_contact_scene(spec);
                    for (const size_t t : type_threads) {
                        measure_scene(
                            "types", scene, spec, DofMapMode::NONE, t, rows);
                    }
                }
            }
        }
    }

    // --- DOF-map side sweep: collision DOF vs. folded vs. mapped ------------
    if (sweep_enabled("dofmap")) {
        const std::vector<size_t> num_collisions = { 5'000, 100'000 };

        for (const size_t nc : num_collisions) {
            SyntheticContactSpec spec;
            spec.num_collisions = nc;
            spec.target_ndof = 300'000;
            spec.dim = 3;
            spec.type = SyntheticContactType::FACE_VERTEX;
            spec.seed = 20260807;

            const AssemblyScene selection_scene =
                build_synthetic_contact_scene(spec);
            measure_scene(
                "dofmap", selection_scene, spec, DofMapMode::NONE,
                threads.back(), rows);
            measure_scene(
                "dofmap", selection_scene, spec, DofMapMode::FOLD,
                threads.back(), rows);

            spec.with_displacement_map = true;
            const AssemblyScene mapped_scene =
                build_synthetic_contact_scene(spec);
            measure_scene(
                "dofmap", mapped_scene, spec, DofMapMode::MAP, threads.back(),
                rows);
        }
    }

    write_csv(rows);
}

TEST_CASE(
    "Benchmark gradient assembly (real scenes)", "[!benchmark][grad-assembly]")
{
    // Spot-check the synthetic-sweep conclusions on realistic active sets and
    // cache state. Scene selection mirrors "Benchmark contact gradient
    // assembly" in benchmark_assembly.cpp.
    const auto spec = GENERATE(from_range(assembly_scene_specs()));

    const std::optional<AssemblyScene> maybe_scene = build_assembly_scene(spec);
    if (!maybe_scene.has_value()) {
        SKIP(
            fmt::format(
                "Scene '{}' is unavailable (missing mesh or no collisions).",
                spec.mesh_name));
    }
    const AssemblyScene& scene = maybe_scene.value();

    fmt::print(
        "\n{} | assembly: {}\n", scene.stats(),
        assembly_is_two_pass(scene.num_collisions(), scene.ndof())
            ? "two-pass"
            : "thread-locals");

    BENCHMARK(fmt::format("{}: gradient", scene.label()))
    {
        return scene.potential().gradient(
            scene.collisions(), scene.mesh(), scene.vertices());
    };

    BENCHMARK(fmt::format("{}: gradient (full DOF)", scene.label()))
    {
        return scene.potential().gradient(
            scene.collisions(), scene.mesh(), scene.vertices(),
            /*in_full_dof=*/true);
    };
}
