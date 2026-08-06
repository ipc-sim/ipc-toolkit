// Baseline measurements for the cost of assembling contact gradients and
// Hessians. These exist to quantify, separately:
//
//   1. per-collision (local) derivative evaluation,
//   2. global assembly (triplets + setFromTriplets),
//   3. the reduced-DOF map (`CollisionMesh::to_full_dof`).
//
// (2) is not measured directly: it is the difference between the full call and
// (1), because the two are interleaved inside a single `tbb::parallel_for`. The
// "Assembly cost breakdown" test computes that difference and prints a table.
//
// Run with:
//   ./ipc_toolkit_tests "[assembly]" --benchmark-samples 20
//
// The scenes are shared with the correctness tests via `assembly_scene.hpp`.

#include "assembly_scene.hpp"

#include <tests/utils.hpp>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/hessian_assembler.hpp>
#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE
#include <ipc/utils/meshfem_hessian_assembler.hpp>
#endif

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>

using namespace ipc;

namespace {

/// @brief Number of assemblies per "solve" in the amortized benchmark.
///
/// A Newton solve performs one Hessian assembly per iteration while the
/// collision set is unchanged, which is exactly the situation a persistent
/// sparsity pattern exploits. Assembling repeatedly in a single benchmark
/// sample keeps that opportunity visible in the baseline numbers.
constexpr int ASSEMBLIES_PER_SOLVE = 10;

/// @brief Evaluate every local Hessian without assembling anything.
///
/// Mirrors the per-collision work inside `Potential::hessian` exactly: the DOF
/// gather plus the local Hessian. Returns a value derived from every result so
/// the computation cannot be optimized away.
double local_hessians_only(
    const ipc::tests::AssemblyScene& scene,
    const PSDProjectionMethod psd = PSDProjectionMethod::NONE)
{
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const Eigen::MatrixXi& edges = scene.mesh().edges();
    const Eigen::MatrixXi& faces = scene.mesh().faces();
    const Eigen::MatrixXd& X = scene.vertices();

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double partial) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const MatrixMax12d local_hess = potential.hessian(
                    collisions[i], collisions[i].dof(X, edges, faces), psd);
                partial += local_hess(0, 0);
            }
            return partial;
        },
        std::plus<double>());
}

/// @brief Evaluate every local gradient without assembling anything.
double local_gradients_only(const ipc::tests::AssemblyScene& scene)
{
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const Eigen::MatrixXi& edges = scene.mesh().edges();
    const Eigen::MatrixXi& faces = scene.mesh().faces();
    const Eigen::MatrixXd& X = scene.vertices();

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double partial) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const VectorMax12d local_grad = potential.gradient(
                    collisions[i], collisions[i].dof(X, edges, faces));
                partial += local_grad(0);
            }
            return partial;
        },
        std::plus<double>());
}

/// @brief Median wall-clock time of `f` over `num_samples` runs, in seconds.
///
/// Used by the breakdown table. Catch2's `BENCHMARK` reports richer statistics
/// but does not expose its measurements programmatically, so the table below
/// does its own (much simpler) timing. Treat the `BENCHMARK` output as
/// authoritative and the table as a summary of the same ratios.
template <typename F> double median_seconds(F&& f, const int num_samples = 5)
{
    std::vector<double> samples;
    samples.reserve(num_samples);
    for (int i = 0; i < num_samples; i++) {
        const auto start = std::chrono::steady_clock::now();
        f();
        const auto end = std::chrono::steady_clock::now();
        samples.push_back(std::chrono::duration<double>(end - start).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

} // namespace

TEST_CASE("Assembly scene statistics", "[!benchmark][assembly]")
{
    fmt::print("\n=== Assembly benchmark scenes ===\n");
    for (const auto& spec : ipc::tests::assembly_scene_specs()) {
        const std::optional<ipc::tests::AssemblyScene> scene =
            ipc::tests::build_assembly_scene(spec);
        if (!scene.has_value()) {
            fmt::print(
                "  {} (dhat={}): UNAVAILABLE (missing mesh or no collisions)\n",
                spec.mesh_name, spec.dhat);
        } else {
            fmt::print("  {}\n", scene->stats());
        }
        // Flush eagerly: stdout is fully buffered when redirected, and an
        // oversized scene can exhaust memory before the buffer is drained.
        std::fflush(stdout);
    }
    fmt::print("\n");
}

// Safely probe a candidate scene for inclusion in `assembly_scene_specs()`.
//
// Building a collision set with a too-large dhat can exhaust host memory, so
// this test (a) sizes dhat relative to the mesh's bounding-box diagonal and
// (b) counts broad-phase candidates first, refusing to build the collision set
// if there are too many. Run one scene per process:
//
//   IPC_ASSEMBLY_PROBE_MESH=puffer-ball/20.ply \
//   IPC_ASSEMBLY_PROBE_DHAT_REL=1e-3 \
//   ./ipc_toolkit_tests "[assembly-probe]"
TEST_CASE("Assembly scene probe", "[.][assembly-probe]")
{
    const char* mesh_name = std::getenv("IPC_ASSEMBLY_PROBE_MESH");
    if (mesh_name == nullptr) {
        SKIP("Set IPC_ASSEMBLY_PROBE_MESH to probe a scene.");
    }
    const char* dhat_rel_str = std::getenv("IPC_ASSEMBLY_PROBE_DHAT_REL");
    const double dhat_rel =
        (dhat_rel_str != nullptr) ? std::atof(dhat_rel_str) : 1e-3;

    // Above this many broad-phase candidates, do not attempt to build the
    // collision set: candidate/collision storage grows superlinearly with dhat
    // and has exhausted memory on a 64 GB host before.
    constexpr size_t MAX_SAFE_CANDIDATES = 10'000'000;

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(ipc::tests::load_mesh(mesh_name, vertices, edges, faces));

    const double bbox_diag =
        (vertices.colwise().maxCoeff() - vertices.colwise().minCoeff()).norm();
    const double dhat = dhat_rel * bbox_diag;

    const CollisionMesh mesh =
        CollisionMesh::build_from_full_mesh(vertices, edges, faces);
    vertices = mesh.vertices(vertices);

    fmt::print(
        "probe {}: V={} E={} F={} bbox_diag={:g} dhat={:g} (rel={:g})\n",
        mesh_name, vertices.rows(), edges.rows(), faces.rows(), bbox_diag, dhat,
        dhat_rel);
    std::fflush(stdout);

    Candidates candidates;
    candidates.build(mesh, vertices, vertices, /*inflation_radius=*/dhat / 2);
    fmt::print("probe {}: candidates={}\n", mesh_name, candidates.size());
    std::fflush(stdout);

    if (candidates.size() > MAX_SAFE_CANDIDATES) {
        fmt::print(
            "probe {}: REFUSING to build collisions (> {} candidates)\n",
            mesh_name, MAX_SAFE_CANDIDATES);
        return;
    }

    NormalCollisions collisions;
    collisions.build(candidates, mesh, vertices, dhat);
    fmt::print("probe {}: collisions={}\n", mesh_name, collisions.size());
    std::fflush(stdout);
}

TEST_CASE("Benchmark contact Hessian assembly", "[!benchmark][assembly]")
{
    const auto spec = GENERATE(from_range(ipc::tests::assembly_scene_specs()));

    const std::optional<ipc::tests::AssemblyScene> maybe_scene =
        ipc::tests::build_assembly_scene(spec);
    if (!maybe_scene.has_value()) {
        SKIP(
            fmt::format(
                "Scene '{}' is unavailable (missing mesh or no collisions).",
                spec.mesh_name));
    }
    const ipc::tests::AssemblyScene& scene = maybe_scene.value();

    const CollisionMesh& mesh = scene.mesh();
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const Eigen::MatrixXd& X = scene.vertices();

    fmt::print("\n{}\n", scene.stats());

    // Precomputed so `to_full_dof` can be timed in isolation.
    const Eigen::SparseMatrix<double> hess =
        potential.hessian(collisions, mesh, X);

    BENCHMARK(fmt::format("{}: local hessians", scene.label()))
    {
        return local_hessians_only(scene);
    };

    BENCHMARK(fmt::format("{}: hessian (collision DOF)", scene.label()))
    {
        return potential.hessian(collisions, mesh, X);
    };

    BENCHMARK(fmt::format("{}: to_full_dof (hessian)", scene.label()))
    {
        return mesh.to_full_dof(hess);
    };

    BENCHMARK(fmt::format("{}: hessian + to_full_dof", scene.label()))
    {
        return mesh.to_full_dof(potential.hessian(collisions, mesh, X));
    };

    // Phase 1: assemble directly in full DOF (folds to_full_dof into the
    // triplet remap).
    BENCHMARK(fmt::format("{}: hessian (full DOF, folded)", scene.label()))
    {
        return potential.hessian(
            collisions, mesh, X, PSDProjectionMethod::NONE,
            /*in_full_dof=*/true);
    };

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE
    // Phase 3: MeshFEMSparse block-CSC backend, cold (pattern built per
    // call, e.g., a fresh assembler every iteration).
    BENCHMARK(fmt::format("{}: hessian (MeshFEM, cold)", scene.label()))
    {
        MeshFEMHessianAssembler assembler;
        potential.assemble_hessian(
            collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
            /*in_full_dof=*/true);
        return &assembler.get_matrix();
    };

    // Without the Eigen conversion: what a caller that consumes the block-CSC
    // format directly (e.g., a block-aware solver) would pay.
    BENCHMARK(fmt::format("{}: hessian (MeshFEM, cold, block)", scene.label()))
    {
        // Pattern + scattered values only. The scatter writes to heap
        // storage through virtual dispatch, so it cannot be optimized away.
        MeshFEMHessianAssembler assembler;
        potential.assemble_hessian(
            collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
            /*in_full_dof=*/true);
    };

    // Phase 4: persistent assembler — the pattern (and the cached Eigen
    // structure) are reused across assemblies, as in a Newton solve.
    {
        MeshFEMHessianAssembler assembler;
        potential.assemble_hessian(
            collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
            /*in_full_dof=*/true); // warm-up: builds the pattern

        BENCHMARK(fmt::format("{}: hessian (MeshFEM, reused)", scene.label()))
        {
            potential.assemble_hessian(
                collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
                /*in_full_dof=*/true);
            return &assembler.get_matrix();
        };

        BENCHMARK(
            fmt::format(
                "{}: {}x (MeshFEM, reused)", scene.label(),
                ASSEMBLIES_PER_SOLVE))
        {
            double checksum = 0;
            for (int i = 0; i < ASSEMBLIES_PER_SOLVE; i++) {
                potential.assemble_hessian(
                    collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
                    /*in_full_dof=*/true);
                checksum += assembler.get_matrix().coeff(0, 0);
            }
            return checksum;
        };

        // Caller-asserted unchanged stencils: change detection skipped too.
        assembler.set_assume_unchanged_stencils(true);
        BENCHMARK(fmt::format("{}: hessian (MeshFEM, assumed)", scene.label()))
        {
            potential.assemble_hessian(
                collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
                /*in_full_dof=*/true);
            return &assembler.get_matrix();
        };
    }
#endif

    BENCHMARK(
        fmt::format(
            "{}: {}x (hessian + to_full_dof)", scene.label(),
            ASSEMBLIES_PER_SOLVE))
    {
        double checksum = 0;
        for (int i = 0; i < ASSEMBLIES_PER_SOLVE; i++) {
            checksum += mesh.to_full_dof(potential.hessian(collisions, mesh, X))
                            .coeff(0, 0);
        }
        return checksum;
    };
}

TEST_CASE("Benchmark contact gradient assembly", "[!benchmark][assembly]")
{
    const auto spec = GENERATE(from_range(ipc::tests::assembly_scene_specs()));

    const std::optional<ipc::tests::AssemblyScene> maybe_scene =
        ipc::tests::build_assembly_scene(spec);
    if (!maybe_scene.has_value()) {
        SKIP(
            fmt::format(
                "Scene '{}' is unavailable (missing mesh or no collisions).",
                spec.mesh_name));
    }
    const ipc::tests::AssemblyScene& scene = maybe_scene.value();

    const CollisionMesh& mesh = scene.mesh();
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const Eigen::MatrixXd& X = scene.vertices();

    const Eigen::VectorXd grad = potential.gradient(collisions, mesh, X);

    BENCHMARK(fmt::format("{}: local gradients", scene.label()))
    {
        return local_gradients_only(scene);
    };

    BENCHMARK(fmt::format("{}: gradient (collision DOF)", scene.label()))
    {
        return potential.gradient(collisions, mesh, X);
    };

    BENCHMARK(fmt::format("{}: to_full_dof (gradient)", scene.label()))
    {
        return mesh.to_full_dof(grad);
    };

    // Phase 1: assemble directly in full DOF.
    BENCHMARK(fmt::format("{}: gradient (full DOF, folded)", scene.label()))
    {
        return potential.gradient(collisions, mesh, X, /*in_full_dof=*/true);
    };
}

TEST_CASE("Assembly cost breakdown", "[!benchmark][assembly]")
{
    // Emits the table that goes into the performance report. Every row is
    // measured on the same process/thread configuration so the shares are
    // directly comparable.
    fmt::print("\n=== Hessian assembly cost breakdown ===\n");
    fmt::print(
        "{:<20} {:>9} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} "
        "{:>10} {:>10} {:>8} {:>8} {:>8} {:>8}\n",
        "scene", "#collis", "local(ms)", "triplet(ms)", "asm(ms)", "full(ms)",
        "fold(ms)", "mfem(ms)", "mfblk(ms)", "mfr(ms)", "mfa(ms)", "local%",
        "asm%", "full%", "speedup");

    for (const auto& spec : ipc::tests::assembly_scene_specs()) {
        const std::optional<ipc::tests::AssemblyScene> maybe_scene =
            ipc::tests::build_assembly_scene(spec);
        if (!maybe_scene.has_value()) {
            continue;
        }
        const ipc::tests::AssemblyScene& scene = maybe_scene.value();

        const CollisionMesh& mesh = scene.mesh();
        const NormalCollisions& collisions = scene.collisions();
        const BarrierPotential& potential = scene.potential();
        const Eigen::MatrixXd& X = scene.vertices();

        const Eigen::SparseMatrix<double> hess =
            potential.hessian(collisions, mesh, X);

        const double t_local =
            median_seconds([&] { (void)local_hessians_only(scene); });
        // The triplet baseline, measured explicitly: hessian() routes through
        // whichever backend is compiled in, so it is no longer the baseline.
        const double t_total = median_seconds([&] {
            TripletHessianAssembler assembler;
            potential.assemble_hessian(collisions, mesh, X, assembler);
            (void)assembler.get_matrix();
        });
        const double t_full =
            median_seconds([&] { (void)mesh.to_full_dof(hess); });
        // Phase 1: fold to_full_dof into assembly.
        const double t_folded = median_seconds([&] {
            (void)potential.hessian(
                collisions, mesh, X, PSDProjectionMethod::NONE,
                /*in_full_dof=*/true);
        });

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE
        // Phase 3: MeshFEM block-CSC backend (pattern rebuilt per call).
        const double t_meshfem = median_seconds([&] {
            MeshFEMHessianAssembler assembler;
            potential.assemble_hessian(
                collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
                /*in_full_dof=*/true);
            (void)assembler.get_matrix();
        });
        // Same, but skipping the block-CSC -> Eigen conversion.
        const double t_meshfem_blk = median_seconds([&] {
            MeshFEMHessianAssembler assembler;
            potential.assemble_hessian(
                collisions, mesh, X, assembler, PSDProjectionMethod::NONE,
                /*in_full_dof=*/true);
        });
        // Phase 4: persistent assembler (pattern + Eigen structure reused).
        MeshFEMHessianAssembler persistent_assembler;
        potential.assemble_hessian(
            collisions, mesh, X, persistent_assembler,
            PSDProjectionMethod::NONE, /*in_full_dof=*/true); // warm-up
        const double t_meshfem_reused = median_seconds([&] {
            potential.assemble_hessian(
                collisions, mesh, X, persistent_assembler,
                PSDProjectionMethod::NONE, /*in_full_dof=*/true);
            (void)persistent_assembler.get_matrix();
        });
        // Same, with caller-asserted unchanged stencils (no detection).
        persistent_assembler.set_assume_unchanged_stencils(true);
        const double t_meshfem_assumed = median_seconds([&] {
            potential.assemble_hessian(
                collisions, mesh, X, persistent_assembler,
                PSDProjectionMethod::NONE, /*in_full_dof=*/true);
            (void)persistent_assembler.get_matrix();
        });
#else
        const double t_meshfem = std::numeric_limits<double>::quiet_NaN();
        const double t_meshfem_blk = std::numeric_limits<double>::quiet_NaN();
        const double t_meshfem_reused =
            std::numeric_limits<double>::quiet_NaN();
        const double t_meshfem_assumed =
            std::numeric_limits<double>::quiet_NaN();
#endif

        // Assembly is what the full call does beyond the local derivatives.
        const double t_asm = t_total - t_local;
        const double t_end_to_end = t_total + t_full;

        constexpr double MS = 1e3;
        fmt::print(
            "{:<20} {:>9} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f} "
            "{:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f} {:>7.1f}% {:>7.1f}% "
            "{:>7.1f}% {:>7.2f}x\n",
            scene.label(), scene.num_collisions(), t_local * MS, t_total * MS,
            t_asm * MS, t_full * MS, t_folded * MS, t_meshfem * MS,
            t_meshfem_blk * MS, t_meshfem_reused * MS, t_meshfem_assumed * MS,
            100.0 * t_local / t_end_to_end, 100.0 * t_asm / t_end_to_end,
            100.0 * t_full / t_end_to_end, t_end_to_end / t_folded);
        std::fflush(stdout);
    }
    fmt::print("\n");
}
