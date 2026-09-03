#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/esp/esp_potential.hpp>
#include <ipc/esp/adaptive_support.hpp>
#include <ipc/distance/line_line.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>
#include <igl/readCSV.h>
#include <ipc/ipc.hpp>

#include "igl/read_triangle_mesh.h"
#include "igl/write_triangle_mesh.h"
#include "ipc/distance/edge_edge.hpp"

#include "ipc/esp/quadrature_potential.hpp"
#include <chrono>

#include <cmath>

using namespace ipc;

namespace {

struct TriMeshData {
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    CollisionMesh mesh;
};

TriMeshData load_triangle_mesh(const std::string& path)
{
    TriMeshData data;
    igl::read_triangle_mesh(path, data.V, data.F);
    igl::edges(data.F, data.E);
    data.mesh = CollisionMesh(data.V, data.E, data.F);
    return data;
}

TriMeshData load_wrapped_sphere()
{
    return load_triangle_mesh(
        (tests::DATA_DIR / "../src/tests/potential/wrapped_sphere.obj")
            .string());
}

CollisionMesh
make_2d_collision_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E)
{
    Eigen::MatrixXi F;
    return CollisionMesh(
        std::vector<bool>(V.rows(), true), std::vector<bool>(V.rows(), false),
        V, E, F);
}

inline std::shared_ptr<Barrier> make_inverse_quadratic_barrier()
{
    return std::make_shared<InversePowerBarrier>(2.0);
}
inline std::shared_ptr<Barrier> make_linear_inverse_barrier()
{
    return std::make_shared<InversePowerBarrier>(1.0);
}

struct EeLimitSweepStats {
    double max_abs_P = 0;
    double max_abs_g = 0;
    double max_dP = 0;         // max |P(eps_{i+1}) - P(eps_i)|
    double max_dg = 0;         // max ||g|(eps_{i+1}) - |g|(eps_i)|
    double max_fd_slope_P = 0; // max |dP / d eps|
    double max_fd_slope_g = 0; // max |d|g| / d eps|
    double max_shift_P = 0;    // max |P(eps) - P(eps_min)|
    double max_shift_g = 0;    // max ||g|(eps) - |g|(eps_min)|
    bool all_finite = true;
};

// Build the EA_EB tet-tet geometry used by the three edge-edge limit tests,
// parametrised by the horizontal offset epsilon.
inline void build_ee_limit_geometry(
    const double epsilon,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
    V.resize(8, 3);
    V << 0, 0, 0, 1, 0, 0, 0.5, -0.5, 1, 0.5, 0.5, 1, epsilon, 0.5, -0.01,
        epsilon, -0.5, -0.01, epsilon + 0.5, 0, -1.01, epsilon - 0.5, 0, -1.01;
    F.resize(8, 3);
    F << 0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3, 4, 5, 6, 4, 6, 7, 4, 5, 7, 5, 6, 7;
    E.resize(12, 2);
    E << 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, 4, 5, 4, 6, 4, 7, 5, 6, 5, 7, 6, 7;
}

// Sweep the EA_EB limit over logspaced epsilons in (eps_min, eps_max], compute
// the potential and gradient norm at each sample, and return max absolute
// finite-difference between consecutive samples (both raw ΔP, Δ|g|, and the
// slope ΔP/Δeps, Δ|g|/Δeps).
inline EeLimitSweepStats ee_limit_fd_sweep(
    std::shared_ptr<Barrier> barrier,
    int n_samples = 25,
    double eps_min = 1e-16,
    double eps_max = 1e-5)
{
    EeLimitSweepStats stats;
    std::vector<double> eps_vec, P_vec, g_vec;
    eps_vec.reserve(n_samples);
    P_vec.reserve(n_samples);
    g_vec.reserve(n_samples);

    const double log_lo = std::log10(eps_min);
    const double log_hi = std::log10(eps_max);

    for (int i = 0; i < n_samples; ++i) {
        const double t = double(i) / double(n_samples - 1);
        const double eps = std::pow(10.0, log_lo + t * (log_hi - log_lo));

        Eigen::MatrixXd V;
        Eigen::MatrixXi E, F;
        build_ee_limit_geometry(eps, V, E, F);

        CollisionMesh mesh(V, E, F);
        const double dhat = 0.1;
        EspParameters params(dhat, 1., 0);
        params.barrier = barrier;

        EspCollisions collisions;
        collisions.build(mesh, V, params);
        EspPotential potential(params);

        const double x = potential(collisions, mesh, V);
        const double gn = potential.gradient(collisions, mesh, V).norm();

        if (!std::isfinite(x) || !std::isfinite(gn))
            stats.all_finite = false;

        stats.max_abs_P = std::max(stats.max_abs_P, std::abs(x));
        stats.max_abs_g = std::max(stats.max_abs_g, std::abs(gn));

        eps_vec.push_back(eps);
        P_vec.push_back(x);
        g_vec.push_back(gn);
    }

    for (size_t i = 1; i < eps_vec.size(); ++i) {
        const double dP = std::abs(P_vec[i] - P_vec[i - 1]);
        const double dg = std::abs(g_vec[i] - g_vec[i - 1]);
        const double deps = std::abs(eps_vec[i] - eps_vec[i - 1]);
        stats.max_dP = std::max(stats.max_dP, dP);
        stats.max_dg = std::max(stats.max_dg, dg);
        if (deps > 0) {
            stats.max_fd_slope_P = std::max(stats.max_fd_slope_P, dP / deps);
            stats.max_fd_slope_g = std::max(stats.max_fd_slope_g, dg / deps);
        }
    }

    // Shift relative to the sample at the smallest eps (first sample).
    if (!eps_vec.empty()) {
        const double P0 = P_vec.front();
        const double g0 = g_vec.front();
        for (size_t i = 0; i < eps_vec.size(); ++i) {
            stats.max_shift_P =
                std::max(stats.max_shift_P, std::abs(P_vec[i] - P0));
            stats.max_shift_g =
                std::max(stats.max_shift_g, std::abs(g_vec[i] - g0));
        }
    }
    return stats;
}

} // anonymous namespace

// When the edge-edge closest point approaches the end points of the edge, the
// potential should converge to a finite number
TEST_CASE(
    "Convergent Quadrature Edge Edge Limit",
    "[esp_potential], [esp_potential_3d]")
{
    auto stats =
        ee_limit_fd_sweep(std::make_shared<NormalizedClampedLogBarrier<>>());
    CHECK(stats.all_finite);
    REQUIRE(stats.max_abs_P < 2);
    REQUIRE(stats.max_abs_g < 200);
    // 10x current measured shift from the smallest-eps sample
    // (current: max|ΔP|≈2.4e-5, max|Δ|g||≈1.64).
    CHECK(stats.max_shift_P < 2.4e-4);
    CHECK(stats.max_shift_g < 16.5);
}

// Same configuration as above, but uses an inverse-quadratic barrier to probe
// whether the ESP potential stays finite under a stronger barrier.
TEST_CASE(
    "Convergent Quadrature Edge Edge Limit (Inverse Quadratic Barrier)",
    "[esp_potential], [esp_potential_3d]")
{
    auto stats = ee_limit_fd_sweep(make_inverse_quadratic_barrier());
    CHECK(stats.all_finite);
    CHECK(stats.max_abs_P < 1e8);
    CHECK(stats.max_abs_g < 1e10);
    // 10x current measured shift from the smallest-eps sample
    // (current: max|ΔP|≈1.4e-3, max|Δ|g||≈1.55).
    CHECK(stats.max_shift_P < 1.4e-2);
    CHECK(stats.max_shift_g < 15.5);
}

// Same configuration but with a linear-inverse barrier (1/d divergence).
TEST_CASE(
    "Convergent Quadrature Edge Edge Limit (Linear Inverse Barrier)",
    "[esp_potential], [esp_potential_3d]")
{
    auto stats = ee_limit_fd_sweep(make_linear_inverse_barrier());
    CHECK(stats.all_finite);
    CHECK(stats.max_abs_P < 1e8);
    CHECK(stats.max_abs_g < 1e10);
    // 10x current measured shift from the smallest-eps sample
    // (current: max|ΔP|≈1.85e-5, max|Δ|g||≈0.044).
    CHECK(stats.max_shift_P < 1.85e-4);
    CHECK(stats.max_shift_g < 0.44);
}

TEST_CASE(
    "Convergent Quadrature Gradient and Hessian",
    "[esp_potential], [esp_potential_3d]")
{
    auto [V, E, F, mesh] = load_wrapped_sphere();

    const double dbar_factor = GENERATE(1.0, 0.7, 0.3);
    // Keep dbar = dhat * dbar_factor ≈ 0.15 so the active contact set is
    // comparable across dbar_factor values.
    const double dhat = 0.15 / dbar_factor;
    CAPTURE(dbar_factor);
    EspParameters params(dhat, dbar_factor, 0);

    const bool use_near_far = GENERATE(true, false);
    const bool use_adaptive = GENERATE(true, false);
    CAPTURE(use_near_far, use_adaptive, dbar_factor);
    EspPotential potential(params, use_near_far);

    // Compute adaptive support once so every FD step uses identical dhat
    // values.
    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;

    EspCollisions collisions;
    collisions.build(mesh, V, params, adaptive.get());
    REQUIRE(!collisions.empty());

    // full finite difference is too expensive, verify directional derivative
    // only
    Eigen::VectorXd test_dir(V.size(), 1);
    for (int i = 0; i < test_dir.size(); i++) {
        test_dir(i) = i;
    }
    test_dir.normalize();

    SECTION("gradient")
    {
        Eigen::VectorXd g = potential.gradient(collisions, mesh, V);

        Eigen::VectorXd fg;
        fd::finite_gradient(
            Eigen::VectorXd::Zero(1),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = V + fd::unflatten(test_dir, 3) * y(0);
                EspCollisions collisions_;
                collisions_.build(mesh, V_, params, adaptive.get());
                return potential(collisions_, mesh, V_);
            },
            fg, fd::AccuracyOrder::FOURTH, 1e-5);

        REQUIRE(
            abs(fg(0) - g.dot(test_dir)) < std::max(fg.norm() * 1e-5, 1e-9));
    }

    SECTION("hessian")
    {
        Eigen::MatrixXd h = potential.hessian(collisions, mesh, V);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            Eigen::VectorXd::Zero(1),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = V + fd::unflatten(test_dir, 3) * y(0);
                EspCollisions collisions_;
                collisions_.build(mesh, V_, params, adaptive.get());
                return potential.gradient(collisions_, mesh, V_);
            },
            fh, fd::AccuracyOrder::FOURTH, 1e-5);

        REQUIRE(
            (fh.col(0) - h * test_dir).norm()
            < std::max(fh.norm() * 1e-6, 1e-9));
    }
}

#if defined(NDEBUG) && !defined(WIN32)
static std::string tagsopt =
    "[esp_potential], [esp_potential_3d]";
#else
static std::string tagsopt =
    "[.][esp_potential], [.][esp_potential_3d]";
#endif

TEST_CASE("Convergent Quadrature Gradient and Hessian Expensive", tagsopt)
{
    auto [V, E, F, mesh] = load_wrapped_sphere();

    const double dhat = 0.1;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    EspParameters params(dhat, dbar_factor, 0);

    const bool use_adaptive = GENERATE(true, false);
    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;

    EspCollisions collisions;
    collisions.build(mesh, V, params, adaptive.get());

    const bool normalize_weights = GENERATE(true, false);
    EspPotential potential(params, normalize_weights);

    SECTION("gradient")
    {
        Eigen::VectorXd g = potential.gradient(collisions, mesh, V);

        Eigen::VectorXd fg;
        fd::finite_gradient(
            fd::flatten(V),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = fd::unflatten(y, 3);
                EspCollisions collisions_;
                collisions_.build(mesh, V_, params, adaptive.get());
                return potential(collisions_, mesh, V_);
            },
            fg, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE((fg - g).norm() < std::max(1e-8, fg.norm()) * 1e-6);
    }

    SECTION("hessian")
    {
        Eigen::MatrixXd h = potential.hessian(collisions, mesh, V);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            fd::flatten(V),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = fd::unflatten(y, 3);
                EspCollisions collisions_;
                collisions_.build(mesh, V_, params, adaptive.get());
                return potential.gradient(collisions_, mesh, V_);
            },
            fh, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE((fh - h).norm() < std::max(1e-8, fh.norm()) * 1e-6);
    }
}

TEST_CASE(
    "Convergent Quadrature Zero on Sphere",
    "[esp_potential], [esp_potential_3d]")
{
    auto [V, E, F, mesh] = load_triangle_mesh(
        (tests::DATA_DIR / "../src/tests/potential/sphere.obj").string());

    const double dhat = 0.2;
    const double dbar_factor = GENERATE(1.0, 0.9);
    EspParameters params(dhat, dbar_factor, 0);

    const bool adaptive_dhat = GENERATE(true, false);
    auto adaptive = adaptive_dhat
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;
    EspCollisions collisions;
    collisions.build(mesh, V, params, adaptive.get());

    EspPotential potential(params);
    double val = potential(collisions, mesh, V);
    REQUIRE(val == 0);

    auto g = potential.gradient(collisions, mesh, V);
    REQUIRE(g.norm() == 0);

    auto H = potential.hessian(collisions, mesh, V);
    REQUIRE(H.norm() == 0);
}

TEST_CASE(
    "Number of Pairs", "[esp_potential], [esp_potential_3d]")
{
    double dhat = -1;
    std::string mesh_name;
    // SECTION("mesh1")
    // {
    //     dhat = 1e-2;
    //     mesh_name = "bunny.ply";
    // }
    SECTION("mesh2")
    {
        dhat = 1e-2;
        mesh_name = "armadillo-rollers/327.ply";
    }

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    REQUIRE(success);
    CAPTURE(mesh_name);

    CollisionMesh mesh;
    mesh = CollisionMesh(
        std::vector<bool>(vertices.rows(), true),
        std::vector<bool>(vertices.rows(), false), vertices, edges, faces);

    {
        EspCollisions collisions;
        EspParameters params(dhat, 1., 0);
        collisions.build(mesh, vertices, params);

        std::cout << "ESP collision size " << collisions.size()
                  << std::endl;
    }

    {
        NormalCollisions collisions;
        collisions.build(mesh, vertices, dhat);

        std::cout << "normal collision size " << collisions.size() << std::endl;
    }

    {
        EspCollisions collisions;
        EspParameters params(dhat, 1., 0);
        collisions.build(mesh, vertices, params);

        std::cout << "ESP collision pairs (before cancellation) "
                  << collisions.num_quadrature_collision_pairs << std::endl;

        const auto dist = collisions.edge_id_count_distribution();
        std::cout << "edge id count distribution (count: num_edges):"
                  << std::endl;
        for (const auto& [count, num_edges] : dist) {
            std::cout << "  " << count << ": " << num_edges << ", ";
        }
    }
}

TEST_CASE(
    "Convergent Quadrature Vertex Hessian",
    "[esp_potential], [esp_potential_3d]")
{
    const auto method = make_default_broad_phase();
    auto [V, E, F, mesh] = load_wrapped_sphere();

    const double dhat = 0.15;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    EspParameters params(dhat, dbar_factor, 0);

    const bool use_adaptive = GENERATE(true, false);
    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;

    Candidates candidates;
    candidates.build(mesh, V, dhat / 2, method.get(), true);
    candidates.convert_candidates_to_sets();
    PointPotential point_potential(mesh, candidates, params);

    for (int vid = 0; vid < V.rows(); ++vid) {
        size_t num_collision_pairs = 0;
        const auto collisions = point_potential.build_collisions_at_vertex(
            V, vid, num_collision_pairs);

        if (collisions->size() == 0) {
            continue;
        }

        std::vector<int> indices;
        {
            Eigen::VectorXd local_grad = PointPotentialHelper::
                evaluate_potential_gradient_at_vertex_with_cached_collisions(
                    V, *collisions, params, adaptive.get());
            indices = collisions->dofs();

            if (local_grad.norm() < 1e-10) {
                continue;
            }
        }

        Eigen::MatrixXd h = PointPotentialHelper::
            evaluate_potential_hessian_at_vertex_with_cached_collisions(
                V, *collisions, params, adaptive.get(),
                PSDProjectionMethod::NONE);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            fd::flatten(V)(indices),
            [&](const Eigen::VectorXd& y) {
                Eigen::VectorXd y_ = fd::flatten(V);
                y_(indices) = y;
                Eigen::MatrixXd V_fd = fd::unflatten(y_, 3);

                return PointPotentialHelper::
                    evaluate_potential_gradient_at_vertex_with_cached_collisions(
                        V_fd, *collisions, params, adaptive.get());
            },
            fh, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE(
            (h - fh).norm() < 1e-6 * std::max({ h.norm(), fh.norm(), 1e-8 }));
    }
}

TEST_CASE(
    "Convergent Quadrature Face Hessian",
    "[esp_potential], [esp_potential_3d]")
{
    const auto method = make_default_broad_phase();
    auto [V, E, F, mesh] = load_wrapped_sphere();

    const double dhat = 0.15;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    EspParameters params(dhat, dbar_factor, 0);

    const bool use_adaptive = GENERATE(true, false);
    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;

    Candidates candidates;
    candidates.build(mesh, V, dhat / 2, method.get(), true);
    candidates.convert_candidates_to_sets();
    PointPotential point_potential(mesh, candidates, params);

    for (int fid = 0; fid < F.rows(); ++fid) {

        size_t num_collision_pairs = 0;
        const auto collisions = point_potential.build_collisions_at_face_center(
            V, fid, num_collision_pairs);

        if (collisions->size() == 0) {
            continue;
        }

        Eigen::Vector3<index_t> vids;
        vids << F(fid, 0), F(fid, 1), F(fid, 2);

        Eigen::RowVector3d face_center =
            (V.row(vids[0]) + V.row(vids[1]) + V.row(vids[2])) / 3.;

        VertexMatrixView<3> V_extended(V, face_center);

        std::vector<int> indices;
        {
            Eigen::VectorXd local_grad = PointPotentialHelper::
                evaluate_potential_gradient_at_face_center_with_cached_collisions(
                    V_extended, *collisions, params, adaptive.get());
            indices = collisions->dofs();

            if (local_grad.norm() < 1e-10) {
                continue;
            }
        }

        Eigen::MatrixXd h = PointPotentialHelper::
            evaluate_potential_hessian_at_face_center_with_cached_collisions(
                V_extended, *collisions, params, adaptive.get(),
                PSDProjectionMethod::NONE);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            fd::flatten(V)(indices),
            [&](const Eigen::VectorXd& y) {
                Eigen::VectorXd y_ = fd::flatten(V);
                y_(indices) = y;
                Eigen::MatrixXd V_fd = fd::unflatten(y_, 3);
                Eigen::RowVector3d face_center_fd =
                    (V_fd.row(vids[0]) + V_fd.row(vids[1]) + V_fd.row(vids[2]))
                    / 3.;
                VertexMatrixView<3> V_fd_extended(V_fd, face_center_fd);

                return PointPotentialHelper::
                    evaluate_potential_gradient_at_face_center_with_cached_collisions(
                        V_fd_extended, *collisions, params, adaptive.get());
            },
            fh, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE(
            (h - fh).norm() < 1e-6 * std::max({ h.norm(), fh.norm(), 1e-8 }));
    }
}

// Test FV-3D mollification: two aligned cubes with vertices approaching face
// edges This configuration makes the mollification issue critical: vertices of
// one cube approach the faces of another cube, with closest points near
// triangle edges
TEST_CASE(
    "ESP potential 3D finite differences (FV mollification)",
    "[esp_potential], [esp_potential_3d]")
{
    const auto method = make_default_broad_phase();

    // Load cube mesh
    std::string cube_path = (tests::DATA_DIR / "cube.obj").string();
    Eigen::MatrixXd V_single;
    Eigen::MatrixXi F_single;
    if (!igl::read_triangle_mesh(cube_path, V_single, F_single)) {
        SKIP("Could not load cube.obj");
    }

    // Create two cubes: one fixed, one translated slightly
    Eigen::MatrixXd V(V_single.rows() * 2, 3);
    V.topRows(V_single.rows()) = V_single;
    // Second cube: translate along x-axis to create face-vertex collisions
    // with aligned vertices approaching the faces of the first cube
    V.bottomRows(V_single.rows()) =
        V_single.rowwise() + Eigen::RowVector3d(1.001, 0, 0);

    Eigen::MatrixXi F(F_single.rows() * 2, 3);
    F.topRows(F_single.rows()) = F_single;
    F.bottomRows(F_single.rows()) =
        F_single.array() + static_cast<int>(V_single.rows());

    Eigen::MatrixXi E;
    igl::edges(F, E);

    CollisionMesh mesh(V, E, F);

    const double dhat = .5;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    EspParameters params(dhat, dbar_factor, 0);

    const bool use_adaptive = GENERATE(true, false);
    CAPTURE(use_adaptive);

    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;
    if (adaptive) {
        adaptive->scale(
            1.2); // manually scale adaptive dhat so energy is not zero
    }

    Candidates candidates;
    candidates.build(mesh, V, dhat / 2, method.get(), true);
    candidates.convert_candidates_to_sets();

    EspCollisions collisions;
    collisions.build(candidates, mesh, V, params, adaptive.get());
    std::cerr << "EspCollisions after build: " << collisions.size()
              << "\n";

    REQUIRE(!collisions.empty());
    REQUIRE(!has_intersections(mesh, V));

    EspPotential potential(params);
    double energy = potential(collisions, mesh, V);
    CAPTURE(energy);
    CHECK(energy > 0);
    CHECK(std::isfinite(energy));

    // Test gradient accuracy: without mollification of (u,v),
    // this will fail when closest point is near face edges where adaptive dhat
    // varies
    Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        fd::flatten(V),
        [&](const Eigen::VectorXd& x) {
            return potential(collisions, mesh, fd::unflatten(x, V.cols()));
        },
        fgrad, fd::AccuracyOrder::SECOND, 1e-8);

    CAPTURE(grad.norm());
    CAPTURE(fgrad.norm());
    const double error = (grad - fgrad).norm();
    const double threshold =
        1e-3 * std::max({ grad.norm(), fgrad.norm(), 1e-8 });
    CAPTURE(error);
    CAPTURE(threshold);
    // Without mollification: FD will mismatch analytical gradient near face
    // edges With mollification: both should agree
    CHECK(error < threshold);
}

// 2D TESTS //

TEST_CASE(
    "ESP potential codim",
    "[esp_potential], [esp_potential_2d]")
{
    const auto method = make_default_broad_phase();
    double dhat = 2;
    const int quadrature_order = 2;
    EspParameters params(dhat, 1., quadrature_order);

    Eigen::MatrixXd vertices(4, 2);
    Eigen::MatrixXi edges(2, 2);

    vertices << -1, 0, 0, 0, 1, 0, 1.5, 0.2;
    edges << 0, 1, 1, 2;

    CollisionMesh mesh = make_2d_collision_mesh(vertices, edges);

    EspCollisions collisions;
    collisions.build(mesh, vertices, params, nullptr, method.get());
    CAPTURE(dhat, method);
    CHECK(!collisions.empty());
    CHECK(!has_intersections(mesh, vertices));

    EspPotential potential(params);
    double energy = potential(collisions, mesh, vertices);
    CHECK(energy != 0);

    // Gradient
    const Eigen::VectorXd grad = potential.gradient(collisions, mesh, vertices);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        fd::flatten(vertices),
        [&](const Eigen::VectorXd& x) {
            return potential(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        },
        fgrad, fd::AccuracyOrder::SECOND, 1e-8);

    REQUIRE(grad.squaredNorm() > 1e-8);
    CHECK((grad - fgrad).norm() / grad.norm() < 1e-4);

    // Hessian
    Eigen::MatrixXd hess = potential.hessian(collisions, mesh, vertices);

    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        fd::flatten(vertices),
        [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, fd::unflatten(x, vertices.cols()));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-8);

    REQUIRE(hess.squaredNorm() > 1e-8);
    CHECK((hess - fhess).norm() / hess.norm() < 1e-3);
}

TEST_CASE(
    "ESP potential 2D no forces",
    "[esp_potential], [esp_potential_2d]")
{
    const auto method = make_default_broad_phase();
    Eigen::MatrixXd V;
    Eigen::MatrixXi E;
    double dhat = 1.;
    const int quadrature_order = GENERATE(1, 2, 7, 10, 14);
    EspParameters params(dhat, 1., quadrature_order);

    const bool use_adaptive = GENERATE(true, false);
    std::string name;
    SECTION("square_1")
    {
        name = "square_1";
        V.resize(4, 2);
        E.resize(4, 2);
        V << -1., -1., 1., -1., 1., 1., -1., 1.;
        E << 0, 1, 1, 2, 2, 3, 3, 0;
    }
    SECTION("square_2")
    {
        name = "square_2";
        V.resize(8, 2);
        E.resize(8, 2);
        V << -1., -1., 0., -1., 1., -1., 1., 0., 1., 1., 0., 1., -1., 1., -1.,
            0.;
        E << 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 0;
    }
    SECTION("circle")
    {
        const int n = GENERATE(5, 6, 7, 8, 9, 10, 50, 100, 200, 2000);
        name = "circle" + std::to_string(n);
        V.resize(n, 2);
        E.resize(n, 2);
        for (int i = 0; i < n; i++) {
            V(i, 0) = std::cos(2 * M_PI * i / n);
            V(i, 1) = std::sin(2 * M_PI * i / n);
        }
        for (int i = 0; i < n; i++) {
            E(i, 0) = i;
            E(i, 1) = (i + 1) % n;
        }
    }

    CollisionMesh mesh = make_2d_collision_mesh(V, E);

    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;
    EspCollisions collisions;
    collisions.build(mesh, V, params, adaptive.get(), method.get());

    REQUIRE(!has_intersections(mesh, V));

    EspPotential potential(params);
    double energy = potential(collisions, mesh, V);
    CAPTURE(name);
    CAPTURE(quadrature_order);
    CHECK(energy == 0);

    Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
    CHECK(grad.squaredNorm() == 0);

    Eigen::MatrixXd hess = potential.hessian(collisions, mesh, V);
    CHECK(hess.squaredNorm() == 0);
}

TEST_CASE(
    "ESP potential 2D finite differences",
    "[esp_potential], [esp_potential_2d]")
{
    const auto method = make_default_broad_phase();
    Eigen::MatrixXd V;
    Eigen::MatrixXi E;
    double dhat = 0.6;
    constexpr double BA = 0; // a small constant to break perfect alignments
    const int quadrature_order = GENERATE(1, 2, 7, 14);
    EspParameters params(dhat, 1., quadrature_order);
    const bool adaptive_dhat = GENERATE(true, false);
    CAPTURE(quadrature_order);
    CAPTURE(adaptive_dhat);

    auto run_checks = [&]() {
        CollisionMesh mesh = make_2d_collision_mesh(V, E);

        auto adaptive = adaptive_dhat
            ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
            : nullptr;

        EspCollisions collisions;
        collisions.build(mesh, V, params, adaptive.get(), method.get());

        REQUIRE(!collisions.empty());
        REQUIRE(!has_intersections(mesh, V));

        EspPotential potential(params);
        double energy = potential(collisions, mesh, V);
        if (!adaptive_dhat)
            CHECK(energy > 0);

        Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
        if (!adaptive_dhat)
            REQUIRE(grad.squaredNorm() > 1e-8);
        Eigen::VectorXd fgrad;
        fd::finite_gradient(
            fd::flatten(V),
            [&](const Eigen::VectorXd& x) {
                return potential(collisions, mesh, fd::unflatten(x, V.cols()));
            },
            fgrad, fd::AccuracyOrder::SECOND, 1e-8);
        CAPTURE(grad.norm());
        CAPTURE(fgrad.norm());
        CHECK(
            (grad - fgrad).norm() < std::max(
                1e-4 * std::max({ grad.norm(), fgrad.norm(), 1e-8 }), 1e-9));

        Eigen::MatrixXd hess = potential.hessian(collisions, mesh, V);
        if (!adaptive_dhat)
            REQUIRE(hess.squaredNorm() > 1e-3);
        Eigen::MatrixXd fhess;
        fd::finite_jacobian(
            fd::flatten(V),
            [&](const Eigen::VectorXd& x) {
                return potential.gradient(
                    collisions, mesh, fd::unflatten(x, V.cols()));
            },
            fhess, fd::AccuracyOrder::SECOND, 1e-12);
        CAPTURE(hess.norm());
        CAPTURE(fhess.norm());
        CHECK(
            (hess - fhess).norm() < std::max(
                3e-3 * std::max({ hess.norm(), fhess.norm(), 1e-8 }), 1e-9));
    };

    SECTION("Corners")
    {
        const double P0x = GENERATE(.01, -.01, 0.0);
        const double P1y = GENERATE(.49, .5, .51);
        CAPTURE(P0x, P1y);
        V.resize(8, 2);
        E.resize(8, 2);
        V << -1., 1., -1., 0., 0., 0., P0x, .5 + BA, 0., 1., 1., 0., 1., 1.,
            .02, P1y;
        E << 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 5, 6, 6, 7, 7, 5;
        run_checks();
    }

    SECTION("squares")
    {
        V.resize(8, 2);
        E.resize(8, 2);
        E << 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4;
        SECTION("horizontal_squares")
        {
            INFO("horizontal_squares");
            V << -1., 1. + BA, -1., 0. + BA, -.1, 0. + BA, -.1, 1. + BA, .1, 1.,
                .1, 0., 1., 0., 1., 1.;
            run_checks();
        }
        SECTION("vertical_squares")
        {
            INFO("vertical_squares");
            V << 0. + BA, -1., 1. + BA, -1., 1. + BA, -.1, 0. + BA, -.1, 0., .1,
                1., .1, 1., 1., 0., 1.;
            run_checks();
        }
    }

    SECTION("mesh_1")
    {
        INFO("mesh 1");
        std::string mesh_name =
            (tests::DATA_DIR / "gcp" / "nonlinear_solve_iter020.obj").string();
        bool success = igl::readCSV(mesh_name + "-v.csv", V);
        success = success && igl::readCSV(mesh_name + "-e.csv", E);
        REQUIRE(success);
        V.col(0) += Eigen::VectorXd::Random(V.rows()) * BA;
        run_checks();
    }

    SECTION("mesh_2")
    {
        INFO("mesh 2");
        std::string mesh_name =
            (tests::DATA_DIR / "gcp" / "simple_2d.obj").string();
        bool success = igl::readCSV(mesh_name + "-v.csv", V);
        success = success && igl::readCSV(mesh_name + "-e.csv", E);
        REQUIRE(success);
        V.col(0) += Eigen::VectorXd::Random(V.rows()) * BA;
        run_checks();
    }
}

// 3D FACE QUADRATURE TESTS //

// Verify that face quadrature gives gradient/hessian consistent with finite
// differences on the wrapped-sphere geometry, for several quadrature orders.
TEST_CASE(
    "Face Quadrature Gradient and Hessian",
    "[esp_potential], [esp_potential_3d]")
{
    auto [V, E, F, mesh] = load_wrapped_sphere();

    const double dhat = 0.15;
    const int quad_order = GENERATE(
        0, 3, 6); // Using fekete rules, orders 1-2-3 and 4-5-6 are the same
    EspParameters params(dhat, 1., quad_order);

    const bool use_adaptive = GENERATE(true, false);
    const bool normalize_weights = GENERATE(true, false);
    EspPotential potential(params, normalize_weights);

    // Compute once so every FD step uses identical dhat values.
    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;
    if (adaptive) {
        adaptive->scale(
            1.2); // manually scale adaptive dhat so energy is not zero
    }

    EspCollisions collisions;
    collisions.build(mesh, V, params, adaptive.get());

    REQUIRE(potential(collisions, mesh, V) != 0);

    // Directional finite-difference to keep the test inexpensive
    Eigen::VectorXd test_dir(V.size());
    for (int i = 0; i < test_dir.size(); i++) {
        test_dir(i) = i;
    }
    test_dir.normalize();

    SECTION("gradient")
    {
        Eigen::VectorXd g = potential.gradient(collisions, mesh, V);

        Eigen::VectorXd fg;
        fd::finite_gradient(
            Eigen::VectorXd::Zero(1),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = V + fd::unflatten(test_dir, 3) * y(0);
                EspCollisions c;
                c.build(mesh, V_, params, adaptive.get());
                return potential(c, mesh, V_);
            },
            fg, fd::AccuracyOrder::SECOND, 1e-7);

        REQUIRE(abs(fg(0) - g.dot(test_dir)) < fg.norm() * 1e-5);
    }

    SECTION("hessian")
    {
        Eigen::MatrixXd h = potential.hessian(collisions, mesh, V);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            Eigen::VectorXd::Zero(1),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = V + fd::unflatten(test_dir, 3) * y(0);
                EspCollisions c;
                c.build(mesh, V_, params, adaptive.get());
                return potential.gradient(c, mesh, V_);
            },
            fh, fd::AccuracyOrder::SECOND, 1e-6);

        REQUIRE((fh.col(0) - h * test_dir).norm() < fh.norm() * 1e-4);
    }
}

// Verify that with normalize_weights = true, the global hessian is PSD whenever
// project_hessian_to_psd is set. The non-normalized branch uses local PSD
// projection which trivially yields a PSD assembly.
TEST_CASE(
    "Convergent Quadrature Hessian PSD",
    "[esp_potential], [esp_potential_3d]")
{
    auto [V, E, F, mesh] = load_wrapped_sphere();

    const double dhat = 0.15;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    EspParameters params(dhat, dbar_factor, 0);

    const bool use_adaptive = GENERATE(true, false);
    const bool normalize_weights = GENERATE(true, false);
    const PSDProjectionMethod psd_method =
        GENERATE(PSDProjectionMethod::CLAMP, PSDProjectionMethod::ABS);

    EspPotential potential(params, normalize_weights);

    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;
    EspCollisions collisions;
    collisions.build(mesh, V, params, adaptive.get());

    Eigen::SparseMatrix<double> H =
        potential.hessian(collisions, mesh, V, psd_method);
    Eigen::MatrixXd Hd(H);
    // Symmetrize numerically to remove tiny asymmetry from triplet ordering.
    Hd = 0.5 * (Hd + Hd.transpose()).eval();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
        Hd, Eigen::EigenvaluesOnly);
    REQUIRE(es.info() == Eigen::Success);
    const double lambda_min = es.eigenvalues().minCoeff();
    const double lambda_max = es.eigenvalues().maxCoeff();
    const double tol = std::max(1e-10, 1e-10 * std::abs(lambda_max));

    INFO(
        "normalize_weights="
        << normalize_weights << " method=" << static_cast<int>(psd_method)
        << " lambda_min=" << lambda_min << " lambda_max=" << lambda_max);
    REQUIRE(lambda_min >= -tol);
}

TEST_CASE("NearFarBarrier decomposition", "[esp_potential][barrier]")
{
    const double dhat = 0.1;
    const double alpha = GENERATE(0.01, 0.25, 0.5, 0.75, 0.99);

    enum class BarrierType {
        ClampedLog,
        ClampedLogSq,
        Cubic,
        TwoStage,
        InversePower1,
        InversePower2
    };

    const BarrierType type = GENERATE(
        BarrierType::ClampedLog, BarrierType::ClampedLogSq, BarrierType::Cubic,
        BarrierType::TwoStage, BarrierType::InversePower1,
        BarrierType::InversePower2);

    auto run_test = [&](const auto& base) {
        NearFarBarrier nf(&base, alpha);

        constexpr int N = 200;
        for (int i = 0; i < N; ++i) {
            const double d = dhat * (i + 1.0) / N;

            const double b = base(d, dhat);
            CHECK(nf.near(d, dhat) + nf.far(d, dhat) == Catch::Approx(b));

            const double db = base.first_derivative(d, dhat);
            CHECK(
                nf.first_derivative_near(d, dhat)
                    + nf.first_derivative_far(d, dhat)
                == Catch::Approx(db));

            const double ddb = base.second_derivative(d, dhat);
            CHECK(
                nf.second_derivative_near(d, dhat)
                    + nf.second_derivative_far(d, dhat)
                == Catch::Approx(ddb));

            CAPTURE(d);
            CAPTURE(alpha);
            CAPTURE(dhat);
            CAPTURE(alpha * dhat);
            CAPTURE(alpha * dhat / 2);
            // Check that near barrier is 0 above alpha*dhat and non-zero below
            constexpr double eps_tol = 1e-9;
            if (d >= alpha * dhat) {
                CHECK(nf.near(d, dhat) == 0.0);
            } else if (d < alpha * dhat - eps_tol) {
                CHECK(nf.near(d, dhat) > 0.0);
            }

            // Check that far barrier is 0 below dhat*alpha/2 and non-zero above
            if (d <= dhat * alpha / 2 || d >= dhat) {
                CHECK(nf.far(d, dhat) == 0.0);
            } else if (d > dhat * alpha / 2 + eps_tol) {
                CHECK(nf.far(d, dhat) > 0.0);
            }
        }
    };

    switch (type) {
    case BarrierType::ClampedLog:
        run_test(ClampedLogBarrier<>());
        break;
    case BarrierType::ClampedLogSq:
        run_test(ClampedLogSqBarrier<>());
        break;
    case BarrierType::Cubic:
        run_test(CubicBarrier<>());
        break;
    case BarrierType::TwoStage:
        run_test(TwoStageBarrier<>());
        break;
    case BarrierType::InversePower1:
        run_test(InversePowerBarrier(1.0));
        break;
    case BarrierType::InversePower2:
        run_test(InversePowerBarrier(2.0));
        break;
    }
}

// ---------------------------------------------------------------------------
// Adaptive support tests
// ---------------------------------------------------------------------------

// With a large dhat, many pairs are in contact without adaptive support.
// After computing the adaptive support (which iteratively shrinks per-vertex
// dhat until every primitive is beyond its own dhat), the potential must be
// exactly zero.
TEST_CASE(
    "Adaptive Support Reduces Potential to Zero (3D)",
    "[adaptive_support], [esp_potential_3d]")
{
    auto [V, E, F, mesh] = load_wrapped_sphere();

    // dhat = 0.3 is intentionally large: many vertex pairs are within range,
    // so the potential without adaptive support is clearly non-zero.
    const double dhat = 10;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    EspParameters params(dhat, dbar_factor, 0);
    EspPotential potential(params);

    // Baseline: without adaptive, potential must be non-zero.
    {
        EspCollisions collisions;
        collisions.build(mesh, V, params);
        const double energy = potential(collisions, mesh, V);
        REQUIRE(energy > 0);
    }

    // With adaptive support the per-primitive dhat values are reduced until
    // no collision pair contributes, so the evaluated potential is exactly 0.
    auto adaptive = EspCollisions::compute_adaptive_dhat(mesh, V, params);
    REQUIRE(adaptive != nullptr);

    // All vertex dhat values must be in (0, params.dhat] after reduction.
    bool any_reduced = false;
    for (int i = 0; i < mesh.num_vertices(); i++) {
        CHECK(adaptive->vertex(i) > 0.0);
        CHECK(adaptive->vertex(i) <= dhat);
        if (adaptive->vertex(i) < dhat)
            any_reduced = true;
    }
    CHECK(any_reduced);

    {
        EspCollisions collisions;
        collisions.build(mesh, V, params, adaptive.get());
        const double energy = potential(collisions, mesh, V);
        CHECK(energy == 0.0);
    }
}

TEST_CASE(
    "Adaptive Support Reduces Potential to Zero (2D)",
    "[adaptive_support], [esp_potential_2d]")
{
    const auto method = make_default_broad_phase();
    Eigen::MatrixXd V;
    Eigen::MatrixXi E;

    SECTION("Corners")
    {
        const double P0x = GENERATE(.01, -.01, 0.0);
        const double P1y = GENERATE(.49, .5, .51);
        CAPTURE(P0x, P1y);
        V.resize(8, 2);
        E.resize(8, 2);
        V << -1., 1., -1., 0., 0., 0., P0x, .5, 0., 1., 1., 0., 1., 1., .02,
            P1y;
        E << 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 5, 6, 6, 7, 7, 5;
    }

    SECTION("squares")
    {
        V.resize(8, 2);
        E.resize(8, 2);
        E << 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4;
        SECTION("horizontal_squares")
        {
            INFO("horizontal_squares");
            V << -1., 1., -1., 0., -.1, 0., -.1, 1., .1, 1., .1, 0., 1., 0., 1.,
                1.;
        }
        SECTION("vertical_squares")
        {
            INFO("vertical_squares");
            V << 0., -1., 1., -1., 1., -.1, 0., -.1, 0., .1, 1., .1, 1., 1., 0.,
                1.;
        }
    }

    CollisionMesh mesh = make_2d_collision_mesh(V, E);
    REQUIRE(!has_intersections(mesh, V));

    const double dhat = 10.;
    const int quad_order = 14;
    EspParameters params(dhat, 1.0, quad_order);
    EspPotential potential(params);

    // Baseline: without adaptive, potential must be non-zero.
    {
        EspCollisions collisions;
        collisions.build(mesh, V, params, nullptr, method.get());
        const double energy = potential(collisions, mesh, V);
        REQUIRE(energy != 0);
    }

    // With adaptive support: per-primitive dhat values fall below 0.2
    // so the barrier is exactly zero for every pair.
    auto adaptive = EspCollisions::compute_adaptive_dhat(mesh, V, params);
    REQUIRE(adaptive != nullptr);

    // Primitive vertices (those on the far edge) must have reduced dhat.
    bool any_reduced = false;
    for (int i = 0; i < mesh.num_vertices(); i++) {
        CHECK(adaptive->vertex(i) > 0.0);
        CHECK(adaptive->vertex(i) <= dhat);
        if (adaptive->vertex(i) < dhat)
            any_reduced = true;
    }
    REQUIRE(any_reduced);

    {
        EspCollisions collisions;
        collisions.build(mesh, V, params, adaptive.get(), method.get());
        const double energy = potential(collisions, mesh, V);
        CHECK(energy == 0.0);
    }
}

// Same check for the 3D face-quadrature variant: ESP quadrature points
// inside each face must also yield a PSD assembly under combined projection.
TEST_CASE(
    "Face Quadrature Hessian PSD",
    "[esp_potential], [esp_potential_3d]")
{
    auto [V, E, F, mesh] = load_wrapped_sphere();

    const double dhat = 0.15;
    const int quad_order = GENERATE(0, 3, 6);
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    EspParameters params(dhat, dbar_factor, quad_order);

    const bool use_adaptive = GENERATE(true, false);
    const bool normalize_weights = GENERATE(true, false);
    const PSDProjectionMethod psd_method =
        GENERATE(PSDProjectionMethod::CLAMP, PSDProjectionMethod::ABS);

    EspPotential potential(params, normalize_weights);

    auto adaptive = use_adaptive
        ? EspCollisions::compute_adaptive_dhat(mesh, V, params)
        : nullptr;
    EspCollisions collisions;
    collisions.build(mesh, V, params, adaptive.get());

    Eigen::SparseMatrix<double> H =
        potential.hessian(collisions, mesh, V, psd_method);
    Eigen::MatrixXd Hd(H);
    Hd = 0.5 * (Hd + Hd.transpose()).eval();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
        Hd, Eigen::EigenvaluesOnly);
    REQUIRE(es.info() == Eigen::Success);
    const double lambda_min = es.eigenvalues().minCoeff();
    const double lambda_max = es.eigenvalues().maxCoeff();
    const double tol = std::max(1e-10, 1e-10 * std::abs(lambda_max));

    INFO(
        "normalize_weights="
        << normalize_weights << " quad_order=" << quad_order
        << " method=" << static_cast<int>(psd_method)
        << " lambda_min=" << lambda_min << " lambda_max=" << lambda_max);
    REQUIRE(lambda_min >= -tol);
}
