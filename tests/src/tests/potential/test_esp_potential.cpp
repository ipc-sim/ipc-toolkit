#include <tests/config.hpp>
#include <tests/dof_layout.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/esp/arbitrary_point_esp.hpp>
#include <ipc/esp/esp_potential.hpp>
#include <ipc/distance/line_line.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>
#include <igl/readCSV.h>
#include <ipc/ipc.hpp>

#include "igl/read_triangle_mesh.h"
#include "ipc/distance/edge_edge.hpp"

#include "ipc/esp/quadrature_potential.hpp"

#include <ipc/distance/distance_type_exact.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/esp/collisions/esp_collision_template.hpp>
#include <ipc/esp/collisions/esp_quadrature.hpp>
#include <ipc/esp/collisions/vertex_matrix_view.hpp>

#include <tbb/global_control.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#ifndef _WIN32
#include <pthread.h>
#endif

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
    // Fail loudly if the mesh is missing (e.g., an outdated test-data
    // checkout) instead of running the tests on an empty mesh.
    REQUIRE(igl::read_triangle_mesh(path, data.V, data.F));
    igl::edges(data.F, data.E);
    data.mesh = CollisionMesh(data.V, data.E, data.F);
    return data;
}

TriMeshData load_wrapped_sphere()
{
    return load_triangle_mesh(
        (tests::DATA_DIR / "wrapped_sphere.obj").string());
}

/// Smallest and largest eigenvalues of a symmetric sparse matrix. Only the
/// rows and columns with a nonzero entry enter the dense solve: every other
/// one contributes an exact zero eigenvalue, so the result is unchanged and
/// the solve is much smaller when few DOFs are in contact.
std::pair<double, double> eigenvalue_range(const Eigen::SparseMatrix<double>& H)
{
    std::vector<int> local(H.rows(), -1);
    int n_active = 0;
    for (int k = 0; k < H.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
            if (it.value() != 0) {
                for (const auto i : { it.row(), it.col() }) {
                    if (local[i] < 0) {
                        local[i] = n_active++;
                    }
                }
            }
        }
    }
    if (n_active == 0) {
        return { 0.0, 0.0 };
    }

    Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(n_active, n_active);
    for (int k = 0; k < H.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
            if (it.value() != 0) {
                Hd(local[it.row()], local[it.col()]) += it.value();
            }
        }
    }
    // Symmetrize numerically to remove tiny asymmetry from triplet ordering.
    Hd = 0.5 * (Hd + Hd.transpose()).eval();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
        Hd, Eigen::EigenvaluesOnly);
    REQUIRE(es.info() == Eigen::Success);
    double lambda_min = es.eigenvalues().minCoeff();
    double lambda_max = es.eigenvalues().maxCoeff();
    if (n_active < H.rows()) { // the inactive DOFs' zero eigenvalues
        lambda_min = std::min(lambda_min, 0.0);
        lambda_max = std::max(lambda_max, 0.0);
    }
    return { lambda_min, lambda_max };
}

/// Positive-weight face quadrature rules for the 3D face-quadrature tests:
/// "vertices+centroid" (the three vertices and the centroid) or
/// "edge-midpoints" (one point on each edge of the face).
FaceQuadRule make_face_quad_rule(const std::string& name)
{
    if (name == "vertices+centroid") {
        return {
            { { { 1., 0., 0. } }, 0.25 },
            { { { 0., 1., 0. } }, 0.25 },
            { { { 0., 0., 1. } }, 0.25 },
            { { { 1. / 3, 1. / 3, 1. / 3 } }, 0.25 },
        };
    }
    REQUIRE(name == "edge-midpoints");
    return {
        { { { 0., .5, .5 } }, 1. / 3 },
        { { { .5, 0., .5 } }, 1. / 3 },
        { { { .5, .5, 0. } }, 1. / 3 },
    };
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
    return std::make_shared<InversePowerBarrier<>>(2.0);
}
inline std::shared_ptr<Barrier> make_linear_inverse_barrier()
{
    return std::make_shared<InversePowerBarrier<>>(1.0);
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
        ESPParameters params(dhat, 1., 0);
        params.barrier = barrier;

        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        ESPPotential potential(params);

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

// Two closed tetrahedra, apex to apex along z and separated by a gap. Each apex
// has valence 3, so ESP has to resolve the vertex-vertex contact through the
// vertex-vertex and edge-vertex candidates.
inline void build_apex_to_apex_tets(
    const double gap,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
    const double c = std::cos(2 * M_PI / 3), s = std::sin(2 * M_PI / 3);
    const double r = std::cos(M_PI / 3), q = std::sin(M_PI / 3);
    V.resize(8, 3);
    V << 0, 0, gap / 2, 1, 0, 1, c, s, 1, c, -s, 1, //
        0, 0, -gap / 2, r, q, -1, -1, 0, -1, r, -q, -1;
    F.resize(8, 3);
    F << 0, 2, 1, 0, 3, 2, 0, 1, 3, 1, 2, 3, //
        4, 5, 6, 4, 6, 7, 4, 7, 5, 5, 7, 6;
    igl::edges(F, E);
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
    // (current: max|ΔP|≈6.4e-4, max|Δ|g||≈94; faces weighted by their area,
    // rule weights summing to 1).
    CHECK(stats.max_shift_P < 6.4e-3);
    CHECK(stats.max_shift_g < 940);
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
    // (current: max|ΔP|≈2.3e-2, max|Δ|g||≈57; faces weighted by their area,
    // rule weights summing to 1).
    CHECK(stats.max_shift_P < 0.23);
    CHECK(stats.max_shift_g < 575);
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
    // (current: max|ΔP|≈4.2e-4, max|Δ|g||≈1.9; faces weighted by their area,
    // rule weights summing to 1).
    CHECK(stats.max_shift_P < 4.2e-3);
    CHECK(stats.max_shift_g < 19);
}

TEST_CASE(
    "Convergent Quadrature Gradient and Hessian",
    "[esp_potential], [esp_potential_3d]")
{
    TriMeshData data = load_wrapped_sphere();
    Eigen::MatrixXd& V = data.V;
    CollisionMesh& mesh = data.mesh;

    const double dbar_factor = GENERATE(0.3, 1.0, 2.0);
    // Keep the edge-edge support min(dbar, dhat) at 0.15 so the active contact
    // set is comparable across dbar_factor values.
    const double dhat = 0.15 / std::min(dbar_factor, 1.0);
    CAPTURE(dbar_factor);
    ESPParameters params(dhat, dbar_factor, 0);

    const bool use_near_far = GENERATE(true, false);
    CAPTURE(use_near_far, dbar_factor);
    ESPPotential potential(params, use_near_far);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);
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
                Eigen::MatrixXd V_ = V + tests::unflatten(test_dir, 3) * y(0);
                ESPCollisions collisions_;
                collisions_.build(mesh, V_, params);
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
                Eigen::MatrixXd V_ = V + tests::unflatten(test_dir, 3) * y(0);
                ESPCollisions collisions_;
                collisions_.build(mesh, V_, params);
                return potential.gradient(collisions_, mesh, V_);
            },
            fh, fd::AccuracyOrder::FOURTH, 1e-5);

        REQUIRE(
            (fh.col(0) - h * test_dir).norm()
            < std::max(fh.norm() * 1e-6, 1e-9));
    }
}

TEST_CASE(
    "Convergent Quadrature Gradient and Hessian Expensive",
    "[.][esp_potential][esp_potential_3d]")
{
    TriMeshData data = load_wrapped_sphere();
    Eigen::MatrixXd& V = data.V;
    CollisionMesh& mesh = data.mesh;

    const double dhat = 0.1;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    ESPParameters params(dhat, dbar_factor, 0);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);

    const bool normalize_weights = GENERATE(true, false);
    ESPPotential potential(params, normalize_weights);

    SECTION("gradient")
    {
        Eigen::VectorXd g = potential.gradient(collisions, mesh, V);

        Eigen::VectorXd fg;
        fd::finite_gradient(
            tests::flatten(V),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = tests::unflatten(y, 3);
                ESPCollisions collisions_;
                collisions_.build(mesh, V_, params);
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
            tests::flatten(V),
            [&](const Eigen::VectorXd& y) {
                Eigen::MatrixXd V_ = tests::unflatten(y, 3);
                ESPCollisions collisions_;
                collisions_.build(mesh, V_, params);
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
    auto [V, E, F, mesh] =
        load_triangle_mesh((tests::DATA_DIR / "sphere.obj").string());

    const double dhat = 0.2;
    const double dbar_factor = GENERATE(1.0, 0.9);
    ESPParameters params(dhat, dbar_factor, 0);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);

    ESPPotential potential(params);
    double val = potential(collisions, mesh, V);
    REQUIRE(val == 0);

    auto g = potential.gradient(collisions, mesh, V);
    REQUIRE(g.norm() == 0);

    auto H = potential.hessian(collisions, mesh, V);
    REQUIRE(H.norm() == 0);
}

TEST_CASE(
    "Convergent Quadrature Vertex Hessian",
    "[esp_potential], [esp_potential_3d]")
{
    const auto method = make_default_broad_phase();
    TriMeshData data = load_wrapped_sphere();
    Eigen::MatrixXd& V = data.V;
    CollisionMesh& mesh = data.mesh;

    const double dhat = 0.15;
    const double dbar_factor = GENERATE(1.0, 0.7, 0.4, 0.1);
    ESPParameters params(dhat, dbar_factor, 0);

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
                    V, *collisions, params);
            indices = collisions->dofs();

            if (local_grad.norm() < 1e-10) {
                continue;
            }
        }

        Eigen::MatrixXd h = PointPotentialHelper::
            evaluate_potential_hessian_at_vertex_with_cached_collisions(
                V, *collisions, params, PSDProjectionMethod::NONE);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            Eigen::VectorXd(tests::flatten(V))(indices),
            [&](const Eigen::VectorXd& y) {
                Eigen::VectorXd y_ = tests::flatten(V);
                y_(indices) = y;
                Eigen::MatrixXd V_fd = tests::unflatten(y_, 3);

                return PointPotentialHelper::
                    evaluate_potential_gradient_at_vertex_with_cached_collisions(
                        V_fd, *collisions, params);
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
    Eigen::MatrixXd V_single;
    Eigen::MatrixXi E_single, F_single;
    REQUIRE(tests::load_mesh("cube.ply", V_single, E_single, F_single));

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
    ESPParameters params(dhat, dbar_factor, 0);

    Candidates candidates;
    candidates.build(mesh, V, dhat / 2, method.get(), true);
    candidates.convert_candidates_to_sets();

    ESPCollisions collisions;
    collisions.build(candidates, mesh, V, params);

    REQUIRE(!collisions.empty());
    REQUIRE(!has_intersections(mesh, V));

    ESPPotential potential(params);
    double energy = potential(collisions, mesh, V);
    CAPTURE(energy);
    CHECK(energy > 0);
    CHECK(std::isfinite(energy));

    // Test gradient accuracy
    Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        tests::flatten(V),
        [&](const Eigen::VectorXd& x) {
            return potential(collisions, mesh, tests::unflatten(x, V.cols()));
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
    "ESP vertex quadrature ignores quad_order in 3D",
    "[esp_potential], [esp_potential_3d]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("two-cubes-close.ply", V, E, F));
    const CollisionMesh mesh(V, E, F);
    const double dhat = 0.1;

    // Without a face quadrature rule, 3D integrates over the face vertices;
    // quad_order only sets the 2D edge quadrature.
    const auto potential_value = [&](const int quad_order) {
        const ESPParameters params(dhat, 1.0, quad_order);
        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        return ESPPotential(params)(collisions, mesh, V);
    };
    const double expected = potential_value(0);
    REQUIRE(expected > 0);
    CHECK(potential_value(1) == Catch::Approx(expected).epsilon(1e-12));
    CHECK(potential_value(2) == Catch::Approx(expected).epsilon(1e-12));
}

TEST_CASE(
    "ESP potential without edge-edge terms",
    "[esp_potential], [esp_potential_3d]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("two-cubes-close.ply", V, E, F));
    const CollisionMesh mesh(V, E, F);

    // dbar_factor = 0 drops the edge-edge terms, so every face's total weight
    // is that of its rule, which sums to 1, and the normalized potential
    // equals the unnormalized one.
    const ESPParameters params(0.15, 0.0);
    ESPCollisions collisions;
    collisions.build(mesh, V, params);

    const ESPPotential normalized(params, /*use_near_far=*/true);
    const ESPPotential unnormalized(params, /*use_near_far=*/false);

    const double energy = unnormalized(collisions, mesh, V);
    REQUIRE(energy > 0);
    CHECK(
        normalized(collisions, mesh, V)
        == Catch::Approx(energy).epsilon(1e-12));

    const Eigen::VectorXd grad = unnormalized.gradient(collisions, mesh, V);
    CHECK(
        (normalized.gradient(collisions, mesh, V) - grad).norm()
        <= 1e-12 * grad.norm());

    const Eigen::SparseMatrix<double> hess =
        unnormalized.hessian(collisions, mesh, V);
    CHECK(
        (normalized.hessian(collisions, mesh, V) - hess).norm()
        <= 1e-12 * hess.norm());
}

TEST_CASE(
    "ESP potential on a 3D mesh without faces",
    "[esp_potential], [esp_potential_3d]")
{
    // Two crossing segments 0.05 apart: an edge-edge pair but no faces. The
    // segments are rods, so they must repel each other even though no face
    // carries the integral.
    Eigen::MatrixXd V(4, 3);
    V << -1, 0, 0, 1, 0, 0, 0, -1, 0.05, 0, 1, 0.05;
    Eigen::MatrixXi E(2, 2);
    E << 0, 1, 2, 3;
    const CollisionMesh mesh(V, E, Eigen::MatrixXi(0, 3));

    const ESPParameters params(0.1, 1.0);
    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    REQUIRE(!collisions.empty());

    const bool use_near_far = GENERATE(true, false);
    CAPTURE(use_near_far);
    const ESPPotential potential(params, use_near_far);
    CHECK(potential(collisions, mesh, V) > 0);

    // Pulling the rods apart lowers the energy, so the force pushes them
    // apart.
    const Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
    REQUIRE(grad.size() == V.size());
    const Eigen::MatrixXd grad_V = tests::unflatten(grad, 3);
    CHECK(grad_V(0, 2) + grad_V(1, 2) > 0);
    CHECK(grad_V(2, 2) + grad_V(3, 2) < 0);

    // The collision set changes with the positions here, so the finite
    // differences rebuild it.
    const auto energy = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd V_ = tests::unflatten(x, 3);
        ESPCollisions collisions_;
        collisions_.build(mesh, V_, params);
        return potential(collisions_, mesh, V_);
    };
    const auto gradient = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd V_ = tests::unflatten(x, 3);
        ESPCollisions collisions_;
        collisions_.build(mesh, V_, params);
        return potential.gradient(collisions_, mesh, V_);
    };

    Eigen::VectorXd fgrad;
    fd::finite_gradient(tests::flatten(V), energy, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));

    const Eigen::MatrixXd hess = potential.hessian(collisions, mesh, V);
    REQUIRE(hess.rows() == V.size());
    REQUIRE(hess.cols() == V.size());
    Eigen::MatrixXd fhess;
    fd::finite_jacobian(tests::flatten(V), gradient, fhess);
    CHECK(fd::compare_hessian(hess, fhess, 1e-3));
}

// A vertex beside the boundary edge of an open sheet sees the sheet's face
// reduced to that edge, and nothing else: a boundary edge has weight zero.
// With the closed-mesh signs, the edge's -1 cancels the face's +1 and the
// vertex feels nothing at any distance.
TEST_CASE(
    "ESP potential beside an open sheet's boundary",
    "[esp_potential], [esp_potential_3d]")
{
    // A unit-square sheet in z = 0 and a tet whose apex approaches the
    // sheet's boundary edge x = 0 from outside, in the sheet's plane.
    Eigen::MatrixXd V(8, 3);
    V << 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,                   // sheet
        0, 0.5, 0, -1, 0.2, -0.3, -1, 0.8, -0.3, -1, 0.5, 0.4; // tet
    Eigen::MatrixXi F(6, 3);
    F << 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 4, 7, 5, 5, 7, 6;
    Eigen::MatrixXi E;
    igl::edges(F, E);
    const CollisionMesh mesh(V, E, F);

    const ESPParameters params(0.1);
    const ESPPotential potential(params);
    const auto energy_at = [&](const double h) {
        Eigen::MatrixXd V_ = V;
        V_(4, 0) = -h;
        ESPCollisions collisions;
        collisions.build(mesh, V_, params);
        return potential(collisions, mesh, V_);
    };

    const double h = GENERATE(0.05, 0.01, 0.001);
    CAPTURE(h);
    const double e = energy_at(h);
    CHECK(e > 0);
    CHECK(energy_at(h / 2) > e); // repulsive: rises as the apex approaches
}

// A vertex resting exactly on another face's interior is an exact contact, so
// the energy must be +inf (which a line search rejects), never NaN (which
// slips past an E_new > E_old check).
TEST_CASE(
    "ESP potential is infinite at an exact contact",
    "[esp_potential], [esp_potential_3d]")
{
    // A unit-square sheet in z = 0 and a tet whose apex touches the interior
    // of the sheet's first face.
    Eigen::MatrixXd V(8, 3);
    V << 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,            // sheet
        0.6, 0.3, 0, 0.3, 0, 1, 0.9, 0, 1, 0.6, 0.6, 1; // tet
    Eigen::MatrixXi F(6, 3);
    F << 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 4, 7, 5, 5, 7, 6;
    Eigen::MatrixXi E;
    igl::edges(F, E);
    const CollisionMesh mesh(V, E, F);

    const ESPParameters params(0.1);
    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    REQUIRE(!collisions.empty());

    const bool use_near_far = GENERATE(true, false);
    CAPTURE(use_near_far);
    const double e = ESPPotential(params, use_near_far)(collisions, mesh, V);
    CHECK(std::isinf(e));
    CHECK(e > 0);
}

// An isolated vertex is its own integration element with a single sample, so
// with unit weights its ESP is exactly the pointwise potential of the rest of
// the mesh at its position. We cross-check that against ArbitraryPointESP on a
// mesh of the rest alone.
TEST_CASE(
    "ESP potential of an isolated vertex matches ArbitraryPointESP",
    "[esp_potential], [esp_potential_3d]")
{
    const std::string target_name =
        GENERATE("triangle interior", "open sheet boundary", "rod");
    const std::string quad_rule = GENERATE("none", "vertices+centroid");
    const bool use_near_far = GENERATE(true, false);
    CAPTURE(target_name, quad_rule, use_near_far);

    // Each target keeps its vertices and, under a face quadrature rule, its
    // face quadrature points more than dhat from q, so only q's own sample
    // sees contact.
    Eigen::MatrixXd V_target;
    Eigen::MatrixXi E_target, F_target(0, 3);
    Eigen::RowVector3d q;
    if (target_name == "triangle interior") {
        V_target.resize(3, 3);
        V_target << 0, 0, 0, 1, 0, 0, 0, 1, 0;
        F_target.resize(1, 3);
        F_target << 0, 1, 2;
        q << 0.15, 0.6, 0.05;
    } else if (target_name == "open sheet boundary") {
        V_target.resize(4, 3);
        V_target << 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0;
        F_target.resize(2, 3);
        F_target << 0, 1, 2, 0, 2, 3;
        q << -0.05, 0.5, 0;
    } else { // a bent rod: two edges in no face
        V_target.resize(3, 3);
        V_target << 0, 0, 0, 1, 0, 0, 2, 0.2, 0;
        E_target.resize(2, 2);
        E_target << 0, 1, 1, 2;
        q << 0.5, 0.05, 0;
    }
    if (F_target.rows() > 0) {
        igl::edges(F_target, E_target);
    }

    ESPParameters params(0.1, 0.2, 1, /*area_weights=*/false);
    if (quad_rule != "none") {
        params.set_quad_rule(make_face_quad_rule(quad_rule));
    }

    // The full mesh: the target plus the isolated vertex q.
    const index_t qi = V_target.rows();
    Eigen::MatrixXd V(qi + 1, 3);
    V << V_target, q;
    const CollisionMesh mesh(V, E_target, F_target);
    REQUIRE(mesh.is_codim_vertex(qi));

    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    const ESPPotential potential(params, use_near_far);

    const CollisionMesh target(V_target, E_target, F_target);
    ArbitraryPointESP<3> pointwise(target, params);
    pointwise.update(V_target);
    const auto [expected, expected_grad, expected_hess] =
        pointwise.evaluate(V_target, q);
    REQUIRE(expected > 0);

    CHECK(
        potential(collisions, mesh, V)
        == Catch::Approx(expected).epsilon(1e-12));

    const Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
    CHECK(
        (tests::unflatten(grad, 3).row(qi).transpose() - expected_grad).norm()
        <= 1e-10 * expected_grad.norm());

    const Eigen::MatrixXd hess = potential.hessian(collisions, mesh, V);
    const std::array<int, 3> q_dofs = tests::vertex_dofs<3>(qi, V.rows());
    CHECK(
        (hess(q_dofs, q_dofs) - expected_hess).norm()
        <= 1e-10 * expected_hess.norm());
}

TEST_CASE("ESP collisions accessors", "[esp_potential], [esp_potential_3d]")
{
    SECTION("edge-edge")
    {
        // Edges (0,1) and (4,5) cross 0.01 apart.
        Eigen::MatrixXd V;
        Eigen::MatrixXi E, F;
        build_ee_limit_geometry(0.25, V, E, F);
        const CollisionMesh mesh(V, E, F);
        const ESPParameters params(0.1, 1.0, 0);

        ESPCollisions collisions;
        CHECK(collisions.empty());
        CHECK(collisions.size() == 0);
        CHECK(
            collisions.compute_minimum_distance(mesh, V)
            == std::numeric_limits<double>::infinity());
        CHECK(collisions.to_string(mesh, V, params).empty());

        collisions.build(mesh, V, params);
        REQUIRE(!collisions.empty());
        CHECK(collisions.size() > 0);
        CHECK(collisions.n_candidates() > 0);
        // compute_minimum_distance returns the squared distance.
        CHECK(
            collisions.compute_minimum_distance(mesh, V)
            == Catch::Approx(1e-4));
        CHECK(
            collisions.to_string(mesh, V, params).find("edge [")
            != std::string::npos);

        collisions.clear();
        CHECK(collisions.empty());
        CHECK(collisions.size() == 0);
    }

    SECTION("vertex and face")
    {
        // Two unit cubes 0.001 apart with aligned vertices.
        Eigen::MatrixXd V_cube;
        Eigen::MatrixXi E_cube, F_cube;
        REQUIRE(tests::load_mesh("cube.ply", V_cube, E_cube, F_cube));
        const int n = V_cube.rows();
        Eigen::MatrixXd V(2 * n, 3);
        V << V_cube, V_cube.rowwise() + Eigen::RowVector3d(1.001, 0, 0);
        Eigen::MatrixXi F(2 * F_cube.rows(), 3);
        F << F_cube, (F_cube.array() + n).matrix();
        Eigen::MatrixXi E;
        igl::edges(F, E);
        const CollisionMesh mesh(V, E, F);

        ESPParameters params(0.5, 1.0, 0);
        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        REQUIRE(!collisions.empty());
        CHECK(
            collisions.to_string(mesh, V, params).find("vert [")
            != std::string::npos);
        CHECK(
            collisions.compute_minimum_distance(mesh, V)
            == Catch::Approx(1e-6));

        params.set_quad_rule(make_face_quad_rule("vertices+centroid"));
        collisions.build(mesh, V, params);
        REQUIRE(!collisions.empty());
        CHECK(
            collisions.to_string(mesh, V, params).find("face [")
            != std::string::npos);
    }
}

TEST_CASE(
    "ESP collisions skip adjacent edge-edge candidates",
    "[esp_potential], [esp_potential_3d]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    build_ee_limit_geometry(0.25, V, E, F);
    const CollisionMesh mesh(V, E, F);
    const ESPParameters params(0.1, 1.0, 0);

    // Edges 0 and 1 share vertex 0; broad phases never emit such a pair.
    REQUIRE(E(0, 0) == E(1, 0));
    Candidates candidates;
    candidates.ee_candidates.emplace_back(0, 1);
    candidates.has_all_types = true; // We fill every type we need by hand.
    ESPCollisions collisions;
    collisions.build(candidates, mesh, V, params);
    CHECK(collisions.empty());

    // Edges (0,1) and (4,5) cross, so their pair is kept.
    REQUIRE(E(6, 0) == 4);
    REQUIRE(E(6, 1) == 5);
    candidates.ee_candidates = { EdgeEdgeCandidate(0, 6) };
    collisions.build(candidates, mesh, V, params);
    CHECK(!collisions.empty());
}

TEST_CASE(
    "ESP contact is unaffected by distant codim. primitives",
    "[esp_potential], [esp_potential_3d]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    build_apex_to_apex_tets(0.05, V, E, F);
    const ESPParameters params(0.1);
    const ESPPotential potential(params);

    const CollisionMesh clean_mesh(V, E, F);
    ESPCollisions collisions;
    collisions.build(clean_mesh, V, params);
    const double expected_energy = potential(collisions, clean_mesh, V);
    const Eigen::VectorXd expected_grad =
        potential.gradient(collisions, clean_mesh, V);
    REQUIRE(expected_energy > 0);

    // Three stray vertices far from the tets, optionally with two of them
    // joined into a dangling edge.
    const bool add_codim_edge = GENERATE(false, true);
    CAPTURE(add_codim_edge);

    Eigen::MatrixXd V_codim(V.rows() + 3, 3);
    V_codim << V, 50, 50, 50, -50, -50, -50, -49, -50, -50;
    Eigen::MatrixXi E_codim = E;
    if (add_codim_edge) {
        E_codim.conservativeResize(E.rows() + 1, Eigen::NoChange);
        E_codim.bottomRows<1>() << 9, 10;
    }
    const CollisionMesh mesh(V_codim, E_codim, F);
    REQUIRE(mesh.num_codim_vertices() > 0);
    REQUIRE(mesh.num_codim_edges() == (add_codim_edge ? 1 : 0));

    collisions.build(mesh, V_codim, params);
    CHECK(
        potential(collisions, mesh, V_codim)
        == Catch::Approx(expected_energy).epsilon(1e-12));

    // We compare per-vertex rows so the check holds in either DOF layout.
    const Eigen::MatrixXd grad =
        tests::unflatten(potential.gradient(collisions, mesh, V_codim), 3);
    CHECK(
        (grad.topRows(V.rows()) - tests::unflatten(expected_grad, 3)).norm()
        <= 1e-12 * expected_grad.norm());
    CHECK(grad.bottomRows(3).isZero());
}

TEST_CASE(
    "ESP collisions reject 3D candidates without all types",
    "[esp_potential], [esp_potential_3d]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    build_apex_to_apex_tets(0.05, V, E, F);
    const CollisionMesh mesh(V, E, F);
    const ESPParameters params(0.1);
    const double inflation_radius = params.dhat / 2;

    ESPCollisions collisions;

    SECTION("IPC candidates")
    {
        Candidates candidates;
        candidates.build(mesh, V, inflation_radius);
        CHECK_THROWS_AS(
            collisions.build(candidates, mesh, V, params),
            std::invalid_argument);
    }

    SECTION("Trajectory candidates")
    {
        Candidates candidates;
        candidates.build(mesh, V, V, inflation_radius);
        CHECK_THROWS_AS(
            collisions.build(candidates, mesh, V, params),
            std::invalid_argument);
    }

    SECTION("All-types candidates")
    {
        Candidates candidates;
        candidates.build(mesh, V, inflation_radius, nullptr, true);
        collisions.build(candidates, mesh, V, params);

        ESPCollisions expected;
        expected.build(mesh, V, params);
        const ESPPotential potential(params);
        CHECK(
            potential(collisions, mesh, V)
            == Catch::Approx(potential(expected, mesh, V)).epsilon(1e-12));
    }
}

TEST_CASE(
    "ESP collisions describe and count their 2D edge quadrature contacts",
    "[esp_potential][esp_potential_2d]")
{
    // Two unit squares 0.2 apart, well within dhat of each other.
    Eigen::MatrixXd V(8, 2);
    V << -1., 1., -1., 0., -.1, 0., -.1, 1., .1, 1., .1, 0., 1., 0., 1., 1.;
    Eigen::MatrixXi E(8, 2), F;
    E << 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4;
    const CollisionMesh mesh(V, E, F);
    const ESPParameters params(0.6, 1., 2);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    REQUIRE(!collisions.empty());

    const std::string description = collisions.to_string(mesh, V, params);
    // One line per collision pair, each starting on a new line.
    CHECK(
        static_cast<size_t>(
            std::count(description.begin(), description.end(), '\n'))
        == collisions.size());
    CHECK(description.find("edge [") != std::string::npos);

    // Every stored pair was counted when it was built.
    CHECK(collisions.num_quadrature_collision_pairs >= collisions.size());
}

TEST_CASE("ESP collision dict index out of range", "[esp_potential]")
{
    ESPCollisionDict<PointType::VERTEX> dict;
    CHECK(dict.size() == 0);
    CHECK_THROWS_AS(dict[0], std::runtime_error);
    const auto& const_dict = dict;
    CHECK_THROWS_AS(const_dict[0], std::runtime_error);
}

TEST_CASE("ESP collision stencils", "[esp_potential]")
{
    constexpr double NONE = std::numeric_limits<double>::max();

    SECTION("3D")
    {
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.2, 0.2, 0.5;
        Eigen::MatrixXi F(1, 3);
        F << 0, 1, 2;
        Eigen::MatrixXi E;
        igl::edges(F, E);
        const CollisionMesh mesh(V, E, F);

        const ESPCollisionTemplate<Vertex3, Vertex3> vv(3, 0, mesh);
        const ESPCollisionTemplate<Edge3P1, Vertex3> ev(0, 3, mesh);
        const ESPCollisionTemplate<Face3P1, Vertex3> fv(0, 3, mesh);
        CHECK(vv.name() == "vv_3d");
        CHECK(ev.name() == "ev_3d");
        CHECK(fv.name() == "fv_3d");
        CHECK(ev.n_vertices_a() == 2);

        CHECK(
            vv.compute_distance(V)
            == Catch::Approx((V.row(3) - V.row(0)).squaredNorm()));
        CHECK(
            ev.compute_distance(V)
            == Catch::Approx(point_edge_distance(
                V.row(3).transpose(), V.row(E(0, 0)).transpose(),
                V.row(E(0, 1)).transpose())));
        CHECK(fv.compute_distance(V) == Catch::Approx(0.25));
        // Without vertex 3 there is no distance to report.
        for (const ESPCollision* c :
             std::initializer_list<const ESPCollision*> { &vv, &ev, &fv }) {
            CHECK(c->compute_distance(V.topRows(3)) == NONE);
        }

        CHECK_THROWS_AS(vv[2], std::runtime_error);
        CHECK_THROWS_AS(
            vv.dof(Eigen::MatrixXd::Zero(4, 4)), std::runtime_error);
        const Eigen::MatrixXd Z22 = Eigen::MatrixXd::Zero(2, 2);
        const Eigen::MatrixXd Z23 = Eigen::MatrixXd::Zero(2, 3);
        CHECK_THROWS_AS(VertexMatrixView<3>(Z22), std::runtime_error);
        CHECK_THROWS_AS(
            VertexMatrixView<3>(Z23, Z22.topRows(1)), std::runtime_error);
    }

    SECTION("2D")
    {
        Eigen::MatrixXd V(3, 2);
        V << 0, 0, 1, 0, 0.5, 0.3;
        Eigen::MatrixXi E(1, 2);
        E << 0, 1;
        const CollisionMesh mesh(V, E, Eigen::MatrixXi());

        const ESPCollisionTemplate<Vertex2, Vertex2> vv(2, 0, mesh);
        const ESPCollisionTemplate<Vertex2, Edge2P1> ve(2, 0, mesh);
        CHECK(vv.name() == "vv_2d_pt");
        CHECK(ve.name() == "ev_2d_pt");
        CHECK(vv.compute_distance(V) == Catch::Approx(0.34));
        CHECK(ve.compute_distance(V) == Catch::Approx(0.09));
        CHECK(vv.compute_distance(V.topRows(2)) == NONE);
        CHECK(ve.compute_distance(V.topRows(2)) == NONE);

        Eigen::VectorXd expected(6);
        expected << 0.5, 0.3, 0, 0, 1, 0;
        CHECK(ve.dof(V) == expected);
    }
}

TEST_CASE("ESP edge-edge distance helpers", "[esp_potential]")
{
    const Eigen::Vector3d ea0(0, 0, 0), ea1(1, 0, 0);

    // Parallel edges 0.1 apart, where the line-line distance divides by zero.
    const Eigen::Vector3d eb0(0.25, 0.1, 0), eb1(0.75, 0.1, 0);
    CHECK(
        edge_edge_distance_parallel_safe(
            ea0, ea1, eb0, eb1, EdgeEdgeDistanceType::EA_EB)
        == Catch::Approx(0.01));

    // Nearly parallel edges that the exact classifier calls EA_EB.
    const Eigen::Vector3d fb0(0.25, -1e-12, 0.1), fb1(0.75, 1e-12, 0.1);
    const EdgeEdgeDistanceType dtype =
        edge_edge_distance_type_exact(ea0, ea1, fb0, fb1);
    CHECK(dtype == EdgeEdgeDistanceType::EA_EB);
    CHECK(
        edge_edge_distance_parallel_safe(ea0, ea1, fb0, fb1, dtype)
        == Catch::Approx(0.01));

    // The closest point on edge A is clamped to the edge.
    const Eigen::Vector3d gb0(2, -1, 1), gb1(2, 1, 1);
    const Eigen::Vector3d hb0(-1, -1, 1), hb1(-1, 1, 1);
    CHECK(closest_point_uv<double>(ea0, ea1, gb0, gb1) == 1.0);
    CHECK(closest_point_uv<double>(ea0, ea1, hb0, hb1) == 0.0);
}

// Edges just above the parallel threshold (sin²θ ≈ 2.5e-16) are still EA_EB
// pairs, so the builder places a virtual vertex at their closest point. These
// two pairs once put it 2.2 edge lengths off the edge (with FMA) or made
// ESPCollisions::build throw (without FMA). The exact closest-point
// parameters come from rational arithmetic on the same doubles.
TEST_CASE(
    "ESP closest point of nearly parallel edges stays on both edges",
    "[esp_potential], [esp_potential_3d]")
{
    struct NearParallelPair {
        Eigen::Vector3d ea0, ea1, eb0, eb1;
        double s, t; // exact closest-point parameters on edges a and b
    };
    const NearParallelPair p = GENERATE(
        values<NearParallelPair>({
            { Eigen::Vector3d(
                  -0.2270162012615336, -0.037097538696173923,
                  0.24897999838085089),
              Eigen::Vector3d(
                  -0.95435862566408902, -0.54571710280362784,
                  1.0400682050760843),
              Eigen::Vector3d(
                  -0.35558625789861842, -0.12718485650007386,
                  0.38892518786460573),
              Eigen::Vector3d(
                  -0.94859663947219408, -0.54186805484382239,
                  1.033908195533616),
              0.42851486912107406, 0.3086223761078891 },
            { Eigen::Vector3d(
                  -0.97676019036967865, 0.37790242165616783,
                  -3.4459492857995926),
              Eigen::Vector3d(
                  -1.3277219282265249, 0.58419418412471269,
                  -2.4629010958870214),
              Eigen::Vector3d(
                  -1.0428448425578856, 0.41555483262155285,
                  -3.2623999701070976),
              Eigen::Vector3d(
                  -1.3887808876763006, 0.61889254451729769,
                  -2.2934288586893583),
              0.37398578721715675, 0.1899773688708917 },
        }));

    REQUIRE(!is_parallel_edge_edge(p.ea0, p.ea1, p.eb0, p.eb1));
    REQUIRE(
        edge_edge_distance_type_exact(p.ea0, p.ea1, p.eb0, p.eb1)
        == EdgeEdgeDistanceType::EA_EB);

    const Eigen::Vector2d st =
        line_line_closest_point_pairs_uv<double>(p.ea0, p.ea1, p.eb0, p.eb1);
    CHECK(st(0) == Catch::Approx(p.s).margin(1e-3));
    CHECK(st(1) == Catch::Approx(p.t).margin(1e-3));
    CHECK(
        closest_point_uv<double>(p.ea0, p.ea1, p.eb0, p.eb1)
        == Catch::Approx(p.s).margin(1e-3));

    // Both pairs are within the edge-edge support, so the builder places a
    // virtual vertex at their closest point.
    Eigen::MatrixXd V(4, 3);
    V << p.ea0.transpose(), p.ea1.transpose(), p.eb0.transpose(),
        p.eb1.transpose();
    Eigen::MatrixXi E(2, 2);
    E << 0, 1, 2, 3;
    const CollisionMesh mesh(V, E, Eigen::MatrixXi());

    const ESPParameters params(0.01);
    ESPCollisions collisions;
    REQUIRE_NOTHROW(collisions.build(mesh, V, params));
    CHECK(!collisions.empty());

    const ESPPotential potential(params);
    const double energy = potential(collisions, mesh, V);
    CHECK(std::isfinite(energy));
    CHECK(energy > 0);
    CHECK(potential.gradient(collisions, mesh, V).allFinite());
}

TEST_CASE(
    "ESP closest point of random nearly parallel edges is accurate",
    "[esp_potential], [esp_potential_3d]")
{
    const double scale = GENERATE(0.01, 1.0, 100.0);
    CAPTURE(scale);

    // We build each pair around a known closest pair: P on edge a at s and
    // Q on edge b at t, separated along the common normal by less than the
    // edge-edge support.
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> uniform(0, 1);
    std::normal_distribution<double> gaussian(0, 1);
    const auto random_vector = [&] {
        return Eigen::Vector3d(gaussian(rng), gaussian(rng), gaussian(rng));
    };
    const double ee_support = ESPParameters(0.01 * scale).ee_support();

    int num_tested = 0;
    for (int i = 0; i < 2000; i++) {
        const Eigen::Vector3d dir_a = random_vector().normalized();
        Eigen::Vector3d perp = random_vector();
        perp = (perp - perp.dot(dir_a) * dir_a).normalized();
        const double sin2 = 2.5e-16 + 5.5e-16 * uniform(rng);
        const double theta = std::asin(std::sqrt(sin2));
        const Eigen::Vector3d dir_b =
            std::cos(theta) * dir_a + std::sin(theta) * perp;

        const Eigen::Vector3d u = scale * (0.5 + uniform(rng)) * dir_a;
        const Eigen::Vector3d v = scale * (0.5 + uniform(rng)) * dir_b;
        const double s = 0.1 + 0.8 * uniform(rng);
        const double t = 0.1 + 0.8 * uniform(rng);
        const double gap = ee_support * (0.05 + 0.9 * uniform(rng));

        const Eigen::Vector3d ea0 = scale * random_vector();
        const Eigen::Vector3d ea1 = ea0 + u;
        const Eigen::Vector3d eb0 =
            ea0 + s * u + gap * dir_a.cross(dir_b).normalized() - t * v;
        const Eigen::Vector3d eb1 = eb0 + v;

        // Only pairs the builder keeps reach the closest-point solve.
        if (is_parallel_edge_edge(ea0, ea1, eb0, eb1)
            || edge_edge_distance_type_exact(ea0, ea1, eb0, eb1)
                != EdgeEdgeDistanceType::EA_EB) {
            continue;
        }
        num_tested++;

        CAPTURE(i, sin2, s, t);
        const Eigen::Vector2d st =
            line_line_closest_point_pairs_uv<double>(ea0, ea1, eb0, eb1);
        REQUIRE(st.allFinite());
        CHECK(st(0) == Catch::Approx(s).margin(0.05));
        CHECK(st(1) == Catch::Approx(t).margin(0.05));
    }
    CHECK(num_tested > 1000);
}

TEST_CASE("ESP parameters dbar_factor range", "[esp_potential]")
{
    CHECK_NOTHROW(ESPParameters(0.1, 0.0));
    CHECK_NOTHROW(ESPParameters(0.1, 0.5));
    CHECK_NOTHROW(ESPParameters(0.1, 1.0));
    CHECK_NOTHROW(ESPParameters(0.1, 2.0));
    CHECK_THROWS_AS(ESPParameters(0.1, -0.1), std::invalid_argument);
    CHECK_THROWS_AS(ESPParameters(0.1, 2.5), std::invalid_argument);
    CHECK_THROWS_AS(ESPParameters(0.1, std::nan("")), std::invalid_argument);
}

TEST_CASE("ESP parameters reject an invalid dhat", "[esp_potential]")
{
    // Zero dhat would silently disable contact, and a negative or non-finite
    // one would give the broad phase a negative or NaN inflation radius.
    const double dhat = GENERATE(
        0.0, -0.1, std::nan(""), std::numeric_limits<double>::infinity());
    CAPTURE(dhat);
    CHECK_THROWS_AS(ESPParameters(dhat), std::invalid_argument);
}

TEST_CASE("ESP parameters reject a negative quad_order", "[esp_potential]")
{
    // 3D ignores quad_order, so zero is accepted here; 2D builds require at
    // least 1 (see the 2D quadrature mismatch test).
    CHECK_NOTHROW(ESPParameters(0.1, 0.2, 0));
    CHECK_THROWS_AS(ESPParameters(0.1, 0.2, -1), std::invalid_argument);
}

TEST_CASE(
    "ESP 2D evaluation requires the build's quadrature order",
    "[esp_potential], [esp_potential_2d]")
{
    // Two parallel segments 0.05 apart.
    Eigen::MatrixXd V(4, 2);
    V << -1, 0, 1, 0, -1, 0.05, 1, 0.05;
    Eigen::MatrixXi E(2, 2);
    E << 0, 1, 2, 3;
    const CollisionMesh mesh = make_2d_collision_mesh(V, E);

    const ESPParameters params(0.1, 1.0, /*quad_order=*/2);
    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    REQUIRE(!collisions.empty());
    CHECK(ESPPotential(params)(collisions, mesh, V) > 0);

    // A different order has a different number of points per edge than the
    // collision sets store.
    const int other_order = GENERATE(1, 3);
    CAPTURE(other_order);
    const ESPPotential other(ESPParameters(0.1, 1.0, other_order));
    CHECK_THROWS_AS(other(collisions, mesh, V), std::invalid_argument);
    CHECK_THROWS_AS(other.gradient(collisions, mesh, V), std::invalid_argument);
    CHECK_THROWS_AS(other.hessian(collisions, mesh, V), std::invalid_argument);

    // 2D edge quadrature needs at least one point per edge.
    CHECK_THROWS_AS(
        collisions.build(mesh, V, ESPParameters(0.1, 1.0, 0)),
        std::invalid_argument);
}

TEST_CASE(
    "ESP 3D evaluation requires the build's face quadrature rule",
    "[esp_potential], [esp_potential_3d]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("two-cubes-close.ply", V, E, F));
    const CollisionMesh mesh(V, E, F);

    ESPParameters params(0.1, 1.0);
    params.set_quad_rule(make_face_quad_rule("edge-midpoints"));
    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    REQUIRE(!collisions.empty());
    CHECK(ESPPotential(params)(collisions, mesh, V) > 0);

    // The default vertex rule is what an ESPPotential constructed before
    // set_quad_rule() holds, since it copies the parameters. The other rules
    // have fewer points, more points, or the same points with other weights.
    const std::string other_rule =
        GENERATE("default", "centroid", "vertices+centroid", "reweighted");
    CAPTURE(other_rule);
    ESPParameters other_params(0.1, 1.0);
    if (other_rule == "centroid") {
        other_params.set_quad_rule({ { { { 1. / 3, 1. / 3, 1. / 3 } }, 1. } });
    } else if (other_rule == "vertices+centroid") {
        other_params.set_quad_rule(make_face_quad_rule(other_rule));
    } else if (other_rule == "reweighted") {
        FaceQuadRule rule = make_face_quad_rule("edge-midpoints");
        rule[0].weight *= 2;
        other_params.set_quad_rule(rule);
    }

    const ESPPotential other(other_params);
    CHECK_THROWS_AS(other(collisions, mesh, V), std::invalid_argument);
    CHECK_THROWS_AS(other.gradient(collisions, mesh, V), std::invalid_argument);
    CHECK_THROWS_AS(other.hessian(collisions, mesh, V), std::invalid_argument);

    // Clearing the set forgets the build's rule along with its collisions.
    collisions.clear();
    CHECK(other(collisions, mesh, V) == 0);
}

TEST_CASE("ESP face quadrature rule validation", "[esp_potential]")
{
    ESPParameters params(0.1);
    // The default is the vertex rule: the corners, weight 1/3 each, all
    // integrated with the shared vertex samples.
    CHECK(params.quad_rule().size() == 3);
    for (const double w : params.corner_weights()) {
        CHECK(w == Catch::Approx(1. / 3));
    }
    CHECK(params.non_corner_rule().empty());

    // Weights are normalized to sum to 1; corner points are split from the
    // others.
    params.set_quad_rule({
        { { { 1., 0., 0. } }, 2. },
        { { { 0., 1., 0. } }, 2. },
        { { { 0., 0., 1. } }, 2. },
        { { { 1. / 3, 1. / 3, 1. / 3 } }, 2. },
    });
    for (const double w : params.corner_weights()) {
        CHECK(w == Catch::Approx(0.25));
    }
    REQUIRE(params.non_corner_rule().size() == 1);
    CHECK(params.non_corner_rule()[0].weight == Catch::Approx(0.25));

    const FaceQuadRule centroid = { { { { 1. / 3, 1. / 3, 1. / 3 } }, 1.0 } };
    CHECK_NOTHROW(params.set_quad_rule(centroid));
    REQUIRE(params.quad_rule().size() == 1);

    // A zero weight is allowed as long as the weights sum to a positive value.
    CHECK_NOTHROW(params.set_quad_rule(
        {
            { { { 1., 0., 0. } }, 0.5 },
            { { { 0., 1., 0. } }, 0.5 },
            { { { 0., 0., 1. } }, 0.0 },
        }));
    CHECK(params.quad_rule().size() == 3);

    // Degree-3 rule with a negative centroid weight (listed last, so the check
    // must reach it): its barrier term would be unbounded below.
    const FaceQuadRule negative = {
        { { { 0.6, 0.2, 0.2 } }, 25. / 48 },
        { { { 0.2, 0.6, 0.2 } }, 25. / 48 },
        { { { 0.2, 0.2, 0.6 } }, 25. / 48 },
        { { { 1. / 3, 1. / 3, 1. / 3 } }, -27. / 48 },
    };
    CHECK_THROWS_AS(params.set_quad_rule(negative), std::invalid_argument);
    CHECK(params.quad_rule().size() == 3); // unchanged on failure

    CHECK_THROWS_AS(
        params.set_quad_rule({ { { { 1. / 3, 1. / 3, 1. / 3 } }, 0.0 } }),
        std::invalid_argument);
    for (const double bad :
         { std::nan(""), std::numeric_limits<double>::infinity() }) {
        CAPTURE(bad);
        CHECK_THROWS_AS(
            params.set_quad_rule(
                {
                    { { { 1., 0., 0. } }, 0.5 },
                    { { { 0., 1., 0. } }, bad },
                }),
            std::invalid_argument);
    }
    CHECK(params.quad_rule().size() == 3); // unchanged on failure

    // A point outside the face (negative barycentric coordinate) or with a
    // non-finite coordinate is rejected.
    for (const double bad : { -0.1, std::nan("") }) {
        CAPTURE(bad);
        CHECK_THROWS_AS(
            params.set_quad_rule(
                {
                    { { { 1., 0., 0. } }, 0.5 },
                    { { { 0.5, 0.5, bad } }, 0.5 }, // the only bad coordinate
                }),
            std::invalid_argument);
    }
    CHECK(params.quad_rule().size() == 3); // unchanged on failure

    // An empty rule is rejected; there is no implicit vertex integration.
    CHECK_THROWS_AS(params.set_quad_rule({}), std::invalid_argument);
    CHECK(params.quad_rule().size() == 3); // unchanged on failure
    CHECK_NOTHROW(params.set_quad_rule(ESPParameters::vertex_rule()));
    CHECK(params.non_corner_rule().empty());
}

TEST_CASE("ESP 2D Gauss-Lobatto weights are positive", "[esp_potential]")
{
    // quad_order n selects the (n + 1)-point Gauss-Lobatto rule, tabulated up
    // to 20 points and computed beyond, so ESPParameters accepts any order.
    for (int order = 1; order <= 25; order++) {
        CAPTURE(order);
        CHECK_NOTHROW(ESPParameters(0.1, 0.2, order));
        const GaussLobatto::Rule& rule = GaussLobatto::get_rule(order);
        REQUIRE(rule.size() == static_cast<size_t>(order + 1));
        double sum = 0;
        for (const EdgeQuadPoint& qp : rule) {
            CHECK(qp.weight > 0);
            CHECK(qp.xi >= 0);
            CHECK(qp.xi <= 1);
            sum += qp.weight;
        }
        CHECK(sum == Catch::Approx(1.0).epsilon(1e-12));
    }

    // A rule needs at least two points.
    CHECK_THROWS_AS(GaussLobatto::get_rule(0), std::runtime_error);
    std::vector<double> x, w;
    CHECK_THROWS_AS(lobatto_compute(0, x, w), std::runtime_error);
    CHECK_THROWS_AS(lobatto_set(1, x, w), std::domain_error);
}

TEST_CASE(
    "NearFarBarrier far part vanishes at alpha = 2", "[esp_potential][barrier]")
{
    const double dhat = 0.1;
    const NormalizedClampedLogBarrier<> base;
    const NearFarBarrier<> nf1(&base, 1.0), nf2(&base, 2.0);
    for (int i = 1; i < 100; i++) {
        const double d = dhat * i / 100.0;
        CAPTURE(d);
        CHECK(nf2.far_value(d, dhat) == 0.0);
        CHECK(nf2.near_value(d, dhat) == base(d, dhat));
        if (d > dhat / 2) {
            CHECK(nf1.far_value(d, dhat) > 0.0); // a real split at alpha = 1
        }
    }
}

// The near/far split is used for every dbar_factor in (0, 2] and its far part
// vanishes at 2, so the potential and its derivatives are continuous in
// dbar_factor, including at 1, where the split used to switch off.
TEST_CASE(
    "ESP potential is continuous in dbar_factor",
    "[esp_potential], [esp_potential_3d]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("two-cubes-close.ply", V, E, F));
    const CollisionMesh mesh(V, E, F);
    const double dhat = 0.15;
    const double alpha = GENERATE(1.0, 2.0);
    const double eps = 1e-7;
    CAPTURE(alpha);

    struct Eval {
        double e;
        Eigen::VectorXd g;
        Eigen::MatrixXd h;
    };
    const auto eval = [&](const double dbar_factor) {
        const ESPParameters params(dhat, dbar_factor, 0);
        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        const ESPPotential potential(params);
        return Eval { potential(collisions, mesh, V),
                      potential.gradient(collisions, mesh, V),
                      Eigen::MatrixXd(potential.hessian(collisions, mesh, V)) };
    };
    const auto check_close = [](const Eval& a, const Eval& b) {
        CHECK(a.e == Catch::Approx(b.e).epsilon(1e-5));
        CHECK((a.g - b.g).norm() <= 1e-5 * b.g.norm());
        CHECK((a.h - b.h).norm() <= 1e-5 * b.h.norm());
    };

    const Eval at = eval(alpha);
    REQUIRE(at.e > 0);
    check_close(eval(alpha - eps), at);
    if (alpha < 2) {
        check_close(eval(alpha + eps), at);
    }
}

// An edge-edge quadrature point whose collision set cancels to empty still has
// a positive weight in its faces' weighted average. Here edges e0 (on T0) and
// e1 (on T1) form an EA_EB pair whose edge-point sets cancel to empty; vertex u
// gives T0 a nonzero vertex potential. Moving vertex v across dhat of the
// closest point on e0 makes that set non-empty; the potential must not jump.
TEST_CASE(
    "ESP edge-edge weights count when their collisions cancel",
    "[esp_potential], [esp_potential_3d]")
{
    const double dhat = 0.5, d0 = 0.2;
    Eigen::MatrixXd V(8, 3);
    V << -1, 0, 0, 1, 0, 0, 0, -1, 0,      // T0, in z = 0
        0, -1, d0, 0, 1, d0, 0, 0, d0 + 1, // T1, in x = 0
        1, 0, -0.2,                        // u
        0, 0, -dhat;                       // v
    Eigen::MatrixXi F(2, 3);
    F << 0, 1, 2, 3, 4, 5;
    Eigen::MatrixXi E;
    igl::edges(F, E);
    const CollisionMesh mesh(V, E, F);

    const ESPParameters params(dhat, 1.0, 0);
    const ESPPotential potential(params);
    const auto energy_at = [&](const double z) {
        V(7, 2) = -z;
        ESPCollisions collisions;
        collisions.build(mesh, V, params);
        return potential(collisions, mesh, V);
    };

    const double eta = 1e-6;
    const double e_out = energy_at(dhat + eta);
    const double e_in = energy_at(dhat - eta);
    CAPTURE(e_out, e_in);
    REQUIRE(e_out != 0);
    CHECK(e_in == Catch::Approx(e_out).epsilon(1e-6));
}

TEST_CASE("ESP potential codim", "[esp_potential], [esp_potential_2d]")
{
    const auto method = make_default_broad_phase();
    double dhat = 2;
    const int quadrature_order = 2;
    ESPParameters params(dhat, 1., quadrature_order);

    Eigen::MatrixXd vertices(4, 2);
    Eigen::MatrixXi edges(2, 2);

    vertices << -1, 0, 0, 0, 1, 0, 1.5, 0.2;
    edges << 0, 1, 1, 2;

    CollisionMesh mesh = make_2d_collision_mesh(vertices, edges);

    ESPCollisions collisions;
    collisions.build(mesh, vertices, params, method.get());
    CAPTURE(dhat, method);
    CHECK(!collisions.empty());
    CHECK(!has_intersections(mesh, vertices));

    ESPPotential potential(params);
    double energy = potential(collisions, mesh, vertices);
    CHECK(energy > 0);

    // Gradient
    const Eigen::VectorXd grad = potential.gradient(collisions, mesh, vertices);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        tests::flatten(vertices),
        [&](const Eigen::VectorXd& x) {
            return potential(
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        },
        fgrad, fd::AccuracyOrder::SECOND, 1e-8);

    REQUIRE(grad.squaredNorm() > 1e-8);
    CHECK((grad - fgrad).norm() / grad.norm() < 1e-4);

    // Hessian
    Eigen::MatrixXd hess = potential.hessian(collisions, mesh, vertices);

    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        tests::flatten(vertices),
        [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, tests::unflatten(x, vertices.cols()));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-8);

    REQUIRE(hess.squaredNorm() > 1e-8);
    CHECK((hess - fhess).norm() / hess.norm() < 1e-3);
}

// The ends of an open polyline have weight zero and an isolated vertex has
// weight one, so a free edge alone feels nothing and an isolated vertex is
// repelled by an edge. With the closed-curve signs (every vertex -1) the free
// edge pulls itself shorter and the vertex is attracted to the edge, with an
// energy unbounded below.
// Even quadrature orders put a quadrature point at the edge's midpoint, right
// under the isolated vertex in the second test below.
TEST_CASE(
    "ESP potential 2D on open and codimensional pieces: free edge",
    "[esp_potential], [esp_potential_2d]")
{
    const int quadrature_order = GENERATE(2, 4);
    CAPTURE(quadrature_order);

    Eigen::MatrixXd V(2, 2);
    V << 0, 0, 1, 0;
    Eigen::MatrixXi E(1, 2);
    E << 0, 1;
    const CollisionMesh mesh = make_2d_collision_mesh(V, E);

    const ESPParameters params(2, 1., quadrature_order);
    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    const ESPPotential potential(params);
    CHECK(potential(collisions, mesh, V) == 0);
    CHECK(potential.gradient(collisions, mesh, V).squaredNorm() == 0);
}

TEST_CASE(
    "ESP potential 2D on open and codimensional pieces: isolated vertex above "
    "an edge",
    "[esp_potential], [esp_potential_2d]")
{
    const int quadrature_order = GENERATE(2, 4);
    const double h = GENERATE(0.1, 1e-2, 1e-5);
    CAPTURE(quadrature_order, h);

    Eigen::MatrixXd V(3, 2);
    V << 0, 0, 1, 0, 0.5, h;
    Eigen::MatrixXi E(1, 2);
    E << 0, 1;
    const CollisionMesh mesh = make_2d_collision_mesh(V, E);

    const ESPParameters params(0.5, 1., quadrature_order);
    const ESPPotential potential(params);
    const auto energy_at = [&](const double y) {
        Eigen::MatrixXd V_ = V;
        V_(2, 1) = y;
        ESPCollisions collisions;
        collisions.build(mesh, V_, params);
        return potential(collisions, mesh, V_);
    };

    const double e = energy_at(h);
    CHECK(e > 0);
    CHECK(energy_at(h / 2) > e);

    // The force on the vertex points away from the edge.
    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    CHECK(potential.gradient(collisions, mesh, V)(5) < 0);
}

namespace {

/// Check that a closed 2D polygon whose edges are all farther than dhat apart
/// gets no energy, forces, or stiffness.
void check_esp_2d_no_forces(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const int quadrature_order)
{
    CAPTURE(quadrature_order);
    const auto method = make_default_broad_phase();
    const ESPParameters params(/*dhat=*/1., 1., quadrature_order);

    CollisionMesh mesh = make_2d_collision_mesh(V, E);

    ESPCollisions collisions;
    collisions.build(mesh, V, params, method.get());

    REQUIRE(!has_intersections(mesh, V));

    ESPPotential potential(params);
    double energy = potential(collisions, mesh, V);
    CHECK(energy == 0);

    Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
    CHECK(grad.squaredNorm() == 0);

    Eigen::MatrixXd hess = potential.hessian(collisions, mesh, V);
    CHECK(hess.squaredNorm() == 0);
}

} // namespace

TEST_CASE("ESP potential 2D no forces", "[esp_potential], [esp_potential_2d]")
{
    const int quadrature_order = GENERATE(1, 14);

    Eigen::MatrixXd V;
    Eigen::MatrixXi E;
    SECTION("square_1")
    {
        V.resize(4, 2);
        E.resize(4, 2);
        V << -1., -1., 1., -1., 1., 1., -1., 1.;
        E << 0, 1, 1, 2, 2, 3, 3, 0;
    }
    SECTION("square_2")
    {
        V.resize(8, 2);
        E.resize(8, 2);
        V << -1., -1., 0., -1., 1., -1., 1., 0., 1., 1., 0., 1., -1., 1., -1.,
            0.;
        E << 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 0;
    }

    check_esp_2d_no_forces(V, E, quadrature_order);
}

TEST_CASE(
    "ESP potential 2D no forces on a circle",
    "[esp_potential], [esp_potential_2d]")
{
    const int quadrature_order = GENERATE(1, 14);
    const int n = GENERATE(5, 6, 7, 8, 9, 10, 50, 100, 200);
    CAPTURE(n);

    // A regular n-gon inscribed in the unit circle.
    Eigen::MatrixXd V(n, 2);
    Eigen::MatrixXi E(n, 2);
    for (int i = 0; i < n; i++) {
        V(i, 0) = std::cos(2 * M_PI * i / n);
        V(i, 1) = std::sin(2 * M_PI * i / n);
        E(i, 0) = i;
        E(i, 1) = (i + 1) % n;
    }

    check_esp_2d_no_forces(V, E, quadrature_order);
}

namespace {

constexpr double BA = 0; // a small constant to break perfect alignments

/// Check the 2D ESP gradient (and optionally the Hessian) against finite
/// differences.
///
/// We check the Hessian against an FD of the analytic gradient on the set
/// built at V: the potential is only C1 across a distance-type switch
/// (point-line vs. point-point distance), so rebuilding at each step is not
/// an option. The fixed set is valid as long as no stored Vertex2-Edge2P1
/// pair (P_E) is within the FD step of a switch, since that collision
/// asserts P_E. Exact ties are classified to endpoint types and are
/// harmless, but rounding can leave P_E pairs ~1e-17 from a switch; as in
/// the GCP tests, such configurations get the gradient check only.
void check_esp_2d_finite_differences(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const int quadrature_order,
    const bool check_hessian)
{
    CAPTURE(quadrature_order);
    const auto method = make_default_broad_phase();
    const ESPParameters params(/*dhat=*/0.6, 1., quadrature_order);

    CollisionMesh mesh = make_2d_collision_mesh(V, E);

    ESPCollisions collisions;
    collisions.build(mesh, V, params, method.get());

    REQUIRE(!collisions.empty());
    REQUIRE(!has_intersections(mesh, V));

    ESPPotential potential(params);
    double energy = potential(collisions, mesh, V);
    CHECK(energy > 0);

    Eigen::VectorXd grad = potential.gradient(collisions, mesh, V);
    REQUIRE(grad.squaredNorm() > 1e-8);

    // Rebuild the collisions at each perturbed configuration: a
    // collision's type encodes its distance type (e.g., Vertex2-Edge2P1
    // only at P_E) and asserts it, so the set built at V cannot be
    // evaluated across a switch, which the energy FD may cross.
    const auto energy_at = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd V_ = tests::unflatten(x, V.cols());
        ESPCollisions collisions_;
        collisions_.build(mesh, V_, params, method.get());
        return potential(collisions_, mesh, V_);
    };

    CAPTURE(grad.norm());
    if (V.size() <= 32) {
        Eigen::VectorXd fgrad;
        fd::finite_gradient(
            tests::flatten(V), energy_at, fgrad, fd::AccuracyOrder::SECOND,
            1e-8);
        CAPTURE(fgrad.norm());
        CHECK(
            (grad - fgrad).norm() < std::max(
                1e-4 * std::max({ grad.norm(), fgrad.norm(), 1e-8 }), 1e-9));
    } else {
        // A full finite difference with a rebuild per evaluation is too
        // expensive on mesh_1 (880 DOFs), verify directional derivative
        // only
        Eigen::VectorXd test_dir(V.size());
        for (int i = 0; i < test_dir.size(); i++) {
            test_dir(i) = i;
        }
        test_dir.normalize();

        Eigen::VectorXd fg;
        fd::finite_gradient(
            Eigen::VectorXd::Zero(1),
            [&](const Eigen::VectorXd& y) {
                return energy_at(tests::flatten(V) + test_dir * y(0));
            },
            fg, fd::AccuracyOrder::SECOND, 1e-8);
        const double g_dir = grad.dot(test_dir);
        CAPTURE(fg(0), g_dir);
        CHECK(
            std::abs(fg(0) - g_dir) < std::max(
                1e-4 * std::max({ std::abs(fg(0)), std::abs(g_dir), 1e-8 }),
                1e-9));
    }

    if (!check_hessian) {
        return;
    }

    Eigen::MatrixXd hess = potential.hessian(collisions, mesh, V);
    REQUIRE(hess.squaredNorm() > 1e-3);
    Eigen::MatrixXd fhess;
    fd::finite_jacobian(
        tests::flatten(V),
        [&](const Eigen::VectorXd& x) {
            return potential.gradient(
                collisions, mesh, tests::unflatten(x, V.cols()));
        },
        fhess, fd::AccuracyOrder::SECOND, 1e-12);
    CAPTURE(hess.norm());
    CAPTURE(fhess.norm());
    CHECK(
        (hess - fhess).norm()
        < std::max(3e-3 * std::max({ hess.norm(), fhess.norm(), 1e-8 }), 1e-9));
}

} // namespace

TEST_CASE(
    "ESP potential 2D finite differences",
    "[esp_potential], [esp_potential_2d]")
{
    const int quadrature_order = GENERATE(1, 14);

    Eigen::MatrixXd V;
    Eigen::MatrixXi E;

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
            // Rounding leaves P_E pairs ~1e-17 from a switch at order 14.
            check_esp_2d_finite_differences(V, E, quadrature_order, false);
        }
        SECTION("vertical_squares")
        {
            INFO("vertical_squares");
            V << 0. + BA, -1., 1. + BA, -1., 1. + BA, -.1, 0. + BA, -.1, 0., .1,
                1., .1, 1., 1., 0., 1.;
            // Rounding leaves P_E pairs ~1e-17 from a switch at order 14.
            check_esp_2d_finite_differences(V, E, quadrature_order, false);
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
        // The closest P_E pair is ~1e-11 from a switch, above h = 1e-12.
        check_esp_2d_finite_differences(V, E, quadrature_order, true);
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
        // Rounding leaves P_E pairs ~1e-17 from a switch at order 14.
        check_esp_2d_finite_differences(V, E, quadrature_order, false);
    }
}

TEST_CASE(
    "ESP potential 2D finite differences at corners",
    "[esp_potential], [esp_potential_2d]")
{
    const int quadrature_order = GENERATE(1, 14);
    const double P0x = GENERATE(.01, -.01, 0.0);
    const double P1y = GENERATE(.49, .5, .51);
    CAPTURE(P0x, P1y);

    Eigen::MatrixXd V(8, 2);
    Eigen::MatrixXi E(8, 2);
    V << -1., 1., -1., 0., 0., 0., P0x, .5 + BA, 0., 1., 1., 0., 1., 1., .02,
        P1y;
    E << 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 5, 6, 6, 7, 7, 5;
    check_esp_2d_finite_differences(V, E, quadrature_order, true);
}

// Each 2D quadrature point's Hessian is projected, so their sum is PSD.
TEST_CASE("ESP potential 2D Hessian PSD", "[esp_potential], [esp_potential_2d]")
{
    const int quadrature_order = GENERATE(1, 2, 7);
    const PSDProjectionMethod psd_method =
        GENERATE(PSDProjectionMethod::CLAMP, PSDProjectionMethod::ABS);
    const ESPParameters params(0.6, 1., quadrature_order);

    Eigen::MatrixXd V(8, 2);
    Eigen::MatrixXi E(8, 2);
    SECTION("corners")
    {
        V << -1., 1., -1., 0., 0., 0., .01, .5, 0., 1., 1., 0., 1., 1., .02,
            .49;
        E << 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 5, 6, 6, 7, 7, 5;
    }
    SECTION("squares")
    {
        V << -1., 1., -1., 0., -.1, 0., -.1, 1., .1, 1., .1, 0., 1., 0., 1., 1.;
        E << 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4;
    }
    const CollisionMesh mesh = make_2d_collision_mesh(V, E);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    REQUIRE(!collisions.empty());

    const auto [lambda_min, lambda_max] = eigenvalue_range(
        ESPPotential(params).hessian(collisions, mesh, V, psd_method));
    CHECK(lambda_min >= -std::max(1e-10, 1e-10 * std::abs(lambda_max)));
}

// 3D FACE QUADRATURE TESTS //

// Verify that face quadrature gives gradient/hessian consistent with finite
// differences on the wrapped-sphere geometry, for two positive-weight face
// quadrature rules (params.set_quad_rule()). The edge-midpoint rule puts each
// point on an edge, which skips the face across that edge.
TEST_CASE(
    "Face Quadrature Gradient and Hessian",
    "[esp_potential], [esp_potential_3d]")
{
    TriMeshData data = load_wrapped_sphere();
    Eigen::MatrixXd& V = data.V;
    CollisionMesh& mesh = data.mesh;

    const double dhat = 0.15;
    const std::string rule_name =
        GENERATE(as<std::string>(), "vertices+centroid", "edge-midpoints");
    CAPTURE(rule_name);
    ESPParameters params(dhat, 1.);
    params.set_quad_rule(make_face_quad_rule(rule_name));

    const bool normalize_weights = GENERATE(true, false);
    ESPPotential potential(params, normalize_weights);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);

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
                Eigen::MatrixXd V_ = V + tests::unflatten(test_dir, 3) * y(0);
                ESPCollisions c;
                c.build(mesh, V_, params);
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
                Eigen::MatrixXd V_ = V + tests::unflatten(test_dir, 3) * y(0);
                ESPCollisions c;
                c.build(mesh, V_, params);
                return potential.gradient(c, mesh, V_);
            },
            fh, fd::AccuracyOrder::SECOND, 1e-6);

        REQUIRE((fh.col(0) - h * test_dir).norm() < fh.norm() * 1e-4);
    }
}

// Verify that the global hessian is PSD whenever project_hessian_to_psd is set.
// With normalized weights (the second ESPPotential argument, use_near_far) and
// 0 < dbar_factor <= 2, the edge-edge weights add indefinite terms, so each
// face's block is projected as a whole; every dbar_factor in that range takes
// this path. The non-normalized branch projects each stencil's block. two-
// cubes-close has many edge-edge pairs and vertex potentials; the wrapped
// sphere has few.
TEST_CASE(
    "Convergent Quadrature Hessian PSD", "[esp_potential], [esp_potential_3d]")
{
    const std::string mesh_name =
        GENERATE(as<std::string>(), "wrapped_sphere", "two-cubes-close");
    auto [V, E, F, mesh] = mesh_name == "wrapped_sphere"
        ? load_wrapped_sphere()
        : load_triangle_mesh(
              (tests::DATA_DIR / "two-cubes-close.ply").string());

    const double dhat = 0.15;
    const double dbar_factor = 0.4;
    ESPParameters params(dhat, dbar_factor, 0);

    const bool normalize_weights = GENERATE(true, false);
    const PSDProjectionMethod psd_method =
        GENERATE(PSDProjectionMethod::CLAMP, PSDProjectionMethod::ABS);

    ESPPotential potential(params, normalize_weights);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);

    const auto [lambda_min, lambda_max] =
        eigenvalue_range(potential.hessian(collisions, mesh, V, psd_method));
    const double tol = std::max(1e-10, 1e-10 * std::abs(lambda_max));

    INFO(
        "mesh=" << mesh_name << " dbar_factor=" << dbar_factor
                << " normalize_weights=" << normalize_weights << " method="
                << static_cast<int>(psd_method) << " lambda_min=" << lambda_min
                << " lambda_max=" << lambda_max);
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
        NearFarBarrier<> nf(&base, alpha);

        constexpr int N = 200;
        for (int i = 0; i < N; ++i) {
            const double d = dhat * (i + 1.0) / N;

            const double b = base(d, dhat);
            CHECK(
                nf.near_value(d, dhat) + nf.far_value(d, dhat)
                == Catch::Approx(b));

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

            // The unsplit evaluation forwards to the base barrier.
            CHECK(nf(d, dhat) == b);
            CHECK(nf.first_derivative(d, dhat) == db);
            CHECK(nf.second_derivative(d, dhat) == ddb);

            CAPTURE(d);
            CAPTURE(alpha);
            CAPTURE(dhat);
            CAPTURE(alpha * dhat);
            CAPTURE(alpha * dhat / 2);
            // Check that near barrier is 0 above alpha*dhat and non-zero below
            constexpr double eps_tol = 1e-9;
            if (d >= alpha * dhat) {
                CHECK(nf.near_value(d, dhat) == 0.0);
            } else if (d < alpha * dhat - eps_tol) {
                CHECK(nf.near_value(d, dhat) > 0.0);
            }

            // Check that far barrier is 0 below dhat*alpha/2 and non-zero above
            if (d <= dhat * alpha / 2 || d >= dhat) {
                CHECK(nf.far_value(d, dhat) == 0.0);
            } else if (d > dhat * alpha / 2 + eps_tol) {
                CHECK(nf.far_value(d, dhat) > 0.0);
            }
        }
        CHECK(nf.units(dhat) == base.units(dhat));
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
        run_test(InversePowerBarrier<>(1.0));
        break;
    case BarrierType::InversePower2:
        run_test(InversePowerBarrier<>(2.0));
        break;
    default:
        FAIL("Unknown barrier type");
    }
}

// At and below alpha * dhat / 2 the far weight vanishes, so the near part is
// the whole barrier and the far part is exactly zero. This must hold even
// where the base barrier is +inf (contact, penetration, or InversePowerBarrier
// overflowing for tiny d), where multiplying by a zero weight would give NaN.
TEST_CASE(
    "NearFarBarrier keeps the whole barrier near at contact",
    "[esp_potential][barrier]")
{
    const double dhat = 1e-3;
    const double alpha = GENERATE(0.2, 2.0);
    const double d = GENERATE(0.0, -1e-4, 1e-200);

    const std::vector<std::shared_ptr<const Barrier>> bases {
        std::make_shared<NormalizedClampedLogBarrier<>>(),
        std::make_shared<ClampedLogBarrier<>>(),
        std::make_shared<InversePowerBarrier<>>(2.0),
        std::make_shared<CubicBarrier<>>(),
    };
    const size_t i = GENERATE(range<size_t>(0, 4));
    const Barrier& base = *bases[i];
    CAPTURE(alpha, d, i);

    const NearFarBarrier<> nf(bases[i], alpha);

    CHECK(nf.near_value(d, dhat) == base(d, dhat));
    CHECK(nf.far_value(d, dhat) == 0.0);
    CHECK(nf.first_derivative_near(d, dhat) == base.first_derivative(d, dhat));
    CHECK(nf.first_derivative_far(d, dhat) == 0.0);
    CHECK(
        nf.second_derivative_near(d, dhat) == base.second_derivative(d, dhat));
    CHECK(nf.second_derivative_far(d, dhat) == 0.0);

    // The near + far = b identity holds at contact too, including b = +inf.
    CHECK(nf.near_value(d, dhat) + nf.far_value(d, dhat) == base(d, dhat));
}

TEST_CASE(
    "NearFarBarrier requires a positive split parameter",
    "[esp_potential][barrier]")
{
    const auto base = std::make_shared<ClampedLogBarrier<>>();
    CHECK_NOTHROW(NearFarBarrier<>(base, 0.2));

    // At alpha = 0 the step collapses to a jump with NaN derivatives.
    const double alpha = GENERATE(
        0.0, -0.5, std::nan(""), std::numeric_limits<double>::infinity());
    CAPTURE(alpha);
    CHECK_THROWS_AS(NearFarBarrier<>(base, alpha), std::invalid_argument);
    CHECK_THROWS_AS(NearFarBarrier<>(base.get(), alpha), std::invalid_argument);
}

// Same check with a face quadrature rule (params.set_quad_rule()), which adds
// face-point stencils to each face's block. dbar_factor = 0 covers the
// per-stencil projection branch of the normalized potential.
TEST_CASE("Face Quadrature Hessian PSD", "[esp_potential], [esp_potential_3d]")
{
    const std::string mesh_name =
        GENERATE(as<std::string>(), "wrapped_sphere", "two-cubes-close");
    auto [V, E, F, mesh] = mesh_name == "wrapped_sphere"
        ? load_wrapped_sphere()
        : load_triangle_mesh(
              (tests::DATA_DIR / "two-cubes-close.ply").string());

    const double dhat = 0.15;
    const std::string rule_name = "vertices+centroid";
    const double dbar_factor = GENERATE(0.4, 0.0);
    ESPParameters params(dhat, dbar_factor);
    params.set_quad_rule(make_face_quad_rule(rule_name));

    const bool normalize_weights = GENERATE(true, false);
    const PSDProjectionMethod psd_method =
        GENERATE(PSDProjectionMethod::CLAMP, PSDProjectionMethod::ABS);

    ESPPotential potential(params, normalize_weights);

    ESPCollisions collisions;
    collisions.build(mesh, V, params);

    const auto [lambda_min, lambda_max] =
        eigenvalue_range(potential.hessian(collisions, mesh, V, psd_method));
    const double tol = std::max(1e-10, 1e-10 * std::abs(lambda_max));

    INFO(
        "mesh=" << mesh_name << " dbar_factor=" << dbar_factor
                << " normalize_weights=" << normalize_weights << " rule="
                << rule_name << " method=" << static_cast<int>(psd_method)
                << " lambda_min=" << lambda_min
                << " lambda_max=" << lambda_max);
    REQUIRE(lambda_min >= -tol);
}

// A block of a taller matrix, such as `V_big.topRows(n)`, binds to the
// potential's vertex argument without a copy, so its columns are further apart
// than its row count. The potential must read it exactly like the equivalent
// contiguous matrix, both when building the collisions and when evaluating.
TEST_CASE(
    "ESP potential reads a strided vertex matrix like a contiguous one",
    "[esp_potential]")
{
    const bool use_near_far = GENERATE(true, false);

    Eigen::MatrixXd V;
    CollisionMesh mesh;
    double dhat = 0.15;
    int quad_order = 0;
    FaceQuadRule face_quad_rule;
    SECTION("3D vertex quadrature")
    {
        TriMeshData data = load_wrapped_sphere();
        V = data.V;
        mesh = data.mesh;
    }
    SECTION("3D face quadrature")
    {
        TriMeshData data = load_wrapped_sphere();
        V = data.V;
        mesh = data.mesh;
        face_quad_rule = make_face_quad_rule("vertices+centroid");
    }
    SECTION("2D")
    {
        // Two rectangles 0.2 apart.
        V.resize(8, 2);
        V << -1, 1, -1, 0, -.1, 0, -.1, 1, .1, 1, .1, 0, 1, 0, 1, 1;
        Eigen::MatrixXi E(8, 2);
        E << 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4;
        mesh = make_2d_collision_mesh(V, E);
        dhat = 0.6;
        quad_order = 2;
    }
    ESPParameters params(dhat, 1., quad_order);
    if (!face_quad_rule.empty()) { // else the default vertex rule
        params.set_quad_rule(face_quad_rule);
    }

    // The rows below V hold values far from the mesh, so reading them instead
    // of V changes the result.
    Eigen::MatrixXd V_big =
        Eigen::MatrixXd::Constant(V.rows() + 7, V.cols(), 1e3);
    V_big.topRows(V.rows()) = V;
    const auto V_strided = V_big.topRows(V.rows());

    ESPCollisions collisions, strided_collisions;
    collisions.build(mesh, V, params);
    strided_collisions.build(mesh, V_strided, params);
    REQUIRE(!collisions.empty());
    CHECK(strided_collisions.size() == collisions.size());

    const ESPPotential potential(params, use_near_far);
    const double P = potential(collisions, mesh, V);
    REQUIRE(P != 0);
    const Eigen::VectorXd g = potential.gradient(collisions, mesh, V);
    const Eigen::SparseMatrix<double> H =
        potential.hessian(collisions, mesh, V);

    for (const ESPCollisions* c : { &collisions, &strided_collisions }) {
        CHECK(
            potential(*c, mesh, V_strided) == Catch::Approx(P).epsilon(1e-12));
        CHECK(
            (potential.gradient(*c, mesh, V_strided) - g).norm()
            <= 1e-12 * g.norm());
        CHECK(
            (potential.hessian(*c, mesh, V_strided) - H).norm()
            <= 1e-12 * H.norm());
    }
}

// The ESP Hessian keeps its per-collision blocks on the stack, so their size
// decides whether a caller's thread survives. macOS gives secondary threads
// (std::thread included) a 512 KiB stack, and tbb::parallel_for also runs work
// on the calling thread, so the Hessian has to fit there. We pin TBB to one
// thread to make the calling thread do all of the work, then run the Hessian
// on a pthread with exactly that stack size. Overflowing it crashes the test.
TEST_CASE(
    "ESP Hessian fits in a 512 KiB thread stack",
    "[esp_potential], [esp_potential_3d]")
{
#ifdef _WIN32
    SKIP("Uses pthreads to choose the thread's stack size");
#else
    // Two triangles 0.01 apart: vertex 3 sits over the interior of edge (0,1)
    // and edge (3,4) crosses it, so we get vertex and edge-edge stencils.
    Eigen::MatrixXd V(6, 3);
    V << 0, 0, 0, 1, 0, 0, 0, 1, 0, //
        0.5, -0.05, 0.01, 0.6, 0.6, 0.01, 0.2, 0.5, 0.01;
    Eigen::MatrixXi F(2, 3);
    F << 0, 1, 2, 3, 4, 5;
    Eigen::MatrixXi E;
    igl::edges(F, E);
    const CollisionMesh mesh(V, E, F);

    const ESPParameters params(0.1); // default 3D parameters
    ESPCollisions collisions;
    collisions.build(mesh, V, params);
    REQUIRE(!collisions.empty());
    const ESPPotential potential(params);

    struct Job {
        const ESPPotential& potential;
        const ESPCollisions& collisions;
        const CollisionMesh& mesh;
        const Eigen::MatrixXd& V;
        Eigen::SparseMatrix<double> hess;
    } job { potential, collisions, mesh, V, {} };

    const tbb::global_control single_thread(
        tbb::global_control::max_allowed_parallelism, 1);

    pthread_attr_t attr;
    REQUIRE(pthread_attr_init(&attr) == 0);
    REQUIRE(pthread_attr_setstacksize(&attr, 512 * 1024) == 0);
    pthread_t thread;
    REQUIRE(
        pthread_create(
            &thread, &attr,
            [](void* arg) -> void* {
                auto& j = *static_cast<Job*>(arg);
                j.hess = j.potential.hessian(j.collisions, j.mesh, j.V);
                return nullptr;
            },
            &job)
        == 0);
    REQUIRE(pthread_join(thread, nullptr) == 0);
    pthread_attr_destroy(&attr);

    CHECK(job.hess.nonZeros() > 0);
#endif
}
