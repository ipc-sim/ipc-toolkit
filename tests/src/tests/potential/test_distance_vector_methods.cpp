#include <ipc/config.hpp>
#include <tests/config.hpp>
#include <tests/dof_layout.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/collisions/normal/vertex_vertex.hpp>
#include <ipc/collisions/normal/edge_vertex.hpp>
#include <ipc/collisions/normal/edge_edge.hpp>
#include <ipc/collisions/normal/face_vertex.hpp>

#include <igl/edges.h>
#include <memory>

using namespace ipc;

// =============================================================================
// Helper: build the explicit distance-vector Hessian for a single unmollified
// collision so we can compare the efficient methods against it.
//
// H_t = w * [4 * f''(d²) * (J*t)(J*t)ᵀ  +  2 * f'(d²) * J*Jᵀ]
//
// where J = ∂t/∂x = [c₀I, c₁I, ..., cₙI]ᵀ  (ndof × dim)
//       t = ∑ cᵢ xᵢ   (distance vector)
//       d² = tᵀt
// =============================================================================

static MatrixMax12d build_gauss_newton_hessian(
    const NormalCollision& collision,
    Eigen::ConstRef<VectorMax12d> positions,
    const BarrierPotential& barrier_potential)
{
    const int n = collision.num_vertices();
    const int d = collision.dim(positions.size());
    const int ndof = n * d;

    VectorMax4d coeffs;
    VectorMax3d t = collision.compute_distance_vector(positions, coeffs);

    // Construct J (ndof × dim)
    MatrixMax12d J;
    J.setZero(ndof, d);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < d; j++) {
            J(i * d + j, j) = coeffs[i];
        }
    }

    // J*t ∈ ℝ^{ndof}
    VectorMax12d Jt = J * t;

    // grad = w * f'(d²) * ∇d(x) = w * f'(d²) * 2 * Jᵀt
    VectorMax12d grad = barrier_potential.gradient(collision, positions);

    // Find a nonzero index in Jt
    double w_grad_f = 0; // = w * f'(d²)
    for (int i = 0; i < ndof; i++) {
        if (std::abs(Jt[i]) > 1e-15) {
            w_grad_f = grad[i] / (2.0 * Jt[i]);
            break;
        }
    }

    // For f"(d²), use the diagonal:
    //   diag[i] = w * 4 * f"(d²) * Jt[i]² + w * 2 * f'(d²) * J_diag[i]
    // where J_diag[i] = cₖ² for the vertex k that owns index i.
    //
    // So: w * 4 * f''(d²) = (diag[i] - 2 * w_grad_f * cₖ²) / Jt[i]²
    //   for some i where Jt[i] ≠ 0.
    VectorMax12d diag_vec =
        barrier_potential.gauss_newton_hessian_diagonal(collision, positions);

    double w_hess_f = 0; // = w * f"(d²)
    for (int i = 0; i < ndof; i++) {
        if (std::abs(Jt[i]) > 1e-15) {
            int vi = i / d;
            double ci2 = coeffs[vi] * coeffs[vi];
            w_hess_f =
                (diag_vec[i] - 2.0 * w_grad_f * ci2) / (4.0 * Jt[i] * Jt[i]);
            break;
        }
    }

    MatrixMax12d H_t(ndof, ndof);
    H_t = 4.0 * w_hess_f * Jt * Jt.transpose()
        + 2.0 * w_grad_f * J * J.transpose();

    return H_t;
}

// =============================================================================
// Tests for CollisionStencil::compute_distance_vector
// =============================================================================

TEST_CASE(
    "distance vector methods match distance and gradient",
    "[stencil][distance_vector]")
{
    std::unique_ptr<CollisionStencil> collision;
    VectorMax12d positions;

    SECTION("Vertex-vertex")
    {
        Eigen::MatrixXd V(2, 3);
        V << -1, 0, 0, /**/ 1, 0, 0;
        Eigen::MatrixXi E, F;
        collision = std::make_unique<VertexVertexCandidate>(0, 1);
        positions = collision->dof(V, E, F);
    }
    SECTION("Edge-vertex")
    {
        Eigen::MatrixXd V(3, 3);
        V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0;
        Eigen::MatrixXi E(1, 2), F;
        E << 0, 1;
        collision = std::make_unique<EdgeVertexCandidate>(0, 2);
        positions = collision->dof(V, E, F);
    }
    SECTION("Edge-edge")
    {
        Eigen::MatrixXd V(4, 3);
        V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, -1, 1, /**/ 0, 1, 1;
        Eigen::MatrixXi E(2, 2), F;
        E << 0, 1, /**/ 2, 3;
        collision = std::make_unique<EdgeEdgeCandidate>(0, 1);
        positions = collision->dof(V, E, F);
    }
    SECTION("Face-vertex")
    {
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0, /**/ 0.25, 0.25, 1;
        Eigen::MatrixXi F(1, 3);
        F << 0, 1, 2;
        Eigen::MatrixXi E;
        igl::edges(F, E);
        collision = std::make_unique<FaceVertexCandidate>(0, 3);
        positions = collision->dof(V, E, F);
    }

    VectorMax3d t = collision->compute_distance_vector(positions);
    CHECK(
        t.squaredNorm()
        == Catch::Approx(collision->compute_distance(positions)));

    MatrixMax<double, 12, 3> dt_dx =
        collision->compute_distance_vector_jacobian(positions);
    VectorMax12d grad = collision->compute_distance_gradient(positions);
    CAPTURE(dt_dx * t, grad);
    CHECK((dt_dx * (2 * t)).isApprox(grad));
}

TEST_CASE(
    "compute_distance_vector with coeffs output matches separate calls",
    "[stencil][distance_vector]")
{
    Eigen::MatrixXd V(4, 3);
    V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, -1, 1, /**/ 0, 1, 1;
    Eigen::MatrixXi E(2, 2), F;
    E << 0, 1, /**/ 2, 3;

    EdgeEdgeCandidate ee(0, 1);
    VectorMax12d positions = ee.dof(V, E, F);

    VectorMax4d coeffs_separate = ee.compute_coefficients(positions);
    VectorMax3d t_separate = ee.compute_distance_vector(positions);

    VectorMax4d coeffs_joint;
    VectorMax3d t_joint = ee.compute_distance_vector(positions, coeffs_joint);

    CHECK(t_joint.isApprox(t_separate));
    CHECK(coeffs_joint.isApprox(coeffs_separate));
}

// =============================================================================
// Tests for CollisionStencil static helper methods
// =============================================================================

TEST_CASE(
    "Collision stencil distance vector diagonals match explicit computation",
    "[stencil][distance_vector]")
{
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0, /**/ 0.25, 0.25, 1;
    Eigen::MatrixXi F(1, 3);
    F << 0, 1, 2;
    Eigen::MatrixXi E;
    igl::edges(F, E);

    FaceVertexCandidate fv(0, 3);
    VectorMax12d positions = fv.dof(V, E, F);
    const int dim = fv.dim(positions.size());

    VectorMax4d coeffs;
    VectorMax3d t = fv.compute_distance_vector(positions, coeffs);

    {
        // Efficient method
        VectorMax12d diag_efficient =
            CollisionStencil::diag_distance_vector_outer(coeffs, dim);

        // Explicit: form J = ∂t/∂x, compute diag(JJᵀ)
        MatrixMax12d J = fv.compute_distance_vector_jacobian(positions);
        MatrixMax12d JJT = J * J.transpose();
        VectorMax12d diag_explicit = JJT.diagonal();

        CHECK(diag_efficient.isApprox(diag_explicit));
    }

    {
        // Efficient method
        VectorMax12d diag_efficient =
            CollisionStencil::diag_distance_vector_t_outer(coeffs, t);

        // Explicit: form J·t, compute diag((J·t)(J·t)ᵀ)
        MatrixMax12d J = fv.compute_distance_vector_jacobian(positions);
        VectorMax12d Jt = J * t;
        VectorMax12d diag_explicit = (Jt.array() * Jt.array()).matrix();

        CHECK(diag_efficient.isApprox(diag_explicit));
    }

    {
        // Some arbitrary direction vector
        VectorMax12d p(12);
        p << 0.1, -0.2, 0.3, 0.4, -0.5, 0.6, -0.7, 0.8, -0.9, 1.0, -1.1, 1.2;

        // Efficient method
        VectorMax3d pTJ_efficient =
            CollisionStencil::contract_distance_vector_jacobian(coeffs, p, dim);

        // Explicit: Jᵀp  (J is ndof×dim, so Jᵀ is dim×ndof, Jᵀp is dim×1)
        MatrixMax12d J = fv.compute_distance_vector_jacobian(positions);
        Eigen::VectorXd pTJ_explicit = J.transpose() * p;

        CHECK(pTJ_efficient.isApprox(pTJ_explicit));
    }
}

// =============================================================================
// Tests for NormalPotential::gauss_newton_hessian_diagonal (single collision)
//
// The distance-vector Hessian is a Gauss-Newton-like approximation:
//   H_t = w * [4·f''(d²)·(J·t)(J·t)ᵀ + 2·f'(d²)·JJᵀ]
// It omits terms from the derivative of the closest-point parameters.
// For VV the coefficients are constant so H_t == H_full exactly.
// For other types we verify self-consistency: the efficient diagonal must
// match the diagonal of the explicitly-constructed distance-vector Hessian.
// =============================================================================

TEST_CASE(
    "gauss_newton_hessian_diagonal matches diagonal of full hessian for VV (exact case)",
    "[potential][barrier_potential][distance_vector][gauss_newton_hessian_diagonal]")
{
    const double dhat = 1.0;
    const double kappa = 1.0;
    BarrierPotential barrier_potential(dhat, kappa);

    Eigen::MatrixXd V(2, 3);
    V << 0, 0, 0, /**/ 0.5, 0, 0;
    Eigen::MatrixXi E, F;

    VertexVertexNormalCollision collision(0, 1);
    collision.dmin = 0;
    collision.weight = 1.0;

    VectorMax12d positions = collision.dof(V, E, F);

    MatrixMax12d hess_full = barrier_potential.hessian(
        collision, positions, PSDProjectionMethod::NONE);
    VectorMax12d diag_standard = hess_full.diagonal();
    VectorMax12d diag_t =
        barrier_potential.gauss_newton_hessian_diagonal(collision, positions);

    // VV coefficients are constant → distance-vector Hessian is exact
    CHECK(diag_t.isApprox(diag_standard));
}

TEST_CASE(
    "gauss_newton_hessian_diagonal is self-consistent with explicit distance-vector Hessian",
    "[potential][barrier_potential][distance_vector][gauss_newton_hessian_diagonal]")
{
    const double dhat = 1.0;
    const double kappa = 1.0;
    BarrierPotential barrier_potential(dhat, kappa);

    std::unique_ptr<NormalCollision> collision;
    VectorMax12d positions;

    SECTION("Edge-vertex")
    {
        Eigen::MatrixXd V(3, 3);
        V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, 0.5, 0;
        Eigen::MatrixXi E(1, 2), F;
        E << 0, 1;

        collision = std::make_unique<EdgeVertexNormalCollision>(0, 2);
        collision->dmin = 0;
        collision->weight = 1.0;

        positions = collision->dof(V, E, F);
    }
    SECTION("Face-vertex")
    {
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0, /**/ 0.25, 0.25, 0.5;
        Eigen::MatrixXi F(1, 3);
        F << 0, 1, 2;
        Eigen::MatrixXi E;
        igl::edges(F, E);

        collision = std::make_unique<FaceVertexNormalCollision>(0, 3);
        collision->dmin = 0;
        collision->weight = 1.0;

        positions = collision->dof(V, E, F);
    }
    SECTION("Face-vertex with weight and dmin")
    {
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0, /**/ 0.25, 0.25, 0.5;
        Eigen::MatrixXi F(1, 3);
        F << 0, 1, 2;
        Eigen::MatrixXi E;
        igl::edges(F, E);

        collision = std::make_unique<FaceVertexNormalCollision>(0, 3);
        collision->dmin = 0.01;
        collision->weight = 2.5;

        positions = collision->dof(V, E, F);
    }
    SECTION("Edge-edge")
    {
        Eigen::MatrixXd V(4, 3);
        V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, -1, 0.5, /**/ 0, 1, 0.5;
        Eigen::MatrixXi E(2, 2), F;
        E << 0, 1, /**/ 2, 3;

        collision = std::make_unique<EdgeEdgeNormalCollision>(0, 1, 0.0);
        collision->dmin = 0;
        collision->weight = 1.0;

        positions = collision->dof(V, E, F);
    }

    // Build the explicit distance-vector Hessian matrix
    MatrixMax12d H_t =
        build_gauss_newton_hessian(*collision, positions, barrier_potential);
    VectorMax12d diag_explicit = H_t.diagonal();

    VectorMax12d diag_efficient =
        barrier_potential.gauss_newton_hessian_diagonal(*collision, positions);

    CHECK(diag_efficient.isApprox(diag_explicit));
}

// =============================================================================
// Tests for NormalPotential::gauss_newton_hessian_quadratic_form (single
// collision)
//
// Same logic: exact for VV, self-consistent with the explicit distance-vector
// Hessian matrix for other types.
// =============================================================================

TEST_CASE(
    "gauss_newton_hessian_quadratic_form matches pTHp from full hessian for VV (exact case)",
    "[potential][barrier_potential][distance_vector][gauss_newton_hessian_quadratic_form]")
{
    const double dhat = 1.0;
    const double kappa = 1.0;
    BarrierPotential barrier_potential(dhat, kappa);

    Eigen::MatrixXd V(2, 3);
    V << 0, 0, 0, /**/ 0.5, 0, 0;
    Eigen::MatrixXi E, F;

    VertexVertexNormalCollision collision(0, 1);
    collision.dmin = 0;
    collision.weight = 1.0;

    VectorMax12d positions = collision.dof(V, E, F);
    VectorMax12d p(6);
    p << 0.1, -0.2, 0.3, -0.4, 0.5, -0.6;

    MatrixMax12d hess_full = barrier_potential.hessian(
        collision, positions, PSDProjectionMethod::NONE);
    double pTHp_standard = p.dot(hess_full * p);
    double pTHp_t = barrier_potential.gauss_newton_hessian_quadratic_form(
        collision, positions, p);

    // VV: exact match
    CHECK(pTHp_t == Catch::Approx(pTHp_standard).epsilon(1e-10));
}

TEST_CASE(
    "gauss_newton_hessian_quadratic_form is self-consistent with explicit distance-vector Hessian",
    "[potential][barrier_potential][distance_vector][gauss_newton_hessian_quadratic_form]")
{
    const double dhat = 1.0;
    const double kappa = 1.0;
    BarrierPotential barrier_potential(dhat, kappa);

    std::unique_ptr<NormalCollision> collision;
    VectorMax12d positions;

    SECTION("Edge-vertex")
    {
        Eigen::MatrixXd V(3, 3);
        V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, 0.5, 0;
        Eigen::MatrixXi E(1, 2), F;
        E << 0, 1;

        collision = std::make_unique<EdgeVertexNormalCollision>(0, 2);

        positions = collision->dof(V, E, F);
    }
    SECTION("Face-vertex")
    {
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, /**/ 1, 0, 0, /**/ 0, 1, 0, /**/ 0.25, 0.25, 0.5;
        Eigen::MatrixXi F(1, 3);
        F << 0, 1, 2;
        Eigen::MatrixXi E;
        igl::edges(F, E);

        collision = std::make_unique<FaceVertexNormalCollision>(0, 3);

        positions = collision->dof(V, E, F);
    }
    SECTION("Edge-edge")
    {
        Eigen::MatrixXd V(4, 3);
        V << -1, 0, 0, /**/ 1, 0, 0, /**/ 0, -1, 0.5, /**/ 0, 1, 0.5;
        Eigen::MatrixXi E(2, 2), F;
        E << 0, 1, /**/ 2, 3;

        collision = std::make_unique<EdgeEdgeNormalCollision>(0, 1, 0.0);

        positions = collision->dof(V, E, F);
    }

    collision->dmin = 0;
    collision->weight = 1.0;

    VectorMax12d p = VectorMax12d::Random(positions.size());

    MatrixMax12d H_t =
        build_gauss_newton_hessian(*collision, positions, barrier_potential);
    double pTHp_explicit = p.dot(H_t * p);

    double pTHp_efficient =
        barrier_potential.gauss_newton_hessian_quadratic_form(
            *collision, positions, p);
    CHECK(pTHp_efficient == Catch::Approx(pTHp_explicit).epsilon(1e-10));
}

// =============================================================================
// Tests for cumulative methods on a mesh
//
// The cumulative methods assemble contributions from all collisions (including
// mollified EE). We test self-consistency by computing per-collision
// contributions and summing them up, then comparing diagonal/quadratic-form.
// =============================================================================

TEST_CASE(
    "Cumulative gauss_newton_hessian_diagonal is self-consistent",
    "[potential][barrier_potential][distance_vector][cumulative]")
{
    const double dhat = sqrt(2.0);
    const double kappa = 1.0;
    const std::string mesh_name = "cube.ply";

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions collisions;
    collisions.set_use_area_weighting(false);
    collisions.build(mesh, vertices, dhat);
    REQUIRE(!collisions.empty());

    BarrierPotential barrier_potential(dhat, kappa);

    // Compute diagonal via the efficient cumulative method
    Eigen::VectorXd diag_efficient =
        barrier_potential.gauss_newton_hessian_diagonal(
            collisions, mesh, vertices);

    // Compute diagonal by summing up individual collision contributions
    const int dim = vertices.cols();
    Eigen::VectorXd diag_reference = Eigen::VectorXd::Zero(vertices.size());

    for (size_t i = 0; i < collisions.size(); i++) {
        const NormalCollision& collision = collisions[i];
        const VectorMax12d positions = collision.dof(vertices, edges, faces);

        VectorMax12d local_diag =
            barrier_potential.gauss_newton_hessian_diagonal(
                collision, positions);

        const auto vids = collision.vertex_ids(edges, faces);
        for (int vi = 0; vi < collision.num_vertices(); vi++) {
            if (vids[vi] >= 0) {
                for (int d = 0; d < dim; d++) {
                    if constexpr (VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor) {
                        diag_reference[vids[vi] * dim + d] +=
                            local_diag[vi * dim + d];
                    } else {
                        diag_reference[vids[vi] + d * vertices.rows()] +=
                            local_diag[vi + d * collision.num_vertices()];
                    }
                }
            }
        }
    }

    REQUIRE(diag_efficient.size() == diag_reference.size());
    CHECK(
        (diag_efficient - diag_reference).norm()
        < 1e-10 * (1.0 + diag_reference.norm()));
}

TEST_CASE(
    "Cumulative gauss_newton_hessian_quadratic_form is self-consistent",
    "[potential][barrier_potential][distance_vector][cumulative]")
{
    const double dhat = sqrt(2.0);
    const double kappa = 1.0;
    const std::string mesh_name = "cube.ply";

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
    CAPTURE(mesh_name);
    REQUIRE(success);

    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions collisions;
    collisions.set_use_area_weighting(false);
    collisions.build(mesh, vertices, dhat);
    REQUIRE(!collisions.empty());

    BarrierPotential barrier_potential(dhat, kappa);

    // Random direction vector
    Eigen::VectorXd p = Eigen::VectorXd::Random(vertices.size());
    p *= 0.1; // keep small for numerical stability

    // Compute pᵀHp via the efficient cumulative method
    double pTHp_efficient =
        barrier_potential.gauss_newton_hessian_quadratic_form(
            collisions, mesh, vertices, p);

    // Compute pᵀHp by summing up individual collision contributions
    const int dim = vertices.cols();
    const auto p_mat =
        p.reshaped<VERTEX_DERIVATIVE_LAYOUT>(vertices.rows(), dim);

    double pTHp_reference = 0.0;
    for (size_t i = 0; i < collisions.size(); i++) {
        const NormalCollision& collision = collisions[i];
        const VectorMax12d positions = collision.dof(vertices, edges, faces);

        const auto vids = collision.vertex_ids(edges, faces);
        VectorMax12d local_p(collision.num_vertices() * dim);
        for (int vi = 0; vi < collision.num_vertices(); vi++) {
            local_p.segment(vi * dim, dim) = p_mat.row(vids[vi]);
        }

        pTHp_reference += barrier_potential.gauss_newton_hessian_quadratic_form(
            collision, positions, local_p);
    }

    CHECK(
        pTHp_efficient
        == Catch::Approx(pTHp_reference).epsilon(1e-10).margin(1e-14));
}