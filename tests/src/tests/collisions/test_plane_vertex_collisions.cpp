#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/candidates/plane_vertex.hpp>
#include <ipc/ccd/additive_ccd.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/normal/plane_vertex.hpp>
#include <ipc/collisions/tangential/plane_vertex.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>

using namespace ipc;

// =============================================================================
// PlaneVertexCandidate tests
// =============================================================================

TEST_CASE("PlaneVertexCandidate", "[candidates][plane-vertex]")
{
    const Eigen::Vector3d n = Eigen::Vector3d::UnitY();
    const Eigen::Vector3d o = Eigen::Vector3d::Zero();
    const Eigen::Hyperplane<double, 3> plane(n, o);

    PlaneVertexCandidate pvc(plane, 0);

    CHECK(pvc.num_vertices() == 1);

    Eigen::MatrixXi edges, faces;
    auto ids = pvc.vertex_ids(edges, faces);
    CHECK(ids == std::array<index_t, 4> { { 0, -1, -1, -1 } });

    SECTION("Distance computation")
    {
        Eigen::Vector3d point(1.0, 2.0, 3.0);
        double dist = pvc.compute_distance(point);
        CHECK(dist == Catch::Approx(4.0)); // squared distance = (2)^2

        Eigen::VectorXd grad = pvc.compute_distance_gradient(point);
        CHECK(grad.size() == 3);
        CHECK(grad.isApprox(Eigen::Vector3d(0, 4, 0)));

        Eigen::MatrixXd hess = pvc.compute_distance_hessian(point);
        CHECK(hess.isApprox(2 * n * n.transpose()));
    }

    SECTION("Coefficients")
    {
        Eigen::Vector3d point(0.5, 1.0, 0.5);
        VectorMax4d coeffs = pvc.compute_coefficients(point);
        CHECK(coeffs.size() == 1);
        CHECK(coeffs(0) == Catch::Approx(1.0));
    }

    SECTION("Normal (exercises unnormalized normal internally)")
    {
        Eigen::Vector3d point(0, 1, 0);
        // compute_normal is public and calls compute_unnormalized_normal
        VectorMax3d normal = pvc.compute_normal(point);
        CHECK(normal.isApprox(n));
    }

    SECTION("CCD")
    {
        AdditiveCCD ccd;
        ccd.conservative_rescaling = 1.0; // No rescaling for testing

        // Vertex starts above the plane and moves below
        Eigen::Vector3d v_t0(0, 1, 0);
        Eigen::Vector3d v_t1(0, -1, 0);
        double toi;
        bool hit =
            pvc.ccd(v_t0, v_t1, toi, /*min_distance=*/0, /*tmax=*/1.0, ccd);
        CHECK(hit);
        CHECK(toi == Catch::Approx(0.5));

        // Vertex starts and stays above the plane
        Eigen::Vector3d v_t0_above(0, 2, 0);
        Eigen::Vector3d v_t1_above(0, 1, 0);
        bool no_hit = pvc.ccd(
            v_t0_above, v_t1_above, toi, /*min_distance=*/0, /*tmax=*/1.0, ccd);
        CHECK_FALSE(no_hit);

        // Vertex starts below and moves further below — already past plane
        Eigen::Vector3d v_t0_below(0, -0.5, 0);
        Eigen::Vector3d v_t1_below(0, -2, 0);
        bool below_hit = pvc.ccd(
            v_t0_below, v_t1_below, toi, /*min_distance=*/0, /*tmax=*/1.0, ccd);
        CHECK_FALSE(below_hit);
    }
}

// =============================================================================
// PlaneVertexNormalCollision constructor tests
// =============================================================================

TEST_CASE(
    "PlaneVertexNormalCollision constructors", "[collision][plane-vertex]")
{
    const Eigen::Vector3d n = Eigen::Vector3d::UnitY();
    const Eigen::Vector3d o = Eigen::Vector3d::Zero();
    const Eigen::Hyperplane<double, 3> plane(n, o);

    SECTION("From PlaneVertexCandidate")
    {
        PlaneVertexCandidate candidate(plane, 5);
        PlaneVertexNormalCollision collision(candidate);
        CHECK(collision.vertex_id == 5);
        CHECK(collision.plane.normal().isApprox(n));
        CHECK(collision.num_vertices() == 1);
    }

    SECTION("From plane, vertex_id, weight, weight_gradient")
    {
        const double weight = 2.5;
        Eigen::SparseVector<double> weight_gradient(9);
        weight_gradient.insert(3) = 1.0;

        PlaneVertexNormalCollision collision(plane, 3, weight, weight_gradient);
        CHECK(collision.vertex_id == 3);
        CHECK(collision.plane.normal().isApprox(n));
        CHECK(collision.weight == Catch::Approx(weight));
        CHECK(collision.weight_gradient.nonZeros() == 1);
    }
}

// =============================================================================
// PlaneVertexTangentialCollision tests
// =============================================================================

TEST_CASE(
    "PlaneVertexTangentialCollision", "[collision][tangential][plane-vertex]")
{
    const Eigen::Vector3d n = Eigen::Vector3d::UnitY();
    const Eigen::Vector3d o = Eigen::Vector3d::Zero();
    const Eigen::Hyperplane<double, 3> plane(n, o);

    PlaneVertexNormalCollision normal_collision(
        plane, 0, 1.0, Eigen::SparseVector<double>(9));

    SECTION("Constructor from NormalCollision")
    {
        PlaneVertexTangentialCollision tc(normal_collision);
        CHECK(tc.vertex_id == 0);
        CHECK(tc.weight == Catch::Approx(1.0));
    }

    SECTION("Constructor with positions and normal_force")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        CHECK(tc.vertex_id == 0);
        CHECK(tc.normal_force_magnitude == Catch::Approx(normal_force));
    }

    SECTION("Constructor with positions and NormalPotential")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        BarrierPotential bp(1.0, 1.0);
        PlaneVertexTangentialCollision tc(normal_collision, pos, bp);
        CHECK(tc.vertex_id == 0);
    }

    SECTION("Tangent basis")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        // Call through base class reference (public in TangentialCollision)
        TangentialCollision& tc_base = tc;
        auto basis = tc_base.compute_tangent_basis(pos);
        // For a Y-normal plane, tangent basis should be in XZ plane
        CHECK(basis.rows() == 3);
        CHECK(basis.cols() == 2);
        // The tangent basis vectors should be orthogonal to the normal
        CHECK(std::abs(n.dot(basis.col(0))) < 1e-10);
        CHECK(std::abs(n.dot(basis.col(1))) < 1e-10);
    }

    SECTION("Tangent basis jacobian")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        TangentialCollision& tc_base = tc;
        auto jac = tc_base.compute_tangent_basis_jacobian(pos);
        CHECK(jac.isZero());
    }

    SECTION("Closest point")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        TangentialCollision& tc_base = tc;
        auto cp = tc_base.compute_closest_point(pos);
        CHECK(cp.size() == 0);
    }

    SECTION("Closest point jacobian")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        TangentialCollision& tc_base = tc;
        auto jac = tc_base.compute_closest_point_jacobian(pos);
        CHECK(jac.rows() == 0);
    }

    SECTION("Relative velocity")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        Eigen::Vector3d velocity(1.0, 2.0, 3.0);
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        TangentialCollision& tc_base = tc;
        auto rv = tc_base.relative_velocity(velocity);
        CHECK(rv.isApprox(velocity));
    }

    SECTION("Relative velocity matrix")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        VectorMax2d cp; // empty closest point for pv
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        TangentialCollision& tc_base = tc;
        auto mat = tc_base.relative_velocity_matrix(cp);
        CHECK(mat.isApprox(MatrixMax<double, 3, 12>::Identity(3, 3)));
    }

    SECTION("Relative velocity matrix jacobian")
    {
        Eigen::Vector3d pos(0, 0.5, 0);
        const double normal_force = 10.0;
        VectorMax2d cp; // empty closest point for pv
        PlaneVertexTangentialCollision tc(normal_collision, pos, normal_force);
        TangentialCollision& tc_base = tc;
        auto jac = tc_base.relative_velocity_matrix_jacobian(cp);
        CHECK(jac.isZero());
    }
}

// =============================================================================
// Integration: Candidates::build with planes
// =============================================================================

TEST_CASE(
    "Candidates build with ground plane", "[candidates][plane-vertex][build]")
{
    // 4 vertices: two above plane, two below
    Eigen::MatrixXd vertices(4, 3);
    vertices << 0, 0.5, 0, // above, close
        1, 2.0, 0,         // above, far
        -1, 0.1, 0,        // above, very close
        0, -0.5, 0;        // below

    CollisionMesh mesh(vertices);
    mesh.planes.emplace_back(Eigen::Vector3d::UnitY(), Eigen::Vector3d::Zero());

    const double inflation_radius = 1.0;

    SECTION("Static build")
    {
        Candidates candidates;
        candidates.build(mesh, vertices, inflation_radius);

        // Vertices 0 (d=0.5), 2 (d=0.1), 3 (d=-0.5) are within
        // inflation_radius of the plane
        CHECK(candidates.pv_candidates.size() >= 2);
        CHECK(!candidates.empty());

        // Verify pv_candidates are included in size()
        size_t expected_size = candidates.vv_candidates.size()
            + candidates.ev_candidates.size() + candidates.ee_candidates.size()
            + candidates.fv_candidates.size() + candidates.pv_candidates.size();
        CHECK(candidates.size() == expected_size);

        // Verify operator[] can reach pv_candidates
        if (!candidates.pv_candidates.empty()) {
            size_t pv_start = candidates.vv_candidates.size()
                + candidates.ev_candidates.size()
                + candidates.ee_candidates.size()
                + candidates.fv_candidates.size();
            auto& pv_ref = candidates[pv_start];
            CHECK(&pv_ref == &candidates.pv_candidates[0]);

            // const version
            const Candidates& const_cand = candidates;
            auto& pv_const_ref = const_cand[pv_start];
            CHECK(&pv_const_ref == &candidates.pv_candidates[0]);
        }
    }

    SECTION("CCD build")
    {
        // Move vertices downward
        Eigen::MatrixXd vertices_t1 = vertices;
        vertices_t1.col(1).array() -= 1.0;

        Candidates candidates;
        candidates.build(mesh, vertices, vertices_t1, inflation_radius);

        // All vertices moving toward or past the plane should be candidates
        CHECK(candidates.pv_candidates.size() >= 2);
    }

    SECTION("Clear resets pv_candidates")
    {
        Candidates candidates;
        candidates.build(mesh, vertices, inflation_radius);
        CHECK(!candidates.pv_candidates.empty());
        candidates.clear();
        CHECK(candidates.pv_candidates.empty());
        CHECK(candidates.empty());
    }
}

// =============================================================================
// Integration: NormalCollisions::build with planes
// =============================================================================

TEST_CASE(
    "NormalCollisions build with ground plane",
    "[collisions][normal][plane-vertex][build]")
{
    // Vertex sitting just above the ground plane
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0, 0.05, 0, // close to ground
        5, 5, 5;            // far away

    CollisionMesh mesh(vertices);
    mesh.init_area_jacobians();
    mesh.planes.emplace_back(Eigen::Vector3d::UnitY(), Eigen::Vector3d::Zero());

    const double dhat = 0.2;

    SECTION("Basic build")
    {
        NormalCollisions collisions;
        collisions.build(mesh, vertices, dhat);

        // Vertex 0 is at distance 0.05 from the plane (< dhat)
        CHECK(!collisions.pv_collisions.empty());
        CHECK(collisions.pv_collisions.size() >= 1);
        CHECK(!collisions.empty());
        CHECK(collisions.size() >= 1);

        // Verify is_plane_vertex
        size_t pv_idx = collisions.vv_collisions.size()
            + collisions.ev_collisions.size() + collisions.ee_collisions.size()
            + collisions.fv_collisions.size();
        CHECK(collisions.is_plane_vertex(pv_idx));
        CHECK(!collisions.is_vertex_vertex(pv_idx));
        CHECK(!collisions.is_edge_vertex(pv_idx));
        CHECK(!collisions.is_edge_edge(pv_idx));
        CHECK(!collisions.is_face_vertex(pv_idx));

        // Verify operator[] reaches pv_collisions
        auto& collision_ref = collisions[pv_idx];
        CHECK(
            &collision_ref
            == static_cast<NormalCollision*>(&collisions.pv_collisions[0]));

        const NormalCollisions& const_collisions = collisions;
        auto& const_ref = const_collisions[pv_idx];
        CHECK(
            &const_ref
            == static_cast<const NormalCollision*>(
                &collisions.pv_collisions[0]));
    }

    SECTION("With area weighting")
    {
        NormalCollisions collisions;
        collisions.set_use_area_weighting(true);
        collisions.build(mesh, vertices, dhat);

        if (!collisions.pv_collisions.empty()) {
            CHECK(collisions.pv_collisions[0].weight > 0);
        }
    }

    SECTION("With shape derivatives")
    {
        NormalCollisions collisions;
        collisions.set_use_area_weighting(true);
        collisions.set_enable_shape_derivatives(true);
        collisions.build(mesh, vertices, dhat);

        if (!collisions.pv_collisions.empty()) {
            CHECK(collisions.pv_collisions[0].weight > 0);
            // weight_gradient should be set
        }
    }

    SECTION("Clear resets pv_collisions")
    {
        NormalCollisions collisions;
        collisions.build(mesh, vertices, dhat);
        CHECK(!collisions.pv_collisions.empty());
        collisions.clear();
        CHECK(collisions.pv_collisions.empty());
        CHECK(collisions.empty());
    }

    SECTION("Barrier potential with plane-vertex")
    {
        NormalCollisions collisions;
        collisions.build(mesh, vertices, dhat);
        REQUIRE(!collisions.pv_collisions.empty());

        BarrierPotential bp(dhat, 1.0);
        double val = bp(collisions, mesh, vertices);
        CHECK(val > 0.0);
    }
}

// =============================================================================
// Integration: TangentialCollisions::build with planes
// =============================================================================

TEST_CASE(
    "TangentialCollisions build with ground plane",
    "[collisions][tangential][plane-vertex][build]")
{
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0, 0.05, 0, // close to ground
        5, 5, 5;            // far away

    CollisionMesh mesh(vertices);
    mesh.init_area_jacobians();
    mesh.planes.emplace_back(Eigen::Vector3d::UnitY(), Eigen::Vector3d::Zero());

    const double dhat = 0.2;
    const double mu = 0.5;

    NormalCollisions normal_collisions;
    normal_collisions.build(mesh, vertices, dhat);
    REQUIRE(!normal_collisions.pv_collisions.empty());

    BarrierPotential bp(dhat, 1.0);

    SECTION("Basic build")
    {
        TangentialCollisions tangential_collisions;
        tangential_collisions.build(mesh, vertices, normal_collisions, bp, mu);

        CHECK(!tangential_collisions.pv_collisions.empty());
        CHECK(tangential_collisions.pv_collisions.size() >= 1);
        CHECK(!tangential_collisions.empty());
        CHECK(tangential_collisions.size() >= 1);

        // Verify friction coefficient is set
        CHECK(tangential_collisions.pv_collisions[0].mu_s == Catch::Approx(mu));
        CHECK(tangential_collisions.pv_collisions[0].mu_k == Catch::Approx(mu));
    }

    SECTION("Per-vertex friction")
    {
        Eigen::VectorXd mu_s = Eigen::VectorXd::Constant(vertices.rows(), 0.3);
        Eigen::VectorXd mu_k = Eigen::VectorXd::Constant(vertices.rows(), 0.2);
        mu_s(0) = 0.8; // vertex 0 has different friction
        mu_k(0) = 0.6;

        TangentialCollisions tangential_collisions;
        tangential_collisions.build(
            mesh, vertices, normal_collisions, bp, mu_s, mu_k);

        CHECK(!tangential_collisions.pv_collisions.empty());
        // For pv collision, mu_s/mu_k should come from just the vertex
        CHECK(
            tangential_collisions.pv_collisions[0].mu_s == Catch::Approx(0.8));
        CHECK(
            tangential_collisions.pv_collisions[0].mu_k == Catch::Approx(0.6));
    }

    SECTION("operator[] reaches pv_collisions")
    {
        TangentialCollisions tangential_collisions;
        tangential_collisions.build(mesh, vertices, normal_collisions, bp, mu);

        size_t pv_idx = tangential_collisions.vv_collisions.size()
            + tangential_collisions.ev_collisions.size()
            + tangential_collisions.ee_collisions.size()
            + tangential_collisions.fv_collisions.size();

        CHECK(
            &tangential_collisions[pv_idx]
            == static_cast<TangentialCollision*>(
                &tangential_collisions.pv_collisions[0]));

        const TangentialCollisions& const_tc = tangential_collisions;
        CHECK(
            &const_tc[pv_idx]
            == static_cast<const TangentialCollision*>(
                &tangential_collisions.pv_collisions[0]));
    }

    SECTION("Clear resets pv_collisions")
    {
        TangentialCollisions tangential_collisions;
        tangential_collisions.build(mesh, vertices, normal_collisions, bp, mu);
        CHECK(!tangential_collisions.pv_collisions.empty());
        tangential_collisions.clear();
        CHECK(tangential_collisions.pv_collisions.empty());
        CHECK(tangential_collisions.empty());
    }

    SECTION("Friction potential with plane-vertex")
    {
        TangentialCollisions tangential_collisions;
        tangential_collisions.build(mesh, vertices, normal_collisions, bp, mu);
        REQUIRE(!tangential_collisions.pv_collisions.empty());

        Eigen::MatrixXd V1 = vertices;
        V1(0, 0) += 0.01; // Small lateral displacement
        Eigen::MatrixXd U = V1 - vertices;

        FrictionPotential fp(1e-3);
        double val = fp(tangential_collisions, mesh, U);
        CHECK(val >= 0.0);
    }
}

// =============================================================================
// Out-of-range tests for operator[]
// =============================================================================

TEST_CASE(
    "NormalCollisions operator[] out of range with pv",
    "[collisions][normal][plane-vertex]")
{
    NormalCollisions collisions;
    collisions.pv_collisions.emplace_back(
        Eigen::Hyperplane<double, 3>(Eigen::Vector3d::UnitY(), 0), 0);

    CHECK(collisions.size() == 1);
    CHECK(
        &collisions[0]
        == static_cast<NormalCollision*>(&collisions.pv_collisions[0]));

    try {
        collisions[1];
        FAIL("Should have thrown an exception");
    } catch (const std::out_of_range& e) {
        SUCCEED("Exception thrown");
    }
}

TEST_CASE(
    "TangentialCollisions operator[] out of range with pv",
    "[collisions][tangential][plane-vertex]")
{
    const Eigen::Vector3d n = Eigen::Vector3d::UnitY();
    const Eigen::Hyperplane<double, 3> plane(n, Eigen::Vector3d::Zero());

    PlaneVertexNormalCollision nc(
        plane, 0, 1.0, Eigen::SparseVector<double>(9));
    TangentialCollisions tangential_collisions;
    tangential_collisions.pv_collisions.emplace_back(nc);

    CHECK(tangential_collisions.size() == 1);
    CHECK(!tangential_collisions.empty());

    CHECK(
        &tangential_collisions[0]
        == static_cast<TangentialCollision*>(
            &tangential_collisions.pv_collisions[0]));

    try {
        tangential_collisions[1];
        FAIL("Should have thrown an exception");
    } catch (const std::out_of_range& e) {
        SUCCEED("Exception thrown");
    }
}
