#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/collisions/normal/ogc.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/face_vertex.hpp>

#include <Eigen/Dense>

using namespace ipc;

TEST_CASE("Check vertex feasible region", "[ogc]")
{
    // Initialize vertices
    Eigen::MatrixXd vertices(3, 3);
    vertices << -1, 0, 0, /**/ 0, 0, 0, /**/ 1, 0, 0;

    // Initialize edges to encode vertex adjacencies: (0-1), (0-2), (1-3), (2-3)
    Eigen::MatrixXi edges(2, 2);
    edges << 0, 1, /**/ 1, 2;

    // Initialize faces (empty for this test)
    Eigen::MatrixXi faces;

    // Create and initialize the CollisionMesh
    CollisionMesh mesh(vertices, edges, faces);

    Eigen::Vector3d x(-1.1, 0.1, 0); // Test point
    CHECK(ogc::check_vertex_feasible_region(mesh, vertices, x, 0));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 1));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 2));

    x.x() = -1e-8;
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 0));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 1));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 2));

    x.x() = 0.0;
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 0));
    CHECK(ogc::check_vertex_feasible_region(mesh, vertices, x, 1));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 2));

    x.x() = 1e-8;
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 0));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 1));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 2));

    x.x() = 1.1;
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 0));
    CHECK(!ogc::check_vertex_feasible_region(mesh, vertices, x, 1));
    CHECK(ogc::check_vertex_feasible_region(mesh, vertices, x, 2));
}

TEST_CASE("Check edge feasible region", "[ogc]")
{
    // Initialize vertices
    Eigen::MatrixXd vertices(5, 3);
    vertices << 0, 0, 1, /**/ 1, 0, 0, /**/ 0, 0, -1, /**/ -1, 0, 0, //
        1e-8, 1, 1e-8;

    // Initialize edges
    Eigen::MatrixXi edges(5, 2);
    edges << 0, 2, /**/ 0, 1, /**/ 1, 2, /**/ 2, 3, /**/ 3, 0;

    // Initialize faces so edge (0,1) is adjacent to vertices 2 and 3
    Eigen::MatrixXi faces(2, 3);
    faces << 0, 1, 2, /**/ 0, 2, 3;

    // Create and initialize the CollisionMesh
    CollisionMesh mesh(vertices, edges, faces);

    CHECK(!ogc::check_edge_feasible_region(mesh, vertices, 4, 0));

    vertices(0, 1) = vertices(2, 1) = 0.25;

    CHECK(ogc::check_edge_feasible_region(mesh, vertices, 4, 0));

    vertices(4, 2) = -1.1;

    CHECK(!ogc::check_edge_feasible_region(mesh, vertices, 4, 0));
}

TEST_CASE("Check edge-vertex feasibility", "[ogc]")
{
    // Initialize vertices
    Eigen::MatrixXd vertices(5, 3);
    vertices << 0, 0, 1, /**/ 1, 0, 0, /**/ 0, 0, -1, /**/ -1, 0, 0, //
        1e-8, 1, 1e-8;

    // Initialize edges
    Eigen::MatrixXi edges(5, 2);
    edges << 0, 2, /**/ 0, 1, /**/ 1, 2, /**/ 2, 3, /**/ 3, 0;

    // Initialize faces so edge (0,1) is adjacent to vertices 2 and 3
    Eigen::MatrixXi faces(2, 3);
    faces << 0, 1, 2, /**/ 0, 2, 3;

    // Create and initialize the CollisionMesh
    CollisionMesh mesh(vertices, edges, faces);

    EdgeVertexCandidate candidate(0, 4);

    CHECK(ogc::is_edge_vertex_feasible(mesh, vertices, candidate));

    vertices(0, 1) = vertices(2, 1) = 0.25;

    CHECK(ogc::is_edge_vertex_feasible(mesh, vertices, candidate));

    vertices(4, 2) = -1.1;

    CHECK(ogc::is_edge_vertex_feasible(mesh, vertices, candidate));

    vertices(4, 2) = 1.1;

    CHECK(ogc::is_edge_vertex_feasible(mesh, vertices, candidate));

    CHECK_THROWS(
        ogc::is_edge_vertex_feasible(
            mesh, vertices, EdgeVertexCandidate(0, 0),
            PointEdgeDistanceType(255)),
        "Invalid PointEdgeDistanceType: 255");
}

TEST_CASE("Check edge-edge feasibility", "[ogc]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("pentagon-dome.obj", vertices, edges, faces));

    vertices.conservativeResize(vertices.rows() * 2, vertices.cols());
    vertices.bottomRows(vertices.rows() / 2) =
        vertices.topRows(vertices.rows() / 2);
    vertices.col(1).array() += 1.5;

    edges.conservativeResize(edges.rows() * 2, edges.cols());
    edges.bottomRows(edges.rows() / 2) =
        edges.topRows(edges.rows() / 2).array() + vertices.rows() / 2;

    faces.conservativeResize(faces.rows() * 2, faces.cols());
    faces.bottomRows(faces.rows() / 2) =
        faces.topRows(faces.rows() / 2).array() + vertices.rows() / 2;

    // Create and initialize the CollisionMesh
    CollisionMesh mesh(vertices, edges, faces);

    CHECK(ogc::is_edge_edge_feasible(mesh, vertices, EdgeEdgeCandidate(0, 5)));

    CHECK_THROWS(
        ogc::is_edge_edge_feasible(
            mesh, vertices, EdgeEdgeCandidate(0, 0), EdgeEdgeDistanceType(255)),
        "Invalid EdgeEdgeDistanceType: 255");
}

TEST_CASE("Check face-vertex feasibility", "[ogc]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(tests::load_mesh("pentagon-dome.obj", vertices, edges, faces));

    vertices.conservativeResize(vertices.rows() + 1, vertices.cols());
    vertices.row(vertices.rows() - 1) = Eigen::RowVector3d(0, 1, 0);

    const int peak_vi = 5;

    // Create and initialize the CollisionMesh
    CollisionMesh mesh(vertices, edges, faces);

    for (int fi = 0; fi < faces.rows(); ++fi) {
        FaceVertexCandidate candidate(fi, peak_vi);

        CHECK(ogc::is_face_vertex_feasible(mesh, vertices, candidate));
    }

    CHECK_THROWS(
        ogc::is_face_vertex_feasible(
            mesh, vertices, FaceVertexCandidate(0, 0),
            PointTriangleDistanceType(255)),
        "Invalid PointTriangleDistanceType: 255");
}