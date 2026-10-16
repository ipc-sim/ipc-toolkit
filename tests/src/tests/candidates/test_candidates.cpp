#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/candidates/edge_face.hpp>
#include <ipc/candidates/face_face.hpp>
#include <ipc/candidates/plane_vertex.hpp>

#include <algorithm>

using namespace ipc;

TEST_CASE("Candidates", "[candidates]")
{
    Candidates candidates;
    candidates.vv_candidates = { { 0, 1 }, { 2, 3 } };
    candidates.ev_candidates = { { 3, 4 }, { 5, 6 } };
    candidates.ee_candidates = { { 6, 7 }, { 8, 9 } };
    candidates.fv_candidates = { { 10, 11 }, { 12, 13 } };
    candidates.pv_candidates.emplace_back(
        Eigen::Hyperplane<double, 3>(Eigen::Vector3d::UnitY(), 0), 14);
    candidates.pv_candidates.emplace_back(
        Eigen::Hyperplane<double, 3>(Eigen::Vector3d::UnitY(), 0), 15);

    CHECK(candidates.size() == 10);

    CHECK(
        dynamic_cast<VertexVertexCandidate&>(candidates[0])
        == candidates.vv_candidates[0]);
    CHECK(
        dynamic_cast<EdgeVertexCandidate&>(candidates[2])
        == candidates.ev_candidates[0]);
    CHECK(
        dynamic_cast<EdgeEdgeCandidate&>(candidates[4])
        == candidates.ee_candidates[0]);
    CHECK(
        dynamic_cast<FaceVertexCandidate&>(candidates[6])
        == candidates.fv_candidates[0]);
    CHECK(&(candidates[8]) == &(candidates.pv_candidates[0]));

    CHECK(&(candidates[0]) == &(candidates.vv_candidates[0]));
    CHECK(&(candidates[2]) == &(candidates.ev_candidates[0]));
    CHECK(&(candidates[4]) == &(candidates.ee_candidates[0]));
    CHECK(&(candidates[6]) == &(candidates.fv_candidates[0]));
    CHECK(&(candidates[8]) == &(candidates.pv_candidates[0]));

    try {
        candidates[candidates.size()];
        FAIL("Should have thrown an exception");
    } catch (const std::out_of_range& e) {
        SUCCEED("Exception thrown");
        CHECK(e.what() == std::string("Candidate index is out of range!"));
    } catch (...) {
        FAIL("Uknown exception thrown");
    }

    const Candidates& const_candidates = candidates;
    CHECK(
        const_candidates[0].as<VertexVertexCandidate>()
        == candidates.vv_candidates[0]);
    CHECK(
        const_candidates[2].as<EdgeVertexCandidate>()
        == candidates.ev_candidates[0]);
    CHECK(
        const_candidates[4].as<EdgeEdgeCandidate>()
        == candidates.ee_candidates[0]);
    CHECK(
        const_candidates[6].as<FaceVertexCandidate>()
        == candidates.fv_candidates[0]);
    CHECK(
        const_candidates[8].as<PlaneVertexCandidate>()
        == candidates.pv_candidates[0]);

    CHECK(&(const_candidates[0]) == &(candidates.vv_candidates[0]));
    CHECK(&(const_candidates[2]) == &(candidates.ev_candidates[0]));
    CHECK(&(const_candidates[4]) == &(candidates.ee_candidates[0]));
    CHECK(&(const_candidates[6]) == &(candidates.fv_candidates[0]));
    CHECK(&(const_candidates[8]) == &(candidates.pv_candidates[0]));

    try {
        const_candidates[candidates.size()];
        FAIL("Should have thrown an exception");
    } catch (const std::out_of_range& e) {
        SUCCEED("Exception thrown");
        CHECK(e.what() == std::string("Candidate index is out of range!"));
    } catch (...) {
        FAIL("Uknown exception thrown");
    }
}

TEST_CASE("Vertex-Vertex Candidate", "[candidates][vertex-vertex]")
{
    CHECK(VertexVertexCandidate(0, 1) == VertexVertexCandidate(0, 1));
    CHECK(VertexVertexCandidate(0, 1) == VertexVertexCandidate(1, 0));
    CHECK(VertexVertexCandidate(0, 1) != VertexVertexCandidate(0, 2));
    CHECK(VertexVertexCandidate(0, 1) != VertexVertexCandidate(2, 0));
    CHECK(VertexVertexCandidate(0, 1) < VertexVertexCandidate(0, 2));
    CHECK(VertexVertexCandidate(0, 1) < VertexVertexCandidate(2, 0));
    CHECK(!(VertexVertexCandidate(1, 1) < VertexVertexCandidate(0, 2)));
}

TEST_CASE("Edge-Vertex Candidate", "[candidates][edge-vertex]")
{
    CHECK(EdgeVertexCandidate(0, 1) == EdgeVertexCandidate(0, 1));
    CHECK(EdgeVertexCandidate(0, 1) != EdgeVertexCandidate(1, 0));
    CHECK(EdgeVertexCandidate(0, 1) != EdgeVertexCandidate(0, 2));
    CHECK(EdgeVertexCandidate(0, 1) != EdgeVertexCandidate(2, 0));
    CHECK(EdgeVertexCandidate(0, 1) < EdgeVertexCandidate(0, 2));
    CHECK(EdgeVertexCandidate(0, 1) < EdgeVertexCandidate(2, 0));
    CHECK(!(EdgeVertexCandidate(1, 1) < EdgeVertexCandidate(0, 2)));
}

TEST_CASE("Edge-Edge Candidate", "[candidates][edge-edge]")
{
    CHECK(EdgeEdgeCandidate(0, 1) == EdgeEdgeCandidate(0, 1));
    CHECK(EdgeEdgeCandidate(0, 1) == EdgeEdgeCandidate(1, 0));
    CHECK(EdgeEdgeCandidate(0, 1) != EdgeEdgeCandidate(0, 2));
    CHECK(EdgeEdgeCandidate(0, 1) != EdgeEdgeCandidate(2, 0));
    CHECK(EdgeEdgeCandidate(0, 1) < EdgeEdgeCandidate(0, 2));
    CHECK(EdgeEdgeCandidate(0, 1) < EdgeEdgeCandidate(2, 0));
    CHECK(!(EdgeEdgeCandidate(1, 1) < EdgeEdgeCandidate(0, 2)));
}

TEST_CASE("Face-Vertex Candidate", "[candidates][face-vertex]")
{
    CHECK(FaceVertexCandidate(0, 1) == FaceVertexCandidate(0, 1));
    CHECK(FaceVertexCandidate(0, 1) != FaceVertexCandidate(1, 0));
    CHECK(FaceVertexCandidate(0, 1) != FaceVertexCandidate(0, 2));
    CHECK(FaceVertexCandidate(0, 1) != FaceVertexCandidate(2, 0));
    CHECK(FaceVertexCandidate(0, 1) < FaceVertexCandidate(0, 2));
    CHECK(FaceVertexCandidate(0, 1) < FaceVertexCandidate(2, 0));
    CHECK(!(FaceVertexCandidate(1, 1) < FaceVertexCandidate(0, 2)));
}

TEST_CASE("Edge-Face Candidate", "[candidates][edge-face]")
{
    CHECK(EdgeFaceCandidate(0, 1) == EdgeFaceCandidate(0, 1));
    CHECK(EdgeFaceCandidate(0, 1) != EdgeFaceCandidate(1, 0));
    CHECK(EdgeFaceCandidate(0, 1) != EdgeFaceCandidate(0, 2));
    CHECK(EdgeFaceCandidate(0, 1) != EdgeFaceCandidate(2, 0));
    CHECK(EdgeFaceCandidate(0, 1) < EdgeFaceCandidate(0, 2));
    CHECK(EdgeFaceCandidate(0, 1) < EdgeFaceCandidate(2, 0));
    CHECK(!(EdgeFaceCandidate(1, 1) < EdgeFaceCandidate(0, 2)));
}

TEST_CASE("Face-Face Candidate", "[candidates][face-face]")
{
    CHECK(FaceFaceCandidate(0, 1) == FaceFaceCandidate(0, 1));
    CHECK(FaceFaceCandidate(0, 1) == FaceFaceCandidate(1, 0));
    CHECK(FaceFaceCandidate(0, 1) != FaceFaceCandidate(0, 2));
    CHECK(FaceFaceCandidate(0, 1) != FaceFaceCandidate(2, 0));
    CHECK(FaceFaceCandidate(0, 1) < FaceFaceCandidate(0, 2));
    CHECK(FaceFaceCandidate(0, 1) < FaceFaceCandidate(2, 0));
    CHECK(!(FaceFaceCandidate(1, 1) < FaceFaceCandidate(0, 2)));
}

TEST_CASE("Plane-Vertex Candidate", "[candidates][plane-vertex]")
{
    Eigen::Hyperplane<double, 3> plane(Eigen::Vector3d::UnitY(), 0);
    CHECK(PlaneVertexCandidate(plane, 1) == PlaneVertexCandidate(plane, 1));
    CHECK(PlaneVertexCandidate(plane, 1) != PlaneVertexCandidate(plane, 2));
    CHECK(PlaneVertexCandidate(plane, 1) < PlaneVertexCandidate(plane, 2));
    CHECK(!(PlaneVertexCandidate(plane, 2) < PlaneVertexCandidate(plane, 1)));
}

TEST_CASE("Candidates adjacency sets", "[candidates]")
{
    // Vertex 2 is 0.05 above the middle of edge 0.
    Eigen::MatrixXd V(4, 2);
    V << 0, 0, 1, 0, 0.5, 0.05, 3, 3;
    Eigen::MatrixXi E(2, 2);
    E << 0, 1, 2, 3;
    const CollisionMesh mesh(V, E, Eigen::MatrixXi());

    // Before convert_candidates_to_sets(), every query is empty.
    Candidates candidates;
    CHECK(candidates.vv_set(mesh, 0).empty());
    CHECK(candidates.ve_set(0).empty());
    CHECK(candidates.vf_set(0).empty());
    CHECK(candidates.ev_set(mesh, 0).empty());
    CHECK(candidates.ee_set(mesh, 0).empty());
    CHECK(candidates.ef_set(mesh, 0).empty());
    CHECK(candidates.fv_set(mesh, 0).empty());
    CHECK(candidates.fe_set(mesh, 0).empty());
    CHECK(candidates.ff_set(mesh, 0).empty());

    candidates.build(mesh, V, /*inflation_radius=*/0.1);
    REQUIRE(!candidates.ev_candidates.empty());
    candidates.convert_candidates_to_sets();
    // In 2D, vv_set also holds the endpoints of the edges near the vertex.
    CHECK(candidates.vv_set(mesh, 2) == std::set<index_t> { 0, 1 });
    CHECK(candidates.ve_set(2).count(0) == 1);
    CHECK(candidates.vf_set(2).empty());
}

TEST_CASE(
    "All-types candidates keep whole-mesh pairs alongside codim. primitives",
    "[candidates]")
{
    // Two triangles meet tip to tip (vertices 0 and 3). Far away, a codim.
    // vertex (9) sits next to a dangling edge (7, 8), and vertex 6 is stray.
    Eigen::MatrixXd V(10, 3);
    V << 0, 0, 0.025, 1, 0, 1, -1, 0, 1,   //
        0, 0, -0.025, 0, 1, -1, 0, -1, -1, //
        50, 50, 50,                        //
        -50, -50, -50, -49, -50, -50, -49.5, -50, -49.99;
    Eigen::MatrixXi F(2, 3);
    F << 0, 1, 2, 3, 4, 5;
    Eigen::MatrixXi E(7, 2);
    E << 0, 1, 1, 2, 2, 0, 3, 4, 4, 5, 5, 3, 7, 8;
    const CollisionMesh mesh(V, E, F);
    REQUIRE(mesh.num_codim_vertices() == 2);
    REQUIRE(mesh.num_codim_edges() == 1);

    const auto contains = [](const auto& list, const auto& candidate) {
        return std::find(list.begin(), list.end(), candidate) != list.end();
    };

    const bool all_types = GENERATE(false, true);
    CAPTURE(all_types);

    Candidates candidates;
    candidates.build(mesh, V, /*inflation_radius=*/0.05, nullptr, all_types);
    CHECK(candidates.has_all_types == all_types);

    // Both builds find the codim. edge next to the codim. vertex.
    CHECK(contains(candidates.ev_candidates, EdgeVertexCandidate(6, 9)));

    if (all_types) {
        // The tips are a surface vertex-vertex pair, which only all_types
        // detects. The codim. primitives must not replace it.
        CHECK(contains(candidates.vv_candidates, VertexVertexCandidate(0, 3)));
    }

    candidates.clear();
    CHECK(!candidates.has_all_types);

    // Trajectory candidates never include the ESP-only types.
    candidates.build(mesh, V, V, /*inflation_radius=*/0.05);
    CHECK(!candidates.has_all_types);
}
