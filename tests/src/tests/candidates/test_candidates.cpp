#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/candidates/candidates.hpp>
#include <ipc/candidates/edge_face.hpp>
#include <ipc/candidates/face_face.hpp>
#include <ipc/candidates/plane_vertex.hpp>

#include <igl/PI.h>
#include <igl/edges.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <string>
#include <vector>

using namespace ipc;

namespace {
std::vector<index_t> ids(const AdjacencyList& list)
{
    return std::vector<index_t>(list.begin(), list.end());
}

std::vector<index_t> ids(const std::set<index_t>& set)
{
    return std::vector<index_t>(set.begin(), set.end());
}

// Whether the list is strictly ascending, i.e. sorted and free of duplicates.
bool is_sorted_unique(const AdjacencyList& list)
{
    return std::adjacent_find(list.begin(), list.end(), std::greater_equal<>())
        == list.end();
}
} // namespace

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
        CHECK(
            e.what()
            == fmt::format(
                "Candidate index {} is out of range! Total candidates: {}",
                candidates.size(), candidates.size()));
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
        CHECK(
            e.what()
            == fmt::format(
                "Candidate index {} is out of range! Total candidates: {}",
                candidates.size(), candidates.size()));
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
    CHECK(ids(candidates.vv_set(mesh, 2)) == std::vector<index_t> { 0, 1 });
    CHECK(candidates.ve_set(2).contains(0));
    CHECK(candidates.vf_set(2).empty());
}

// Two tetrahedra apex to apex, so every primitive pairing has candidates.
TEST_CASE(
    "Candidate adjacency lists join candidate partners with mesh neighbors",
    "[candidates]")
{
    const double c = std::cos(2 * igl::PI / 3), s = std::sin(2 * igl::PI / 3);
    const double r = std::cos(igl::PI / 3), q = std::sin(igl::PI / 3);
    Eigen::MatrixXd V(8, 3);
    V << 0, 0, 0.05, 1, 0, 1, c, s, 1, c, -s, 1, //
        0, 0, -0.05, r, q, -1, -1, 0, -1, r, -q, -1;
    Eigen::MatrixXi F(8, 3);
    F << 0, 2, 1, 0, 3, 2, 0, 1, 3, 1, 2, 3, //
        4, 5, 6, 4, 6, 7, 4, 7, 5, 5, 7, 6;
    Eigen::MatrixXi E;
    igl::edges(F, E);
    const CollisionMesh mesh(V, E, F);

    Candidates candidates;
    candidates.build(
        mesh, V, /*inflation_radius=*/1, /*broad_phase=*/nullptr,
        /*all_types=*/true);
    REQUIRE(!candidates.vv_candidates.empty());
    REQUIRE(!candidates.ev_candidates.empty());
    REQUIRE(!candidates.ee_candidates.empty());
    REQUIRE(!candidates.fv_candidates.empty());
    REQUIRE(!candidates.ef_candidates.empty());
    REQUIRE(!candidates.ff_candidates.empty());
    candidates.convert_candidates_to_sets();

    // The expected neighbors, computed straight from the candidate lists and
    // the mesh's vertex ids without using any of the adjacency structure.
    const size_t nv = mesh.num_vertices(), ne = mesh.num_edges(),
                 nf = mesh.num_faces();
    std::vector<std::set<index_t>> vv(nv), ve(nv), vf(nv), ev(ne), ee(ne),
        ef(ne), fv(nf), fe(nf), ff(nf);
    for (const auto& [a, b] : candidates.vv_candidates) {
        vv[a].insert(b);
        vv[b].insert(a);
    }
    for (const auto& [e, v] : candidates.ev_candidates) {
        ev[e].insert(v);
        ve[v].insert(e);
    }
    for (const auto& [a, b] : candidates.ee_candidates) {
        ee[a].insert(b);
        ee[b].insert(a);
    }
    for (const auto& [f, v] : candidates.fv_candidates) {
        fv[f].insert(v);
        vf[v].insert(f);
    }
    for (const auto& [e, f] : candidates.ef_candidates) {
        ef[e].insert(f);
        fe[f].insert(e);
    }
    for (const auto& [a, b] : candidates.ff_candidates) {
        ff[a].insert(b);
        ff[b].insert(a);
    }

    const auto has_vertex = [&](const Eigen::MatrixXi& M, index_t i,
                                index_t v) {
        return (M.row(i).array() == v).any();
    };
    for (index_t e = 0; e < index_t(ne); ++e) {
        for (index_t v : { mesh.edges()(e, 0), mesh.edges()(e, 1) }) {
            ev[e].insert(v);
            for (index_t other = 0; other < index_t(ne); ++other) {
                if (has_vertex(mesh.edges(), other, v)) {
                    ee[e].insert(other);
                }
            }
            for (index_t f = 0; f < index_t(nf); ++f) {
                if (has_vertex(mesh.faces(), f, v)) {
                    ef[e].insert(f);
                }
            }
        }
        ee[e].erase(e);
        for (index_t f = 0; f < index_t(nf); ++f) {
            if (has_vertex(mesh.faces(), f, mesh.edges()(e, 0))
                && has_vertex(mesh.faces(), f, mesh.edges()(e, 1))) {
                ef[e].erase(f); // A face containing the edge is not near it.
            }
        }
    }
    for (index_t f = 0; f < index_t(nf); ++f) {
        for (index_t lv = 0; lv < 3; ++lv) {
            const index_t v = mesh.faces()(f, lv);
            fv[f].insert(v);
            for (index_t e = 0; e < index_t(ne); ++e) {
                if (has_vertex(mesh.edges(), e, v)) {
                    fe[f].insert(e);
                }
            }
            for (index_t other = 0; other < index_t(nf); ++other) {
                if (has_vertex(mesh.faces(), other, v)) {
                    ff[f].insert(other);
                }
            }
        }
        ff[f].erase(f);
    }
    for (index_t v = 0; v < index_t(nv); ++v) {
        vv[v].erase(v);
    }

    const auto check = [&](const AdjacencyList& list,
                           const std::set<index_t>& want) {
        CHECK(is_sorted_unique(list));
        CHECK(ids(list) == ids(want));
    };
    for (index_t v = 0; v < index_t(nv); ++v) {
        CAPTURE(v);
        check(candidates.vv_set(mesh, v), vv[v]);
        check(candidates.ve_set(v), ve[v]);
        check(candidates.vf_set(v), vf[v]);
    }
    for (index_t e = 0; e < index_t(ne); ++e) {
        CAPTURE(e);
        check(candidates.ev_set(mesh, e), ev[e]);
        check(candidates.ee_set(mesh, e), ee[e]);
        check(candidates.ef_set(mesh, e), ef[e]);
    }
    for (index_t f = 0; f < index_t(nf); ++f) {
        CAPTURE(f);
        check(candidates.fv_set(mesh, f), fv[f]);
        check(candidates.fe_set(mesh, f), fe[f]);
        check(candidates.ff_set(mesh, f), ff[f]);
    }

    SECTION("Ids with no candidates have empty vertex lists")
    {
        CHECK(candidates.ve_set(1000).empty());
        CHECK(candidates.vf_set(-1).empty());
    }

    SECTION("A copy shares the lists")
    {
        const Candidates copy = candidates;
        for (index_t v = 0; v < index_t(nv); ++v) {
            CHECK(ids(copy.ve_set(v)) == ids(candidates.ve_set(v)));
            CHECK(ids(copy.ee_set(mesh, v)) == ids(candidates.ee_set(mesh, v)));
        }
    }

    SECTION("Clearing the candidates empties every query")
    {
        candidates.clear();
        CHECK(candidates.ve_set(0).empty());
        CHECK(candidates.ee_set(mesh, 0).empty());
    }
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
