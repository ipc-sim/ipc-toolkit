#include <catch2/catch_test_macros.hpp>

#include <ipc/candidates/candidates.hpp>

using namespace ipc;

TEST_CASE("Candidates", "[candidates]")
{
    Candidates candidates;
    candidates.vv_candidates = { { 0, 1 }, { 2, 3 } };
    candidates.ev_candidates = { { 3, 4 }, { 5, 6 } };
    candidates.ee_candidates = { { 6, 7 }, { 8, 9 } };
    candidates.fv_candidates = { { 10, 11 }, { 12, 13 } };

    CHECK(candidates.size() == 8);

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

    CHECK(&(candidates[0]) == &(candidates.vv_candidates[0]));
    CHECK(&(candidates[2]) == &(candidates.ev_candidates[0]));
    CHECK(&(candidates[4]) == &(candidates.ee_candidates[0]));
    CHECK(&(candidates[6]) == &(candidates.fv_candidates[0]));

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
        dynamic_cast<const VertexVertexCandidate&>(const_candidates[0])
        == candidates.vv_candidates[0]);
    CHECK(
        dynamic_cast<const EdgeVertexCandidate&>(const_candidates[2])
        == candidates.ev_candidates[0]);
    CHECK(
        dynamic_cast<const EdgeEdgeCandidate&>(const_candidates[4])
        == candidates.ee_candidates[0]);
    CHECK(
        dynamic_cast<const FaceVertexCandidate&>(const_candidates[6])
        == candidates.fv_candidates[0]);

    CHECK(&(const_candidates[0]) == &(candidates.vv_candidates[0]));
    CHECK(&(const_candidates[2]) == &(candidates.ev_candidates[0]));
    CHECK(&(const_candidates[4]) == &(candidates.ee_candidates[0]));
    CHECK(&(const_candidates[6]) == &(candidates.fv_candidates[0]));

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
