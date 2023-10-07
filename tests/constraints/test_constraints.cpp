#include <catch2/catch_all.hpp>

#include <ipc/collisions/collision_constraints.hpp>

using namespace ipc;

TEST_CASE("Vertex-Vertex Constraint", "[constraint][vertex-vertex]")
{
    CHECK(VertexVertexConstraint(0, 1) == VertexVertexConstraint(0, 1));
    CHECK(VertexVertexConstraint(0, 1) == VertexVertexConstraint(1, 0));
    CHECK(VertexVertexConstraint(0, 1) != VertexVertexConstraint(0, 2));
    CHECK(VertexVertexConstraint(0, 1) != VertexVertexConstraint(2, 0));
    CHECK(VertexVertexConstraint(0, 1) < VertexVertexConstraint(0, 2));
    CHECK(VertexVertexConstraint(0, 1) < VertexVertexConstraint(2, 0));

    CHECK(
        VertexVertexConstraint(VertexVertexCandidate(0, 1))
        == VertexVertexConstraint(0, 1));
}

TEST_CASE("Edge-Vertex Constraint", "[constraint][edge-vertex]")
{
    CHECK(EdgeVertexConstraint(0, 1) == EdgeVertexConstraint(0, 1));
    CHECK(EdgeVertexConstraint(0, 1) != EdgeVertexConstraint(1, 0));
    CHECK(EdgeVertexConstraint(0, 1) != EdgeVertexConstraint(0, 2));
    CHECK(EdgeVertexConstraint(0, 1) != EdgeVertexConstraint(2, 0));
    CHECK(EdgeVertexConstraint(0, 1) < EdgeVertexConstraint(0, 2));
    CHECK(!(EdgeVertexConstraint(1, 1) < EdgeVertexConstraint(0, 2)));
    CHECK(EdgeVertexConstraint(0, 1) < EdgeVertexConstraint(2, 0));

    CHECK(
        EdgeVertexConstraint(EdgeVertexCandidate(0, 1))
        == EdgeVertexConstraint(0, 1));
}

TEST_CASE("Edge-Edge Constraint", "[constraint][edge-edge]")
{
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) == EdgeEdgeConstraint(0, 1, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) == EdgeEdgeConstraint(1, 0, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) != EdgeEdgeConstraint(0, 2, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) != EdgeEdgeConstraint(2, 0, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) < EdgeEdgeConstraint(0, 2, 1.0));
    CHECK(EdgeEdgeConstraint(0, 1, 0.0) < EdgeEdgeConstraint(2, 0, 1.0));

    CHECK(
        EdgeEdgeConstraint(EdgeEdgeCandidate(0, 1), 0.0)
        == EdgeEdgeConstraint(0, 1, 0.0));
}

TEST_CASE("Face-Vertex Constraint", "[constraint][face-vertex]")
{
    CHECK(FaceVertexConstraint(0, 1) == FaceVertexConstraint(0, 1));
    CHECK(FaceVertexConstraint(0, 1) != FaceVertexConstraint(1, 0));
    CHECK(FaceVertexConstraint(0, 1) != FaceVertexConstraint(0, 2));
    CHECK(FaceVertexConstraint(0, 1) != FaceVertexConstraint(2, 0));
    CHECK(FaceVertexConstraint(0, 1) < FaceVertexConstraint(0, 2));
    CHECK(!(FaceVertexConstraint(1, 1) < FaceVertexConstraint(0, 2)));
    CHECK(FaceVertexConstraint(0, 1) < FaceVertexConstraint(2, 0));

    CHECK(
        FaceVertexConstraint(FaceVertexCandidate(0, 1))
        == FaceVertexConstraint(0, 1));
}

TEST_CASE("Plane-Vertex Constraint", "[constraint][plane-vertex]")
{
    Eigen::MatrixXi E, F;
    const Eigen::Vector3d n(0, 1, 0), o(0, 0, 0);
    const PlaneVertexConstraint c(o, n, 0);
    CHECK(c.num_vertices() == 1);
    CHECK(c.vertex_ids(E, F) == std::array<long, 4> { { 0, -1, -1, -1 } });
    CHECK(c.plane_origin == o);
    CHECK(c.plane_normal == n);
    CHECK(c.vertex_id == 0);

    CHECK(c.compute_distance(Eigen::RowVector3d(0, -2, 0), E, F) == 4.0);
    CHECK(c.compute_distance(Eigen::RowVector3d(0, 2, 0), E, F) == 4.0);
    CHECK(
        c.compute_distance_gradient(Eigen::RowVector3d(0, 2, 0), E, F)
        == Eigen::Vector3d(0, 4, 0));
    CHECK(
        c.compute_distance_hessian(Eigen::RowVector3d(0, 2, 0), E, F)
        == 2 * n * n.transpose());
}