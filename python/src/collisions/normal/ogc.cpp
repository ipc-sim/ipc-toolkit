#include <common.hpp>

#include <ipc/collisions/normal/ogc.hpp>

using namespace ipc;

void define_ogc(py::module_& m)
{
    py::module_ ogc =
        m.def_submodule("ogc", "Offset Geometric Contact (OGC) helpers");

    ogc.def(
        "check_vertex_feasible_region",
        py::overload_cast<
            const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
            Eigen::ConstRef<VectorMax3d>, const index_t>(
            &ipc::ogc::check_vertex_feasible_region),
        R"ipc_Qu8mg5v7(
        Check if point `x` is in the feasible region of vertex `vi`.

        Parameters:
            mesh: Collision mesh containing the vertex adjacencies
            vertices: Matrix of current vertex positions (rowwise)
            point: Position of the point to check
            vi: Index of the vertex for which to check the feasible region

        Returns:
            True if the point is in the feasible region, false otherwise
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "point"_a, "vi"_a);

    ogc.def(
        "check_vertex_feasible_region",
        py::overload_cast<
            const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
            const index_t, const index_t>(
            &ipc::ogc::check_vertex_feasible_region),
        R"ipc_Qu8mg5v7(
        Check if vertex `xi` is in the feasible region of vertex `vi`.

        Parameters:
            mesh: Collision mesh containing the vertex adjacencies
            vertices: Matrix of current vertex positions (rowwise)
            xi: Index of the vertex to check
            vi: Index of the vertex for which to check the feasible region

        Returns:
            True if the vertex is in the feasible region, false otherwise
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "xi"_a, "vi"_a);

    ogc.def(
        "check_edge_feasible_region", &ogc::check_edge_feasible_region,
        R"ipc_Qu8mg5v7(
        Check if vertex `xi` is in the feasible region of edge `ei`.

        Parameters:
            mesh: Collision mesh containing the edge adjacencies
            vertices: Matrix of current vertex positions (rowwise)
            xi: Index of the vertex to check
            ei: Index of the edge for which to check the feasible region

        Returns:
            True if the vertex is in the feasible region, false otherwise
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "xi"_a, "ei"_a);

    ogc.def(
        "is_edge_vertex_feasible", &ogc::is_edge_vertex_feasible,
        R"ipc_Qu8mg5v7(
        Check if the edge-vertex candidate is feasible.

        Parameters:
            mesh: Collision mesh containing the edge adjacencies
            vertices: Matrix of current vertex positions (rowwise)
            candidate: Edge-vertex candidate to check
            dtype: Edge-vertex distance type to use for the check.

        Returns:
            True if the edge-vertex candidate is in the feasible region, false otherwise
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "candidate"_a,
        "dtype"_a = PointEdgeDistanceType::AUTO);

    ogc.def(
        "is_edge_edge_feasible", &ogc::is_edge_edge_feasible,
        R"ipc_Qu8mg5v7(
        Check if the edge-edge candidate is feasible.

        Parameters:
            mesh: Collision mesh containing the edge adjacencies
            vertices: Matrix of current vertex positions (rowwise)
            candidate: Edge-edge candidate to check
            dtype: Edge-edge distance type to use for the check.

        Returns:
            True if the edge-edge candidate is in the feasible region, false otherwise
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "candidate"_a,
        "dtype"_a = EdgeEdgeDistanceType::AUTO);

    ogc.def(
        "is_face_vertex_feasible", &ogc::is_face_vertex_feasible,
        R"ipc_Qu8mg5v7(
        Check if the face-vertex candidate is feasible.

        Parameters:
            mesh: Collision mesh containing the face adjacencies
            vertices: Matrix of current vertex positions (rowwise)
            candidate: Face-vertex candidate to check
            dtype: Point-triangle distance type to use for the check.

        Returns:
            True if the face-vertex candidate is in the feasible region, false otherwise
        )ipc_Qu8mg5v7",
        "mesh"_a, "vertices"_a, "candidate"_a,
        "dtype"_a = PointTriangleDistanceType::AUTO);
}
