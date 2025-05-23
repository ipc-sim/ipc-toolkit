#include <common.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/logger.hpp>

namespace py = pybind11;
using namespace ipc;

#ifdef IPC_TOOLKIT_WITH_ABSEIL
using MapCanCollide = std::unordered_map<
    std::pair<size_t, size_t>,
    bool,
    absl::Hash<std::pair<size_t, size_t>>>;
#else
using MapCanCollide = std::unordered_map<
    std::pair<size_t, size_t>,
    bool,
    Hash<std::pair<size_t, size_t>>>;
#endif

/// @brief A functor which the value of the pair if it is in the map, otherwise a default value.
class SparseCanCollide {
public:
    /// @brief Construct a new Sparse Can Collide object.
    /// @param explicit_values A map from vertex pairs to whether they can collide. Only the upper triangle is used. The map is assumed to be symmetric.
    /// @param default_value The default value to return if the pair is not in the map.
    SparseCanCollide(const MapCanCollide& explicit_values, bool default_value)
        : m_explicit_values(explicit_values)
        , m_default_value(default_value)
    {
    }

    /// @brief Can two vertices collide?
    /// @param i Index of the first vertex.
    /// @param j Index of the second vertex.
    /// @return The value of the pair if it is in the map, otherwise the default value.
    bool operator()(size_t i, size_t j) const
    {
        auto it = m_explicit_values.find({ std::min(i, j), std::max(i, j) });

        assert(
            m_explicit_values.find({ std::max(i, j), std::min(i, j) })
            == m_explicit_values.end());

        if (it != m_explicit_values.end()) {
            return it->second;
        }
        return m_default_value;
    }

private:
    const MapCanCollide m_explicit_values;
    const bool m_default_value;
};

class VertexPatchesCanCollide {
public:
    /// @brief Construct a new Vertex Patches Can Collide object.
    /// @param vertex_patches Vector of patches labels for each vertex.
    VertexPatchesCanCollide(Eigen::ConstRef<Eigen::VectorXi> vertex_patches)
        : m_vertex_patches(vertex_patches)
    {
    }

    /// @brief Can two vertices collide?
    /// @param i Index of the first vertex.
    /// @param j Index of the second vertex.
    /// @return true if the vertices are in different patches
    bool operator()(size_t i, size_t j) const
    {
        assert(i < m_vertex_patches.size());
        assert(j < m_vertex_patches.size());
        return m_vertex_patches[i] != m_vertex_patches[j];
    }

    size_t num_vertices() const { return m_vertex_patches.size(); }

private:
    const Eigen::VectorXi m_vertex_patches;
};

void define_collision_mesh(py::module_& m)
{
    py::class_<SparseCanCollide>(
        m, "SparseCanCollide",
        "A functor which the value of the pair if it is in the map, otherwise a default value.")
        .def(
            py::init<const MapCanCollide&, bool>(),
            R"ipc_Qu8mg5v7(
            Construct a new Sparse Can Collide object.

            Parameters:
                explicit_values: A map from vertex pairs to whether they can collide. Only the upper triangle is used. The map is assumed to be symmetric.
                default_value: The default value to return if the pair is not in the map.
            )ipc_Qu8mg5v7",
            py::arg("explicit_values"), py::arg("default_value"))
        .def(
            "__call__", &SparseCanCollide::operator(), R"ipc_Qu8mg5v7(
            Can two vertices collide?

            Parameters:
                i: Index of the first vertex.
                j: Index of the second vertex.

            Returns:
                The value of the pair if it is in the map, otherwise the default value.
            )ipc_Qu8mg5v7",
            py::arg("i"), py::arg("j"));

    py::class_<VertexPatchesCanCollide>(
        m, "VertexPatchesCanCollide",
        "A functor which returns true if the vertices are in different patches.")
        .def(
            py::init<Eigen::ConstRef<Eigen::VectorXi>>(),
            py::arg("vertex_patches"),
            R"ipc_Qu8mg5v7(
            Construct a new Vertex Patches Can Collide object.

            Parameters:
                vertex_patches: Vector of patches labels for each vertex.
            )ipc_Qu8mg5v7")
        .def(
            "__call__", &VertexPatchesCanCollide::operator(), R"ipc_Qu8mg5v7(
            Can two vertices collide?

            Parameters:
                i: Index of the first vertex.
                j: Index of the second vertex.

            Returns:
                True if the vertices are in different patches.
            )ipc_Qu8mg5v7",
            py::arg("i"), py::arg("j"));

    py::class_<CollisionMesh>(m, "CollisionMesh")
        .def(
            py::init<
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                const Eigen::SparseMatrix<double>&>(),
            R"ipc_Qu8mg5v7(
            Construct a new Collision Mesh object directly from the collision mesh vertices.

            Parameters:
                rest_positions: The vertices of the collision mesh at rest (#V × dim).
                edges: The edges of the collision mesh (#E × 2).
                faces: The faces of the collision mesh (#F × 3).
                displacement_map: The displacement mapping from displacements on the full mesh to the collision mesh.
            )ipc_Qu8mg5v7",
            py::arg("rest_positions"), py::arg("edges") = Eigen::MatrixXi(),
            py::arg("faces") = Eigen::MatrixXi(),
            py::arg("displacement_map") = Eigen::SparseMatrix<double>())
        .def(
            py::init<
                const std::vector<bool>&, Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                const Eigen::SparseMatrix<double>&>(),
            R"ipc_Qu8mg5v7(
            Construct a new Collision Mesh object from a full mesh vertices.

            Parameters:
                include_vertex: Vector of bools indicating whether each vertex should be included in the collision mesh.
                full_rest_positions: The vertices of the full mesh at rest (#V × dim).
                edges: The edges of the collision mesh indexed into the full mesh vertices (#E × 2).
                faces: The faces of the collision mesh indexed into the full mesh vertices (#F × 3).
                displacement_map: The displacement mapping from displacements on the full mesh to the collision mesh.
            )ipc_Qu8mg5v7",
            py::arg("include_vertex"), py::arg("full_rest_positions"),
            py::arg("edges") = Eigen::MatrixXi(),
            py::arg("faces") = Eigen::MatrixXi(),
            py::arg("displacement_map") = Eigen::SparseMatrix<double>())
        .def_static(
            "build_from_full_mesh", &CollisionMesh::build_from_full_mesh,
            R"ipc_Qu8mg5v7(
            Helper function that automatically builds include_vertex using construct_is_on_surface.

            Parameters:
                full_rest_positions: The full vertices at rest (#FV × dim).
                edges: The edge matrix of mesh (#E × 2).
                faces: The face matrix of mesh (#F × 3).

            Returns:
                Constructed CollisionMesh.
            )ipc_Qu8mg5v7",
            py::arg("full_rest_positions"), py::arg("edges"),
            py::arg("faces") = Eigen::MatrixXi())
        .def(
            "init_adjacencies", &CollisionMesh::init_adjacencies,
            "Initialize vertex-vertex and edge-vertex adjacencies.")
        .def(
            "init_area_jacobians", &CollisionMesh::init_area_jacobians,
            "Initialize vertex and edge areas.")
        .def_property_readonly(
            "num_vertices", &CollisionMesh::num_vertices,
            "Get the number of vertices in the collision mesh.")
        .def_property_readonly(
            "num_codim_vertices", &CollisionMesh::num_codim_vertices,
            "Get the number of codimensional vertices in the collision mesh.")
        .def_property_readonly(
            "num_codim_edges", &CollisionMesh::num_codim_edges,
            "Get the number of codimensional edges in the collision mesh.")
        .def_property_readonly(
            "num_edges", &CollisionMesh::num_edges,
            "Get the number of edges in the collision mesh.")
        .def_property_readonly(
            "num_faces", &CollisionMesh::num_faces,
            "Get the number of faces in the collision mesh.")
        .def_property_readonly(
            "dim", &CollisionMesh::dim, "Get the dimension of the mesh.")
        .def_property_readonly(
            "ndof", &CollisionMesh::ndof,
            "Get the number of degrees of freedom in the collision mesh.")
        .def_property_readonly(
            "full_num_vertices", &CollisionMesh::full_num_vertices,
            "Get the number of vertices in the full mesh.")
        .def_property_readonly(
            "full_ndof", &CollisionMesh::full_ndof,
            "Get the number of degrees of freedom in the full mesh.")
        .def_property_readonly(
            "rest_positions", &CollisionMesh::rest_positions,
            "Get the vertices of the collision mesh at rest (#V × dim).")
        .def_property_readonly(
            "codim_vertices", &CollisionMesh::codim_vertices,
            "Get the indices of codimensional vertices of the collision mesh (#CV x 1).")
        .def_property_readonly(
            "codim_edges", &CollisionMesh::codim_edges,
            "Get the indices of codimensional edges of the collision mesh (#CE x 1).")
        .def_property_readonly(
            "edges", &CollisionMesh::edges,
            "Get the edges of the collision mesh (#E × 2).")
        .def_property_readonly(
            "faces", &CollisionMesh::faces,
            "Get the faces of the collision mesh (#F × 3).")
        .def_property_readonly(
            "faces_to_edges", &CollisionMesh::faces_to_edges,
            "Get the mapping from faces to edges of the collision mesh (#F × 3).")
        .def(
            "vertices", &CollisionMesh::vertices,
            R"ipc_Qu8mg5v7(
            Compute the vertex positions from the positions of the full mesh.

            Parameters:
                full_positions: The vertex positions of the full mesh (#FV × dim).

            Returns:
                The vertex positions of the collision mesh (#V × dim).
            )ipc_Qu8mg5v7",
            py::arg("full_positions"))
        .def(
            "displace_vertices", &CollisionMesh::displace_vertices,
            R"ipc_Qu8mg5v7(
            Compute the vertex positions from vertex displacements on the full mesh.

            Parameters:
                full_displacements: The vertex displacements on the full mesh (#FV × dim).

            Returns:
                The vertex positions of the collision mesh (#V × dim).
            )ipc_Qu8mg5v7",
            py::arg("full_displacements"))
        .def(
            "map_displacements", &CollisionMesh::map_displacements,
            R"ipc_Qu8mg5v7(
            Map vertex displacements on the full mesh to vertex displacements on the collision mesh.

            Parameters:
                full_displacements: The vertex displacements on the full mesh (#FV × dim).

            Returns:
                The vertex displacements on the collision mesh (#V × dim).
            )ipc_Qu8mg5v7",
            py::arg("full_displacements"))
        .def(
            "to_full_vertex_id", &CollisionMesh::to_full_vertex_id,
            R"ipc_Qu8mg5v7(
            Map a vertex ID to the corresponding vertex ID in the full mesh.

            Parameters:
                id: Vertex ID in the collision mesh.

            Returns:
                Vertex ID in the full mesh.
            )ipc_Qu8mg5v7",
            py::arg("id"))
        .def(
            "to_full_dof",
            py::overload_cast<Eigen::ConstRef<Eigen::VectorXd>>(
                &CollisionMesh::to_full_dof, py::const_),
            R"ipc_Qu8mg5v7(
            Map a vector quantity on the collision mesh to the full mesh.

            This is useful for mapping gradients from the collision mesh to the full mesh (i.e., applies the chain-rule).

            Parameters:
                x: Vector quantity on the collision mesh with size equal to ndof().

            Returns:
                Vector quantity on the full mesh with size equal to full_ndof().
            )ipc_Qu8mg5v7",
            py::arg("x"))
        .def(
            "to_full_dof",
            py::overload_cast<const Eigen::SparseMatrix<double>&>(
                &CollisionMesh::to_full_dof, py::const_),
            R"ipc_Qu8mg5v7(
            Map a matrix quantity on the collision mesh to the full mesh.

            This is useful for mapping Hessians from the collision mesh to the full mesh (i.e., applies the chain-rule).

            Parameters:
                X: Matrix quantity on the collision mesh with size equal to ndof() × ndof().

            Returns:
                Matrix quantity on the full mesh with size equal to full_ndof() × full_ndof().
            )ipc_Qu8mg5v7",
            py::arg("X"))
        .def_property_readonly(
            "vertex_vertex_adjacencies",
            &CollisionMesh::vertex_vertex_adjacencies,
            "Get the vertex-vertex adjacency matrix.")
        .def_property_readonly(
            "vertex_edge_adjacencies", &CollisionMesh::vertex_edge_adjacencies,
            "Get the vertex-edge adjacency matrix.")
        .def_property_readonly(
            "edge_vertex_adjacencies", &CollisionMesh::edge_vertex_adjacencies,
            "Get the edge-vertex adjacency matrix.")
        .def(
            "are_adjacencies_initialized",
            &CollisionMesh::are_adjacencies_initialized,
            "Determine if the adjacencies have been initialized by calling init_adjacencies().")
        .def(
            "is_vertex_on_boundary", &CollisionMesh::is_vertex_on_boundary,
            R"ipc_Qu8mg5v7(
            Is a vertex on the boundary of the collision mesh?

            Parameters:
                vi: Vertex ID.

            Returns:
                True if the vertex is on the boundary of the collision mesh.
            )ipc_Qu8mg5v7",
            py::arg("vi"))
        .def(
            "vertex_area", &CollisionMesh::vertex_area,
            R"ipc_Qu8mg5v7(
            Get the barycentric area of a vertex.

            Parameters:
                vi: Vertex ID.

            Returns:
                Barycentric area of vertex vi.
            )ipc_Qu8mg5v7",
            py::arg("vi"))
        .def_property_readonly(
            "vertex_areas", &CollisionMesh::vertex_areas,
            "Get the barycentric area of the vertices.")
        .def(
            "vertex_area_gradient",
            [](const CollisionMesh& self,
               const size_t vi) -> Eigen::SparseMatrix<double> {
                return self.vertex_area_gradient(vi);
            },
            R"ipc_Qu8mg5v7(
            Get the gradient of the barycentric area of a vertex wrt the rest positions of all points.

            Parameters:
                vi: Vertex ID.

            Returns:
                Gradient of the barycentric area of vertex vi wrt the rest positions of all points.
            )ipc_Qu8mg5v7",
            py::arg("vi"))
        .def(
            "edge_area", &CollisionMesh::edge_area,
            R"ipc_Qu8mg5v7(
            Get the barycentric area of an edge.

            Parameters:
                ei: Edge ID.

            Returns:
                Barycentric area of edge ei.
            )ipc_Qu8mg5v7",
            py::arg("ei"))
        .def(
            "edge_areas", &CollisionMesh::edge_areas,
            "Get the barycentric area of the edges.")
        .def(
            "edge_area_gradient",
            [](const CollisionMesh& self,
               const size_t ei) -> Eigen::SparseMatrix<double> {
                return self.edge_area_gradient(ei);
            },
            R"ipc_Qu8mg5v7(
            Get the gradient of the barycentric area of an edge wrt the rest positions of all points.

            Parameters:
                ei: Edge ID.

            Returns:
                Gradient of the barycentric area of edge ei wrt the rest positions of all points.
            )ipc_Qu8mg5v7",
            py::arg("ei"))
        .def(
            "are_area_jacobians_initialized",
            &CollisionMesh::are_area_jacobians_initialized,
            "Determine if the area Jacobians have been initialized by calling init_area_jacobians().")
        .def_static(
            "construct_is_on_surface", &CollisionMesh::construct_is_on_surface,
            R"ipc_Qu8mg5v7(
            Construct a vector of bools indicating whether each vertex is on the surface.

            Parameters:
                num_vertices: The number of vertices in the mesh.
                edges: The surface edges of the mesh (#E × 2).
                codim_vertices: The indices of codimensional vertices (#CV x 1).

            Returns:
                A vector of bools indicating whether each vertex is on the surface.
            )ipc_Qu8mg5v7",
            py::arg("num_vertices"), py::arg("edges"),
            py::arg("codim_vertices") = Eigen::VectorXi())
        .def_static(
            "construct_faces_to_edges",
            &CollisionMesh::construct_faces_to_edges,
            R"ipc_Qu8mg5v7(
            Construct a matrix that maps from the faces' edges to rows in the edges matrix.

            Parameters:
                faces: The face matrix of mesh (#F × 3).
                edges: The edge matrix of mesh (#E × 2).

            Returns:
                Matrix that maps from the faces' edges to rows in the edges matrix.
            )ipc_Qu8mg5v7",
            py::arg("faces"), py::arg("edges"))
        .def_property(
            "can_collide", [](CollisionMesh& self) { return self.can_collide; },
            [](CollisionMesh& self, const py::object& can_collide) {
                if (py::isinstance<SparseCanCollide>(can_collide)) {

                    self.can_collide = py::cast<SparseCanCollide>(can_collide);

                } else if (py::isinstance<VertexPatchesCanCollide>(
                               can_collide)) {

                    const VertexPatchesCanCollide& vertex_patches_can_collide =
                        py::cast<VertexPatchesCanCollide>(can_collide);

                    if (self.num_vertices()
                        != vertex_patches_can_collide.num_vertices()) {
                        throw py::value_error(
                            "The number of vertices in the VertexPatchesCanCollide object must match the number of vertices in the CollisionMesh.");
                    }

                    self.can_collide = vertex_patches_can_collide;

                } else if (py::isinstance<py::function>(can_collide)) {

                    logger().warn(
                        "Using a custom function for can_collide is deprecated because it is slow. "
                        "Please use a SparseCanCollide or VertexPatchesCanCollide object.");
                    self.can_collide = [can_collide](size_t i, size_t j) {
                        return py::cast<bool>(can_collide(i, j));
                    };

                } else {
                    throw py::value_error(
                        "Unknown type for can_collide. Must be a SparseCanCollide, VertexPatchesCanCollide, or a function.");
                }
            },
            R"ipc_Qu8mg5v7(
            A function that takes two vertex IDs and returns true if the vertices (and faces or edges containing the vertices) can collide.

            By default all primitives can collide with all other primitives.
            )ipc_Qu8mg5v7")
        .def_property(
            "vertex_materials",
            [](const CollisionMesh& self) {
                Eigen::VectorXi materials(self.num_vertices());
                for (int i = 0; i < self.num_vertices(); i++) {
                    materials(i) = self.vertex_material(i);
                }
                return materials;
            },
            [](CollisionMesh& self, const Eigen::VectorXi& materials) {
                self.set_vertex_materials(materials);
            },
            "Material IDs for vertices in the collision mesh")
        .def_property(
            "edge_materials",
            [](const CollisionMesh& self) {
                Eigen::VectorXi materials(self.num_edges());
                for (int i = 0; i < self.num_edges(); i++) {
                    materials(i) = self.edge_material(i);
                }
                return materials;
            },
            [](CollisionMesh& self, const Eigen::VectorXi& materials) {
                self.set_edge_materials(materials);
            },
            "Material IDs for edges in the collision mesh")
        .def_property(
            "face_materials",
            [](const CollisionMesh& self) {
                Eigen::VectorXi materials(self.num_faces());
                for (int i = 0; i < self.num_faces(); i++) {
                    materials(i) = self.face_material(i);
                }
                return materials;
            },
            [](CollisionMesh& self, const Eigen::VectorXi& materials) {
                self.set_face_materials(materials);
            },
            "Material IDs for faces in the collision mesh")
        .def(
            "set_single_material_id", &CollisionMesh::set_single_material_id,
            "Set a single material ID for the entire mesh",
            py::arg("material_id"))
        .def(
            "has_material_ids", &CollisionMesh::has_material_ids,
            "Check if material IDs are being used for this mesh");
}
