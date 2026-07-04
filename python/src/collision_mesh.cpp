#include <common.hpp>

#include <ipc/collision_filter.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/utils/logger.hpp>

using namespace ipc;

namespace {

struct PairHash {
    std::size_t operator()(const std::pair<size_t, size_t>& p) const noexcept
    {
        const size_t h1 = std::hash<size_t>()(p.first);
        const size_t h2 = std::hash<size_t>()(p.second);
        // A simplified version of boost::hash_combine
        return h1 ^ (h2 + 0x9e3779b97f4a7c15 + (h1 << 6) + (h1 >> 2));
    }
};

using MapCanCollide =
    std::unordered_map<std::pair<size_t, size_t>, bool, PairHash>;

CollisionFilter
make_sparse_filter(MapCanCollide explicit_values, bool default_value)
{
    return [m_explicit_values = std::move(explicit_values),
            m_default_value = default_value](size_t i, size_t j) {
        auto it = m_explicit_values.find({ std::min(i, j), std::max(i, j) });
        assert(
            m_explicit_values.find({ std::max(i, j), std::min(i, j) })
            == m_explicit_values.end());
        if (it != m_explicit_values.end()) {
            return it->second;
        }
        return m_default_value;
    };
}

} // namespace

void define_collision_mesh(py::module_& m)
{
    py::class_<CollisionFilter>(
        m, "CollisionFilter",
        R"ipc_Qu8mg5v7(
        A composable, type-erased collision filter.

        Wraps any callable ``bool(int, int)`` and supports logical composition
        via ``|`` (union), ``&`` (intersection), and ``~`` (negation) operators.
        The default-constructed filter accepts all pairs.

        Example:
            .. code-block:: python

                patches = CollisionFilter(lambda i, j: patches[i] != patches[j])
                static  = make_static_obstacle_filter(n_dynamic)
                active  = patches & static
                if active(i, j):
                    ...
        )ipc_Qu8mg5v7")
        .def(py::init<>(), "Default filter: accept all pairs.")
        .def(
            py::init([](const py::function& fn) {
                return CollisionFilter([fn](size_t i, size_t j) -> bool {
                    py::gil_scoped_acquire gil;
                    return py::cast<bool>(fn(i, j));
                });
            }),
            "Construct from a Python callable ``bool(int, int)``. "
            "Python-backed filters acquire the GIL on each call and may not "
            "be safe or performant for parallel broad-phase use.",
            "fn"_a)
        .def(
            "__call__", &CollisionFilter::operator(),
            "Test whether two vertices may collide.", "vi"_a, "vj"_a)
        .def(
            "__or__",
            [](const CollisionFilter& a, const CollisionFilter& b) {
                return a | b;
            },
            "Union: accept if EITHER filter passes.")
        .def(
            "__and__",
            [](const CollisionFilter& a, const CollisionFilter& b) {
                return a & b;
            },
            "Intersection: accept only if BOTH filters pass.")
        .def(
            "__invert__", [](const CollisionFilter& f) { return !f; },
            "Negation: accept only if this filter rejects.")
        .def(
            "__ior__",
            [](CollisionFilter& f,
               const CollisionFilter& b) -> CollisionFilter& { return f |= b; })
        .def(
            "__iand__",
            [](CollisionFilter& f, const CollisionFilter& b)
                -> CollisionFilter& { return f &= b; });

    m.def(
        "make_sparse_filter", &make_sparse_filter,
        R"ipc_Qu8mg5v7(
        Create a filter from a sparse map of explicit vertex-pair values.

        Pairs present in ``explicit_values`` use the stored boolean; all
        other pairs fall back to ``default_value``.  Only the upper triangle
        of the pair space is used — keys must satisfy ``i < j``.

        Parameters:
            explicit_values: Dict mapping ``(i, j)`` pairs (``i < j``) to
                whether those two vertices can collide.
            default_value: Value returned for pairs not in the map.

        Returns:
            A CollisionFilter backed by the sparse map.
        )ipc_Qu8mg5v7",
        "explicit_values"_a, "default_value"_a);

    m.def(
        "make_vertex_patches_filter", &make_vertex_patches_filter,
        R"ipc_Qu8mg5v7(
        Create a filter that only allows collisions between vertices in different patches (e.g., different garment panels or bodies).

        Parameters:
            patch_ids: Per-vertex patch label vector (one entry per vertex).

        Returns:
            A CollisionFilter that blocks same-patch pairs.
        )ipc_Qu8mg5v7",
        "patch_ids"_a);

    m.def(
        "make_static_obstacle_filter", &make_static_obstacle_filter,
        R"ipc_Qu8mg5v7(
        Create a filter that prevents static obstacles from colliding with each other.
        A vertex is considered "static" if its index is >= n_dynamic.
        Pairs where both vertices are static are rejected.

        Parameters:
            n_dynamic: Number of dynamic (simulated) vertices; static vertices occupy indices [n_dynamic, n_verts).

        Returns:
            A CollisionFilter that blocks static-static pairs.
        )ipc_Qu8mg5v7",
        "n_dynamic"_a);

    m.def(
        "make_connected_components_filter", &make_connected_components_filter,
        R"ipc_Qu8mg5v7(
        Create a filter that prevents self-collisions within a connected
        component of the face mesh. Two vertices in the same connected
        component are blocked; cross-component pairs are allowed.

        Parameters:
            faces: Face index matrix (#F × 3).

        Returns:
            A CollisionFilter that blocks intra-component pairs.
        )ipc_Qu8mg5v7",
        "faces"_a);

    py::class_<Eigen::Hyperplane<double, 3>>(m, "Hyperplane")
        .def(py::init<>())
        .def(py::init<Eigen::Vector3d, Eigen::Vector3d>())
        .def(
            "normal",
            [](const Eigen::Hyperplane<double, 3>& self) {
                return self.normal();
            })
        .def(
            "offset",
            [](const Eigen::Hyperplane<double, 3>& self) {
                return self.offset();
            })
        .def("origin", [](const Eigen::Hyperplane<double, 3>& self) {
            return -self.offset() * self.normal();
        });

    py::class_<CollisionMesh, std::shared_ptr<CollisionMesh>>(
        m, "CollisionMesh")
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
            "rest_positions"_a, "edges"_a = Eigen::MatrixXi(),
            "faces"_a = Eigen::MatrixXi(),
            "displacement_map"_a = Eigen::SparseMatrix<double>())
        .def(
            py::init<
                const std::vector<bool>&, const std::vector<bool>&,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>,
                const Eigen::SparseMatrix<double>&>(),
            R"ipc_Qu8mg5v7(
            Construct a new Collision Mesh object from a full mesh vertices.

            Parameters:
                include_vertex: Vector of bools indicating whether each vertex should be included in the collision mesh.
                orient_vertex: Vector of bools indicating whether each vertex is orientable.
                full_rest_positions: The vertices of the full mesh at rest (#V × dim).
                edges: The edges of the collision mesh indexed into the full mesh vertices (#E × 2).
                faces: The faces of the collision mesh indexed into the full mesh vertices (#F × 3).
                displacement_map: The displacement mapping from displacements on the full mesh to the collision mesh.
            )ipc_Qu8mg5v7",
            "include_vertex"_a, "orient_vertex"_a, "full_rest_positions"_a,
            "edges"_a = Eigen::MatrixXi(), "faces"_a = Eigen::MatrixXi(),
            "displacement_map"_a = Eigen::SparseMatrix<double>())
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
            "full_rest_positions"_a, "edges"_a, "faces"_a = Eigen::MatrixXi())
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
            "full_positions"_a)
        .def(
            "displace_vertices", &CollisionMesh::displace_vertices,
            R"ipc_Qu8mg5v7(
            Compute the vertex positions from vertex displacements on the full mesh.

            Parameters:
                full_displacements: The vertex displacements on the full mesh (#FV × dim).

            Returns:
                The vertex positions of the collision mesh (#V × dim).
            )ipc_Qu8mg5v7",
            "full_displacements"_a)
        .def(
            "map_displacements", &CollisionMesh::map_displacements,
            R"ipc_Qu8mg5v7(
            Map vertex displacements on the full mesh to vertex displacements on the collision mesh.

            Parameters:
                full_displacements: The vertex displacements on the full mesh (#FV × dim).

            Returns:
                The vertex displacements on the collision mesh (#V × dim).
            )ipc_Qu8mg5v7",
            "full_displacements"_a)
        .def(
            "to_full_vertex_id",
            py::overload_cast<const index_t>(
                &CollisionMesh::to_full_vertex_id, py::const_),
            R"ipc_Qu8mg5v7(
            Map a vertex ID to the corresponding vertex ID in the full mesh.

            Parameters:
                id: Vertex ID in the collision mesh.

            Returns:
                Vertex ID in the full mesh.
            )ipc_Qu8mg5v7",
            "id"_a)
        .def(
            "to_full_vertex_id",
            py::overload_cast<>(&CollisionMesh::to_full_vertex_id, py::const_),
            R"ipc_Qu8mg5v7(
            Get the complete mapping of vertex IDs to their corresponding vertex IDs in the full mesh.

            Returns:
                Vector of size num_vertices() where each entry is the full vertex ID corresponding to the collision mesh vertex ID.
            )ipc_Qu8mg5v7")
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
            "x"_a)
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
            "X"_a)
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
            "vi"_a)
        .def(
            "vertex_area", &CollisionMesh::vertex_area,
            R"ipc_Qu8mg5v7(
            Get the barycentric area of a vertex.

            Parameters:
                vi: Vertex ID.

            Returns:
                Barycentric area of vertex vi.
            )ipc_Qu8mg5v7",
            "vi"_a)
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
            "vi"_a)
        .def(
            "edge_area", &CollisionMesh::edge_area,
            R"ipc_Qu8mg5v7(
            Get the barycentric area of an edge.

            Parameters:
                ei: Edge ID.

            Returns:
                Barycentric area of edge ei.
            )ipc_Qu8mg5v7",
            "ei"_a)
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
            "ei"_a)
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
            "num_vertices"_a, "edges"_a, "codim_vertices"_a = Eigen::VectorXi())
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
            "faces"_a, "edges"_a)
        .def_property(
            "can_collide", [](CollisionMesh& self) { return self.can_collide; },
            [](CollisionMesh& self, const py::object& can_collide) {
                if (py::isinstance<CollisionFilter>(can_collide)) {
                    self.can_collide = py::cast<CollisionFilter>(can_collide);
                } else if (py::isinstance<py::function>(can_collide)) {
                    logger().warn(
                        "Using a custom Python function for can_collide is deprecated. Please use a CollisionFilter object.");
                    self.can_collide =
                        CollisionFilter([can_collide](size_t i, size_t j) {
                            py::gil_scoped_acquire gil;
                            return py::cast<bool>(can_collide(i, j));
                        });
                } else {
                    throw py::value_error(
                        "can_collide must be a CollisionFilter or a callable "
                        "bool(int, int). Use make_sparse_filter(), "
                        "make_vertex_patches_filter(), or similar factory "
                        "functions to create a CollisionFilter.");
                }
            },
            R"ipc_Qu8mg5v7(
            A function that takes two vertex IDs and returns true if the vertices (and faces or edges containing the vertices) can collide.

            By default all primitives can collide with all other primitives.
            )ipc_Qu8mg5v7")
        .def_readwrite(
            "planes", &CollisionMesh::planes,
            R"ipc_Qu8mg5v7(
            A vector of planes in the collision mesh.

            Each plane is represented as a `Hyperplane` object (wrapping `Eigen::Hyperplane<double, 3>`).
            In Python, a `Hyperplane` can be constructed from either:

            * a normal and a point on the plane: ``Hyperplane(normal, point)``, or
            * a normal and an offset: ``Hyperplane(normal, offset)``,

            where ``normal`` and ``point`` are 3D vectors and ``offset`` is a scalar.
            )ipc_Qu8mg5v7");
}
