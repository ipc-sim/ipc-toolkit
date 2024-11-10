#pragma once

#include <ipc/utils/unordered_map_and_set.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

/// @brief A class for encapsolating the transformation/selections needed to go from a volumetric FE mesh to a surface collision mesh.
class CollisionMesh {
public:
    /// @brief Construct a new Collision Mesh object.
    /// Collision Mesh objects are immutable, so use the other constructors.
    CollisionMesh() = default;

    /// @brief Construct a new Collision Mesh object directly from the collision mesh vertices.
    /// @param rest_positions The vertices of the collision mesh at rest (#V × dim).
    /// @param edges The edges of the collision mesh (#E × 2).
    /// @param faces The faces of the collision mesh (#F × 3).
    /// @param displacement_map The displacement mapping from displacements on the full mesh to the collision mesh.
    CollisionMesh(
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXi& edges = Eigen::MatrixXi(),
        const Eigen::MatrixXi& faces = Eigen::MatrixXi(),
        const Eigen::SparseMatrix<double>& displacement_map =
            Eigen::SparseMatrix<double>());

    /// @brief Construct a new Collision Mesh object from a full mesh vertices.
    /// @param include_vertex Vector of bools indicating whether each vertex should be included in the collision mesh.
    /// @param full_rest_positions The vertices of the full mesh at rest (#V × dim).
    /// @param edges The edges of the collision mesh indexed into the full mesh vertices (#E × 2).
    /// @param faces The faces of the collision mesh indexed into the full mesh vertices (#F × 3).
    /// @param displacement_map The displacement mapping from displacements on the full mesh to the collision mesh.
    CollisionMesh(
        const std::vector<bool>& include_vertex,
        const Eigen::MatrixXd& full_rest_positions,
        const Eigen::MatrixXi& edges = Eigen::MatrixXi(),
        const Eigen::MatrixXi& faces = Eigen::MatrixXi(),
        const Eigen::SparseMatrix<double>& displacement_map =
            Eigen::SparseMatrix<double>());

    /// @brief Helper function that automatically builds include_vertex using construct_is_on_surface.
    /// @param full_rest_positions The full vertices at rest (#FV × dim).
    /// @param edges The edge matrix of mesh (#E × 2).
    /// @param faces The face matrix of mesh (#F × 3).
    /// @return Constructed CollisionMesh.
    static CollisionMesh build_from_full_mesh(
        const Eigen::MatrixXd& full_rest_positions,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces = Eigen::MatrixXi())
    {
        return CollisionMesh(
            construct_is_on_surface(full_rest_positions.rows(), edges),
            full_rest_positions, edges, faces);
    }

    // The following functions are used to initialize optional data structures.

    /// @brief Initialize vertex-vertex and edge-vertex adjacencies.
    void init_adjacencies();

    /// @brief Initialize vertex and edge areas.
    void init_area_jacobians();

    /// @brief Destroy the Collision Mesh object
    ~CollisionMesh() = default;

    /// @brief Get the number of vertices in the collision mesh.
    size_t num_vertices() const { return m_vertex_to_full_vertex.size(); }

    /// @brief Get the number of codimensional vertices in the collision mesh.
    size_t num_codim_vertices() const { return codim_vertices().size(); }

    /// @brief Get the number of codimensional edges in the collision mesh.
    size_t num_codim_edges() const { return codim_edges().size(); }

    /// @brief Get the number of edges in the collision mesh.
    size_t num_edges() const { return edges().rows(); }

    /// @brief Get the number of faces in the collision mesh.
    size_t num_faces() const { return faces().rows(); }

    /// @brief Get the dimension of the mesh.
    size_t dim() const { return m_rest_positions.cols(); }

    /// @brief Get the number of degrees of freedom in the collision mesh.
    size_t ndof() const { return num_vertices() * dim(); }

    /// @brief Get the number of vertices in the full mesh.
    size_t full_num_vertices() const { return m_full_vertex_to_vertex.size(); }

    /// @brief Get the number of degrees of freedom in the full mesh.
    size_t full_ndof() const { return full_num_vertices() * dim(); }

    /// @brief Get the vertices of the collision mesh at rest (#V × dim).
    const Eigen::MatrixXd& rest_positions() const { return m_rest_positions; }

    /// @brief Get the indices of codimensional vertices of the collision mesh (#CV x 1).
    const Eigen::VectorXi& codim_vertices() const { return m_codim_vertices; }

    /// @brief Get the indices of codimensional edges of the collision mesh (#CE x 1).
    const Eigen::VectorXi& codim_edges() const { return m_codim_edges; }

    /// @brief Get the edges of the collision mesh (#E × 2).
    const Eigen::MatrixXi& edges() const { return m_edges; }

    /// @brief Get the faces of the collision mesh (#F × 3).
    const Eigen::MatrixXi& faces() const { return m_faces; }

    /// @brief Get the mapping from faces to edges of the collision mesh (#F × 3).
    const Eigen::MatrixXi& faces_to_edges() const { return m_faces_to_edges; }

    // const std::vector<std::vector<int>>& vertices_to_edges() const
    // {
    //     return m_vertices_to_edges;
    // }

    // const std::vector<std::vector<int>>& vertices_to_faces() const
    // {
    //     return m_vertices_to_faces;
    // }

    // -----------------------------------------------------------------------

    /// @brief Compute the vertex positions from the positions of the full mesh.
    /// @param full_positions The vertex positions of the full mesh (#FV × dim).
    /// @return The vertex positions of the collision mesh (#V × dim).
    Eigen::MatrixXd vertices(const Eigen::MatrixXd& full_positions) const;

    /// @brief Compute the vertex positions from vertex displacements on the full mesh.
    /// @param full_displacements The vertex displacements on the full mesh (#FV × dim).
    /// @return The vertex positions of the collision mesh (#V × dim).
    Eigen::MatrixXd
    displace_vertices(const Eigen::MatrixXd& full_displacements) const;

    /// @brief Map vertex displacements on the full mesh to vertex displacements on the collision mesh.
    /// @param full_displacements The vertex displacements on the full mesh (#FV × dim).
    /// @return The vertex displacements on the collision mesh (#V × dim).
    Eigen::MatrixXd
    map_displacements(const Eigen::MatrixXd& full_displacements) const;

    /// @brief Map a vertex ID to the corresponding vertex ID in the full mesh.
    /// @param id Vertex ID in the collision mesh.
    /// @return Vertex ID in the full mesh.
    size_t to_full_vertex_id(const size_t id) const
    {
        assert(id < num_vertices());
        return m_vertex_to_full_vertex[id];
    }

    /// @brief Map a vector quantity on the collision mesh to the full mesh.
    /// This is useful for mapping gradients from the collision mesh to the full
    /// mesh (i.e., applies the chain-rule).
    /// @param x Vector quantity on the collision mesh with size equal to ndof().
    /// @return Vector quantity on the full mesh with size equal to full_ndof().
    Eigen::VectorXd to_full_dof(const Eigen::VectorXd& x) const;

    /// @brief Map a matrix quantity on the collision mesh to the full mesh.
    /// This is useful for mapping Hessians from the collision mesh to the full
    /// mesh (i.e., applies the chain-rule).
    /// @param X Matrix quantity on the collision mesh with size equal to ndof() × ndof().
    /// @return Matrix quantity on the full mesh with size equal to full_ndof() × full_ndof().
    Eigen::SparseMatrix<double>
    to_full_dof(const Eigen::SparseMatrix<double>& X) const;

    // -----------------------------------------------------------------------

    /// @brief Get the vertex-vertex adjacency matrix.
    const std::vector<unordered_set<int>>& vertex_vertex_adjacencies() const
    {
        if (!are_adjacencies_initialized()) {
            throw std::runtime_error(
                "Vertex-vertex adjacencies not initialized. Call init_adjacencies() first.");
        }
        return m_vertex_vertex_adjacencies;
    }

    /// @brief Get the vertex-edge adjacency matrix.
    const std::vector<unordered_set<int>>& vertex_edge_adjacencies() const
    {
        if (!are_adjacencies_initialized()) {
            throw std::runtime_error(
                "Vertex-edge adjacencies not initialized. Call init_adjacencies() first.");
        }
        return m_vertex_edge_adjacencies;
    }

    /// @brief Get the edge-vertex adjacency matrix.
    const std::vector<unordered_set<int>>& edge_vertex_adjacencies() const
    {
        if (!are_adjacencies_initialized()) {
            throw std::runtime_error(
                "Edge-vertex adjacencies not initialized. Call init_adjacencies() first.");
        }
        return m_edge_vertex_adjacencies;
    }

    /// @brief Determine if the adjacencies have been initialized by calling init_adjacencies().
    bool are_adjacencies_initialized() const
    {
        return !m_vertex_vertex_adjacencies.empty()
            && !m_vertex_edge_adjacencies.empty()
            && !m_edge_vertex_adjacencies.empty();
    }

    /// @brief Is a vertex on the boundary of the collision mesh?
    /// @param vi Vertex ID.
    /// @return True if the vertex is on the boundary of the collision mesh.
    bool is_vertex_on_boundary(const int vi) const
    {
        return m_is_vertex_on_boundary[vi];
    }

    /// @brief Get the barycentric area of a vertex.
    /// @param vi Vertex ID.
    /// @return Barycentric area of vertex vi.
    double vertex_area(const size_t vi) const { return m_vertex_areas[vi]; }

    /// @brief Get the barycentric area of the vertices.
    const Eigen::VectorXd& vertex_areas() const { return m_vertex_areas; }

    /// @brief Get the gradient of the barycentric area of a vertex wrt the rest positions of all points.
    /// @param vi Vertex ID.
    /// @return Gradient of the barycentric area of vertex vi wrt the rest positions of all points.
    const Eigen::SparseVector<double>&
    vertex_area_gradient(const size_t vi) const
    {
        if (!are_area_jacobians_initialized()) {
            throw std::runtime_error(
                "Vertex area Jacobian not initialized. Call init_area_jacobians() first.");
        }
        return m_vertex_area_jacobian[vi];
    }

    /// @brief Get the barycentric area of an edge.
    /// @param ei Edge ID.
    /// @return Barycentric area of edge ei.
    double edge_area(const size_t ei) const { return m_edge_areas[ei]; }

    /// @brief Get the barycentric area of the edges.
    const Eigen::VectorXd& edge_areas() const { return m_edge_areas; }

    /// @brief Get the gradient of the barycentric area of an edge wrt the rest positions of all points.
    /// @param ei Edge ID.
    /// @return Gradient of the barycentric area of edge ei wrt the rest positions of all points.
    const Eigen::SparseVector<double>& edge_area_gradient(const size_t ei) const
    {
        if (!are_area_jacobians_initialized()) {
            throw std::runtime_error(
                "Edge area Jacobian not initialized. Call init_area_jacobians() first.");
        }
        return m_edge_area_jacobian[ei];
    }

    /// @brief Determine if the area Jacobians have been initialized by calling init_area_jacobians().
    bool are_area_jacobians_initialized() const
    {
        return m_vertex_area_jacobian.size() == num_vertices()
            && m_edge_area_jacobian.size() == num_edges();
    }

    // -----------------------------------------------------------------------

    /// @brief Construct a vector of bools indicating whether each vertex is on the surface.
    /// @param num_vertices The number of vertices in the mesh.
    /// @param edges The surface edges of the mesh (#E × 2).
    /// @param codim_vertices The indices of codimensional vertices (#CV x 1).
    /// @return A vector of bools indicating whether each vertex is on the surface.
    static std::vector<bool> construct_is_on_surface(
        const long num_vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::VectorXi& codim_vertices = Eigen::VectorXi());

    /// @brief Construct a matrix that maps from the faces' edges to rows in the edges matrix.
    /// @param faces The face matrix of mesh (#F × 3).
    /// @param edges The edge matrix of mesh (#E × 2).
    /// @return Matrix that maps from the faces' edges to rows in the edges matrix.
    static Eigen::MatrixXi construct_faces_to_edges(
        const Eigen::MatrixXi& faces, const Eigen::MatrixXi& edges);

    /// A function that takes two vertex IDs and returns true if the vertices
    /// (and faces or edges containing the vertices) can collide. By default all
    /// primitives can collide with all other primitives.
    std::function<bool(size_t, size_t)> can_collide = default_can_collide;

protected:
    // -----------------------------------------------------------------------
    // Helper initialization functions

    void init_codim_vertices();
    void init_codim_edges();

    /// @brief Initialize the selection matrix from full vertices/DOF to collision vertices/DOF.
    void init_selection_matrices(const int dim);

    /// @brief Initialize vertex and edge areas.
    void init_areas();

    /// @brief Convert a matrix meant for M_V * vertices to M_dof * x by duplicating the entries dim times.
    static Eigen::SparseMatrix<double> vertex_matrix_to_dof_matrix(
        const Eigen::SparseMatrix<double>& M_V, int dim);

    // -----------------------------------------------------------------------

    /// @brief The full vertex positions at rest (#FV × dim).
    Eigen::MatrixXd m_full_rest_positions;
    /// @brief The vertex positions at rest (#V × dim).
    Eigen::MatrixXd m_rest_positions;
    /// @brief The indices of codimensional vertices (#CV x 1).
    Eigen::VectorXi m_codim_vertices;
    /// @brief The indices of codimensional edges (#CE x 1).
    Eigen::VectorXi m_codim_edges;
    /// @brief Edges as rows of indicies into vertices (#E × 2).
    Eigen::MatrixXi m_edges;
    /// @brief Triangular faces as rows of indicies into vertices (#F × 3).
    Eigen::MatrixXi m_faces;
    /// @brief Map from faces edges to rows of edges (#F × 3).
    Eigen::MatrixXi m_faces_to_edges;

    /// @brief Map from full vertices to collision vertices.
    /// @note Negative values indicate full vertex is dropped.
    Eigen::VectorXi m_full_vertex_to_vertex;
    /// @brief Map from collision vertices to full vertices.
    Eigen::VectorXi m_vertex_to_full_vertex;

    /// @brief Selection matrix S ∈ ℝ^{collision×full} for vertices
    Eigen::SparseMatrix<double> m_select_vertices;
    /// @brief Selection matrix S ∈ ℝ^{(dim*collision)×(dim*full)} for DOF
    Eigen::SparseMatrix<double> m_select_dof;

    /// @brief Mapping from full displacements to collision displacements
    /// @note this is premultiplied by m_select_vertices
    Eigen::SparseMatrix<double> m_displacement_map;
    /// @brief Mapping from full displacements DOF to collision displacements DOF
    /// @note this is premultiplied by m_select_dof
    Eigen::SparseMatrix<double> m_displacement_dof_map;

    /// @brief Vertices adjacent to vertices
    std::vector<unordered_set<int>> m_vertex_vertex_adjacencies;
    /// @brief Edges adjacent to vertices
    std::vector<unordered_set<int>> m_vertex_edge_adjacencies;
    /// @brief Vertices adjacent to edges
    std::vector<unordered_set<int>> m_edge_vertex_adjacencies;

    // std::vector<std::vector<int>> m_vertices_to_faces;
    // std::vector<std::vector<int>> m_vertices_to_edges;

    /// @brief Is vertex on the boundary of the triangle mesh in 3D or polyline in 2D?
    std::vector<bool> m_is_vertex_on_boundary;

    /// @brief Vertex areas
    /// 2D: 1/2 sum of length of connected edges
    /// 3D: 1/3 sum of area of connected triangles
    Eigen::VectorXd m_vertex_areas;
    /// @brief Edge areas
    /// 3D: 1/3 sum of area of connected triangles
    Eigen::VectorXd m_edge_areas;

    // Stored as a std::vector so it is easier to access the rows directly.
    /// @brief The rows of the Jacobian of the vertex areas vector.
    std::vector<Eigen::SparseVector<double>> m_vertex_area_jacobian;
    /// @brief The rows of the Jacobian of the edge areas vector.
    std::vector<Eigen::SparseVector<double>> m_edge_area_jacobian;

private:
    /// @brief By default all primitives can collide with all other primitives.
    static int default_can_collide(size_t, size_t) { return true; }
};

} // namespace ipc
