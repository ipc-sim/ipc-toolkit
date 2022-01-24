#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

class SurfaceMesh {
public:
    SurfaceMesh() { }

    SurfaceMesh(
        const Eigen::MatrixXd& surface_vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces);

    SurfaceMesh(
        const std::vector<bool>& is_on_surface,
        const Eigen::MatrixXd& full_vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces);

    ~SurfaceMesh() { }

    const Eigen::MatrixXd& vertices_at_rest() const
    {
        return m_vertices_at_rest;
    }
    const Eigen::MatrixXi& edges() const { return m_edges; }
    const Eigen::MatrixXi& faces() const { return m_faces; }
    const Eigen::MatrixXi& faces_to_edges() const { return m_faces_to_edges; }

    size_t surface_size() const { return m_surface_to_full.size(); }
    size_t full_size() const { return m_full_size; }
    const Eigen::VectorXi& surface_to_full() const { return m_surface_to_full; }
    const Eigen::VectorXi& full_to_surface() const { return m_full_to_surface; }

    Eigen::VectorXd map_surface_to_full(const Eigen::VectorXd& x) const;
    Eigen::SparseMatrix<double>
    map_surface_to_full(const Eigen::SparseMatrix<double>& X) const;

    Eigen::MatrixXd surface_vertices(const Eigen::MatrixXd& full_V) const;

    static std::vector<bool> construct_is_on_surface(
        const int num_vertices, const Eigen::MatrixXi& edges);

    static SurfaceMesh build_from_full_mesh(
        const Eigen::MatrixXd& full_vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces)
    {
        return SurfaceMesh(
            construct_is_on_surface(full_vertices_at_rest.rows(), edges),
            full_vertices_at_rest, edges, faces);
    }

    /// Construct a matrix that maps from the faces' edges to rows in the edges
    /// matrix.
    static Eigen::MatrixXi construct_faces_to_edges(
        const Eigen::MatrixXi& faces, const Eigen::MatrixXi& edges);

    /// A function that takes two vertex IDs (row numbers in V) and returns true
    /// if the vertices (and faces or edges containing the vertices) can
    /// collide. By default all primitives can collide with all other
    /// primitives.
    std::function<bool(size_t, size_t)> can_collide = [](size_t, size_t) {
        return true;
    };

protected:
    Eigen::MatrixXd m_vertices_at_rest;
    /// @param[in] E Edges as rows of indicies into V.
    Eigen::MatrixXi m_edges;
    /// @param[in] F Triangular faces as rows of indicies into V.
    Eigen::MatrixXi m_faces;
    /// @param[in] F2E Map from F edges to rows of E.
    Eigen::MatrixXi m_faces_to_edges;

    size_t m_full_size;
    Eigen::VectorXi m_surface_to_full;
    Eigen::VectorXi m_full_to_surface;
};

} // namespace ipc