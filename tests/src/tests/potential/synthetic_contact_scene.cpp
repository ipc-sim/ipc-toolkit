#include "synthetic_contact_scene.hpp"

#include <ipc/distance/edge_edge_mollifier.hpp>

#include <spdlog/fmt/fmt.h>

#include <cassert>
#include <cmath>
#include <random>
#include <stdexcept>

namespace ipc::tests {

namespace {

    /// @brief Multiplicative hash stride used to decorrelate the two sides of
    /// a sampled collision pair (Knuth).
    constexpr size_t HASH_STRIDE = 2654435761u;

} // namespace

AssemblyScene build_synthetic_contact_scene(const SyntheticContactSpec& spec)
{
    const int dim = spec.dim;
    if (dim == 2 && spec.type != SyntheticContactType::EDGE_VERTEX) {
        throw std::invalid_argument("2D synthetic scenes must be EDGE_VERTEX");
    }
    if (dim == 3 && spec.type == SyntheticContactType::EDGE_VERTEX) {
        throw std::invalid_argument("EDGE_VERTEX synthetic scenes must be 2D");
    }

    // Element structure below needs a handful of vertices to exist.
    constexpr size_t MIN_VERTICES = 8;
    const size_t num_verts = std::max(
        MIN_VERTICES, size_t(std::llround(double(spec.target_ndof) / dim)));

    // All vertices inside a unit box, and dhat larger than the box diagonal:
    // every sampled collision is active no matter which vertices it touches.
    const double dhat = 3.0;

    std::mt19937 rng(spec.seed);
    std::uniform_real_distribution<double> unit(-0.5, 0.5);

    Eigen::MatrixXd vertices(num_verts, dim);
    for (Eigen::Index i = 0; i < vertices.size(); i++) {
        vertices.data()[i] = unit(rng);
    }

    // --- Elements and collisions --------------------------------------------
    //
    // Elements (faces/edges) are built from disjoint vertex tuples at the
    // front of the vertex pool; the remaining "free" vertices are used as the
    // point side of point-element collisions. Collisions are sampled by
    // striding through the element/free pools, so vertices are shared between
    // collisions once num_collisions exceeds the pool sizes (the dense
    // self-contact regime).

    Eigen::MatrixXi edges, faces;
    NormalCollisions collisions;

    switch (spec.type) {
    case SyntheticContactType::FACE_VERTEX: {
        const size_t n_faces = std::max(size_t(1), num_verts / 6);
        assert(3 * n_faces < num_verts);
        const size_t free_start = 3 * n_faces;
        const size_t n_free = num_verts - free_start;

        faces.resize(n_faces, 3);
        edges.resize(3 * n_faces, 2);
        for (size_t j = 0; j < n_faces; j++) {
            faces.row(j) << int(3 * j), int(3 * j + 1), int(3 * j + 2);
            edges.row(3 * j) << int(3 * j), int(3 * j + 1);
            edges.row(3 * j + 1) << int(3 * j + 1), int(3 * j + 2);
            edges.row(3 * j + 2) << int(3 * j + 2), int(3 * j);
        }

        collisions.fv_collisions.reserve(spec.num_collisions);
        for (size_t i = 0; i < spec.num_collisions; i++) {
            const index_t face_id = index_t(i % n_faces);
            const index_t vertex_id =
                index_t(free_start + (i * HASH_STRIDE) % n_free);
            collisions.fv_collisions.emplace_back(
                FaceVertexCandidate(face_id, vertex_id));
        }
        break;
    }

    case SyntheticContactType::EDGE_EDGE: {
        const size_t n_edges = std::max(size_t(2), num_verts / 2);
        assert(2 * n_edges <= num_verts);

        edges.resize(n_edges, 2);
        for (size_t j = 0; j < n_edges; j++) {
            edges.row(j) << int(2 * j), int(2 * j + 1);
        }

        collisions.ee_collisions.reserve(spec.num_collisions);
        for (size_t i = 0; i < spec.num_collisions; i++) {
            const size_t ea = i % n_edges;
            const size_t eb =
                (ea + 1 + (i * HASH_STRIDE) % (n_edges - 1)) % n_edges;
            assert(ea != eb);
            const double eps_x = edge_edge_mollifier_threshold(
                vertices.row(edges(ea, 0)), vertices.row(edges(ea, 1)),
                vertices.row(edges(eb, 0)), vertices.row(edges(eb, 1)));
            collisions.ee_collisions.emplace_back(
                index_t(ea), index_t(eb), eps_x);
        }
        break;
    }

    case SyntheticContactType::EDGE_VERTEX: {
        const size_t n_edges = std::max(size_t(1), num_verts / 4);
        assert(2 * n_edges < num_verts);
        const size_t free_start = 2 * n_edges;
        const size_t n_free = num_verts - free_start;

        edges.resize(n_edges, 2);
        for (size_t j = 0; j < n_edges; j++) {
            edges.row(j) << int(2 * j), int(2 * j + 1);
        }

        collisions.ev_collisions.reserve(spec.num_collisions);
        for (size_t i = 0; i < spec.num_collisions; i++) {
            const index_t edge_id = index_t(i % n_edges);
            const index_t vertex_id =
                index_t(free_start + (i * HASH_STRIDE) % n_free);
            collisions.ev_collisions.emplace_back(
                EdgeVertexCandidate(edge_id, vertex_id));
        }
        break;
    }

    default:
        throw std::invalid_argument("unknown SyntheticContactType");
    }

    // --- Collision mesh ------------------------------------------------------
    //
    // Pad the full mesh with interior vertices (see AssemblySceneSpec), and
    // include every non-padding vertex in the collision mesh explicitly (some
    // "free" vertices belong to no element, and construct_is_on_surface would
    // drop them).

    const size_t num_interior =
        size_t(std::llround(spec.interior_vertex_ratio * num_verts));
    Eigen::MatrixXd full_vertices(num_verts + num_interior, dim);
    full_vertices.topRows(num_verts) = vertices;
    if (num_interior > 0) {
        full_vertices.bottomRows(num_interior).rowwise() =
            vertices.colwise().mean();
    }

    std::vector<bool> include_vertex(full_vertices.rows(), false);
    for (size_t i = 0; i < num_verts; i++) {
        include_vertex[i] = true;
    }

    Eigen::SparseMatrix<double> displacement_map;
    if (spec.with_displacement_map) {
        // Identity: a semantic no-op that disables the pure-selection DOF map
        // (see SyntheticContactSpec::with_displacement_map).
        displacement_map.resize(full_vertices.rows(), full_vertices.rows());
        displacement_map.setIdentity();
    }

    const CollisionMesh mesh(
        include_vertex, std::vector<bool>(full_vertices.rows(), false),
        full_vertices, edges, faces, displacement_map);
    assert(mesh.num_vertices() == num_verts);
    assert(mesh.is_selection_dof_map() == !spec.with_displacement_map);

    const BarrierPotential potential(dhat, /*stiffness=*/1.0);

    return AssemblyScene(
        fmt::format(
            "synthetic-{}-{}d-c{}-n{}{}", to_string(spec.type), dim,
            spec.num_collisions, num_verts * dim,
            spec.with_displacement_map ? "-dispmap" : ""),
        mesh, mesh.vertices(full_vertices), std::move(collisions), potential);
}

} // namespace ipc::tests
