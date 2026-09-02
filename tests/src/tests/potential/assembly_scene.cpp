#include "assembly_scene.hpp"

#include <tests/utils.hpp>

#include <spdlog/fmt/fmt.h>

namespace ipc::tests {

AssemblyScene::AssemblyScene(
    std::string label,
    const CollisionMesh& mesh,
    Eigen::MatrixXd vertices,
    NormalCollisions collisions,
    const BarrierPotential& potential)
    : m_label(std::move(label))
    , m_mesh(mesh)
    , m_vertices(std::move(vertices))
    , m_collisions(std::move(collisions))
    , m_potential(potential)
{
}

std::string AssemblyScene::stats() const
{
    return fmt::format(
        "{}: {} collisions | {} collision vertices ({} DOF) | "
        "{} full vertices ({} DOF) | dim={}",
        label(), num_collisions(), num_vertices(), ndof(), full_num_vertices(),
        full_ndof(), m_mesh.dim());
}

std::optional<AssemblyScene> build_assembly_scene(const AssemblySceneSpec& spec)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    if (!load_mesh(spec.mesh_name, vertices, edges, faces)) {
        return std::nullopt;
    }

    // Pad the full mesh with interior vertices that no edge or face
    // references, so that the collision mesh is a strict subset of the full
    // mesh (see AssemblySceneSpec::interior_vertex_ratio).
    const Eigen::Index num_surface_vertices = vertices.rows();
    const Eigen::Index num_interior_vertices = Eigen::Index(
        std::llround(spec.interior_vertex_ratio * num_surface_vertices));

    Eigen::MatrixXd full_vertices(
        num_surface_vertices + num_interior_vertices, vertices.cols());
    full_vertices.topRows(num_surface_vertices) = vertices;
    if (num_interior_vertices > 0) {
        // Position is irrelevant (these are excluded from the collision mesh),
        // but the centroid keeps any bounding-box computation well behaved.
        full_vertices.bottomRows(num_interior_vertices).rowwise() =
            vertices.colwise().mean();
    }

    const CollisionMesh mesh(
        CollisionMesh::construct_is_on_surface(full_vertices.rows(), edges),
        std::vector<bool>(full_vertices.rows(), false), full_vertices, edges,
        faces);

    Eigen::MatrixXd collision_vertices = mesh.vertices(full_vertices);

    NormalCollisions collisions;
    collisions.build(mesh, collision_vertices, spec.dhat);
    if (collisions.empty()) {
        return std::nullopt;
    }

    // A stiffness of 1 keeps the numbers interpretable; assembly cost is
    // independent of its value.
    const BarrierPotential potential(spec.dhat, /*stiffness=*/1.0);

    std::optional<AssemblyScene> scene;
    scene.emplace(
        spec.label.empty() ? spec.mesh_name : spec.label, mesh,
        std::move(collision_vertices), std::move(collisions), potential);
    return scene;
}

const std::vector<AssemblySceneSpec>& assembly_scene_specs()
{
    // dhat values were chosen so each scene has a non-trivial collision set;
    // see the "Assembly scene statistics" test, which prints the actual counts.
    // WARNING: increasing `dhat` grows the collision set superlinearly. The
    // values below are measured (see the "Assembly scene statistics" test) and
    // span 390 -> 512k collisions. Do not raise them without checking the
    // resulting candidate/collision count first: `dhat` an order of magnitude
    // larger on `cloth_ball92.ply` exhausts memory on a 64 GB host.
    //
    // The simulation-frame scenes use dhat = 1e-3 * bbox diagonal.
    static const std::vector<AssemblySceneSpec> specs = {
        { "two-cubes-close.ply", 1e-1, 1.0, "two-cubes" },
        { "bunny.ply", 1e-2, 1.0, "bunny" },
        { "cloth_ball92.ply", 1e-3, 1.0, "cloth-ball" },
        { "cloth-funnel/227.ply", 4.0e-3, 1.0, "cloth-funnel" },
        { "armadillo-rollers/326.ply", 1e-3, 1.0, "armadillo-rollers" },
        { "n-body-simulation/balls16_18.ply", 1.78e-2, 1.0, "n-body" },
        { "rod-twist/3036.ply", 1.09e-3, 1.0, "rod-twist" },
        { "puffer-ball/20.ply", 1.44e-4, 1.0, "puffer-ball" },
    };
    return specs;
}

} // namespace ipc::tests
