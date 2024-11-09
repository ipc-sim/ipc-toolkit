#include "friction_collisions.hpp"

#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <ipc/friction/smooth_friction_mollifier.hpp>

#include <stdexcept> // std::out_of_range
#include <optional>

namespace ipc {

void FrictionCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Collisions& collisions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double, std::optional<BlendType>)>&
        blend_mu)
{
    assert(mus.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    auto& FC_vv = vv_collisions;
    auto& FC_ev = ev_collisions;
    auto& FC_ee = ee_collisions;
    auto& FC_fv = fv_collisions;

    FC_vv.reserve(collisions.vv_collisions.size());
    for (const auto& c_vv : collisions.vv_collisions) {
        FC_vv.emplace_back(
            c_vv, c_vv.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

        FC_vv.back().mu = blend_mu(mus(v0i), mus(v1i), std::nullopt);
    }

    FC_ev.reserve(collisions.ev_collisions.size());
    for (const auto& c_ev : collisions.ev_collisions) {
        FC_ev.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

        const double edge_mu =
            (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
        FC_ev.back().mu = blend_mu(edge_mu, mus(vi), std::nullopt);
    }

    FC_ee.reserve(collisions.ee_collisions.size());
    for (const auto& c_ee : collisions.ee_collisions) {
        const auto& [ea0i, ea1i, eb0i, eb1i] = c_ee.vertex_ids(edges, faces);
        const Eigen::Vector3d ea0 = vertices.row(ea0i);
        const Eigen::Vector3d ea1 = vertices.row(ea1i);
        const Eigen::Vector3d eb0 = vertices.row(eb0i);
        const Eigen::Vector3d eb1 = vertices.row(eb1i);

        // Skip EE collisions that are close to parallel
        if (edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1) < c_ee.eps_x) {
            continue;
        }

        FC_ee.emplace_back(
            c_ee, c_ee.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);

        double ea_mu =
            (mus(ea1i) - mus(ea0i)) * FC_ee.back().closest_point[0] + mus(ea0i);
        double eb_mu =
            (mus(eb1i) - mus(eb0i)) * FC_ee.back().closest_point[1] + mus(eb0i);
        FC_ee.back().mu = blend_mu(ea_mu, eb_mu, std::nullopt);
    }

    FC_fv.reserve(collisions.fv_collisions.size());
    for (const auto& c_fv : collisions.fv_collisions) {
        FC_fv.emplace_back(
            c_fv, c_fv.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [vi, f0i, f1i, f2i] = FC_fv.back().vertex_ids(edges, faces);

        double face_mu = mus(f0i)
            + FC_fv.back().closest_point[0] * (mus(f1i) - mus(f0i))
            + FC_fv.back().closest_point[1] * (mus(f2i) - mus(f0i));
        FC_fv.back().mu = blend_mu(face_mu, mus(vi), std::nullopt);
    }
}

void FrictionCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Collisions& collisions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const double mu,
    const double static_mu,
    const double kinetic_mu,
    const std::map<std::tuple<int, int>, std::pair<double, double>>& pairwise_friction,
    const std::function<double(double, double, std::optional<BlendType>)>& blend_mu)
{
    // Clear any existing collisions
    clear();

    // Default blend_mu to transition behavior if not provided
    auto effective_blend_mu = blend_mu;
    if (!effective_blend_mu) {
        effective_blend_mu = [static_mu, kinetic_mu](double mu1, double mu2, std::optional<BlendType>) {
            return ipc::blend_mu(mu1, mu2, BlendType::TRANSITION);
        };
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    if (mu < 0 || static_mu < 0 || kinetic_mu < 0) {
        throw std::invalid_argument(
            "Friction coefficient must be non-negative.");
    }

    // Handle Vertex-Vertex Collisions
    for (const auto& c_vv : collisions.vv_collisions) {
        int v0i = c_vv.vertex_ids(edges, faces)[0];
        int v1i = c_vv.vertex_ids(edges, faces)[1];
        auto [pair_static_mu, pair_kinetic_mu] =
            get_pairwise_friction_coefficients(
                mesh, v0i, v1i, pairwise_friction, static_mu, kinetic_mu);

        vv_collisions.emplace_back(
            c_vv, c_vv.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness, pair_static_mu, pair_kinetic_mu);
    }

    // Handle Edge-Vertex Collisions
    for (const auto& c_ev : collisions.ev_collisions) {
        int vi = c_ev.vertex_ids(edges, faces)[0];
        int e0i = c_ev.vertex_ids(edges, faces)[1];
        int e1i = c_ev.vertex_ids(edges, faces)[2];

        // Determine material IDs for the edge vertices
        int mat_e0 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[e0i];
        int mat_e1 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[e1i];

        // Choose a representative material ID for the edge (e.g., minimum)
        int mat_edge = std::min(mat_e0, mat_e1);

        // Retrieve friction coefficients between vertex and edge
        auto [pair_static_mu, pair_kinetic_mu] =
            get_pairwise_friction_coefficients(
                mesh, vi, mat_edge, pairwise_friction, static_mu, kinetic_mu);

        ev_collisions.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness, pair_static_mu, pair_kinetic_mu);
    }

    // Handle Edge-Edge Collisions
    for (const auto& c_ee : collisions.ee_collisions) {
        int ea0i = c_ee.vertex_ids(edges, faces)[0];
        int ea1i = c_ee.vertex_ids(edges, faces)[1];
        int eb0i = c_ee.vertex_ids(edges, faces)[2];
        int eb1i = c_ee.vertex_ids(edges, faces)[3];

        // Determine material IDs for both edges
        int mat_ea0 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[ea0i];
        int mat_ea1 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[ea1i];
        int mat_eb0 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[eb0i];
        int mat_eb1 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[eb1i];

        // Choose representative material IDs for each edge (e.g., minimum)
        int mat_edge_a = std::min(mat_ea0, mat_ea1);
        int mat_edge_b = std::min(mat_eb0, mat_eb1);

        // Retrieve friction coefficients between the two edges
        auto [pair_static_mu, pair_kinetic_mu] =
            get_pairwise_friction_coefficients(
                mesh, mat_edge_a, mat_edge_b, pairwise_friction, static_mu,
                kinetic_mu);

        ee_collisions.emplace_back(
            c_ee, c_ee.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness, pair_static_mu, pair_kinetic_mu);
    }

    // Handle Face-Vertex Collisions
    for (const auto& c_fv : collisions.fv_collisions) {
        int vi = c_fv.vertex_ids(edges, faces)[0];
        int f0i = c_fv.vertex_ids(edges, faces)[1];
        int f1i = c_fv.vertex_ids(edges, faces)[2];
        int f2i = c_fv.vertex_ids(edges, faces)[3];

        // Determine material IDs for the face vertices
        int mat_f0 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[f0i];
        int mat_f1 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[f1i];
        int mat_f2 = (mesh.vertex_material_ids().empty())
            ? 0
            : mesh.vertex_material_ids()[f2i];

        // Choose a representative material ID for the face (e.g., minimum)
        int mat_face = std::min({ mat_f0, mat_f1, mat_f2 });

        // Retrieve friction coefficients between vertex and face
        auto [pair_static_mu, pair_kinetic_mu] =
            get_pairwise_friction_coefficients(
                mesh, vi, mat_face, pairwise_friction, static_mu, kinetic_mu);

        fv_collisions.emplace_back(
            c_fv, c_fv.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness, pair_static_mu, pair_kinetic_mu);
    }
}

// ============================================================================

size_t FrictionCollisions::size() const
{
    return vv_collisions.size() + ev_collisions.size() + ee_collisions.size()
        + fv_collisions.size();
}

bool FrictionCollisions::empty() const
{
    return vv_collisions.empty() && ev_collisions.empty()
        && ee_collisions.empty() && fv_collisions.empty();
}

void FrictionCollisions::clear()
{
    vv_collisions.clear();
    ev_collisions.clear();
    ee_collisions.clear();
    fv_collisions.clear();
}

FrictionCollision& FrictionCollisions::operator[](size_t i)
{
    if (i < vv_collisions.size()) {
        return vv_collisions[i];
    }
    i -= vv_collisions.size();
    if (i < ev_collisions.size()) {
        return ev_collisions[i];
    }
    i -= ev_collisions.size();
    if (i < ee_collisions.size()) {
        return ee_collisions[i];
    }
    i -= ee_collisions.size();
    if (i < fv_collisions.size()) {
        return fv_collisions[i];
    }
    throw std::out_of_range("Friction collision index is out of range!");
}

const FrictionCollision& FrictionCollisions::operator[](size_t i) const
{
    if (i < vv_collisions.size()) {
        return vv_collisions[i];
    }
    i -= vv_collisions.size();
    if (i < ev_collisions.size()) {
        return ev_collisions[i];
    }
    i -= ev_collisions.size();
    if (i < ee_collisions.size()) {
        return ee_collisions[i];
    }
    i -= ee_collisions.size();
    if (i < fv_collisions.size()) {
        return fv_collisions[i];
    }
    throw std::out_of_range("Friction collision index is out of range!");
}

std::pair<double, double> FrictionCollisions::retrieve_friction_coefficients(
    int id1,
    int id2,
    const std::map<std::tuple<int, int>, std::pair<double, double>>&
        pairwise_friction,
    std::optional<double> static_mu,
    std::optional<double> kinetic_mu) const
{
    auto it = pairwise_friction.find(std::make_tuple(id1, id2));
    if (it != pairwise_friction.end()) {
        return it->second;
    }
    if (static_mu && kinetic_mu) {
        return { static_mu.value(), kinetic_mu.value() };
    }
    throw std::runtime_error(
        "No friction coefficients available for the given pair.");
}

std::pair<double, double>
FrictionCollisions::get_pairwise_friction_coefficients(
    const CollisionMesh& mesh,
    int vertex1,
    int vertex2,
    const std::map<std::tuple<int, int>, std::pair<double, double>>&
        pairwise_friction,
    double default_static_mu,
    double default_kinetic_mu) const
{
    // Retrieve material IDs for the vertices
    int mat1 = (mesh.vertex_material_ids().empty())
        ? 0
        : mesh.vertex_material_ids()[vertex1];
    int mat2 = (mesh.vertex_material_ids().empty())
        ? 0
        : mesh.vertex_material_ids()[vertex2];

    // Order the material pair consistently
    std::tuple<int, int> ordered_material_pair = { mat1, mat2 };
    if (mat1 < mat2) {
        ordered_material_pair = { mat2, mat1 };
    } else {
        ordered_material_pair = { mat1, mat2 };
    }

    auto it = pairwise_friction.find(ordered_material_pair);
    if (it != pairwise_friction.end()) {
        return it->second;
    } else {
        return { default_static_mu, default_kinetic_mu };
    }
}

} // namespace ipc
