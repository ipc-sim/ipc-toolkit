#include "tangential_collisions.hpp"

#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const NormalCollisions& collisions,
    const NormalPotential& normal_potential,
    const double normal_stiffness,
    Eigen::ConstRef<Eigen::VectorXd> mus,
    const std::function<double(double, double)>& blend_mu)
{
    assert(mus.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    const auto& C_vv = collisions.vv_collisions;
    const auto& C_ev = collisions.ev_collisions;
    const auto& C_ee = collisions.ee_collisions;
    const auto& C_fv = collisions.fv_collisions;

    auto& FC_vv = this->vv_collisions;
    auto& FC_ev = this->ev_collisions;
    auto& FC_ee = this->ee_collisions;
    auto& FC_fv = this->fv_collisions;

    // Check if we need to use material IDs
    bool use_materials = mesh.has_material_ids();

    FC_vv.reserve(C_vv.size());
    for (const auto& c_vv : C_vv) {
        FC_vv.emplace_back(
            c_vv, c_vv.dof(vertices, edges, faces), normal_potential,
            normal_stiffness);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

        // Only set material IDs if needed
        if (use_materials) {
            FC_vv.back().material_id1 = mesh.vertex_material(v0i);
            FC_vv.back().material_id2 = mesh.vertex_material(v1i);
        }

        FC_vv.back().mu = default_blend_mu(mus(v0i), mus(v1i), blend_type);
        FC_vv.back().s_mu = -1;
        FC_vv.back().k_mu = -1;
    }

    FC_ev.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        FC_ev.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), normal_potential,
            normal_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

        // Only set material IDs if needed
        if (use_materials) {
            FC_ev.back().material_id1 = mesh.edge_material(FC_ev.back().edge_id);
            FC_ev.back().material_id2 = mesh.vertex_material(vi);
        }

        const double edge_mu =
            (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
        FC_ev.back().mu = default_blend_mu(edge_mu, mus(vi), blend_type);
        FC_ev.back().s_mu = -1;
        FC_ev.back().k_mu = -1;
    }

    FC_ee.reserve(C_ee.size());
    for (const auto& c_ee : C_ee) {
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
            c_ee, c_ee.dof(vertices, edges, faces), normal_potential,
            normal_stiffness);

        // Only set material IDs if needed
        if (use_materials) {
            FC_ee.back().material_id1 = mesh.edge_material(FC_ee.back().edge0_id);
            FC_ee.back().material_id2 = mesh.edge_material(FC_ee.back().edge1_id);
        }

        double ea_mu =
            (mus(ea1i) - mus(ea0i)) * FC_ee.back().closest_point[0] + mus(ea0i);
        double eb_mu =
            (mus(eb1i) - mus(eb0i)) * FC_ee.back().closest_point[1] + mus(eb0i);
        FC_ee.back().mu = default_blend_mu(ea_mu, eb_mu, blend_type);
        FC_ee.back().s_mu = -1;
        FC_ee.back().k_mu = -1;
    }

    FC_fv.reserve(C_fv.size());
    for (const auto& c_fv : C_fv) {
        FC_fv.emplace_back(
            c_fv, c_fv.dof(vertices, edges, faces), normal_potential,
            normal_stiffness);
        const auto& [vi, f0i, f1i, f2i] = FC_fv.back().vertex_ids(edges, faces);

        // Only set material IDs if needed
        if (use_materials) {
            FC_fv.back().material_id1 = mesh.face_material(FC_fv.back().face_id);
            FC_fv.back().material_id2 = mesh.vertex_material(vi);
        }

        double face_mu = mus(f0i)
            + FC_fv.back().closest_point[0] * (mus(f1i) - mus(f0i))
            + FC_fv.back().closest_point[1] * (mus(f2i) - mus(f0i));
        FC_fv.back().mu = default_blend_mu(face_mu, mus(vi), blend_type);
        FC_fv.back().s_mu = -1;
        FC_fv.back().k_mu = -1;
    }
}

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const NormalCollisions& collisions,
    const NormalPotential& normal_potential,
    double barrier_stiffness,
    double mu,
    double s_mu,
    double k_mu)
{
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();
    clear();

    auto setFrictionParams = [](auto& collision, double g_mu, double g_s_mu, double g_k_mu) {
        collision.mu = g_mu;
        collision.s_mu = g_s_mu;
        collision.k_mu = g_k_mu;
    };

    // Check if we need to use material IDs
    bool use_materials = mesh.has_material_ids();

    const auto& C_vv = collisions.vv_collisions;
    const auto& C_ev = collisions.ev_collisions;
    const auto& C_ee = collisions.ee_collisions;
    const auto& C_fv = collisions.fv_collisions;

    auto& FC_vv = this->vv_collisions;
    auto& FC_ev = this->ev_collisions;
    auto& FC_ee = this->ee_collisions;
    auto& FC_fv = this->fv_collisions;

    FC_vv.reserve(C_vv.size());
    for (const auto& c_vv : C_vv) {
        FC_vv.emplace_back(
            c_vv, c_vv.dof(vertices, edges, faces), normal_potential,
            barrier_stiffness);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);
        
        // Only set material IDs if needed
        if (use_materials) {
            FC_vv.back().material_id1 = mesh.vertex_material(v0i);
            FC_vv.back().material_id2 = mesh.vertex_material(v1i);
        }
        
        setFrictionParams(FC_vv.back(), mu, s_mu, k_mu);
    }

    FC_ev.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        FC_ev.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), normal_potential,
            barrier_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);
        
        // Only set material IDs if needed
        if (use_materials) {
            FC_ev.back().material_id1 = mesh.edge_material(FC_ev.back().edge_id);
            FC_ev.back().material_id2 = mesh.vertex_material(vi);
        }
        
        setFrictionParams(FC_ev.back(), mu, s_mu, k_mu);
    }

    FC_ee.reserve(C_ee.size());
    for (const auto& c_ee : C_ee) {
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
            c_ee, c_ee.dof(vertices, edges, faces), normal_potential,
            barrier_stiffness);

        // Only set material IDs if needed
        if (use_materials) {
            FC_ee.back().material_id1 = mesh.edge_material(FC_ee.back().edge0_id);
            FC_ee.back().material_id2 = mesh.edge_material(FC_ee.back().edge1_id);
        }

        setFrictionParams(FC_ee.back(), mu, s_mu, k_mu);
    }

    FC_fv.reserve(C_fv.size());
    for (const auto& c_fv : C_fv) {
        FC_fv.emplace_back(
            c_fv, c_fv.dof(vertices, edges, faces), normal_potential,
            barrier_stiffness);
        const auto& [vi, f0i, f1i, f2i] = FC_fv.back().vertex_ids(edges, faces);

        // Only set material IDs if needed
        if (use_materials) {
            FC_fv.back().material_id1 = mesh.face_material(FC_fv.back().face_id);
            FC_fv.back().material_id2 = mesh.vertex_material(vi);
        }

        setFrictionParams(FC_fv.back(), mu, s_mu, k_mu);
    }
}

// ============================================================================

size_t TangentialCollisions::size() const
{
    return vv_collisions.size() + ev_collisions.size() + ee_collisions.size()
        + fv_collisions.size();
}

bool TangentialCollisions::empty() const
{
    return vv_collisions.empty() && ev_collisions.empty()
        && ee_collisions.empty() && fv_collisions.empty();
}

void TangentialCollisions::clear()
{
    vv_collisions.clear();
    ev_collisions.clear();
    ee_collisions.clear();
    fv_collisions.clear();
}

TangentialCollision& TangentialCollisions::operator[](size_t i)
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

const TangentialCollision& TangentialCollisions::operator[](size_t i) const
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

} // namespace ipc
