#include "tangential_collisions.hpp"

#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const NormalCollisions& collisions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double, BlendType)>& blend_mu,
    const BlendType blend_type)
{
    assert(mus.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    const auto& C_vv = collisions.vv_collisions;
    const auto& C_ev = collisions.ev_collisions;
    const auto& C_ee = collisions.ee_collisions;
    const auto& C_fv = collisions.fv_collisions;

    vv_collisions.reserve(C_vv.size());
    for (const auto& c_vv : C_vv) {
        vv_collisions.emplace_back(
            c_vv, c_vv.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [v0i, v1i, _, __] = vv_collisions.back().vertex_ids(edges, faces);

        vv_collisions.back().mu = blend_mu(mus(v0i), mus(v1i), blend_type);
        vv_collisions.back().s_mu = -1;
        vv_collisions.back().k_mu = -1;
    }

    ev_collisions.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        ev_collisions.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [vi, e0i, e1i, _] = ev_collisions.back().vertex_ids(edges, faces);

        const double edge_mu =
            (mus(e1i) - mus(e0i)) * ev_collisions.back().closest_point[0] + mus(e0i);
        ev_collisions.back().mu = blend_mu(edge_mu, mus(vi), blend_type);
        ev_collisions.back().s_mu = -1;
        ev_collisions.back().k_mu = -1;
    }

    ee_collisions.reserve(C_ee.size());
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

        ee_collisions.emplace_back(
            c_ee, c_ee.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);

        double ea_mu =
            (mus(ea1i) - mus(ea0i)) * ee_collisions.back().closest_point[0] + mus(ea0i);
        double eb_mu =
            (mus(eb1i) - mus(eb0i)) * ee_collisions.back().closest_point[1] + mus(eb0i);
        ee_collisions.back().mu = blend_mu(ea_mu, eb_mu, blend_type);
        ee_collisions.back().s_mu = -1;
        ee_collisions.back().k_mu = -1;
    }

    fv_collisions.reserve(C_fv.size());
    for (const auto& c_fv : C_fv) {
        fv_collisions.emplace_back(
            c_fv, c_fv.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [vi, f0i, f1i, f2i] = fv_collisions.back().vertex_ids(edges, faces);

        double face_mu = mus(f0i)
            + fv_collisions.back().closest_point[0] * (mus(f1i) - mus(f0i))
            + fv_collisions.back().closest_point[1] * (mus(f2i) - mus(f0i));
        fv_collisions.back().mu = blend_mu(face_mu, mus(vi), blend_type);
        fv_collisions.back().s_mu = -1;
        fv_collisions.back().k_mu = -1;
    }
}

// The other build functions go here as previously defined

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
