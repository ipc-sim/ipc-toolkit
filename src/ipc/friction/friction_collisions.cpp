#include "friction_collisions.hpp"

#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

void FrictionCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Collisions& collisions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu)
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

        FC_vv.back().mu = blend_mu(mus(v0i), mus(v1i));
    }

    FC_ev.reserve(collisions.ev_collisions.size());
    for (const auto& c_ev : collisions.ev_collisions) {
        FC_ev.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

        const double edge_mu =
            (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
        FC_ev.back().mu = blend_mu(edge_mu, mus(vi));
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
        FC_ee.back().mu = blend_mu(ea_mu, eb_mu);
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
        FC_fv.back().mu = blend_mu(face_mu, mus(vi));
    }
}

void FrictionCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Collisions& collisions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const double static_mu,
    const double kinetic_mu,
    const std::map<std::tuple<int, int>, std::pair<double, double>>& pairwise_friction)
{
    clear();  // Clear any existing collisions

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    auto get_pairwise_friction = [&](int mat1, int mat2) -> std::pair<double, double> {
        auto it = pairwise_friction.find(std::make_tuple(mat1, mat2));
        if (it != pairwise_friction.end()) {
            return it->second;
        } else {
            return {static_mu, kinetic_mu};
        }
    };

    // Handle Vertex-Vertex Collisions
    for (const auto& c_vv : collisions.vv_collisions) {
        int v0i = c_vv.vertex_ids(edges, faces)[0];
        int v1i = c_vv.vertex_ids(edges, faces)[1];
        auto [pair_static_mu, pair_kinetic_mu] = get_pairwise_friction(v0i, v1i);

        vv_collisions.emplace_back(
            c_vv,
            c_vv.dof(vertices, edges, faces),
            barrier_potential,
            barrier_stiffness,
            pair_static_mu,
            pair_kinetic_mu
        );
    }

    // Handle Edge-Vertex Collisions
    for (const auto& c_ev : collisions.ev_collisions) {
        int vi = c_ev.vertex_ids(edges, faces)[0];
        int e0i = c_ev.vertex_ids(edges, faces)[1];
        int e1i = c_ev.vertex_ids(edges, faces)[2];

        auto [pair_static_mu, pair_kinetic_mu] = get_pairwise_friction(vi, std::min(e0i, e1i));

        ev_collisions.emplace_back(
            c_ev,
            c_ev.dof(vertices, edges, faces),
            barrier_potential,
            barrier_stiffness,
            pair_static_mu,
            pair_kinetic_mu
        );
    }

    // Handle Edge-Edge Collisions
    for (const auto& c_ee : collisions.ee_collisions) {
        int ea0i = c_ee.vertex_ids(edges, faces)[0];
        int ea1i = c_ee.vertex_ids(edges, faces)[1];
        int eb0i = c_ee.vertex_ids(edges, faces)[2];
        int eb1i = c_ee.vertex_ids(edges, faces)[3];

        auto [pair_static_mu, pair_kinetic_mu] = get_pairwise_friction(std::min(ea0i, ea1i), std::min(eb0i, eb1i));

        ee_collisions.emplace_back(
            c_ee,
            c_ee.dof(vertices, edges, faces),
            barrier_potential,
            barrier_stiffness,
            pair_static_mu,
            pair_kinetic_mu
        );
    }

    // Handle Face-Vertex Collisions
    for (const auto& c_fv : collisions.fv_collisions) {
        int vi = c_fv.vertex_ids(edges, faces)[0];
        int f0i = c_fv.vertex_ids(edges, faces)[1];
        int f1i = c_fv.vertex_ids(edges, faces)[2];
        int f2i = c_fv.vertex_ids(edges, faces)[3];

        auto [pair_static_mu, pair_kinetic_mu] = get_pairwise_friction(vi, std::min({f0i, f1i, f2i}));

        fv_collisions.emplace_back(
            c_fv,
            c_fv.dof(vertices, edges, faces),
            barrier_potential,
            barrier_stiffness,
            pair_static_mu,
            pair_kinetic_mu
        );
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

} // namespace ipc
