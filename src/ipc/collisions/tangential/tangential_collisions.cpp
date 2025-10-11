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
    const double _normal_stiffness,
    Eigen::ConstRef<Eigen::VectorXd> mu_s,
    Eigen::ConstRef<Eigen::VectorXd> mu_k,
    const std::function<double(double, double)>& blend_mu)
{
    normal_stiffness = _normal_stiffness;
    assert(mu_s.size() == vertices.rows());
    assert(mu_k.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    const auto& C_vv = collisions.vv_collisions;
    const auto& C_ev = collisions.ev_collisions;
    const auto& C_ee = collisions.ee_collisions;
    const auto& C_fv = collisions.fv_collisions;
    auto& [FC_vv, FC_ev, FC_ee, FC_fv, kappa] = *this;

    FC_vv.reserve(C_vv.size());
    for (const auto& c_vv : C_vv) {
        FC_vv.emplace_back(
            c_vv, c_vv.dof(vertices, edges, faces), normal_potential,
            normal_stiffness);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

        FC_vv.back().mu_s = blend_mu(mu_s(v0i), mu_s(v1i));
        FC_vv.back().mu_k = blend_mu(mu_k(v0i), mu_k(v1i));
    }

    FC_ev.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        FC_ev.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), normal_potential,
            normal_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

        const double edge_mu_s =
            (mu_s(e1i) - mu_s(e0i)) * FC_ev.back().closest_point[0] + mu_s(e0i);
        FC_ev.back().mu_s = blend_mu(edge_mu_s, mu_s(vi));
        const double edge_mu_k =
            (mu_k(e1i) - mu_k(e0i)) * FC_ev.back().closest_point[0] + mu_k(e0i);
        FC_ev.back().mu_k = blend_mu(edge_mu_k, mu_k(vi));
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

        double ea_mu_s =
            (mu_s(ea1i) - mu_s(ea0i)) * FC_ee.back().closest_point[0]
            + mu_s(ea0i);
        double eb_mu_s =
            (mu_s(eb1i) - mu_s(eb0i)) * FC_ee.back().closest_point[1]
            + mu_s(eb0i);
        FC_ee.back().mu_s = blend_mu(ea_mu_s, eb_mu_s);

        double ea_mu_k =
            (mu_k(ea1i) - mu_k(ea0i)) * FC_ee.back().closest_point[0]
            + mu_k(ea0i);
        double eb_mu_k =
            (mu_k(eb1i) - mu_k(eb0i)) * FC_ee.back().closest_point[1]
            + mu_k(eb0i);
        FC_ee.back().mu_k = blend_mu(ea_mu_k, eb_mu_k);
    }

    FC_fv.reserve(C_fv.size());
    for (const auto& c_fv : C_fv) {
        FC_fv.emplace_back(
            c_fv, c_fv.dof(vertices, edges, faces), normal_potential,
            normal_stiffness);
        const auto& [vi, f0i, f1i, f2i] = FC_fv.back().vertex_ids(edges, faces);

        double face_mu_s = mu_s(f0i)
            + FC_fv.back().closest_point[0] * (mu_s(f1i) - mu_s(f0i))
            + FC_fv.back().closest_point[1] * (mu_s(f2i) - mu_s(f0i));
        FC_fv.back().mu_s = blend_mu(face_mu_s, mu_s(vi));

        double face_mu_k = mu_k(f0i)
            + FC_fv.back().closest_point[0] * (mu_k(f1i) - mu_k(f0i))
            + FC_fv.back().closest_point[1] * (mu_k(f2i) - mu_k(f0i));
        FC_fv.back().mu_k = blend_mu(face_mu_k, mu_k(vi));
    }
}

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const SmoothCollisions& collisions,
    const ParameterType& params,
    const double _normal_stiffness,
    Eigen::ConstRef<Eigen::VectorXd> mu_s,
    Eigen::ConstRef<Eigen::VectorXd> mu_k,
    const std::function<double(double, double)>& blend_mu)
{
    normal_stiffness = _normal_stiffness;
    assert(mu_k.size() == vertices.rows());
    assert(mu_s.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    auto& [FC_vv, FC_ev, FC_ee, FC_fv, kappa] = *this;

    // FC_vv.reserve(C_vv.size());
    for (size_t i = 0; i < collisions.size(); i++) {
        const auto& cc = collisions[i];
        Eigen::VectorXd contact_potential_grad =
            cc.gradient(cc.dof(vertices), params);
        const double contact_force =
            normal_stiffness * contact_potential_grad.norm();

        if (mesh.dim() == 3) {
            TangentialCollision* ptr = nullptr;
            if (const auto* const cvv = dynamic_cast<
                    const SmoothCollisionTemplate<Point3, Point3>*>(&cc)) {
                Eigen::VectorXd collision_points = cvv->core_dof(vertices);
                FC_vv.emplace_back(
                    VertexVertexNormalCollision(
                        cc[0], cc[1], 1., Eigen::SparseVector<double>()),
                    collision_points, contact_force);
                const auto& [v0i, v1i, _, __] =
                    FC_vv.back().vertex_ids(edges, faces);

                FC_vv.back().mu_k = blend_mu(mu_k(v0i), mu_k(v1i));
                FC_vv.back().mu_s = blend_mu(mu_s(v0i), mu_s(v1i));
                ptr = &(FC_vv.back());
            } else if (
                const auto* const cev =
                    dynamic_cast<const SmoothCollisionTemplate<Edge3, Point3>*>(
                        &cc)) {
                Eigen::VectorXd collision_points = cev->core_dof(vertices);
                collision_points =
                    collision_points({ 6, 7, 8, 0, 1, 2, 3, 4, 5 })
                        .eval(); // {edge, point} -> {point, edge}
                FC_ev.emplace_back(
                    EdgeVertexNormalCollision(
                        cc[0], cc[1], 1., Eigen::SparseVector<double>()),
                    collision_points, contact_force);
                const auto& [vi, e0i, e1i, _] =
                    FC_ev.back().vertex_ids(edges, faces);

                const double edge_mu_k =
                    (mu_k(e1i) - mu_k(e0i)) * FC_ev.back().closest_point[0]
                    + mu_k(e0i);
                FC_ev.back().mu_k = blend_mu(edge_mu_k, mu_k(vi));

                const double edge_mu_s =
                    (mu_s(e1i) - mu_s(e0i)) * FC_ev.back().closest_point[0]
                    + mu_s(e0i);
                FC_ev.back().mu_s = blend_mu(edge_mu_s, mu_s(vi));

                ptr = &(FC_ev.back());
            } else if (
                const auto* const cee =
                    dynamic_cast<const SmoothCollisionTemplate<Edge3, Edge3>*>(
                        &cc)) {
                Eigen::VectorXd collision_points = cee->core_dof(vertices);
                const auto vert_ids = cee->core_vertex_ids();
                const Eigen::Vector3d ea0 = vertices.row(vert_ids[0]);
                const Eigen::Vector3d ea1 = vertices.row(vert_ids[1]);
                const Eigen::Vector3d eb0 = vertices.row(vert_ids[2]);
                const Eigen::Vector3d eb1 = vertices.row(vert_ids[3]);

                // Skip EE collisions that are close to parallel
                if (edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1)
                    < edge_edge_mollifier_threshold(ea0, ea1, eb0, eb1)) {
                    continue;
                }

                FC_ee.emplace_back(
                    EdgeEdgeNormalCollision(
                        cc[0], cc[1], 0., EdgeEdgeDistanceType::EA_EB),
                    collision_points, contact_force);

                double ea_mu_k = (mu_k(vert_ids[1]) - mu_k(vert_ids[0]))
                        * FC_ee.back().closest_point[0]
                    + mu_k(vert_ids[0]);
                double eb_mu_k = (mu_k(vert_ids[3]) - mu_k(vert_ids[2]))
                        * FC_ee.back().closest_point[1]
                    + mu_k(vert_ids[2]);
                FC_ee.back().mu_k = blend_mu(ea_mu_k, eb_mu_k);

                double ea_mu_s = (mu_s(vert_ids[1]) - mu_s(vert_ids[0]))
                        * FC_ee.back().closest_point[0]
                    + mu_s(vert_ids[0]);
                double eb_mu_s = (mu_s(vert_ids[3]) - mu_s(vert_ids[2]))
                        * FC_ee.back().closest_point[1]
                    + mu_s(vert_ids[2]);
                FC_ee.back().mu_s = blend_mu(ea_mu_s, eb_mu_s);

                ptr = &(FC_ee.back());
            } else if (
                const auto* const cfv =
                    dynamic_cast<const SmoothCollisionTemplate<Face, Point3>*>(
                        &cc)) {
                Eigen::VectorXd collision_points = cfv->core_dof(vertices);
                collision_points =
                    collision_points({ 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8 })
                        .eval(); // {face, point} -> {point, face}
                FC_fv.emplace_back(
                    FaceVertexNormalCollision(
                        cc[0], cc[1], 1., Eigen::SparseVector<double>()),
                    collision_points, contact_force);
                const auto& [vi, f0i, f1i, f2i] =
                    FC_fv.back().vertex_ids(edges, faces);

                double face_mu_k = mu_k(f0i)
                    + FC_fv.back().closest_point[0] * (mu_k(f1i) - mu_k(f0i))
                    + FC_fv.back().closest_point[1] * (mu_k(f2i) - mu_k(f0i));
                FC_fv.back().mu_k = blend_mu(face_mu_k, mu_k(vi));

                double face_mu_s = mu_s(f0i)
                    + FC_fv.back().closest_point[0] * (mu_s(f1i) - mu_s(f0i))
                    + FC_fv.back().closest_point[1] * (mu_s(f2i) - mu_s(f0i));
                FC_fv.back().mu_s = blend_mu(face_mu_s, mu_s(vi));

                ptr = &(FC_fv.back());
            }
            if (ptr) {
                ptr->smooth_collision = collisions.collisions[i];
            }
        } else {
            TangentialCollision* ptr = nullptr;
            if (const auto* const cvv = dynamic_cast<
                    const SmoothCollisionTemplate<Point2, Point2>*>(&cc)) {
                Eigen::VectorXd collision_points = cvv->core_dof(vertices);
                FC_vv.emplace_back(
                    VertexVertexNormalCollision(
                        cc[0], cc[1], 1., Eigen::SparseVector<double>()),
                    collision_points, contact_force);
                const auto& [v0i, v1i, _, __] =
                    FC_vv.back().vertex_ids(edges, faces);

                FC_vv.back().mu_s = blend_mu(mu_s(v0i), mu_s(v1i));
                FC_vv.back().mu_k = blend_mu(mu_k(v0i), mu_k(v1i));
                ptr = &(FC_vv.back());
            } else if (
                const auto* const cev =
                    dynamic_cast<const SmoothCollisionTemplate<Edge2, Point2>*>(
                        &cc)) {
                Eigen::VectorXd collision_points = cev->core_dof(vertices);
                collision_points =
                    collision_points({ 4, 5, 0, 1, 2, 3 })
                        .eval(); // {edge, point} -> {point, edge}
                FC_ev.emplace_back(
                    EdgeVertexNormalCollision(
                        cc[0], cc[1], 1., Eigen::SparseVector<double>()),
                    collision_points, contact_force);
                const auto& [vi, e0i, e1i, _] =
                    FC_ev.back().vertex_ids(edges, faces);

                const double edge_mu_k =
                    (mu_k(e1i) - mu_k(e0i)) * FC_ev.back().closest_point[0]
                    + mu_k(e0i);
                FC_ev.back().mu_k = blend_mu(edge_mu_k, mu_k(vi));

                const double edge_mu_s =
                    (mu_s(e1i) - mu_s(e0i)) * FC_ev.back().closest_point[0]
                    + mu_s(e0i);
                FC_ev.back().mu_s = blend_mu(edge_mu_s, mu_s(vi));

                ptr = &(FC_ev.back());
            }
            if (ptr) {
                ptr->smooth_collision = collisions.collisions[i];
            }
        }
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
