#include "friction_collisions.hpp"

#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

template <int dim>
void FrictionCollisions::build_for_smooth_contact(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const SmoothCollisions<dim>& collisions,
    const ParameterType &params,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu)
{
    barrier_stiffness_ = barrier_stiffness;
    assert(mus.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    auto& [FC_vv, FC_ev, FC_ee, FC_fv, kappa] = *this;

    // FC_vv.reserve(C_vv.size());
    for (size_t i = 0; i < collisions.size(); i++) {
        const auto& cc = collisions[i];
        Eigen::VectorXd contact_potential_grad = cc.gradient(cc.dof(vertices, edges, faces), params);
        const double contact_force = barrier_stiffness * contact_potential_grad.norm();

        if constexpr (dim == 3)
        {
            FrictionCollision* ptr = nullptr;
            if (const auto cvv = dynamic_cast<const SmoothCollisionTemplate<max_vert_3d, Point3, Point3> *>(&cc))
            {
                Eigen::VectorXd collision_points = cvv->core_dof(vertices, edges, faces);
                FC_vv.emplace_back(
                    VertexVertexCollision(cc[0], cc[1], 1., Eigen::SparseVector<double>()), collision_points, contact_force);
                const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

                FC_vv.back().mu = blend_mu(mus(v0i), mus(v1i));
                ptr = &(FC_vv.back());
            }
            else if (const auto cev = dynamic_cast<const SmoothCollisionTemplate<max_vert_3d, Edge3, Point3> *>(&cc))
            {
                Eigen::VectorXd collision_points = cev->core_dof(vertices, edges, faces);
                collision_points = collision_points({6,7,8,0,1,2,3,4,5}).eval();
                FC_ev.emplace_back(
                    EdgeVertexCollision(cc[0], cc[1], 1., Eigen::SparseVector<double>()), collision_points, contact_force);
                const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

                const double edge_mu =
                    (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
                FC_ev.back().mu = blend_mu(edge_mu, mus(vi));
                ptr = &(FC_ev.back());
            }
            else if (const auto cee = dynamic_cast<const SmoothCollisionTemplate<max_vert_3d, Edge3, Edge3> *>(&cc))
            {
                Eigen::VectorXd collision_points = cee->core_dof(vertices, edges, faces);
                const auto vert_ids = cee->core_vertex_ids(edges, faces);
                const Eigen::Vector3d ea0 = vertices.row(vert_ids[0]);
                const Eigen::Vector3d ea1 = vertices.row(vert_ids[1]);
                const Eigen::Vector3d eb0 = vertices.row(vert_ids[2]);
                const Eigen::Vector3d eb1 = vertices.row(vert_ids[3]);

                // Skip EE collisions that are close to parallel
                if (edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1) < edge_edge_mollifier_threshold(ea0, ea1, eb0, eb1)) {
                    continue;
                }

                FC_ee.emplace_back(
                    EdgeEdgeCollision(cc[0], cc[1], 0., EdgeEdgeDistanceType::EA_EB), collision_points, contact_force);

                double ea_mu =
                    (mus(vert_ids[1]) - mus(vert_ids[0])) * FC_ee.back().closest_point[0] + mus(vert_ids[0]);
                double eb_mu =
                    (mus(vert_ids[3]) - mus(vert_ids[2])) * FC_ee.back().closest_point[1] + mus(vert_ids[2]);
                FC_ee.back().mu = blend_mu(ea_mu, eb_mu);
                ptr = &(FC_ee.back());
            }
            else if (const auto cfv = dynamic_cast<const SmoothCollisionTemplate<max_vert_3d, Face, Point3> *>(&cc))
            {
                Eigen::VectorXd collision_points = cfv->core_dof(vertices, edges, faces);
                collision_points = collision_points({9,10,11,0,1,2,3,4,5,6,7,8}).eval();
                FC_fv.emplace_back(
                    FaceVertexCollision(cc[0], cc[1], 1., Eigen::SparseVector<double>()), collision_points, contact_force);
                const auto& [vi, f0i, f1i, f2i] = FC_fv.back().vertex_ids(edges, faces);

                double face_mu = mus(f0i)
                    + FC_fv.back().closest_point[0] * (mus(f1i) - mus(f0i))
                    + FC_fv.back().closest_point[1] * (mus(f2i) - mus(f0i));
                FC_fv.back().mu = blend_mu(face_mu, mus(vi));
                ptr = &(FC_fv.back());
            }
            if (ptr)
                ptr->smooth_collision_3d = collisions.collisions[i];
        }
        else
        {
            FrictionCollision* ptr = nullptr;
            if (const auto cvv = dynamic_cast<const SmoothCollisionTemplate<max_vert_2d, Point2, Point2> *>(&cc))
            {
                Eigen::VectorXd collision_points = cvv->core_dof(vertices, edges, faces);
                FC_vv.emplace_back(
                    VertexVertexCollision(cc[0], cc[1], 1., Eigen::SparseVector<double>()), collision_points, contact_force);
                const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

                FC_vv.back().mu = blend_mu(mus(v0i), mus(v1i));
                ptr = &(FC_vv.back());
            }
            else if (const auto cev = dynamic_cast<const SmoothCollisionTemplate<max_vert_2d, Edge2, Point2> *>(&cc))
            {
                Eigen::VectorXd collision_points = cev->core_dof(vertices, edges, faces);
                collision_points = collision_points({4,5,0,1,2,3}).eval();
                FC_ev.emplace_back(
                    EdgeVertexCollision(cc[0], cc[1], 1., Eigen::SparseVector<double>()), collision_points, contact_force);
                const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

                const double edge_mu =
                    (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
                FC_ev.back().mu = blend_mu(edge_mu, mus(vi));
                ptr = &(FC_ev.back());
            }
            if (ptr)
                ptr->smooth_collision_2d = collisions.collisions[i];
        }
    }
}

template
void FrictionCollisions::build_for_smooth_contact<2>(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const SmoothCollisions<2>& collisions,
    const ParameterType &params,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu);

template
void FrictionCollisions::build_for_smooth_contact<3>(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const SmoothCollisions<3>& collisions,
    const ParameterType &params,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu);

void FrictionCollisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Collisions& collisions,
    const BarrierPotential& barrier_potential,
    const double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu)
{
    barrier_stiffness_ = barrier_stiffness;
    assert(mus.size() == vertices.rows());

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
            c_vv, c_vv.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

        FC_vv.back().mu = blend_mu(mus(v0i), mus(v1i));
    }

    FC_ev.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        FC_ev.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);
        const auto& [vi, e0i, e1i, _] = FC_ev.back().vertex_ids(edges, faces);

        const double edge_mu =
            (mus(e1i) - mus(e0i)) * FC_ev.back().closest_point[0] + mus(e0i);
        FC_ev.back().mu = blend_mu(edge_mu, mus(vi));
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
            c_ee, c_ee.dof(vertices, edges, faces), barrier_potential,
            barrier_stiffness);

        double ea_mu =
            (mus(ea1i) - mus(ea0i)) * FC_ee.back().closest_point[0] + mus(ea0i);
        double eb_mu =
            (mus(eb1i) - mus(eb0i)) * FC_ee.back().closest_point[1] + mus(eb0i);
        FC_ee.back().mu = blend_mu(ea_mu, eb_mu);
    }

    FC_fv.reserve(C_fv.size());
    for (const auto& c_fv : C_fv) {
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
