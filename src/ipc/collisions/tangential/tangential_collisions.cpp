#include "tangential_collisions.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/friction/smooth_mu.hpp>
#include <ipc/esp/collisions/esp_quadrature.hpp>
#include <ipc/esp/collisions/vertex_matrix_view.hpp>
#include <ipc/esp/esp_parameters.hpp>
#include <ipc/esp/esp_potential.hpp>
#include <ipc/math/math.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/gcp/distance/edge_edge.hpp>
#include <ipc/tangent/closest_point.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <atomic>
#include <stdexcept> // std::out_of_range
#include <tuple>

namespace ipc {

namespace {

    void update_lagged_mu_collision(
        TangentialCollision& collision,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
        Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
        Eigen::ConstRef<Eigen::MatrixXd> velocities)
    {
        const VectorMax12d rest_dof =
            collision.dof(rest_positions, edges, faces);
        const VectorMax12d lag_dof =
            collision.dof(lagged_displacements, edges, faces);
        const VectorMax12d vel_dof = collision.dof(velocities, edges, faces);
        const VectorMax12d lagged_pos = rest_dof + lag_dof;

        const MatrixMax<double, 3, 2> P =
            collision.compute_tangent_basis(lagged_pos);
        const VectorMax2d beta = collision.compute_closest_point(lagged_pos);
        const MatrixMax<double, 3, 12> Gamma =
            collision.relative_velocity_jacobian(beta);
        const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;
        const VectorMax2d tau = T.transpose() * vel_dof;

        const int tangent_dim = tau.size();
        const VectorMax2d tau_aniso =
            collision.mu_aniso.head(tangent_dim).cwiseProduct(tau);

        if (tangent_dim <= 1 || collision.mu_s_aniso.squaredNorm() == 0) {
            collision.mu_s_effective_lagged = collision.mu_s;
            collision.mu_k_effective_lagged = collision.mu_k;
            return;
        }

        const Eigen::Vector2d tau_aniso_2d = tau_aniso;
        std::tie(
            collision.mu_s_effective_lagged, collision.mu_k_effective_lagged) =
            anisotropic_mu_eff_from_tau_aniso(
                tau_aniso_2d, collision.mu_s_aniso, collision.mu_k_aniso,
                collision.mu_s, collision.mu_k, false);
    }

} // namespace

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const NormalCollisions& collisions,
    const NormalPotential& normal_potential,
    Eigen::ConstRef<Eigen::VectorXd> mu_s,
    Eigen::ConstRef<Eigen::VectorXd> mu_k,
    const std::function<double(double, double)>& blend_mu)
{
    assert(mu_s.size() == vertices.rows());
    assert(mu_k.size() == vertices.rows());
    IPC_TOOLKIT_PROFILE_BLOCK("TangentialCollisions::build");

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    const auto& C_vv = collisions.vv_collisions;
    const auto& C_ev = collisions.ev_collisions;
    const auto& C_ee = collisions.ee_collisions;
    const auto& C_fv = collisions.fv_collisions;
    const auto& C_pv = collisions.pv_collisions;
    auto& [FC_vv, FC_ev, FC_ee, FC_fv, FC_pv] = *this;

    FC_vv.reserve(C_vv.size());
    for (const auto& c_vv : C_vv) {
        FC_vv.emplace_back(
            c_vv, c_vv.dof(vertices, edges, faces), normal_potential);
        const auto& [v0i, v1i, _, __] = FC_vv.back().vertex_ids(edges, faces);

        FC_vv.back().mu_s = blend_mu(mu_s(v0i), mu_s(v1i));
        FC_vv.back().mu_k = blend_mu(mu_k(v0i), mu_k(v1i));
    }

    FC_ev.reserve(C_ev.size());
    for (const auto& c_ev : C_ev) {
        FC_ev.emplace_back(
            c_ev, c_ev.dof(vertices, edges, faces), normal_potential);
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
            c_ee, c_ee.dof(vertices, edges, faces), normal_potential);

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
            c_fv, c_fv.dof(vertices, edges, faces), normal_potential);
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

    FC_pv.reserve(C_pv.size());
    for (const auto& c_pv : C_pv) {
        FC_pv.emplace_back(
            c_pv, c_pv.dof(vertices, edges, faces), normal_potential);
        const auto& [vi, _0, _1, _2] = FC_pv.back().vertex_ids(edges, faces);
        FC_pv.back().mu_s = mu_s(vi);
        FC_pv.back().mu_k = mu_k(vi);
    }

    reset_lagged_anisotropic_friction_coefficients();
}

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const GcpCollisions& collisions,
    const GcpParameters& params,
    const double normal_stiffness,
    Eigen::ConstRef<Eigen::VectorXd> mu_s,
    Eigen::ConstRef<Eigen::VectorXd> mu_k,
    const std::function<double(double, double)>& blend_mu)
{
    assert(mu_k.size() == vertices.rows());
    assert(mu_s.size() == vertices.rows());
    IPC_TOOLKIT_PROFILE_BLOCK("TangentialCollisions::build(smooth)");

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    clear();

    auto& [FC_vv, FC_ev, FC_ee, FC_fv, FC_pv] = *this;

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
                    const GcpCollisionTemplate<Point3, Point3>*>(&cc)) {
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
                    dynamic_cast<const GcpCollisionTemplate<Edge3, Point3>*>(
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
                    dynamic_cast<const GcpCollisionTemplate<Edge3, Edge3>*>(
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
                    dynamic_cast<const GcpCollisionTemplate<Face, Point3>*>(
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
                ptr->gcp_collision = collisions.collisions[i];
            }
        } else {
            TangentialCollision* ptr = nullptr;
            if (const auto* const cvv = dynamic_cast<
                    const GcpCollisionTemplate<Point2, Point2>*>(&cc)) {
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
                    dynamic_cast<const GcpCollisionTemplate<Edge2, Point2>*>(
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
                ptr->gcp_collision = collisions.collisions[i];
            }
        }
    }

    reset_lagged_anisotropic_friction_coefficients();
}

void TangentialCollisions::reset_lagged_anisotropic_friction_coefficients()
{
    for (auto& c : vv_collisions) {
        c.mu_s_effective_lagged = c.mu_s;
        c.mu_k_effective_lagged = c.mu_k;
    }
    for (auto& c : ev_collisions) {
        c.mu_s_effective_lagged = c.mu_s;
        c.mu_k_effective_lagged = c.mu_k;
    }
    for (auto& c : ee_collisions) {
        c.mu_s_effective_lagged = c.mu_s;
        c.mu_k_effective_lagged = c.mu_k;
    }
    for (auto& c : fv_collisions) {
        c.mu_s_effective_lagged = c.mu_s;
        c.mu_k_effective_lagged = c.mu_k;
    }
    for (auto& c : pv_collisions) {
        c.mu_s_effective_lagged = c.mu_s;
        c.mu_k_effective_lagged = c.mu_k;
    }
}

void TangentialCollisions::update_lagged_anisotropic_friction_coefficients(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
    Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
    Eigen::ConstRef<Eigen::MatrixXd> velocities)
{
    assert(rest_positions.rows() == lagged_displacements.rows());
    assert(rest_positions.rows() == velocities.rows());
    assert(rest_positions.cols() == lagged_displacements.cols());
    assert(rest_positions.cols() == velocities.cols());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    for (auto& c : vv_collisions) {
        update_lagged_mu_collision(
            c, edges, faces, rest_positions, lagged_displacements, velocities);
    }
    for (auto& c : ev_collisions) {
        update_lagged_mu_collision(
            c, edges, faces, rest_positions, lagged_displacements, velocities);
    }
    for (auto& c : ee_collisions) {
        update_lagged_mu_collision(
            c, edges, faces, rest_positions, lagged_displacements, velocities);
    }
    for (auto& c : fv_collisions) {
        update_lagged_mu_collision(
            c, edges, faces, rest_positions, lagged_displacements, velocities);
    }
    for (auto& c : pv_collisions) {
        update_lagged_mu_collision(
            c, edges, faces, rest_positions, lagged_displacements, velocities);
    }
}

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const EspCollisions& collisions,
    const EspParameters& params,
    const double normal_stiffness,
    Eigen::ConstRef<Eigen::VectorXd> mu_s,
    Eigen::ConstRef<Eigen::VectorXd> mu_k,
    const bool normalize_weights,
    const std::function<double(double, double)>& blend_mu)
{
    assert(mu_s.size() == vertices.rows());
    assert(mu_k.size() == vertices.rows());

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();
    const int dim = mesh.dim();

    clear();

    auto& [FC_vv, FC_ev, FC_ee, FC_fv, FC_pv] = *this;

    auto assign_ev_mu = [&](EdgeVertexTangentialCollision& tc) {
        const auto& [vi, e0i, e1i, _] = tc.vertex_ids(edges, faces);
        const double edge_mu_s =
            (mu_s(e1i) - mu_s(e0i)) * tc.closest_point[0] + mu_s(e0i);
        tc.mu_s = blend_mu(edge_mu_s, mu_s(vi));
        const double edge_mu_k =
            (mu_k(e1i) - mu_k(e0i)) * tc.closest_point[0] + mu_k(e0i);
        tc.mu_k = blend_mu(edge_mu_k, mu_k(vi));
    };

    if (dim == 2) {
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(params.quad_order);
        const index_t n_verts = vertices.rows();

        auto compute_contact_force_2d = [&](const EspCollision& cc,
                                            const VertexMatrixView<2>& V_ext,
                                            const double outer_w) -> double {
            const Eigen::VectorXd positions = cc.dof(V_ext);
            double d2 = 0;
            switch (cc.type()) {
            case EspCollisionType::VERTEX_VERTEX:
                d2 = point_point_distance(
                    Eigen::Vector2d(positions.segment<2>(0)),
                    Eigen::Vector2d(positions.segment<2>(2)));
                break;
            case EspCollisionType::EDGE_VERTEX:
                d2 = point_edge_distance(
                    Eigen::Vector2d(positions.segment<2>(0)),
                    Eigen::Vector2d(positions.segment<2>(2)),
                    Eigen::Vector2d(positions.segment<2>(4)));
                break;
            default:
                return 0;
            }
            const double dist = std::sqrt(d2);
            const double dhat_val = params.dhat;
            return (dist > 0 && dist < dhat_val) ? outer_w * normal_stiffness
                    * std::abs(params.barrier->first_derivative(dist, dhat_val))
                                                 : 0.0;
        };

        for (const auto& [ei, qp_dicts] : collisions.edge_collisions_2d) {
            const index_t e0 = edges(ei, 0);
            const index_t e1 = edges(ei, 1);
            const double L = mesh.edge_length(ei);

            for (size_t qi = 0; qi < qp_dicts.size(); ++qi) {
                const auto& dict_ptr = qp_dicts[qi];
                if (!dict_ptr || dict_ptr->size() == 0)
                    continue;
                const auto& qp = rule[qi];
                const std::array<double, 2> lambda = { { 1.0 - qp.xi, qp.xi } };
                const Eigen::RowVector2d virtual_pos =
                    lambda[0] * vertices.row(e0) + lambda[1] * vertices.row(e1);
                VertexMatrixView<2> V_ext(vertices, virtual_pos);
                const double outer_w = L * qp.weight;

                for (int j = 0; j < dict_ptr->size(); ++j) {
                    const auto& cc = (*dict_ptr)[j];
                    const double contact_force =
                        compute_contact_force_2d(cc, V_ext, outer_w);
                    if (contact_force == 0)
                        continue;

                    switch (cc.type()) {
                    case EspCollisionType::VERTEX_VERTEX: {
                        // One vertex is virtual (n_verts), other is real.
                        // Elevate to EdgeVertex: edge ei vs the real vertex.
                        const index_t v0 = cc.vertex_id(0);
                        const index_t v1 = cc.vertex_id(1);
                        const index_t v_real = (v0 == n_verts) ? v1 : v0;

                        Eigen::Matrix<double, 6, 1> cp;
                        cp.segment<2>(0) = vertices.row(v_real).transpose();
                        cp.segment<2>(2) = vertices.row(e0).transpose();
                        cp.segment<2>(4) = vertices.row(e1).transpose();

                        FC_ev.emplace_back(
                            EdgeVertexNormalCollision(
                                ei, v_real, cc.weight,
                                Eigen::SparseVector<double>()),
                            cp, contact_force);
                        FC_ev.back().weight = cc.weight;
                        assign_ev_mu(FC_ev.back());
                        break;
                    }
                    case EspCollisionType::EDGE_VERTEX: {
                        // Virtual vertex on edge ei vs real edge ej.
                        // Distribute to the two endpoints of ej weighted by
                        // the projection parameter u (parallel to 3D face
                        // dict EV elevation).
                        const index_t ej = cc[1]; // Edge2P1 id
                        const index_t ea = cc.vertex_id(1);
                        const index_t eb = cc.vertex_id(2);
                        const Eigen::Vector2d ea_pos =
                            vertices.row(ea).transpose();
                        const Eigen::Vector2d eb_pos =
                            vertices.row(eb).transpose();
                        const Eigen::Vector2d vp = virtual_pos.transpose();

                        double u = point_edge_closest_point(vp, ea_pos, eb_pos);
                        if (!std::isfinite(u))
                            break;
                        u = std::clamp(u, 0.0, 1.0);

                        auto emit_ev = [&](index_t v_edge, double w) {
                            if (w <= 0)
                                return;
                            Eigen::Matrix<double, 6, 1> cp;
                            cp.segment<2>(0) = vertices.row(v_edge).transpose();
                            cp.segment<2>(2) = vertices.row(e0).transpose();
                            cp.segment<2>(4) = vertices.row(e1).transpose();
                            FC_ev.emplace_back(
                                EdgeVertexNormalCollision(
                                    ei, v_edge, cc.weight,
                                    Eigen::SparseVector<double>()),
                                cp, w * contact_force);
                            FC_ev.back().weight = cc.weight;
                            assign_ev_mu(FC_ev.back());
                        };
                        emit_ev(ea, 1.0 - u);
                        emit_ev(eb, u);
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
        }
    } else {
        // 3D: collisions are stored in per-primitive dicts.
        // VERTEX dicts have only real vertex IDs.
        // EDGE and FACE dicts contain a virtual vertex (ID = n_verts)
        // representing an edge-edge closest point or face quadrature point.
        // We "elevate" sub-collisions with virtual vertices to higher-type
        // tangential collisions that use only real vertex IDs.
        const index_t n_verts = vertices.rows();

        // Helper: compute contact force magnitude for a sub-collision.
        //
        // Uses the scalar derivative of the log-barrier w.r.t. distance,
        // scaled by outer quadrature weight and barrier stiffness.
        auto compute_contact_force = [&](const EspCollision& cc,
                                         const VertexMatrixView<3>& V_ext,
                                         const double outer_w) -> double {
            const Eigen::VectorXd positions = cc.dof(V_ext);
            double d2 = 0;
            switch (cc.type()) {
            case EspCollisionType::VERTEX_VERTEX:
                d2 = point_point_distance(
                    Eigen::Vector3d(positions.segment<3>(0)),
                    Eigen::Vector3d(positions.segment<3>(3)));
                break;
            case EspCollisionType::EDGE_VERTEX:
                d2 = point_edge_distance(
                    Eigen::Vector3d(positions.segment<3>(6)),
                    Eigen::Vector3d(positions.segment<3>(0)),
                    Eigen::Vector3d(positions.segment<3>(3)));
                break;
            case EspCollisionType::FACE_VERTEX:
                d2 = point_triangle_distance(
                    Eigen::Vector3d(positions.segment<3>(9)),
                    Eigen::Vector3d(positions.segment<3>(0)),
                    Eigen::Vector3d(positions.segment<3>(3)),
                    Eigen::Vector3d(positions.segment<3>(6)));
                break;
            default:
                return 0;
            }
            const double dist = sqrt(d2);
            const double dhat_val = params.dhat;
            return (dist > 0 && dist < dhat_val) ? outer_w * normal_stiffness
                    * std::abs(params.barrier->first_derivative(dist, dhat_val))
                                                 : 0.0;
        };

        // Helper: assign EE mu values
        auto assign_ee_mu = [&](EdgeEdgeTangentialCollision& tc) {
            const auto& [ea0i, ea1i, eb0i, eb1i] = tc.vertex_ids(edges, faces);
            double ea_mu_s =
                (mu_s(ea1i) - mu_s(ea0i)) * tc.closest_point[0] + mu_s(ea0i);
            double eb_mu_s =
                (mu_s(eb1i) - mu_s(eb0i)) * tc.closest_point[1] + mu_s(eb0i);
            tc.mu_s = blend_mu(ea_mu_s, eb_mu_s);
            double ea_mu_k =
                (mu_k(ea1i) - mu_k(ea0i)) * tc.closest_point[0] + mu_k(ea0i);
            double eb_mu_k =
                (mu_k(eb1i) - mu_k(eb0i)) * tc.closest_point[1] + mu_k(eb0i);
            tc.mu_k = blend_mu(ea_mu_k, eb_mu_k);
        };

        // Helper: assign FV mu values
        auto assign_fv_mu = [&](FaceVertexTangentialCollision& tc) {
            const auto& [vi, f0i, f1i, f2i] = tc.vertex_ids(edges, faces);
            double face_mu_s = mu_s(f0i)
                + tc.closest_point[0] * (mu_s(f1i) - mu_s(f0i))
                + tc.closest_point[1] * (mu_s(f2i) - mu_s(f0i));
            tc.mu_s = blend_mu(face_mu_s, mu_s(vi));
            double face_mu_k = mu_k(f0i)
                + tc.closest_point[0] * (mu_k(f1i) - mu_k(f0i))
                + tc.closest_point[1] * (mu_k(f2i) - mu_k(f0i));
            tc.mu_k = blend_mu(face_mu_k, mu_k(vi));
        };

        // --- Precompute per-face normalized outer scale ---
        // HOP outer contribution per face f is: (area_f / 9) [* / total_w_f]
        // depending on normalize_weights. total_w_f sums active EE mollifiers
        // and either face quadrature weights (when active) or 3 (for the 3
        // face vertices, when face quadrature is not active).
        const bool has_face_quad = params.quad_order > 0;
        const auto& face_quad_rule = params.get_quad_rule();
        double sum_face_qp_w = 0.0;
        for (const auto& qp : face_quad_rule)
            sum_face_qp_w += qp.weight;

        // When face quadrature is active, vertices are already included
        // in the quadrature rule, so don't count the 3 vertex contributions.
        const double base_w = has_face_quad
            ? (face_quad_rule.empty() ? 3.0 : sum_face_qp_w)
            : 3.0;
        Eigen::VectorXd total_w_per_face =
            Eigen::VectorXd::Constant(faces.rows(), base_w);
        if (normalize_weights) {
            // Add per-face sum of active EE mollifiers (EA_EB only).
            for (const auto& [ei_pair, dict_ptr] :
                 collisions.edge_edge_collisions) {
                if (dict_ptr->ee_dtype() != EdgeEdgeDistanceType::EA_EB)
                    continue;
                const auto [e0, e1] = ei_pair;
                const index_t e00 = edges(e0, 0), e01 = edges(e0, 1);
                const index_t e10 = edges(e1, 0), e11 = edges(e1, 1);
                if (e00 == e10 || e00 == e11 || e01 == e10 || e01 == e11)
                    continue;
                const double dist_sqr = edge_edge_distance(
                    vertices.row(e00), vertices.row(e01), vertices.row(e10),
                    vertices.row(e11), EdgeEdgeDistanceType::EA_EB);
                const double dist = std::sqrt(dist_sqr);
                const auto mtypes = edge_edge_mollifier_type(
                    vertices.row(e00).transpose(),
                    vertices.row(e01).transpose(),
                    vertices.row(e10).transpose(),
                    vertices.row(e11).transpose(), dist_sqr);
                double mol =
                    Math<double>::cubic_spline(dist / params.dbar) * 1.5;
                mol *= edge_edge_mollifier<double>(
                    vertices.row(e00).transpose(),
                    vertices.row(e01).transpose(),
                    vertices.row(e10).transpose(),
                    vertices.row(e11).transpose(), mtypes, dist_sqr);
                // Attribute to every face containing e0 (HOP loop iterates
                // over faces and each face's 3 edges).
                for (index_t f = 0; f < faces.rows(); f++) {
                    for (int le = 0; le < 3; le++) {
                        if (mesh.faces_to_edges()(f, le) == e0) {
                            total_w_per_face(f) += mol;
                            break;
                        }
                    }
                }
            }
        }

        // Per-face normalized scale: (area_f / 9) * [1/total_w_f if normalize]
        Eigen::VectorXd face_scale(faces.rows());
        for (index_t f = 0; f < faces.rows(); f++) {
            const double w_f = mesh.face_areas()(f) / 9.0;
            face_scale(f) =
                normalize_weights ? (w_f / total_w_per_face(f)) : w_f;
        }
        // Precompute per-vertex HOP outer weight = sum_{f ∋ v} face_scale(f).
        Eigen::VectorXd v_outer_w = Eigen::VectorXd::Zero(n_verts);
        for (index_t f = 0; f < faces.rows(); f++) {
            for (int lv = 0; lv < 3; lv++)
                v_outer_w(faces(f, lv)) += face_scale(f);
        }

        // ---- VERTEX dicts: all vertex IDs are real ----
        // Skip when face quadrature is active (quad_order > 0), which
        // already includes vertices, matching the normal potential's behavior.
        if (!has_face_quad)
            for (const auto& [vi, dict_ptr] : collisions.vertex_collisions) {
                VertexMatrixView<3> V_view(vertices);
                const double v_w = v_outer_w(vi);
                for (int j = 0; j < dict_ptr->size(); j++) {
                    const auto& cc = (*dict_ptr)[j];
                    const double contact_force =
                        compute_contact_force(cc, V_view, v_w);
                    if (contact_force == 0)
                        continue;

                    switch (cc.type()) {
                    case EspCollisionType::VERTEX_VERTEX: {
                        const index_t v0 = cc[0];
                        const index_t v1 = cc[1];
                        Vector6d cp;
                        cp.head<3>() = vertices.row(v0);
                        cp.tail<3>() = vertices.row(v1);
                        FC_vv.emplace_back(
                            VertexVertexNormalCollision(
                                v0, v1, cc.weight,
                                Eigen::SparseVector<double>()),
                            cp, contact_force);
                        FC_vv.back().weight = cc.weight;
                        const auto& [v0i, v1i, _, __] =
                            FC_vv.back().vertex_ids(edges, faces);
                        FC_vv.back().mu_s = blend_mu(mu_s(v0i), mu_s(v1i));
                        FC_vv.back().mu_k = blend_mu(mu_k(v0i), mu_k(v1i));
                        break;
                    }
                    case EspCollisionType::EDGE_VERTEX: {
                        const index_t edge_id = cc[0];
                        const index_t vert_id = cc[1];
                        const index_t ea0 = edges(edge_id, 0);
                        const index_t ea1 = edges(edge_id, 1);
                        Vector9d cp;
                        cp.segment<3>(0) = vertices.row(vert_id);
                        cp.segment<3>(3) = vertices.row(ea0);
                        cp.segment<3>(6) = vertices.row(ea1);
                        FC_ev.emplace_back(
                            EdgeVertexNormalCollision(
                                edge_id, vert_id, cc.weight,
                                Eigen::SparseVector<double>()),
                            cp, contact_force);
                        FC_ev.back().weight = cc.weight;
                        assign_ev_mu(FC_ev.back());
                        break;
                    }
                    case EspCollisionType::FACE_VERTEX: {
                        const index_t face_id = cc[0];
                        const index_t vert_id = cc[1];
                        const index_t f0 = faces(face_id, 0);
                        const index_t f1 = faces(face_id, 1);
                        const index_t f2 = faces(face_id, 2);
                        Vector12d cp;
                        cp.segment<3>(0) = vertices.row(vert_id);
                        cp.segment<3>(3) = vertices.row(f0);
                        cp.segment<3>(6) = vertices.row(f1);
                        cp.segment<3>(9) = vertices.row(f2);
                        FC_fv.emplace_back(
                            FaceVertexNormalCollision(
                                face_id, vert_id, cc.weight,
                                Eigen::SparseVector<double>()),
                            cp, contact_force);
                        FC_fv.back().weight = cc.weight;
                        assign_fv_mu(FC_fv.back());
                        break;
                    }
                    default:
                        break;
                    }
                }
            }

        // Precompute per-edge HOP outer weight = sum_{f ∋ e} face_scale(f).
        Eigen::VectorXd e_outer_w = Eigen::VectorXd::Zero(edges.rows());
        for (index_t f = 0; f < faces.rows(); f++) {
            for (int le = 0; le < 3; le++)
                e_outer_w(mesh.faces_to_edges()(f, le)) += face_scale(f);
        }

        // ---- EDGE dicts: virtual vertex at edge-edge closest point ----
        for (const auto& [ei_pair, dict_ptr] :
             collisions.edge_edge_collisions) {
            const auto [e0, e1] = ei_pair;
            const auto dtype = dict_ptr->ee_dtype();
            const index_t e00 = edges(e0, 0), e01 = edges(e0, 1);
            const index_t e10 = edges(e1, 0), e11 = edges(e1, 1);

            // The HO potential only contributes for EA_EB; skip otherwise so
            // friction matches exactly.
            if (dtype != EdgeEdgeDistanceType::EA_EB)
                continue;

            // Compute virtual vertex position on edge e0
            // (same logic as quadrature_potential.cpp)
            double closest_uv = line_line_closest_point_pairs_uv<double>(
                vertices.row(e00).transpose(), vertices.row(e01).transpose(),
                vertices.row(e10).transpose(),
                vertices.row(e11).transpose())(0);
            if (!std::isfinite(closest_uv))
                continue;

            const Eigen::RowVector3d virtual_pos =
                closest_uv * (vertices.row(e01) - vertices.row(e00))
                + vertices.row(e00);
            VertexMatrixView<3> V_ext(vertices, virtual_pos);

            // HO potential outer factor for an edge-edge dict (for each face
            // f containing edge e0): area_f/9 * mollifier. Since the mollifier
            // depends only on the four edge endpoints (not f), it factors out
            // and the per-dict outer weight is mollifier * sum_{f∋e0} area_f/9.
            const double dist_sqr_ee = edge_edge_distance(
                vertices.row(e00), vertices.row(e01), vertices.row(e10),
                vertices.row(e11), dtype);
            const double dist_ee = std::sqrt(dist_sqr_ee);
            const auto mtypes = edge_edge_mollifier_type(
                vertices.row(e00).transpose(), vertices.row(e01).transpose(),
                vertices.row(e10).transpose(), vertices.row(e11).transpose(),
                dist_sqr_ee);
            double mollifier =
                Math<double>::cubic_spline(dist_ee / params.dbar) * 1.5;
            mollifier *= edge_edge_mollifier<double>(
                vertices.row(e00).transpose(), vertices.row(e01).transpose(),
                vertices.row(e10).transpose(), vertices.row(e11).transpose(),
                mtypes, dist_sqr_ee);
            const double edge_outer_w = mollifier * e_outer_w(e0);
            if (edge_outer_w == 0)
                continue;
            for (int j = 0; j < dict_ptr->size(); j++) {
                const auto& cc = (*dict_ptr)[j];
                const double contact_force =
                    compute_contact_force(cc, V_ext, edge_outer_w);
                if (contact_force == 0)
                    continue;

                switch (cc.type()) {
                case EspCollisionType::VERTEX_VERTEX: {
                    // One vertex is virtual (n_verts), the other is real.
                    // Elevate to EdgeVertex: edge e0 contains the virtual
                    // vertex, paired with the real vertex.
                    const index_t v0 = cc[0];
                    const index_t v1 = cc[1];
                    const index_t v_real = (v0 == n_verts) ? v1 : v0;

                    Vector9d cp;
                    // Order: [vertex, edge_v0, edge_v1]
                    cp.segment<3>(0) = vertices.row(v_real);
                    cp.segment<3>(3) = vertices.row(e00);
                    cp.segment<3>(6) = vertices.row(e01);

                    FC_ev.emplace_back(
                        EdgeVertexNormalCollision(
                            e0, v_real, cc.weight,
                            Eigen::SparseVector<double>()),
                        cp, contact_force);
                    FC_ev.back().weight = cc.weight;
                    assign_ev_mu(FC_ev.back());
                    break;
                }
                case EspCollisionType::EDGE_VERTEX: {
                    // Real edge (cc[0]) vs virtual vertex (cc[1]) on e0.
                    // Elevate to EdgeEdge: edge e0 vs the real edge.
                    const index_t other_e = cc[0];
                    const index_t oe0 = edges(other_e, 0);
                    const index_t oe1 = edges(other_e, 1);

                    const Eigen::Vector3d e00_pos = vertices.row(e00);
                    const Eigen::Vector3d e01_pos = vertices.row(e01);
                    const Eigen::Vector3d oe0_pos = vertices.row(oe0);
                    const Eigen::Vector3d oe1_pos = vertices.row(oe1);

                    // Mollifier check: skip near-parallel edges
                    if (edge_edge_cross_squarednorm(
                            e00_pos, e01_pos, oe0_pos, oe1_pos)
                        < edge_edge_mollifier_threshold(
                            e00_pos, e01_pos, oe0_pos, oe1_pos)) {
                        continue;
                    }

                    Vector12d cp;
                    // Order: [ea0, ea1, eb0, eb1]
                    cp.segment<3>(0) = e00_pos;
                    cp.segment<3>(3) = e01_pos;
                    cp.segment<3>(6) = oe0_pos;
                    cp.segment<3>(9) = oe1_pos;

                    FC_ee.emplace_back(
                        EdgeEdgeNormalCollision(
                            e0, other_e, 0., EdgeEdgeDistanceType::EA_EB),
                        cp, contact_force);
                    FC_ee.back().weight = cc.weight;
                    assign_ee_mu(FC_ee.back());
                    break;
                }
                case EspCollisionType::FACE_VERTEX: {
                    // Real face (cc[0]) vs virtual vertex on edge e0.
                    // No EdgeFace tangential type exists; resolve via the
                    // point-triangle distance type and elevate to EV (when
                    // closest to a face vertex) or EE (when closest to a
                    // face edge). Interior cases are skipped.
                    const index_t fi = cc[0];
                    const index_t fa = faces(fi, 0);
                    const index_t fb = faces(fi, 1);
                    const index_t fc = faces(fi, 2);
                    const Eigen::Vector3d fa_pos = vertices.row(fa);
                    const Eigen::Vector3d fb_pos = vertices.row(fb);
                    const Eigen::Vector3d fc_pos = vertices.row(fc);
                    const Eigen::Vector3d vp = virtual_pos.transpose();

                    const auto dt = point_triangle_distance_type(
                        vp, fa_pos, fb_pos, fc_pos);

                    auto elevate_to_ev = [&](index_t v_face) {
                        Vector9d cp;
                        cp.segment<3>(0) = vertices.row(v_face);
                        cp.segment<3>(3) = vertices.row(e00);
                        cp.segment<3>(6) = vertices.row(e01);
                        FC_ev.emplace_back(
                            EdgeVertexNormalCollision(
                                e0, v_face, cc.weight,
                                Eigen::SparseVector<double>()),
                            cp, contact_force);
                        FC_ev.back().weight = cc.weight;
                        assign_ev_mu(FC_ev.back());
                    };
                    auto elevate_to_ee = [&](int face_edge_local) {
                        const index_t face_edge_id =
                            mesh.faces_to_edges()(fi, face_edge_local);
                        const index_t fe0 = edges(face_edge_id, 0);
                        const index_t fe1 = edges(face_edge_id, 1);
                        const Eigen::Vector3d e00_pos = vertices.row(e00);
                        const Eigen::Vector3d e01_pos = vertices.row(e01);
                        const Eigen::Vector3d fe0_pos = vertices.row(fe0);
                        const Eigen::Vector3d fe1_pos = vertices.row(fe1);
                        // Mollifier check: skip near-parallel edges
                        if (edge_edge_cross_squarednorm(
                                e00_pos, e01_pos, fe0_pos, fe1_pos)
                            < edge_edge_mollifier_threshold(
                                e00_pos, e01_pos, fe0_pos, fe1_pos)) {
                            return;
                        }
                        Vector12d cp;
                        cp.segment<3>(0) = e00_pos;
                        cp.segment<3>(3) = e01_pos;
                        cp.segment<3>(6) = fe0_pos;
                        cp.segment<3>(9) = fe1_pos;
                        FC_ee.emplace_back(
                            EdgeEdgeNormalCollision(
                                e0, face_edge_id, 0.,
                                EdgeEdgeDistanceType::EA_EB),
                            cp, contact_force);
                        FC_ee.back().weight = cc.weight;
                        assign_ee_mu(FC_ee.back());
                    };

                    switch (dt) {
                    case PointTriangleDistanceType::P_T0:
                        elevate_to_ev(fa);
                        break;
                    case PointTriangleDistanceType::P_T1:
                        elevate_to_ev(fb);
                        break;
                    case PointTriangleDistanceType::P_T2:
                        elevate_to_ev(fc);
                        break;
                    case PointTriangleDistanceType::P_E0:
                        elevate_to_ee(0);
                        break;
                    case PointTriangleDistanceType::P_E1:
                        elevate_to_ee(1);
                        break;
                    case PointTriangleDistanceType::P_E2:
                        elevate_to_ee(2);
                        break;
                    default:
                        // P_T (interior): no elevation possible without an
                        // EdgeFace tangential type. Skip.
                        break;
                    }
                    break;
                }
                default:
                    break;
                }
            }
        }

        // ---- FACE dicts: virtual vertex at face quadrature point ----
        for (const auto& [fi, dicts] : collisions.face_collisions) {
            const index_t f0 = faces(fi, 0);
            const index_t f1 = faces(fi, 1);
            const index_t f2 = faces(fi, 2);
            // HO potential outer face weight (optionally normalized).
            const double face_w = face_scale(fi);

            for (size_t qi = 0; qi < dicts.size(); qi++) {
                const auto& dict_ptr = dicts[qi];
                const auto& qp = face_quad_rule[qi];

                // Compute virtual vertex position from barycentric coords
                const Eigen::RowVector3d virtual_pos =
                    qp.lambda[0] * vertices.row(f0)
                    + qp.lambda[1] * vertices.row(f1)
                    + qp.lambda[2] * vertices.row(f2);
                VertexMatrixView<3> V_ext(vertices, virtual_pos);

                // HO potential per-quadrature factor: area/9 * qp.weight.
                const double fq_outer_w = face_w * qp.weight;

                for (int j = 0; j < dict_ptr->size(); j++) {
                    const auto& cc = (*dict_ptr)[j];
                    const double contact_force =
                        compute_contact_force(cc, V_ext, fq_outer_w);
                    if (contact_force == 0)
                        continue;

                    switch (cc.type()) {
                    case EspCollisionType::VERTEX_VERTEX: {
                        // One vertex is virtual (on face fi), other is real.
                        // Elevate to FaceVertex: face fi paired with the
                        // real vertex.
                        const index_t v0 = cc[0];
                        const index_t v1 = cc[1];
                        const index_t v_real = (v0 == n_verts) ? v1 : v0;

                        Vector12d cp;
                        // Order: [vertex, face_v0, face_v1, face_v2]
                        cp.segment<3>(0) = vertices.row(v_real);
                        cp.segment<3>(3) = vertices.row(f0);
                        cp.segment<3>(6) = vertices.row(f1);
                        cp.segment<3>(9) = vertices.row(f2);

                        FC_fv.emplace_back(
                            FaceVertexNormalCollision(
                                fi, v_real, cc.weight,
                                Eigen::SparseVector<double>()),
                            cp, contact_force);
                        FC_fv.back().weight = cc.weight;
                        assign_fv_mu(FC_fv.back());
                        break;
                    }
                    case EspCollisionType::EDGE_VERTEX: {
                        // Real edge (cc[0]) vs virtual vertex on face fi.
                        // Distribute the sub-collision's contact force to
                        // the two edge endpoints using the projection
                        // parameter u (clamped to [0, 1]). Each endpoint
                        // gets an FV(tangential) entry with face = fi and
                        // vertex = endpoint, weighted accordingly. This
                        // handles interior (P_E) as well as P_E0 / P_E1.
                        const index_t other_e = cc[0];
                        const index_t oe0 = edges(other_e, 0);
                        const index_t oe1 = edges(other_e, 1);
                        const Eigen::Vector3d oe0_pos = vertices.row(oe0);
                        const Eigen::Vector3d oe1_pos = vertices.row(oe1);
                        const Eigen::Vector3d vp = virtual_pos.transpose();

                        double u =
                            point_edge_closest_point(vp, oe0_pos, oe1_pos);
                        if (!std::isfinite(u))
                            break;
                        u = std::clamp(u, 0.0, 1.0);
                        const double w0 = 1.0 - u;
                        const double w1 = u;

                        auto emit_fv = [&](index_t v_edge, double w) {
                            if (w <= 0)
                                return;
                            Vector12d cp;
                            cp.segment<3>(0) = vertices.row(v_edge);
                            cp.segment<3>(3) = vertices.row(f0);
                            cp.segment<3>(6) = vertices.row(f1);
                            cp.segment<3>(9) = vertices.row(f2);
                            FC_fv.emplace_back(
                                FaceVertexNormalCollision(
                                    fi, v_edge, cc.weight,
                                    Eigen::SparseVector<double>()),
                                cp, w * contact_force);
                            FC_fv.back().weight = cc.weight;
                            assign_fv_mu(FC_fv.back());
                        };
                        emit_fv(oe0, w0);
                        emit_fv(oe1, w1);
                        break;
                    }
                    case EspCollisionType::FACE_VERTEX: {
                        // Real face (cc[0]) vs virtual vertex on face fi.
                        // Distribute the sub-collision's contact force to
                        // the three cube-face vertices using barycentric
                        // coordinates of the projection of the virtual
                        // point onto the cube face. Each vertex gets an
                        // FV(tangential) entry with face = fi (plane face)
                        // and vertex = cube vertex, weighted by its
                        // barycentric coefficient. This handles interior
                        // (P_T) and edge cases as well as vertex cases —
                        // needed for flat-on-flat sliding.
                        const index_t fj = cc[0];
                        const index_t ja = faces(fj, 0);
                        const index_t jb = faces(fj, 1);
                        const index_t jc = faces(fj, 2);
                        const Eigen::Vector3d ja_pos = vertices.row(ja);
                        const Eigen::Vector3d jb_pos = vertices.row(jb);
                        const Eigen::Vector3d jc_pos = vertices.row(jc);
                        const Eigen::Vector3d vp = virtual_pos.transpose();

                        Eigen::Vector2d bary = point_triangle_closest_point(
                            vp, ja_pos, jb_pos, jc_pos);
                        if (!bary.allFinite())
                            break;
                        double beta = std::clamp(bary(0), 0.0, 1.0);
                        double gamma = std::clamp(bary(1), 0.0, 1.0);
                        if (beta + gamma > 1.0) {
                            const double s = beta + gamma;
                            beta /= s;
                            gamma /= s;
                        }
                        const double alpha = 1.0 - beta - gamma;

                        auto emit_fv = [&](index_t v_other, double w) {
                            if (w <= 0)
                                return;
                            Vector12d cp;
                            cp.segment<3>(0) = vertices.row(v_other);
                            cp.segment<3>(3) = vertices.row(f0);
                            cp.segment<3>(6) = vertices.row(f1);
                            cp.segment<3>(9) = vertices.row(f2);
                            FC_fv.emplace_back(
                                FaceVertexNormalCollision(
                                    fi, v_other, cc.weight,
                                    Eigen::SparseVector<double>()),
                                cp, w * contact_force);
                            FC_fv.back().weight = cc.weight;
                            assign_fv_mu(FC_fv.back());
                        };
                        emit_fv(ja, alpha);
                        emit_fv(jb, beta);
                        emit_fv(jc, gamma);
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
        }
    }
}

// ============================================================================

size_t TangentialCollisions::size() const
{
    return vv_collisions.size() + ev_collisions.size() + ee_collisions.size()
        + fv_collisions.size() + pv_collisions.size();
}

bool TangentialCollisions::empty() const
{
    return vv_collisions.empty() && ev_collisions.empty()
        && ee_collisions.empty() && fv_collisions.empty()
        && pv_collisions.empty();
}

void TangentialCollisions::clear()
{
    vv_collisions.clear();
    ev_collisions.clear();
    ee_collisions.clear();
    fv_collisions.clear();
    pv_collisions.clear();
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
    i -= fv_collisions.size();
    if (i < pv_collisions.size()) {
        return pv_collisions[i];
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
    i -= fv_collisions.size();
    if (i < pv_collisions.size()) {
        return pv_collisions[i];
    }
    throw std::out_of_range("Friction collision index is out of range!");
}

} // namespace ipc
