#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/esp/collisions/esp_quadrature.hpp>
#include <ipc/esp/collisions/vertex_matrix_view.hpp>
#include <ipc/esp/esp_collision_maps.hpp>
#include <ipc/esp/esp_distance.hpp>
#include <ipc/esp/esp_parameters.hpp>
#include <ipc/math/math.hpp>
#include <ipc/tangent/closest_point.hpp>

namespace ipc {

void TangentialCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ESPCollisions& collisions,
    const ESPParameters& params,
    const double normal_stiffness,
    Eigen::ConstRef<Eigen::VectorXd> mu_s,
    Eigen::ConstRef<Eigen::VectorXd> mu_k,
    const std::function<double(double, double)>& blend_mu,
    const bool normalize_weights)
{
    assert(mu_s.size() == vertices.rows());
    assert(mu_k.size() == vertices.rows());
    if (params.local_quadrature_transfer && mesh.dim() == 3
        && params.dbar_factor() > 0) {
        // The legacy friction weights do not represent transferred quadrature.
        if (!mu_s.isZero() || !mu_k.isZero()) {
            throw std::invalid_argument(
                "ESP local quadrature transfer currently supports frictionless contact only.");
        }
        clear();
        return;
    }
    // The quadrature loops below index params' rule by the collision sets'
    // per-point entries, so the two must come from the same rule.
    collisions.check_quadrature(params);

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();
    const int dim = mesh.dim();

    clear();

    auto& FC_vv = vv_collisions;
    auto& FC_ev = ev_collisions;
    auto& FC_ee = ee_collisions;
    auto& FC_fv = fv_collisions;

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

        auto compute_contact_force_2d = [&](const ESPCollision& cc,
                                            const VertexMatrixView<2>& V_ext,
                                            const double outer_w) -> double {
            const Eigen::VectorXd positions = cc.dof(V_ext);
            double d2 = 0;
            switch (cc.type()) {
            case ESPCollisionType::VERTEX_VERTEX:
                d2 = point_point_distance(
                    Eigen::Vector2d(positions.segment<2>(0)),
                    Eigen::Vector2d(positions.segment<2>(2)));
                break;
            case ESPCollisionType::EDGE_VERTEX:
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

        for (const auto& edge_collisions :
             collisions.maps().edge_collisions_2d) {
            const index_t ei = edge_collisions.first;
            const auto& qp_dicts = edge_collisions.second;
            const index_t e0 = edges(ei, 0);
            const index_t e1 = edges(ei, 1);
            const double L = mesh.edge_length(ei);

            for (size_t qi = 0; qi < qp_dicts.size(); ++qi) {
                const auto& dict_ptr = qp_dicts[qi];
                if (!dict_ptr || dict_ptr->size() == 0) {
                    continue;
                }
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
                    if (contact_force == 0) {
                        continue;
                    }

                    switch (cc.type()) {
                    case ESPCollisionType::VERTEX_VERTEX: {
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
                    case ESPCollisionType::EDGE_VERTEX: {
                        // Virtual vertex on edge ei vs real edge ej.
                        // Distribute to the two endpoints of ej weighted by
                        // the projection parameter u (parallel to 3D face
                        // dict EV elevation).
                        const index_t ea = cc.vertex_id(1);
                        const index_t eb = cc.vertex_id(2);
                        const Eigen::Vector2d ea_pos =
                            vertices.row(ea).transpose();
                        const Eigen::Vector2d eb_pos =
                            vertices.row(eb).transpose();
                        const Eigen::Vector2d vp = virtual_pos.transpose();

                        double u = point_edge_closest_point(vp, ea_pos, eb_pos);
                        if (!std::isfinite(u)) {
                            break;
                        }
                        u = std::clamp(u, 0.0, 1.0);

                        auto emit_ev = [&](index_t v_edge, double w) {
                            if (w <= 0) {
                                return;
                            }
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
        // Uses the scalar derivative of params.barrier w.r.t. distance,
        // scaled by outer quadrature weight and barrier stiffness.
        auto compute_contact_force = [&](const ESPCollision& cc,
                                         const VertexMatrixView<3>& V_ext,
                                         const double outer_w) -> double {
            const Eigen::VectorXd positions = cc.dof(V_ext);
            double d2 = 0;
            switch (cc.type()) {
            case ESPCollisionType::VERTEX_VERTEX:
                d2 = point_point_distance(
                    Eigen::Vector3d(positions.segment<3>(0)),
                    Eigen::Vector3d(positions.segment<3>(3)));
                break;
            case ESPCollisionType::EDGE_VERTEX:
                d2 = point_edge_distance(
                    Eigen::Vector3d(positions.segment<3>(6)),
                    Eigen::Vector3d(positions.segment<3>(0)),
                    Eigen::Vector3d(positions.segment<3>(3)));
                break;
            case ESPCollisionType::FACE_VERTEX:
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
        // ESP outer contribution per face f is: area_f [* / total_w_f]
        // depending on normalize_weights. total_w_f sums active EE mollifiers
        // and the face rule's weights, which sum to 1. The rule's corner
        // points are the vertex samples (weighted by the face's corner
        // weights); its other points are the per-face samples. Under
        // ESPParameters::set_auto_lattice() each face has its own rule.
        Eigen::VectorXd total_w_per_face =
            Eigen::VectorXd::Constant(faces.rows(), 1.0);
        if (normalize_weights) {
            // Add per-face sum of active EE mollifiers (EA_EB only).
            for (const auto& [ei_pair, dict_ptr] :
                 collisions.maps().edge_edge_collisions) {
                // The builder keeps only EA_EB pairs of edges that share no
                // vertex.
                assert(dict_ptr->ee_dtype() == EdgeEdgeDistanceType::EA_EB);
                const index_t e0 = ei_pair.first;
                const index_t e1 = ei_pair.second;
                const index_t e00 = edges(e0, 0), e01 = edges(e0, 1);
                const index_t e10 = edges(e1, 0), e11 = edges(e1, 1);
                assert(e00 != e10 && e00 != e11 && e01 != e10 && e01 != e11);
                const double dist_sqr = edge_edge_distance_parallel_safe(
                    vertices.row(e00), vertices.row(e01), vertices.row(e10),
                    vertices.row(e11), EdgeEdgeDistanceType::EA_EB);
                const double dist = std::sqrt(dist_sqr);
                const auto mtypes = edge_edge_mollifier_type(
                    vertices.row(e00).transpose(),
                    vertices.row(e01).transpose(),
                    vertices.row(e10).transpose(),
                    vertices.row(e11).transpose(), dist_sqr);
                double mol =
                    Math<double>::cubic_spline(dist / params.ee_support())
                    * ee_mollifier_config().ee_weight /* PROTOTYPE */;
                mol *= edge_edge_mollifier<double>(
                    vertices.row(e00).transpose(),
                    vertices.row(e01).transpose(),
                    vertices.row(e10).transpose(),
                    vertices.row(e11).transpose(), mtypes, dist_sqr);
                // Attribute to every face containing e0 (ESP loop iterates
                // over faces and each face's 3 edges).
                for (const index_t f : mesh.edges_to_faces()[e0]) {
                    total_w_per_face(f) += mol;
                }
            }
        }

        // Per-face normalized scale: area_f * [1/total_w_f if normalize]
        Eigen::VectorXd face_scale(faces.rows());
        for (index_t f = 0; f < faces.rows(); f++) {
            const double w_f = mesh.face_areas()(f);
            face_scale(f) =
                normalize_weights ? (w_f / total_w_per_face(f)) : w_f;
        }
        // Per-vertex ESP outer weight:
        // sum_{f ∋ v} face_scale(f) * corner weight of v in f.
        Eigen::VectorXd v_outer_w = Eigen::VectorXd::Zero(n_verts);
        for (index_t f = 0; f < faces.rows(); f++) {
            for (int lv = 0; lv < 3; lv++) {
                v_outer_w(faces(f, lv)) += face_scale(f)
                    * params.face_rule(mesh.face_areas()(f)).corner_weights[lv];
            }
        }

        // ---- VERTEX dicts: all vertex IDs are real ----
        {
            for (const auto& [vi, dict_ptr] :
                 collisions.maps().vertex_collisions) {
                VertexMatrixView<3> V_view(vertices);
                const double v_w = v_outer_w(vi);
                for (int j = 0; j < dict_ptr->size(); j++) {
                    const auto& cc = (*dict_ptr)[j];
                    const double contact_force =
                        compute_contact_force(cc, V_view, v_w);
                    if (contact_force == 0) {
                        continue;
                    }

                    switch (cc.type()) {
                    case ESPCollisionType::VERTEX_VERTEX: {
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
                    case ESPCollisionType::EDGE_VERTEX: {
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
                    case ESPCollisionType::FACE_VERTEX: {
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
        }

        // Precompute per-edge ESP outer weight = sum_{f ∋ e} face_scale(f).
        Eigen::VectorXd e_outer_w = Eigen::VectorXd::Zero(edges.rows());
        for (index_t f = 0; f < faces.rows(); f++) {
            for (int le = 0; le < 3; le++) {
                e_outer_w(mesh.faces_to_edges()(f, le)) += face_scale(f);
            }
        }

        // ---- EDGE dicts: virtual vertex at edge-edge closest point ----
        for (const auto& [ei_pair, dict_ptr] :
             collisions.maps().edge_edge_collisions) {
            const index_t e0 = ei_pair.first;
            const index_t e1 = ei_pair.second;
            const auto dtype = dict_ptr->ee_dtype();
            const index_t e00 = edges(e0, 0), e01 = edges(e0, 1);
            const index_t e10 = edges(e1, 0), e11 = edges(e1, 1);

            // The builder keeps only EA_EB pairs.
            assert(dtype == EdgeEdgeDistanceType::EA_EB);

            // Compute virtual vertex position on edge e0
            // (same logic as quadrature_potential.cpp)
            const double closest_uv = closest_point_uv<double>(
                vertices.row(e00).transpose(), vertices.row(e01).transpose(),
                vertices.row(e10).transpose(), vertices.row(e11).transpose());
            if (!std::isfinite(closest_uv)) {
                continue;
            }

            const Eigen::RowVector3d virtual_pos =
                closest_uv * (vertices.row(e01) - vertices.row(e00))
                + vertices.row(e00);
            VertexMatrixView<3> V_ext(vertices, virtual_pos);

            // ESP potential outer factor for an edge-edge dict (for each face
            // f containing edge e0): area_f * mollifier. Since the mollifier
            // depends only on the four edge endpoints (not f), it factors out
            // and the per-dict outer weight is mollifier * sum_{f∋e0} area_f.
            const double dist_sqr_ee = edge_edge_distance_parallel_safe(
                vertices.row(e00), vertices.row(e01), vertices.row(e10),
                vertices.row(e11), dtype);
            const double dist_ee = std::sqrt(dist_sqr_ee);
            const auto mtypes = edge_edge_mollifier_type(
                vertices.row(e00).transpose(), vertices.row(e01).transpose(),
                vertices.row(e10).transpose(), vertices.row(e11).transpose(),
                dist_sqr_ee);
            double mollifier =
                Math<double>::cubic_spline(dist_ee / params.ee_support()) * ee_mollifier_config().ee_weight /* PROTOTYPE */;
            mollifier *= edge_edge_mollifier<double>(
                vertices.row(e00).transpose(), vertices.row(e01).transpose(),
                vertices.row(e10).transpose(), vertices.row(e11).transpose(),
                mtypes, dist_sqr_ee);
            const double edge_outer_w = mollifier * e_outer_w(e0);
            if (edge_outer_w == 0) {
                continue;
            }
            for (int j = 0; j < dict_ptr->size(); j++) {
                const auto& cc = (*dict_ptr)[j];
                const double contact_force =
                    compute_contact_force(cc, V_ext, edge_outer_w);
                if (contact_force == 0) {
                    continue;
                }

                switch (cc.type()) {
                case ESPCollisionType::VERTEX_VERTEX: {
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
                case ESPCollisionType::EDGE_VERTEX: {
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
                case ESPCollisionType::FACE_VERTEX: {
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
        for (const auto& face_collisions : collisions.maps().face_collisions) {
            const index_t fi = face_collisions.first;
            const auto& dicts = face_collisions.second;
            const index_t f0 = faces(fi, 0);
            const index_t f1 = faces(fi, 1);
            const index_t f2 = faces(fi, 2);
            // ESP potential outer face weight (optionally normalized).
            const double face_w = face_scale(fi);
            const FaceQuadRule& face_quad_rule =
                *params.face_rule(mesh.face_areas()(fi)).non_corner;

            for (size_t qi = 0; qi < dicts.size(); qi++) {
                const auto& dict_ptr = dicts[qi];
                const auto& qp = face_quad_rule[qi];

                // Compute virtual vertex position from barycentric coords
                const Eigen::RowVector3d virtual_pos =
                    qp.lambda[0] * vertices.row(f0)
                    + qp.lambda[1] * vertices.row(f1)
                    + qp.lambda[2] * vertices.row(f2);
                VertexMatrixView<3> V_ext(vertices, virtual_pos);

                // ESP potential per-quadrature factor: area * qp.weight.
                const double fq_outer_w = face_w * qp.weight;

                for (int j = 0; j < dict_ptr->size(); j++) {
                    const auto& cc = (*dict_ptr)[j];
                    const double contact_force =
                        compute_contact_force(cc, V_ext, fq_outer_w);
                    if (contact_force == 0) {
                        continue;
                    }

                    switch (cc.type()) {
                    case ESPCollisionType::VERTEX_VERTEX: {
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
                    case ESPCollisionType::EDGE_VERTEX: {
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
                        if (!std::isfinite(u)) {
                            break;
                        }
                        u = std::clamp(u, 0.0, 1.0);
                        const double w0 = 1.0 - u;
                        const double w1 = u;

                        auto emit_fv = [&](index_t v_edge, double w) {
                            if (w <= 0) {
                                return;
                            }
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
                    case ESPCollisionType::FACE_VERTEX: {
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
                        if (!bary.allFinite()) {
                            break;
                        }
                        double beta = std::clamp(bary(0), 0.0, 1.0);
                        double gamma = std::clamp(bary(1), 0.0, 1.0);
                        if (beta + gamma > 1.0) {
                            const double s = beta + gamma;
                            beta /= s;
                            gamma /= s;
                        }
                        const double alpha = 1.0 - beta - gamma;

                        auto emit_fv = [&](index_t v_other, double w) {
                            if (w <= 0) {
                                return;
                            }
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

    reset_lagged_anisotropic_friction_coefficients();
}

} // namespace ipc
