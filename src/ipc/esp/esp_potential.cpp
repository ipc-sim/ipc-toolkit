#include "esp_potential.hpp"

#include "ipc/barrier/barrier.hpp"
#include "ipc/distance/edge_edge.hpp"
#include "ipc/distance/edge_edge_mollifier.hpp"
#include "ipc/esp/collisions/esp_quadrature.hpp"
#include "ipc/esp/collisions/vertex_matrix_view.hpp"
#include "ipc/esp/quadrature_potential.hpp"
#include "ipc/gcp/distance/mollifier.hpp"
#include "ipc/gcp/distance/point_face.hpp"

#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/profile_registry.hpp>

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <unordered_map>

namespace ipc {

constexpr double face_quadrature_weight_scale = 1.0;

namespace {
    // Adapt mollifier order to the barrier singularity.
    // - Log / default barriers: order 1 (no extra power).
    // - InversePowerBarrier(p): order = round(p) + 1, so p=1 -> 2, p=2 -> 3.
    int mollifier_order_for_barrier(const std::shared_ptr<Barrier>& barrier)
    {
        if (const auto* ip =
                dynamic_cast<const InversePowerBarrier*>(barrier.get())) {
            const int p = static_cast<int>(std::lround(ip->power()));
            return std::max(1, p + 1);
        }
        return 1;
    }

    template <class T> inline T pow_int(T x, int n)
    {
        T r = T(1);
        for (int i = 0; i < n; ++i) {
            r = r * x;
        }
        return r;
    }
} // namespace

double ESPPotential::operator()(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    IPC_PROFILE_SCOPE("ho.potential_eval");
    assert(X.rows() == mesh.num_vertices());

    m_edge_evaluation_count.clear();

    if (collisions.empty()) {
        return 0;
    }

    double result = 0;

    if (mesh.dim() == 2) {
        tbb::enumerable_thread_specific<double> potential_storage(0.0);
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(params.quad_order);

        // Collect active edge ids into a flat vector for parallel indexing.
        std::vector<index_t> active_edges;
        active_edges.reserve(collisions.edge_collisions_2d.size());
        for (const auto& [ei, _] : collisions.edge_collisions_2d) {
            active_edges.push_back(ei);
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, active_edges.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                double& total = potential_storage.local();
                for (size_t k = r.begin(); k < r.end(); ++k) {
                    const index_t ei = active_edges[k];
                    const auto& qp_dicts = collisions.edge_collisions_2d.at(ei);
                    const double L = mesh.edge_area(ei);
                    const double w_edge = params.area_weights ? L : 1.;
                    const index_t e0 = mesh.edges()(ei, 0);
                    const index_t e1 = mesh.edges()(ei, 1);

                    for (size_t qi = 0; qi < qp_dicts.size(); ++qi) {
                        const auto& dict = *qp_dicts[qi];
                        if (dict.size() == 0) {
                            continue;
                        }
                        const auto& qp = rule[qi];
                        const std::array<double, 2> lambda = { { 1.0 - qp.xi,
                                                                 qp.xi } };
                        const Eigen::RowVector2d q_pos =
                            lambda[0] * X.row(e0) + lambda[1] * X.row(e1);
                        VertexMatrixView<2> X_ext(X, q_pos);
                        total += w_edge * qp.weight
                            * PointPotentialHelper::
                                evaluate_potential_at_edge_qp(
                                     X_ext, dict, params,
                                     collisions.adaptive_dhat.get());
                    }
                }
            });

        for (const double v : potential_storage) {
            result += v;
        }
    } else if (mesh.dim() == 3) {
        {
            tbb::enumerable_thread_specific<double> potential_storage(0.0);
            tbb::enumerable_thread_specific<CountMap> count_storage {
                CountMap()
            };
            tbb::enumerable_thread_specific<size_t> fq_point_storage(size_t(0));

            // Disable near/far splitting for edge cases
            const double dbar_factor = params.dbar_factor();
            const bool use_nf =
                use_near_far && dbar_factor > 0 && dbar_factor < 1;
            const bool skip_ee =
                (dbar_factor == 0); // Skip EE pairs when dbar_factor == 0

            std::unique_ptr<NearFarBarrier> nf_barrier;
            if (use_nf) {
                nf_barrier = std::make_unique<NearFarBarrier>(
                    params.barrier, dbar_factor);
            }

            auto loop_body = [&](const tbb::blocked_range<index_t>& r) {
                double& total = potential_storage.local();
                CountMap& local_counts = count_storage.local();
                size_t& local_fq_points = fq_point_storage.local();
                for (index_t f = r.begin(); f < r.end(); f++) {
                    const double area = mesh.face_areas()(f);
                    const double w = params.area_weights ? (area / 9.) : 1.;

                    double total_w = 0, total_w_near = 0, total_w_far = 0;
                    double total_p = 0, total_p_near = 0, total_p_far = 0;

                    for (index_t le = 0; le < 3; le++) {
                        const index_t edge_id = mesh.faces_to_edges()(f, le);
                        const index_t ea = mesh.edges()(edge_id, 0);
                        const index_t eb = mesh.edges()(edge_id, 1);

                        for (index_t other_edge_id :
                             collisions.m_candidates.ee_set(edge_id)) {
                            const index_t ec = mesh.edges()(other_edge_id, 0);
                            const index_t ed = mesh.edges()(other_edge_id, 1);

                            // Skip adjacent edges
                            if (ea == ec || ea == ed || eb == ec || eb == ed) {
                                continue;
                            }

                            // Skip EE pairs entirely when dbar_factor == 0
                            if (skip_ee) {
                                continue;
                            }

                            if (auto iter =
                                    collisions.edge_edge_collisions.find(
                                        std::make_pair(edge_id, other_edge_id));
                                iter != collisions.edge_edge_collisions.end()) {

                                const auto dtype = iter->second->ee_dtype();

                                // Skip non EA_EB collision types
                                if (dtype != EdgeEdgeDistanceType::EA_EB) {
                                    continue;
                                }

                                const double dist = sqrt(edge_edge_distance(
                                    X.row(ea), X.row(eb), X.row(ec), X.row(ed),
                                    dtype));

                                const double uv = closest_point_uv<double>(
                                    X.row(ea), X.row(eb), X.row(ec), X.row(ed),
                                    dtype);

                                const Eigen::RowVector3d ee_closest_point =
                                    uv * (X.row(eb) - X.row(ea)) + X.row(ea);

                                const double dist_sqr = edge_edge_distance(
                                    X.row(ea), X.row(eb), X.row(ec), X.row(ed),
                                    dtype);
                                const auto mtypes = edge_edge_mollifier_type(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), dist_sqr);

                                double mollifier = Math<double>::cubic_spline(
                                                       dist / params.dbar)
                                    * 1.5;
                                mollifier *= edge_edge_mollifier<double>(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), mtypes, dist_sqr);
                                mollifier = pow_int(
                                    mollifier,
                                    mollifier_order_for_barrier(
                                        params.barrier));

                                if (use_nf) {
                                    const double P_near = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                                            VertexMatrixView<3>(
                                                X, ee_closest_point),
                                            *(iter->second), params,
                                            collisions.adaptive_dhat.get(),
                                            dtype, *nf_barrier);
                                    total_w_near += mollifier;
                                    total_p_near += mollifier * P_near;
                                } else {
                                    const double P_val = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions(
                                            VertexMatrixView<3>(
                                                X, ee_closest_point),
                                            *(iter->second), params,
                                            collisions.adaptive_dhat.get(),
                                            dtype);
                                    total_w += mollifier;
                                    total_p += mollifier * P_val;
                                }
                                local_counts[edge_id]++;
                            }
                        }
                    }

                    // Face-interior quadrature points controlled by
                    // params.quad_order.
                    const auto& face_quad_rule = params.get_quad_rule();
                    {
                        auto iter = collisions.face_collisions.find(f);
                        for (size_t qi = 0; qi < face_quad_rule.size(); qi++) {
                            const auto& qp = face_quad_rule[qi];
                            if (use_nf) {
                                total_w_near +=
                                    face_quadrature_weight_scale * qp.weight;
                                total_w_far +=
                                    face_quadrature_weight_scale * qp.weight;
                            } else {
                                total_w +=
                                    face_quadrature_weight_scale * qp.weight;
                            }
                            if (iter != collisions.face_collisions.end()) {
                                local_fq_points++;
                                const Eigen::RowVector3d q_pos =
                                    qp.lambda[0] * X.row(mesh.faces()(f, 0))
                                    + qp.lambda[1] * X.row(mesh.faces()(f, 1))
                                    + qp.lambda[2] * X.row(mesh.faces()(f, 2));
                                if (use_nf) {
                                    auto [fq_near, fq_far] = PointPotentialHelper::
                                        evaluate_potential_at_face_center_with_cached_collisions_nearfar(
                                            VertexMatrixView<3>(X, q_pos),
                                            *iter->second[qi], params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    total_p_near += face_quadrature_weight_scale
                                        * qp.weight * fq_near;
                                    total_p_far += face_quadrature_weight_scale
                                        * qp.weight * fq_far;
                                } else {
                                    const double fq_val = PointPotentialHelper::
                                        evaluate_potential_at_face_center_with_cached_collisions(
                                            VertexMatrixView<3>(X, q_pos),
                                            *iter->second[qi], params,
                                            collisions.adaptive_dhat.get());
                                    total_p += face_quadrature_weight_scale
                                        * qp.weight * fq_val;
                                }
                            }
                        }
                    }

                    // Only integrate on vertices explicitly if there is no
                    // ESP quadrature, since that includes verts
                    if (face_quad_rule.empty()) {
                        for (index_t lv = 0; lv < 3; lv++) {
                            const index_t v = mesh.faces()(f, lv);
                            if (use_nf) {
                                total_w_near += 1.;
                                total_w_far += 1.;
                            } else {
                                total_w += 1.;
                            }
                            if (auto iter =
                                    collisions.vertex_collisions.find(v);
                                iter != collisions.vertex_collisions.end()) {
                                if (use_nf) {
                                    auto [vt_near, vt_far] = PointPotentialHelper::
                                        evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                                            X, *(iter->second), params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    total_p_near += vt_near;
                                    total_p_far += vt_far;
                                } else {
                                    const double vt_val = PointPotentialHelper::
                                        evaluate_potential_at_vertex_with_cached_collisions(
                                            X, *(iter->second), params,
                                            collisions.adaptive_dhat.get());
                                    total_p += vt_val;
                                }
                            }
                        }
                    }

                    if (use_nf) {
                        assert(total_w_near >= total_w_far - 1e-14);
                        if (total_w_far == 0) {
                            assert(total_p_far == 0);
                        }
                        if (total_w_near > 0 && total_w_far > 0) {
                            total += w
                                * (total_p_near / total_w_near
                                   + total_p_far / total_w_far);
                        } else if (total_w_near > 0) {
                            total += w * (total_p_near / total_w_near);
                        }
                    } else if (use_near_far) {
                        assert(total_w > 0);
                        total += w * (total_p / total_w);
                    } else {
                        total += w * total_p;
                    }
                }
            };

            tbb::parallel_for(
                tbb::blocked_range<index_t>(0, mesh.num_faces()), loop_body);

            for (const auto& local_potential : potential_storage) {
                result += local_potential;
            }
            /*
            size_t total_fq_points = 0;
            for (const auto& n : fq_point_storage) {
                total_fq_points += n;
            }
            logger().debug("[ESPPotential] face quadrature points
            evaluated: {}", total_fq_points);
            */

            for (const auto& local_counts : count_storage) {
                for (const auto& [id, count] : local_counts) {
                    m_edge_evaluation_count[id] += count;
                }
            }
        }
    }
    return result;
}

Eigen::VectorXd ESPPotential::gradient(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    IPC_PROFILE_SCOPE("ho.potential_gradient");
    assert(X.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(X.size());
    }

    const int dim = X.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(X.size()));

    if (mesh.dim() == 2) {
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(params.quad_order);

        std::vector<index_t> active_edges;
        active_edges.reserve(collisions.edge_collisions_2d.size());
        for (const auto& [ei, _] : collisions.edge_collisions_2d) {
            active_edges.push_back(ei);
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, active_edges.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                Eigen::VectorXd& global_grad = storage.local();

                for (size_t k = r.begin(); k < r.end(); ++k) {
                    const index_t ei = active_edges[k];
                    const auto& qp_dicts = collisions.edge_collisions_2d.at(ei);
                    const double L = mesh.edge_area(ei);
                    const double w_edge = params.area_weights ? L : 1.;
                    const index_t e0 = mesh.edges()(ei, 0);
                    const index_t e1 = mesh.edges()(ei, 1);

                    for (size_t qi = 0; qi < qp_dicts.size(); ++qi) {
                        const auto& dict = *qp_dicts[qi];
                        if (dict.size() == 0) {
                            continue;
                        }
                        const auto& qp = rule[qi];
                        const std::array<double, 2> lambda = { { 1.0 - qp.xi,
                                                                 qp.xi } };
                        const Eigen::RowVector2d q_pos =
                            lambda[0] * X.row(e0) + lambda[1] * X.row(e1);
                        VertexMatrixView<2> X_ext(X, q_pos);

                        const Eigen::VectorXd local_grad =
                            w_edge * qp.weight
                            * PointPotentialHelper::
                                evaluate_potential_gradient_at_edge_qp(
                                    X_ext, dict, params,
                                    collisions.adaptive_dhat.get(), lambda);

                        local_gradient_to_global_gradient(
                            local_grad, dict.vertex_ids(), dim, global_grad);
                    }
                }
            });
    } else if (mesh.dim() == 3) {
        {
            using T = ADGrad<12>;

            const double dbar_factor = params.dbar_factor();
            const bool use_nf_grad =
                use_near_far && dbar_factor > 0 && dbar_factor < 1;
            const bool skip_ee_grad = (dbar_factor == 0);

            auto loop_body = [&](const tbb::blocked_range<index_t>& r) {
                Eigen::VectorXd& grad = storage.local();
                for (index_t f = r.begin(); f < r.end(); f++) {
                    const double area = mesh.face_areas()(f);
                    const double w = params.area_weights ? (area / 9.) : 1.;
                    // Pass 1: collect all quadrature contributions for this
                    // face
                    struct EEGradEntry {
                        const ESPCollisionDict<PointType::EDGE>* dict;
                        double mol_val;
                        Eigen::Vector<double, 12> mol_grad;
                        double P;
                        Eigen::VectorXd grad_P;
                    };
                    struct ConstGradEntry {
                        const std::vector<index_t>* dofs;
                        Eigen::VectorXd grad_P_near,
                            grad_P_far;       // near/far or single (for non-nf)
                        double P_near, P_far; // near/far or single (for non-nf)
                    };
                    std::vector<EEGradEntry> ee_cache;
                    std::vector<ConstGradEntry> const_cache;
                    double total_w = 0;
                    double total_p = 0;
                    double total_w_near = 0, total_p_near = 0;
                    double total_w_far = 0, total_p_far = 0;

                    std::unique_ptr<NearFarBarrier> nf_barrier;
                    if (use_nf_grad) {
                        nf_barrier = std::make_unique<NearFarBarrier>(
                            params.barrier, params.dbar_factor());
                    }

                    for (index_t le = 0; le < 3; le++) {
                        const index_t edge_id = mesh.faces_to_edges()(f, le);
                        const index_t ea = mesh.edges()(edge_id, 0);
                        const index_t eb = mesh.edges()(edge_id, 1);

                        for (index_t other_edge_id :
                             collisions.m_candidates.ee_set(edge_id)) {
                            const index_t ec = mesh.edges()(other_edge_id, 0);
                            const index_t ed = mesh.edges()(other_edge_id, 1);

                            // Skip adjacent edges
                            if (ea == ec || ea == ed || eb == ec || eb == ed) {
                                continue;
                            }

                            // Skip EE pairs entirely when dbar_factor == 0
                            if (skip_ee_grad) {
                                continue;
                            }

                            if (auto iter =
                                    collisions.edge_edge_collisions.find(
                                        std::make_pair(edge_id, other_edge_id));
                                iter != collisions.edge_edge_collisions.end()) {

                                const auto dtype = iter->second->ee_dtype();

                                // Skip non EA_EB collision types
                                if (dtype != EdgeEdgeDistanceType::EA_EB) {
                                    continue;
                                }

                                Eigen::Vector<double, 12> positions;
                                positions << X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose();

                                Eigen::Matrix<T, 4, 3> positionsT =
                                    slice_positions<T, 4, 3>(positions);

                                const T dist = sqrt(
                                    edge_edge_sqr_distance<T>(
                                        positionsT.row(0), positionsT.row(1),
                                        positionsT.row(2), positionsT.row(3),
                                        dtype));

                                const T uv = closest_point_uv<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3),
                                    dtype);

                                const Eigen::RowVector3<T> ee_closest_point_T =
                                    uv * (positionsT.row(1) - positionsT.row(0))
                                    + positionsT.row(0);
                                const Eigen::RowVector3<double>
                                    ee_closest_point(
                                        ee_closest_point_T(0).val,
                                        ee_closest_point_T(1).val,
                                        ee_closest_point_T(2).val);

                                const T dist_sqr = edge_edge_sqr_distance<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3),
                                    dtype);
                                const auto mtypes = edge_edge_mollifier_type(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), dist_sqr.val);

                                T mollifier =
                                    Math<T>::cubic_spline(dist / params.dbar)
                                    * 1.5;
                                mollifier *= edge_edge_mollifier<T>(
                                    positionsT.row(0).transpose(),
                                    positionsT.row(1).transpose(),
                                    positionsT.row(2).transpose(),
                                    positionsT.row(3).transpose(), mtypes,
                                    dist_sqr);
                                mollifier = pow_int(
                                    mollifier,
                                    mollifier_order_for_barrier(
                                        params.barrier));

                                const ESPCollisionDict<PointType::EDGE>& dict =
                                    *(iter->second);

                                VertexMatrixView<3> X_extended(
                                    X, ee_closest_point);
                                assert(X_extended.rows() == X.rows() + 1);
                                assert(
                                    X_extended.m_A == X.data()
                                    && "VertexMatrixView has made a deepcopy!");

                                double P;
                                if (use_nf_grad) {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            dtype, *nf_barrier);
                                } else {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            dtype);
                                }
                                Eigen::VectorXd grad_P;
                                if (use_nf_grad) {
                                    grad_P = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                                            T>(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            ee_closest_point_T, *nf_barrier);
                                } else {
                                    grad_P = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions<
                                            T>(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            ee_closest_point_T);
                                }

                                ee_cache.push_back(
                                    { &dict, mollifier.val, mollifier.grad, P,
                                      grad_P });
                                total_w += mollifier.val;
                                total_p += mollifier.val * P;
                                if (use_nf_grad) {
                                    total_w_near += mollifier.val;
                                    total_p_near += mollifier.val * P;
                                }
                            }
                        }
                    }

                    // Face-interior quadrature points
                    const auto& face_quad_rule = params.get_quad_rule();
                    {
                        auto iter = collisions.face_collisions.find(f);
                        for (size_t qi = 0; qi < face_quad_rule.size(); qi++) {
                            const auto& qp = face_quad_rule[qi];
                            const double qp_weight_scale =
                                face_quadrature_weight_scale * qp.weight;
                            total_w += qp_weight_scale;
                            if (use_nf_grad) {
                                total_w_near += qp_weight_scale;
                                total_w_far += qp_weight_scale;
                            }
                            if (iter != collisions.face_collisions.end()
                                && qi < iter->second.size()) {
                                const auto& dict = *iter->second[qi];
                                const Eigen::RowVector3d q_pos =
                                    qp.lambda[0] * X.row(mesh.faces()(f, 0))
                                    + qp.lambda[1] * X.row(mesh.faces()(f, 1))
                                    + qp.lambda[2] * X.row(mesh.faces()(f, 2));
                                VertexMatrixView<3> X_qp(X, q_pos);

                                if (use_nf_grad) {
                                    auto [P_n, P_f] = PointPotentialHelper::
                                        evaluate_potential_at_face_center_with_cached_collisions_nearfar(
                                            X_qp, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    auto [grad_n, grad_f] = PointPotentialHelper::
                                        evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            qp.lambda, *nf_barrier);
                                    const_cache.push_back(
                                        ConstGradEntry {
                                            &dict.dofs(),
                                            qp_weight_scale * grad_n,
                                            qp_weight_scale * grad_f,
                                            qp_weight_scale * P_n,
                                            qp_weight_scale * P_f });
                                    total_p_near += qp_weight_scale * P_n;
                                    total_p_far += qp_weight_scale * P_f;
                                } else {
                                    const double P = PointPotentialHelper::
                                        evaluate_potential_at_face_center_with_cached_collisions(
                                            X_qp, dict, params,
                                            collisions.adaptive_dhat.get());
                                    const Eigen::VectorXd grad_P =
                                        PointPotentialHelper::
                                            evaluate_potential_gradient_at_face_interior_point_with_cached_collisions(
                                                X_qp, dict, params,
                                                collisions.adaptive_dhat.get(),
                                                qp.lambda);
                                    const_cache.push_back(
                                        ConstGradEntry {
                                            &dict.dofs(),
                                            qp_weight_scale * grad_P,
                                            Eigen::VectorXd::Zero(0),
                                            qp_weight_scale * P, 0 });
                                    total_p += qp_weight_scale * P;
                                }
                            }
                        }
                    }

                    if (face_quad_rule.empty()) {
                        for (index_t lv = 0; lv < 3; lv++) {
                            const index_t v = mesh.faces()(f, lv);
                            total_w += 1.;
                            if (use_nf_grad) {
                                total_w_near += 1.0;
                                total_w_far += 1.0;
                            }
                            if (auto iter =
                                    collisions.vertex_collisions.find(v);
                                iter != collisions.vertex_collisions.end()) {
                                if (use_nf_grad) {
                                    auto [P_n, P_f] = PointPotentialHelper::
                                        evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                                            X, (*iter->second), params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    auto [grad_n, grad_f] = PointPotentialHelper::
                                        evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                                            X, (*iter->second), params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    const_cache.push_back(
                                        ConstGradEntry {
                                            &(*iter->second).dofs(), grad_n,
                                            grad_f, P_n, P_f });
                                    total_p_near += P_n;
                                    total_p_far += P_f;
                                } else {
                                    const double P = PointPotentialHelper::
                                        evaluate_potential_at_vertex_with_cached_collisions(
                                            X, (*iter->second), params,
                                            collisions.adaptive_dhat.get());
                                    const Eigen::VectorXd grad_P =
                                        PointPotentialHelper::
                                            evaluate_potential_gradient_at_vertex_with_cached_collisions(
                                                X, (*iter->second), params,
                                                collisions.adaptive_dhat.get());
                                    const_cache.push_back(
                                        ConstGradEntry {
                                            &(*iter->second).dofs(), grad_P,
                                            Eigen::VectorXd::Zero(0), P, 0 });
                                    total_p += P;
                                }
                            }
                        }
                    }

                    // Pass 2: apply gradient
                    if (use_nf_grad) {
                        assert(total_w_near > 0);
                        if (total_w_far == 0) {
                            assert(total_p_far == 0);
                        }
                        const double avg_P_near = total_p_near / total_w_near;
                        for (const auto& e : ee_cache) {
                            grad(e.dict->dofs()) +=
                                (w / total_w_near * e.mol_val) * e.grad_P;
                            grad(e.dict->primary_dofs()) +=
                                (w / total_w_near * (e.P - avg_P_near))
                                * e.mol_grad;
                        }
                        for (const auto& e : const_cache) {
                            if (total_w_far > 0) {
                                grad(*e.dofs) +=
                                    (w / total_w_near) * e.grad_P_near
                                    + (w / total_w_far) * e.grad_P_far;
                            } else {
                                grad(*e.dofs) +=
                                    (w / total_w_near) * e.grad_P_near;
                            }
                        }
                    } else if (use_near_far) {
                        // Normalized but without NearFarBarrier splitting
                        // (dbar_factor not in (0,1))
                        assert(total_w > 0);
                        const double avg_P = total_p / total_w;
                        for (const auto& e : ee_cache) {
                            grad(e.dict->dofs()) +=
                                (w / total_w * e.mol_val) * e.grad_P;
                            grad(e.dict->primary_dofs()) +=
                                (w / total_w * (e.P - avg_P)) * e.mol_grad;
                        }
                        for (const auto& e : const_cache) {
                            grad(*e.dofs) += (w / total_w) * e.grad_P_near;
                        }
                    } else {
                        // Unnormalized
                        for (const auto& e : ee_cache) {
                            grad(e.dict->dofs()) += w * e.mol_val * e.grad_P;
                            grad(e.dict->primary_dofs()) +=
                                w * e.P * e.mol_grad;
                        }
                        for (const auto& e : const_cache) {
                            grad(*e.dofs) += w * e.grad_P_near;
                        }
                    }
                }
            };

            tbb::parallel_for(
                tbb::blocked_range<index_t>(0, mesh.num_faces()), loop_body);
        }
    }

    Eigen::VectorXd grad;
    grad.setZero(X.size());
    for (const auto& local_storage : storage) {
        grad += local_storage;
    }

    return grad;
}

Eigen::SparseMatrix<double> ESPPotential::hessian(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    IPC_PROFILE_SCOPE("ho.potential_hessian");
    assert(X.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(X.size(), X.size());
    }

    const int dim = X.cols();
    const int ndof = X.size();

    const int max_triplets_size = int(1e7);
    const int buffer_size = std::min(max_triplets_size, ndof);
    tbb::enumerable_thread_specific<LocalThreadMatStorage> storage(
        LocalThreadMatStorage(buffer_size, ndof, ndof));

    if (mesh.dim() == 2) {
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(params.quad_order);

        std::vector<index_t> active_edges;
        active_edges.reserve(collisions.edge_collisions_2d.size());
        for (const auto& [ei, _] : collisions.edge_collisions_2d) {
            active_edges.push_back(ei);
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, active_edges.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                auto& hess_triplets = storage.local();

                for (size_t k = r.begin(); k < r.end(); ++k) {
                    const index_t ei = active_edges[k];
                    const auto& qp_dicts = collisions.edge_collisions_2d.at(ei);
                    const double L = mesh.edge_area(ei);
                    const double w_edge = params.area_weights ? L : 1.;
                    const index_t e0 = mesh.edges()(ei, 0);
                    const index_t e1 = mesh.edges()(ei, 1);

                    for (size_t qi = 0; qi < qp_dicts.size(); ++qi) {
                        const auto& dict = *qp_dicts[qi];
                        if (dict.size() == 0) {
                            continue;
                        }
                        const auto& qp = rule[qi];
                        const std::array<double, 2> lambda = { { 1.0 - qp.xi,
                                                                 qp.xi } };
                        const Eigen::RowVector2d q_pos =
                            lambda[0] * X.row(e0) + lambda[1] * X.row(e1);
                        VertexMatrixView<2> X_ext(X, q_pos);

                        const Eigen::MatrixXd local_hess =
                            w_edge * qp.weight
                            * PointPotentialHelper::
                                evaluate_potential_hessian_at_edge_qp(
                                    X_ext, dict, params,
                                    collisions.adaptive_dhat.get(), lambda,
                                    project_hessian_to_psd);

                        ProfileRegistry::instance().add_value(
                            "ho.local_hessian.size", local_hess.rows());
                        local_hessian_to_global_triplets(
                            local_hess, dict.vertex_ids(), dim,
                            *(hess_triplets.cache));
                    }
                }
            });
    } else if (mesh.dim() == 3) {
        // When use_near_far is on, the per-face hessian is assembled as
        //   Term A (sum of per-stencil H(p_i))
        // + Term B (negative weighted sum of H(mol_i))
        // + Term C (sign-indefinite cross terms ~ sym(G⊗∇Z)).
        // PSD-projecting Terms A and B individually is not enough — Term C is
        // never PSD on its own. To still guarantee a PSD per-face contribution
        // (and thus a PSD global hessian), defer all per-stencil projections
        // and project the assembled per-face block once, over the union of
        // involved DOFs. With use_near_far = false the original
        // local-projection path is preserved exactly.
        const double dbar_factor_hess = params.dbar_factor();
        const bool use_nf_hess =
            use_near_far && dbar_factor_hess > 0 && dbar_factor_hess < 1;
        const bool skip_ee_hess = (dbar_factor_hess == 0);
        const bool combined_psd_projection =
            use_nf_hess && project_hessian_to_psd != PSDProjectionMethod::NONE;
        const PSDProjectionMethod inner_psd_method = combined_psd_projection
            ? PSDProjectionMethod::NONE
            : project_hessian_to_psd;
        {
            using T = ADHessian<12>;

            auto loop_body = [&](const tbb::blocked_range<index_t>& r) {
                auto& hess_triplets = storage.local();
                for (index_t f = r.begin(); f < r.end(); f++) {
                    const double area = mesh.face_areas()(f);
                    const double w = params.area_weights ? (area / 9.) : 1.;

                    // Pass 1: collect all quadrature contributions for this
                    // face
                    struct EEHessEntry {
                        const ESPCollisionDict<PointType::EDGE>* dict;
                        double mol_val;
                        Eigen::Vector<double, 12> mol_grad; // on primary_dofs
                        Eigen::Matrix<double, 12, 12>
                            mol_hess; // H(mol) on primary_dofs
                        double P;     // only near component for use_nf_hess
                        Eigen::VectorXd grad_P;     // indexed by dict->dofs()
                        Eigen::MatrixXd local_hess; // H(mol*P), PSD-projected
                    };
                    struct ConstHessEntry {
                        const std::vector<index_t>* vertex_ids;
                        const std::vector<index_t>* dofs;
                        double P_near, P_far;
                        Eigen::VectorXd grad_P_near,
                            grad_P_far; // indexed by dofs
                        Eigen::MatrixXd local_hess_near,
                            local_hess_far; // H(P_near), H(P_far)
                    };
                    std::vector<EEHessEntry> ee_cache;
                    std::vector<ConstHessEntry> const_cache;
                    double total_w = 0;
                    double total_p = 0;
                    double total_w_near = 0, total_p_near = 0;
                    double total_w_far = 0, total_p_far = 0;

                    // Construct NearFarBarrier if needed
                    std::unique_ptr<NearFarBarrier> nf_barrier;
                    if (use_nf_hess) {
                        nf_barrier = std::make_unique<NearFarBarrier>(
                            params.barrier, params.dbar_factor());
                    }

                    for (index_t le = 0; le < 3; le++) {
                        const index_t edge_id = mesh.faces_to_edges()(f, le);
                        const index_t ea = mesh.edges()(edge_id, 0);
                        const index_t eb = mesh.edges()(edge_id, 1);

                        for (index_t other_edge_id :
                             collisions.m_candidates.ee_set(edge_id)) {
                            const index_t ec = mesh.edges()(other_edge_id, 0);
                            const index_t ed = mesh.edges()(other_edge_id, 1);

                            // Skip adjacent edges
                            if (ea == ec || ea == ed || eb == ec || eb == ed) {
                                continue;
                            }

                            // Skip EE pairs entirely when dbar_factor == 0
                            if (skip_ee_hess) {
                                continue;
                            }

                            if (auto iter =
                                    collisions.edge_edge_collisions.find(
                                        std::make_pair(edge_id, other_edge_id));
                                iter != collisions.edge_edge_collisions.end()) {

                                const auto dtype = iter->second->ee_dtype();

                                // Skip non EA_EB collision types
                                if (dtype != EdgeEdgeDistanceType::EA_EB) {
                                    continue;
                                }

                                Eigen::Vector<double, 12> positions;
                                positions << X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose();

                                Eigen::Matrix<T, 4, 3> positionsT =
                                    slice_positions<T, 4, 3>(positions);

                                const T dist = sqrt(
                                    edge_edge_sqr_distance<T>(
                                        positionsT.row(0), positionsT.row(1),
                                        positionsT.row(2), positionsT.row(3),
                                        dtype));

                                const T uv = closest_point_uv<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3),
                                    dtype);

                                const Eigen::RowVector3<T> ee_closest_point_T =
                                    uv * (positionsT.row(1) - positionsT.row(0))
                                    + positionsT.row(0);
                                const Eigen::RowVector3<double>
                                    ee_closest_point(
                                        ee_closest_point_T(0).val,
                                        ee_closest_point_T(1).val,
                                        ee_closest_point_T(2).val);

                                const T dist_sqr = edge_edge_sqr_distance<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3),
                                    dtype);
                                const auto mtypes = edge_edge_mollifier_type(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), dist_sqr.val);

                                T mollifier =
                                    Math<T>::cubic_spline(dist / params.dbar)
                                    * 1.5;
                                mollifier *= edge_edge_mollifier<T>(
                                    positionsT.row(0).transpose(),
                                    positionsT.row(1).transpose(),
                                    positionsT.row(2).transpose(),
                                    positionsT.row(3).transpose(), mtypes,
                                    dist_sqr);
                                mollifier = pow_int(
                                    mollifier,
                                    mollifier_order_for_barrier(
                                        params.barrier));

                                const ESPCollisionDict<PointType::EDGE>& dict =
                                    *(iter->second);

                                VertexMatrixView<3> X_extended(
                                    X, ee_closest_point);
                                assert(
                                    X_extended.m_A == X.data()
                                    && "VertexMatrixView has made a deepcopy!");

                                double P;
                                Eigen::VectorXd grad_P;
                                Eigen::MatrixXd base_hess;
                                if (use_nf_hess) {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            dtype, *nf_barrier);
                                    grad_P = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                                            T>(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            ee_closest_point_T, *nf_barrier);
                                    base_hess = PointPotentialHelper::
                                        evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions_near(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            ee_closest_point_T, *nf_barrier);
                                } else {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            dtype);
                                    grad_P = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions<
                                            T>(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            ee_closest_point_T);
                                    base_hess = PointPotentialHelper::
                                        evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions(
                                            X_extended, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            ee_closest_point_T);
                                }
                                Eigen::MatrixXd local_hess =
                                    base_hess * mollifier.val;

                                for (index_t i = 0; i < 4; i++) {
                                    for (index_t j = 0; j < 4; j++) {
                                        local_hess.block<3, 3>(
                                            dict.primary_local_ids()[i] * 3,
                                            dict.primary_local_ids()[j] * 3) +=
                                            P
                                            * mollifier.Hess.block<3, 3>(
                                                i * 3, j * 3);
                                    }
                                }

                                // GCC false-positive: mollifier.grad is a
                                // fixed-size member, not a pointer.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
                                for (index_t i = 0; i < 4; i++) {
                                    const Eigen::MatrixXd tmp =
                                        mollifier.grad.segment<3>(i * 3)
                                        * grad_P.transpose();
                                    local_hess.middleRows(
                                        dict.primary_local_ids()[i] * 3, 3) +=
                                        tmp;
                                    local_hess.middleCols(
                                        dict.primary_local_ids()[i] * 3, 3) +=
                                        tmp.transpose();
                                }
#pragma GCC diagnostic pop

                                if (project_hessian_to_psd
                                        != PSDProjectionMethod::NONE
                                    && !combined_psd_projection) {
                                    ProfileRegistry::instance().add_value(
                                        "ho.psd_projection.size",
                                        local_hess.rows());
                                    local_hess = project_to_psd(
                                        local_hess, project_hessian_to_psd);
                                }

                                ee_cache.push_back(
                                    { &dict, mollifier.val, mollifier.grad,
                                      mollifier.Hess, P, grad_P,
                                      std::move(local_hess) });
                                total_w += mollifier.val;
                                total_p += mollifier.val * P;
                                if (use_nf_hess) {
                                    total_w_near += mollifier.val;
                                    total_p_near += mollifier.val * P;
                                }
                            }
                        }
                    }

                    // Face-interior quadrature points
                    const auto& face_quad_rule = params.get_quad_rule();
                    {
                        auto iter = collisions.face_collisions.find(f);
                        for (size_t qi = 0; qi < face_quad_rule.size(); qi++) {
                            const auto& qp = face_quad_rule[qi];
                            total_w += face_quadrature_weight_scale * qp.weight;
                            if (use_nf_hess) {
                                total_w_near +=
                                    face_quadrature_weight_scale * qp.weight;
                                total_w_far +=
                                    face_quadrature_weight_scale * qp.weight;
                            }
                            if (iter != collisions.face_collisions.end()
                                && qi < iter->second.size()) {
                                const auto& dict = *iter->second[qi];
                                const Eigen::RowVector3d q_pos =
                                    qp.lambda[0] * X.row(mesh.faces()(f, 0))
                                    + qp.lambda[1] * X.row(mesh.faces()(f, 1))
                                    + qp.lambda[2] * X.row(mesh.faces()(f, 2));
                                VertexMatrixView<3> X_qp(X, q_pos);
                                ConstHessEntry entry;
                                entry.vertex_ids = &dict.vertex_ids();
                                entry.dofs = &dict.dofs();

                                if (use_nf_hess) {
                                    auto [P_n, P_f] = PointPotentialHelper::
                                        evaluate_potential_at_face_center_with_cached_collisions_nearfar(
                                            X_qp, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    auto [grad_n, grad_f] = PointPotentialHelper::
                                        evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            qp.lambda, *nf_barrier);
                                    auto [hess_n, hess_f] = PointPotentialHelper::
                                        evaluate_potential_hessian_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            qp.lambda, inner_psd_method,
                                            *nf_barrier);
                                    entry.P_near = face_quadrature_weight_scale
                                        * qp.weight * P_n;
                                    entry.P_far = face_quadrature_weight_scale
                                        * qp.weight * P_f;
                                    entry.grad_P_near =
                                        face_quadrature_weight_scale * qp.weight
                                        * grad_n;
                                    entry.grad_P_far =
                                        face_quadrature_weight_scale * qp.weight
                                        * grad_f;
                                    entry.local_hess_near =
                                        face_quadrature_weight_scale * qp.weight
                                        * hess_n;
                                    entry.local_hess_far =
                                        face_quadrature_weight_scale * qp.weight
                                        * hess_f;
                                    total_p_near += entry.P_near;
                                    total_p_far += entry.P_far;
                                } else {
                                    entry.P_near =
                                        face_quadrature_weight_scale * qp.weight
                                        * PointPotentialHelper::
                                            evaluate_potential_at_face_center_with_cached_collisions(
                                                X_qp, dict, params,
                                                collisions.adaptive_dhat.get());
                                    entry.grad_P_near =
                                        face_quadrature_weight_scale * qp.weight
                                        * PointPotentialHelper::
                                            evaluate_potential_gradient_at_face_interior_point_with_cached_collisions(
                                                X_qp, dict, params,
                                                collisions.adaptive_dhat.get(),
                                                qp.lambda);
                                    entry.local_hess_near =
                                        face_quadrature_weight_scale * qp.weight
                                        * PointPotentialHelper::
                                            evaluate_potential_hessian_at_face_interior_point_with_cached_collisions(
                                                X_qp, dict, params,
                                                collisions.adaptive_dhat.get(),
                                                qp.lambda, inner_psd_method);
                                    entry.P_far = 0;
                                    entry.grad_P_far = Eigen::VectorXd::Zero(0);
                                    entry.local_hess_far =
                                        Eigen::MatrixXd::Zero(0, 0);
                                    total_p += entry.P_near;
                                }
                                const_cache.push_back(std::move(entry));
                            }
                        }
                    }

                    if (face_quad_rule.empty()) {
                        for (index_t lv = 0; lv < 3; lv++) {
                            const index_t v = mesh.faces()(f, lv);
                            total_w += 1.;
                            if (use_nf_hess) {
                                total_w_near += 1.0;
                                total_w_far += 1.0;
                            }
                            if (auto iter =
                                    collisions.vertex_collisions.find(v);
                                iter != collisions.vertex_collisions.end()) {
                                const auto& dict = *iter->second;
                                ConstHessEntry entry;
                                entry.vertex_ids = &dict.vertex_ids();
                                entry.dofs = &dict.dofs();

                                if (use_nf_hess) {
                                    auto [P_n, P_f] = PointPotentialHelper::
                                        evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                                            X, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    auto [grad_n, grad_f] = PointPotentialHelper::
                                        evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                                            X, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            *nf_barrier);
                                    auto [hess_n, hess_f] = PointPotentialHelper::
                                        evaluate_potential_hessian_at_vertex_with_cached_collisions_nearfar(
                                            X, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            inner_psd_method, *nf_barrier);
                                    entry.P_near = P_n;
                                    entry.P_far = P_f;
                                    entry.grad_P_near = grad_n;
                                    entry.grad_P_far = grad_f;
                                    entry.local_hess_near = hess_n;
                                    entry.local_hess_far = hess_f;
                                    total_p_near += entry.P_near;
                                    total_p_far += entry.P_far;
                                } else {
                                    entry.P_near = PointPotentialHelper::
                                        evaluate_potential_at_vertex_with_cached_collisions(
                                            X, dict, params,
                                            collisions.adaptive_dhat.get());
                                    entry.grad_P_near = PointPotentialHelper::
                                        evaluate_potential_gradient_at_vertex_with_cached_collisions(
                                            X, dict, params,
                                            collisions.adaptive_dhat.get());
                                    entry
                                        .local_hess_near = PointPotentialHelper::
                                        evaluate_potential_hessian_at_vertex_with_cached_collisions(
                                            X, dict, params,
                                            collisions.adaptive_dhat.get(),
                                            inner_psd_method);
                                    entry.P_far = 0;
                                    entry.grad_P_far = Eigen::VectorXd::Zero(0);
                                    entry.local_hess_far =
                                        Eigen::MatrixXd::Zero(0, 0);
                                    total_p += entry.P_near;
                                }
                                const_cache.push_back(std::move(entry));
                            }
                        }
                    }

                    // Pass 2: apply hessian
                    if (use_nf_hess) {
                        assert(total_w_near > 0);
                        if (total_w_far == 0) {
                            assert(total_p_far == 0);
                        }
                        // Apply quotient rule separately for near and far
                        // components
                        const double avg_P_near = total_p_near / total_w_near;
                        const double scale_C_near =
                            -(w / (total_w_near * total_w_near));

                        if (combined_psd_projection) {
                            // Nothing contributes from this face: skip combined
                            // block.
                            if (ee_cache.empty() && const_cache.empty()) {
                                continue;
                            }
                            // Build the union of vertex IDs across every
                            // contributing stencil so the entire face block can
                            // be assembled into one dense matrix and projected
                            // once. This is the only way Term C
                            // (sym(G⊗∇Z), sign-indefinite) can be made PSD.
                            std::vector<index_t> union_vids;
                            union_vids.reserve(
                                4 * ee_cache.size() + 3 * const_cache.size());
                            for (const auto& e : ee_cache) {
                                for (index_t vid : e.dict->vertex_ids()) {
                                    union_vids.push_back(vid);
                                }
                                for (index_t vid :
                                     e.dict->primary_vertex_ids()) {
                                    if (vid >= 0) {
                                        union_vids.push_back(vid);
                                    }
                                }
                            }
                            for (const auto& e : const_cache) {
                                for (index_t vid : *e.vertex_ids) {
                                    union_vids.push_back(vid);
                                }
                            }
                            std::sort(union_vids.begin(), union_vids.end());
                            union_vids.erase(
                                std::unique(
                                    union_vids.begin(), union_vids.end()),
                                union_vids.end());

                            std::unordered_map<index_t, int> vid_to_local;
                            vid_to_local.reserve(union_vids.size());
                            for (int i = 0;
                                 i < static_cast<int>(union_vids.size()); i++) {
                                vid_to_local[union_vids[i]] = i;
                            }

                            const int n_union_dofs =
                                static_cast<int>(union_vids.size()) * dim;
                            Eigen::MatrixXd H_face = Eigen::MatrixXd::Zero(
                                n_union_dofs, n_union_dofs);

                            // Scatter a (n*dim)x(n*dim) block keyed by vertex
                            // ids into H_face using the union vid mapping.
                            auto add_block = [&](const Eigen::MatrixXd& block,
                                                 auto&& vids, double scale) {
                                const int n =
                                    static_cast<int>(block.rows()) / dim;
                                std::vector<int> lvi(n);
                                for (int i = 0; i < n; i++) {
                                    lvi[i] = vid_to_local.at(vids[i]);
                                }
                                for (int i = 0; i < n; i++) {
                                    for (int j = 0; j < n; j++) {
                                        H_face.block(
                                            lvi[i] * dim, lvi[j] * dim, dim,
                                            dim) += scale
                                            * block.block(
                                                i * dim, j * dim, dim, dim);
                                    }
                                }
                            };

                            // Maps a global DOF index (= vid*dim + d) to its
                            // local index inside H_face.
                            auto global_dof_to_local = [&](index_t gd) {
                                return vid_to_local.at(gd / dim) * dim
                                    + static_cast<int>(gd % dim);
                            };

                            // Adds scale * sym(outer(g_vec, gradz_vec)) to
                            // H_face.
                            auto add_sym_correction_dense =
                                [&](const std::vector<index_t>& g_dofs,
                                    const Eigen::Ref<const Eigen::VectorXd>&
                                        g_vec,
                                    const std::vector<index_t>& gradz_dofs,
                                    const Eigen::Ref<const Eigen::VectorXd>&
                                        gradz_vec,
                                    double scale) {
                                    for (int a = 0;
                                         a < static_cast<int>(g_dofs.size());
                                         a++) {
                                        const int la =
                                            global_dof_to_local(g_dofs[a]);
                                        for (int b = 0; b < static_cast<int>(
                                                            gradz_dofs.size());
                                             b++) {
                                            const int lb = global_dof_to_local(
                                                gradz_dofs[b]);
                                            const double v =
                                                scale * g_vec[a] * gradz_vec[b];
                                            H_face(la, lb) += v;
                                            H_face(lb, la) += v;
                                        }
                                    }
                                };

                            // Term A: (w/total_w_near) * H(p_near) +
                            // (w/total_w_far) * H(p_far)
                            for (const auto& e : ee_cache) {
                                add_block(
                                    e.local_hess, e.dict->vertex_ids(),
                                    w / total_w_near);
                            }
                            for (const auto& e : const_cache) {
                                add_block(
                                    e.local_hess_near, *e.vertex_ids,
                                    w / total_w_near);
                                if (total_w_far > 0) {
                                    add_block(
                                        e.local_hess_far, *e.vertex_ids,
                                        w / total_w_far);
                                }
                            }

                            // Term B: -(w*avg_P_near/total_w_near) * Σ_i
                            // H(mol_i) (Z_far has no V-dependence, so no far
                            // Term B.)
                            const double scale_B_near =
                                -(w * avg_P_near / total_w_near);
                            for (const auto& e : ee_cache) {
                                add_block(
                                    e.mol_hess, e.dict->primary_vertex_ids(),
                                    scale_B_near);
                            }

                            // Term C: -(w/total_w_near²) * sym(G_near ⊗
                            // ∇Z_near) (Z_far has no V-dependence, so no far
                            // Term C.)
                            for (const auto& ei : ee_cache) {
                                const auto& prim_dofs_i =
                                    ei.dict->primary_dofs();
                                const Eigen::Vector<double, 12>& mol_grad_i =
                                    ei.mol_grad;
                                // EE-EE near interactions
                                for (const auto& ek : ee_cache) {
                                    add_sym_correction_dense(
                                        ek.dict->primary_dofs(),
                                        (ek.P - avg_P_near) * ek.mol_grad,
                                        prim_dofs_i, mol_grad_i, scale_C_near);
                                    add_sym_correction_dense(
                                        ek.dict->dofs(), ek.mol_val * ek.grad_P,
                                        prim_dofs_i, mol_grad_i, scale_C_near);
                                }
                                // EE-const near interactions
                                for (const auto& ej : const_cache) {
                                    add_sym_correction_dense(
                                        *ej.dofs, ej.grad_P_near, prim_dofs_i,
                                        mol_grad_i, scale_C_near);
                                }
                            }

                            ProfileRegistry::instance().add_value(
                                "ho.psd_projection.size", H_face.rows());
                            H_face =
                                project_to_psd(H_face, project_hessian_to_psd);

                            ProfileRegistry::instance().add_value(
                                "ho.local_hessian.size", H_face.rows());
                            local_hessian_to_global_triplets(
                                H_face, union_vids, dim,
                                *(hess_triplets.cache));
                        } else {
                            // Adds scale * sym(outer(g_vec, gradz_vec)) to
                            // triplets.
                            auto add_sym_correction =
                                [&](const std::vector<index_t>& g_dofs,
                                    const Eigen::Ref<const Eigen::VectorXd>&
                                        g_vec,
                                    const std::vector<index_t>& gradz_dofs,
                                    const Eigen::Ref<const Eigen::VectorXd>&
                                        gradz_vec,
                                    double scale) {
                                    for (int a = 0;
                                         a < static_cast<int>(g_dofs.size());
                                         a++) {
                                        for (int b = 0; b < static_cast<int>(
                                                            gradz_dofs.size());
                                             b++) {
                                            const double v =
                                                scale * g_vec[a] * gradz_vec[b];
                                            hess_triplets.cache->add_value(
                                                0, g_dofs[a], gradz_dofs[b], v);
                                            hess_triplets.cache->add_value(
                                                0, gradz_dofs[b], g_dofs[a], v);
                                        }
                                    }
                                };

                            // Term A: (w/total_w_near) * H(p_near) +
                            // (w/total_w_far) * H(p_far)
                            for (const auto& e : ee_cache) {
                                ProfileRegistry::instance().add_value(
                                    "ho.local_hessian.size",
                                    e.local_hess.rows());
                                local_hessian_to_global_triplets(
                                    (w / total_w_near) * e.local_hess,
                                    e.dict->vertex_ids(), dim,
                                    *(hess_triplets.cache));
                            }
                            for (const auto& e : const_cache) {
                                ProfileRegistry::instance().add_value(
                                    "ho.local_hessian.size",
                                    e.local_hess_near.rows());
                                local_hessian_to_global_triplets(
                                    (w / total_w_near) * e.local_hess_near,
                                    *e.vertex_ids, dim, *(hess_triplets.cache));
                                if (total_w_far > 0) {
                                    ProfileRegistry::instance().add_value(
                                        "ho.local_hessian.size",
                                        e.local_hess_far.rows());
                                    local_hessian_to_global_triplets(
                                        (w / total_w_far) * e.local_hess_far,
                                        *e.vertex_ids, dim,
                                        *(hess_triplets.cache));
                                }
                            }

                            // Term B: -(w*avg_P_near/total_w_near) * Σ_i
                            // H(mol_i) (Z_far has no V-dependence, so no far
                            // Term B.)
                            for (const auto& e : ee_cache) {
                                ProfileRegistry::instance().add_value(
                                    "ho.local_hessian.size", e.mol_hess.rows());
                                local_hessian_to_global_triplets(
                                    -(w * avg_P_near / total_w_near)
                                        * e.mol_hess,
                                    e.dict->primary_vertex_ids(), dim,
                                    *(hess_triplets.cache));
                            }

                            // Term C: -(w/total_w_near²) * sym(G_near ⊗
                            // ∇Z_near) (Z_far has no V-dependence, so no far
                            // Term C.)
                            for (const auto& ei : ee_cache) {
                                const auto& prim_dofs_i =
                                    ei.dict->primary_dofs();
                                const Eigen::Vector<double, 12>& mol_grad_i =
                                    ei.mol_grad;
                                // EE-EE near interactions
                                for (const auto& ek : ee_cache) {
                                    add_sym_correction(
                                        ek.dict->primary_dofs(),
                                        (ek.P - avg_P_near) * ek.mol_grad,
                                        prim_dofs_i, mol_grad_i, scale_C_near);
                                    add_sym_correction(
                                        ek.dict->dofs(), ek.mol_val * ek.grad_P,
                                        prim_dofs_i, mol_grad_i, scale_C_near);
                                }
                                // EE-const near interactions
                                for (const auto& ej : const_cache) {
                                    add_sym_correction(
                                        *ej.dofs, ej.grad_P_near, prim_dofs_i,
                                        mol_grad_i, scale_C_near);
                                }
                            }
                        }
                    } else if (use_near_far) {
                        // Normalized (use_near_far=true) but without
                        // NearFarBarrier splitting
                        assert(total_w > 0);
                        const double avg_P = total_p / total_w;
                        const double scale_C = -(w / (total_w * total_w));

                        auto add_sym_correction_norm =
                            [&](const std::vector<index_t>& g_dofs,
                                const Eigen::Ref<const Eigen::VectorXd>& g_vec,
                                const std::vector<index_t>& gradz_dofs,
                                const Eigen::Ref<const Eigen::VectorXd>&
                                    gradz_vec) {
                                for (int a = 0;
                                     a < static_cast<int>(g_dofs.size()); a++) {
                                    for (int b = 0; b
                                         < static_cast<int>(gradz_dofs.size());
                                         b++) {
                                        const double v =
                                            scale_C * g_vec[a] * gradz_vec[b];
                                        hess_triplets.cache->add_value(
                                            0, g_dofs[a], gradz_dofs[b], v);
                                        hess_triplets.cache->add_value(
                                            0, gradz_dofs[b], g_dofs[a], v);
                                    }
                                }
                            };

                        // Term A: (w/total_w) * H(p_sum)
                        for (const auto& e : ee_cache) {
                            local_hessian_to_global_triplets(
                                (w / total_w) * e.local_hess,
                                e.dict->vertex_ids(), dim,
                                *(hess_triplets.cache));
                        }
                        for (const auto& e : const_cache) {
                            local_hessian_to_global_triplets(
                                (w / total_w) * e.local_hess_near,
                                *e.vertex_ids, dim, *(hess_triplets.cache));
                        }

                        // Term B: -(w*avg_P/total_w) * Σ_i H(mol_i)
                        for (const auto& e : ee_cache) {
                            local_hessian_to_global_triplets(
                                -(w * avg_P / total_w) * e.mol_hess,
                                e.dict->primary_vertex_ids(), dim,
                                *(hess_triplets.cache));
                        }

                        // Term C: -(w/total_w²) * sym(G ⊗ ∇Z)
                        for (const auto& ei : ee_cache) {
                            const auto& prim_dofs_i = ei.dict->primary_dofs();
                            const Eigen::Vector<double, 12>& mol_grad_i =
                                ei.mol_grad;
                            for (const auto& ek : ee_cache) {
                                add_sym_correction_norm(
                                    ek.dict->primary_dofs(),
                                    (ek.P - avg_P) * ek.mol_grad, prim_dofs_i,
                                    mol_grad_i);
                                add_sym_correction_norm(
                                    ek.dict->dofs(), ek.mol_val * ek.grad_P,
                                    prim_dofs_i, mol_grad_i);
                            }
                            for (const auto& ej : const_cache) {
                                add_sym_correction_norm(
                                    *ej.dofs, ej.grad_P_near, prim_dofs_i,
                                    mol_grad_i);
                            }
                        }
                    } else {
                        // Unnormalized: H(w*p_sum) = w * H(p_sum)
                        for (const auto& e : ee_cache) {
                            ProfileRegistry::instance().add_value(
                                "ho.local_hessian.size", e.local_hess.rows());
                            local_hessian_to_global_triplets(
                                w * e.local_hess, e.dict->vertex_ids(), dim,
                                *(hess_triplets.cache));
                        }
                        for (const auto& e : const_cache) {
                            ProfileRegistry::instance().add_value(
                                "ho.local_hessian.size",
                                e.local_hess_near.rows());
                            local_hessian_to_global_triplets(
                                w * e.local_hess_near, *e.vertex_ids, dim,
                                *(hess_triplets.cache));
                        }
                    }
                }
            };

            tbb::parallel_for(
                tbb::blocked_range<index_t>(0, mesh.num_faces()), loop_body);
        }
    }

    Eigen::SparseMatrix<double> hess(ndof, ndof);

    // Assemble the stiffness matrix by concatenating the tuples in each local
    // storage

    // Collect thread storages
    std::vector<LocalThreadMatStorage*> storages(storage.size());
    int index = 0;
    for (auto& local_storage : storage) {
        storages[index++] = &local_storage;
    }

    tbb::parallel_for(size_t(0), storages.size(), [&](size_t i) {
        storages[i]->cache->prune();
    });

    if (storage.empty()) {
        return Eigen::SparseMatrix<double>();
    }

    // Prepares for parallel concatenation
    std::vector<int> offsets(storage.size());

    index = 0;
    int triplet_count = 0;
    for (auto& local_storage : storage) {
        offsets[index++] = triplet_count;
        triplet_count += local_storage.cache->triplet_count();
    }

    std::vector<Eigen::Triplet<double>> triplets;

    assert(!storages.empty());
    if (triplet_count >= triplets.max_size()) {
        // Serial fallback version in case the vector of triplets cannot be
        // allocated

        logger().warn(
            "Cannot allocate space for triplets, switching to serial assembly.");

        // Serially merge local storages
        for (LocalThreadMatStorage& local_storage : storage) {
            hess += local_storage.cache->get_matrix(false); // will also prune
        }
        hess.makeCompressed();
    } else {
        triplets.resize(triplet_count);

        // Parallel copy into triplets
        tbb::parallel_for(size_t(0), storages.size(), [&](size_t i) {
            const SparseMatrixCache& cache =
                dynamic_cast<const SparseMatrixCache&>(*storages[i]->cache);
            int offset = offsets[i];

            std::copy(
                cache.entries().begin(), cache.entries().end(),
                triplets.begin() + offset);
            offset += cache.entries().size();

            if (cache.mat().nonZeros() > 0) {
                int count = 0;
                for (int k = 0; k < cache.mat().outerSize(); ++k) {
                    for (Eigen::SparseMatrix<double>::InnerIterator it(
                             cache.mat(), k);
                         it; ++it) {
                        assert(count < cache.mat().nonZeros());
                        triplets[offset + count++] = Eigen::Triplet<double>(
                            it.row(), it.col(), it.value());
                    }
                }
            }
        });

        // Sort and assemble
        hess.setFromTriplets(triplets.begin(), triplets.end());
    }

    return hess;
}

double ESPPotential::operator()(
    const ESPCollision& collision,
    Eigen::ConstRef<Eigen::VectorXd> positions) const
{
    return collision.weight * collision(positions, params);
}

Eigen::VectorXd ESPPotential::gradient(
    const ESPCollision& collision,
    Eigen::ConstRef<Eigen::VectorXd> positions) const
{
    return collision.weight * collision.gradient(positions, params);
}

Eigen::MatrixXd ESPPotential::hessian(
    const ESPCollision& collision,
    Eigen::ConstRef<Eigen::VectorXd> positions,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    Eigen::MatrixXd hess =
        collision.weight * collision.hessian(positions, params);
    if (project_hessian_to_psd != PSDProjectionMethod::NONE) {
        ProfileRegistry::instance().add_value(
            "ho.psd_projection.size", hess.rows());
    }
    return project_to_psd(hess, project_hessian_to_psd);
}

} // namespace ipc