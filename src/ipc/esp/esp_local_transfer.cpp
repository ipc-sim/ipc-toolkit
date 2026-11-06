#include "esp_local_transfer.hpp"

#include <ipc/esp/esp_collision_maps.hpp>
#include <ipc/esp/esp_distance.hpp>
#include <ipc/esp/quadrature_potential.hpp>
#include <ipc/gcp/distance/mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <type_traits>

namespace ipc::detail {
namespace {
    // A scalar's value and derivatives on one element's union of real vertices.
    // Small-stencil AD handles geometry; these jets combine existing pointwise
    // derivatives, retaining cross-stencil terms in the quotient-rule Hessian.
    template <int Order> struct Jet {
        double value;
        Eigen::VectorXd grad;
        Eigen::MatrixXd hess;

        Jet(int ndof = 0, double v = 0) : value(v) { set_zero(ndof, v); }
        void set_zero(int ndof, double v = 0)
        {
            value = v;
            if constexpr (Order >= 1) {
                grad.setZero(ndof);
            }
            if constexpr (Order == 2) {
                hess.setZero(ndof, ndof);
            }
        }
        Jet& operator+=(const Jet& b)
        {
            value += b.value;
            if constexpr (Order >= 1) {
                grad += b.grad;
            }
            if constexpr (Order == 2) {
                hess += b.hess;
            }
            return *this;
        }
    };

    template <int O> Jet<O> scale(Jet<O> a, double s)
    {
        a.value *= s;
        if constexpr (O >= 1) {
            a.grad *= s;
        }
        if constexpr (O == 2) {
            a.hess *= s;
        }
        return a;
    }
    template <int O> Jet<O> product(const Jet<O>& a, const Jet<O>& b)
    {
        Jet<O> out(a.grad.size(), a.value * b.value);
        if constexpr (O >= 1) {
            out.grad = b.value * a.grad + a.value * b.grad;
        }
        if constexpr (O == 2) {
            out.hess = b.value * a.hess + a.value * b.hess
                + a.grad * b.grad.transpose() + b.grad * a.grad.transpose();
        }
        return out;
    }
    template <int O> Jet<O> quotient(const Jet<O>& a, const Jet<O>& b)
    {
        Jet<O> inv(b.grad.size(), 1.0 / b.value);
        if constexpr (O >= 1) {
            inv.grad = -b.grad / (b.value * b.value);
        }
        if constexpr (O == 2) {
            inv.hess = -b.hess / (b.value * b.value)
                + (2.0 / (b.value * b.value * b.value)) * b.grad
                    * b.grad.transpose();
        }
        return product(a, inv);
    }

    // Rest-shape weights of the integration elements (faces, rods,
    // particles), as in ESPPotential.
    double
    face_weight(const CollisionMesh& mesh, const ESPParameters& p, index_t f)
    {
        return p.area_weights ? mesh.face_areas()(f) * p.element_weight_scale(3)
                              : 1.0;
    }
    double
    rod_weight(const CollisionMesh& mesh, const ESPParameters& p, index_t e)
    {
        return p.area_weights ? mesh.edge_area(e) * p.element_weight_scale(2)
                              : 1.0;
    }

    // Each vertex sample's coefficient: the sum over its elements of the
    // element weight times its corner weight. Under transfer the fixed
    // samples are not normalized, so this is constant.
    Eigen::VectorXd
    vertex_coefficients(const CollisionMesh& mesh, const ESPParameters& p)
    {
        Eigen::VectorXd coef = Eigen::VectorXd::Zero(mesh.num_vertices());
        for (index_t f = 0; f < index_t(mesh.num_faces()); ++f) {
            const auto& cw = p.face_rule(mesh.face_areas()(f)).corner_weights;
            const double w = face_weight(mesh, p, f);
            for (int i = 0; i < 3; ++i) {
                coef[mesh.faces()(f, i)] += w * cw[i];
            }
        }
        for (const index_t e : mesh.codim_edges()) {
            const double w = rod_weight(mesh, p, e);
            coef[mesh.edges()(e, 0)] += 0.5 * w;
            coef[mesh.edges()(e, 1)] += 0.5 * w;
        }
        for (const index_t v : mesh.codim_vertices()) {
            coef[v] += p.area_weights ? mesh.vertex_area(v) : 1.0;
        }
        return coef;
    }

    // One element's share of an edge's EE term. Each face corner splits its
    // weight evenly between the face's two edges at that corner, so the EE
    // samples on an edge take weight from a pool holding half of each
    // endpoint's corner weight, and the pools of different edges are
    // independent. A rod's single edge holds both endpoints' full weight.
    // An EE sample takes at most c = the smaller endpoint share of the pool,
    // so its area scale compares c times the element's measure with the
    // patch (ESPParameters::ee_area_scaling).
    struct Pool {
        double coef;  // w_e c s_e, multiplying sum_j m_j (N(q_j) - R_j)
        double scale; // s_e, multiplying sum_j m_j in the denominator
    };

    void edge_pools(
        const CollisionMesh& mesh,
        const ESPParameters& p,
        index_t e,
        std::vector<Pool>& pools)
    {
        pools.clear();
        const double r = p.ee_area_scaling ? p.ee_area_reference() : 0.0;
        const index_t a = mesh.edges()(e, 0), b = mesh.edges()(e, 1);
        const auto& faces = mesh.edges_to_faces()[e];
        for (const index_t f : faces) {
            const double area = mesh.face_areas()(f);
            const auto& cw = p.face_rule(area).corner_weights;
            double wa = 0, wb = 0;
            for (int i = 0; i < 3; ++i) {
                if (mesh.faces()(f, i) == a) {
                    wa = cw[i];
                } else if (mesh.faces()(f, i) == b) {
                    wb = cw[i];
                }
            }
            const double c = 0.5 * std::min(wa, wb);
            const double s = ee_area_scale_factor(c * area, r, /*face=*/true);
            pools.push_back({ face_weight(mesh, p, f) * c * s, s });
        }
        if (faces.empty()) { // a rod
            const double c = 0.5;
            const double s =
                ee_area_scale_factor(c * mesh.edge_area(e), r, /*face=*/false);
            pools.push_back({ rod_weight(mesh, p, e) * c * s, s });
        }
    }

    template <typename T> double value(const T& x)
    {
        if constexpr (std::is_same_v<T, double>) {
            return x;
        } else {
            return x.val;
        }
    }

    template <int Order>
    using GeometryScalar = std::conditional_t<
        Order == 0,
        double,
        std::conditional_t<Order == 1, ADGrad<12>, ADHessian<12>>>;

    template <int Order> struct VertexSample {
        index_t vertex;
        const ESPCollisionDict<PointType::VERTEX>* dict;
        Jet<Order> near, far;
        VertexSample(index_t v, const ESPCollisionDict<PointType::VERTEX>* d)
            : vertex(v)
            , dict(d)
        {
        }
    };

    template <int Order> struct EdgeSample {
        index_t ea, eb;
        const ESPCollisionDict<PointType::EDGE>* dict;
        GeometryScalar<Order> u, m;
        Jet<Order> potential;
        EdgeSample(
            index_t a, index_t b, const ESPCollisionDict<PointType::EDGE>* d)
            : ea(a)
            , eb(b)
            , dict(d)
        {
        }
    };

    // All entries are rebuilt for this X, collision set, parameters, and
    // derivative order: each vertex sample and each EE sample is evaluated
    // once, and the EE samples are grouped by source edge.
    template <int Order> struct EvaluationCache {
        EvaluationCache() = default;
        EvaluationCache(const EvaluationCache&) = delete;
        EvaluationCache(EvaluationCache&&) = default;
        std::vector<VertexSample<Order>> vertices;
        std::vector<const VertexSample<Order>*> vertex_at;
        std::vector<EdgeSample<Order>> edges;
        std::vector<size_t> edge_offsets;
    };

    template <int Order>
    void evaluate_vertex_sample(
        VertexSample<Order>& sample,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params,
        const NearFarBarrier<>& nf)
    {
        const auto& dict = *sample.dict;
        std::tie(sample.near.value, sample.far.value) = PointPotentialHelper::
            evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                X, dict, params, nf);
        if constexpr (Order >= 1) {
            std::tie(sample.near.grad, sample.far.grad) = PointPotentialHelper::
                evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                    X, dict, params, nf);
        }
        if constexpr (Order == 2) {
            std::tie(sample.near.hess, sample.far.hess) = PointPotentialHelper::
                evaluate_potential_hessian_at_vertex_with_cached_collisions_nearfar(
                    X, dict, params, PSDProjectionMethod::NONE, nf);
        }
    }

    // The EE sample's source-edge parameter u, mollifier weight m (without
    // the element's area scale), and near potential at its closest point.
    template <int Order>
    void evaluate_edge_sample(
        EdgeSample<Order>& sample,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params,
        const NearFarBarrier<>& nf)
    {
        using T = GeometryScalar<Order>;
        const std::array<index_t, 4> prim = { {
            mesh.edges()(sample.ea, 0),
            mesh.edges()(sample.ea, 1),
            mesh.edges()(sample.eb, 0),
            mesh.edges()(sample.eb, 1),
        } };
        Eigen::Matrix<T, 4, 3> positions;
        if constexpr (Order == 0) {
            for (int j = 0; j < 4; ++j) {
                positions.row(j) = X.row(prim[j]);
            }
        } else {
            Eigen::Vector<double, 12> flat;
            for (int j = 0; j < 4; ++j) {
                flat.segment<3>(3 * j) = X.row(prim[j]);
            }
            positions = slice_positions<T, 4, 3>(flat);
        }
        sample.u = closest_point_uv<T>(
            positions.row(0), positions.row(1), positions.row(2),
            positions.row(3));
        const T d2 = line_line_sqr_distance<T>(
            positions.row(0), positions.row(1), positions.row(2),
            positions.row(3));
        const auto types = edge_edge_mollifier_type(
            X.row(prim[0]), X.row(prim[1]), X.row(prim[2]), X.row(prim[3]),
            value(d2), params.ee_mollifier);
        sample.m =
            Math<T>::cubic_spline(sqrt(d2) / params.ee_support())
            * params.ee_weight
            * edge_edge_mollifier<T>(
                positions.row(0).transpose(), positions.row(1).transpose(),
                positions.row(2).transpose(), positions.row(3).transpose(),
                types, d2, params.ee_mollifier);
        const T initial = sample.m;
        for (int j = 1; j < params.ee_mollifier_exponent(); ++j) {
            sample.m = sample.m * initial;
        }
        const Eigen::Vector3<T> q =
            (positions.row(0)
             + sample.u * (positions.row(1) - positions.row(0)))
                .transpose();
        Eigen::RowVector3d qvalue;
        for (int j = 0; j < 3; ++j) {
            qvalue(j) = value(q(j));
        }
        const VertexMatrixView<3> extended(X, qvalue);
        const auto& dict = *sample.dict;
        sample.potential.value = PointPotentialHelper::
            evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                extended, dict, params, dict.ee_dtype(), nf);
        if constexpr (Order >= 1) {
            sample.potential.grad = PointPotentialHelper::
                evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                    T>(extended, dict, params, q, nf);
        }
        if constexpr (Order == 2) {
            sample.potential.hess = PointPotentialHelper::
                evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions_near(
                    extended, dict, params, q, nf);
        }
    }

    template <int Order>
    EvaluationCache<Order> build_cache(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params)
    {
        EvaluationCache<Order> cache;
        const auto& maps = collisions.maps();
        cache.vertex_at.resize(mesh.num_vertices(), nullptr);
        cache.vertices.reserve(maps.vertex_collisions.size());
        for (const auto& [v, dict] : maps.vertex_collisions) {
            cache.vertices.emplace_back(v, dict.get());
            cache.vertex_at[v] = &cache.vertices.back();
        }

        struct Pair {
            index_t ea, eb;
            const ESPCollisionDict<PointType::EDGE>* dict;
        };
        std::vector<Pair> pairs;
        pairs.reserve(maps.edge_edge_collisions.size());
        for (const auto& [edges, dict] : maps.edge_edge_collisions) {
            pairs.push_back({ edges.first, edges.second, dict.get() });
        }
        std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
            return std::tie(a.ea, a.eb) < std::tie(b.ea, b.eb);
        });
        cache.edges.reserve(pairs.size());
        cache.edge_offsets.resize(mesh.num_edges() + 1, 0);
        for (size_t begin = 0; begin < pairs.size();) {
            const auto ea = pairs[begin].ea;
            // The old loop traverses this sorted list for every incident
            // face. Query it once per active edge, retaining the same pair
            // membership and target-edge order (including empty EE dicts).
            const auto neighbors = collisions.m_candidates.ee_set(mesh, ea);
            size_t end = begin;
            for (; end < pairs.size() && pairs[end].ea == ea; ++end) {
                const auto& pair = pairs[end];
                if (neighbors.contains(pair.eb)) {
                    cache.edges.emplace_back(ea, pair.eb, pair.dict);
                    ++cache.edge_offsets[ea + 1];
                }
            }
            begin = end;
        }
        std::partial_sum(
            cache.edge_offsets.begin(), cache.edge_offsets.end(),
            cache.edge_offsets.begin());

        const NearFarBarrier<> nf(params.barrier, params.dbar_factor());
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, cache.vertices.size()),
            [&](const auto& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    evaluate_vertex_sample(cache.vertices[i], X, params, nf);
                }
            });
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, cache.edges.size()),
            [&](const auto& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    evaluate_edge_sample(cache.edges[i], mesh, X, params, nf);
                }
            });
        return cache;
    }

    // Per-thread output of one evaluation.
    template <int Order> struct Result {
        double energy = 0;
        Eigen::VectorXd gradient;
        std::vector<Eigen::Triplet<double>> hessian;
        std::vector<index_t> ids; // scratch: a block's vertex ids
        std::vector<Pool> pools;  // scratch: an edge's pools
        explicit Result(int n)
        {
            if constexpr (Order == 1) {
                gradient = Eigen::VectorXd::Zero(n);
            }
        }

        // Add c * (value, gradient, Hessian) of a term on vertices @p vids.
        template <typename IDs>
        void
        add(const IDs& vids,
            double v,
            const Eigen::VectorXd& g,
            Eigen::MatrixXd h,
            PSDProjectionMethod projection,
            int num_vertices)
        {
            if constexpr (Order == 0) {
                energy += v;
            }
            if constexpr (Order == 1) {
                local_gradient_to_global_gradient(g, vids, 3, gradient);
            }
            if constexpr (Order == 2) {
                if (projection != PSDProjectionMethod::NONE) {
                    h = project_to_psd(h, projection);
                }
                local_hessian_to_global_triplets(
                    h, vids, 3, hessian, num_vertices);
            }
        }
    };

    // E = sum_e w_e [ sum_i w_i (N_i + F_i) + sum_q w_q (N_q + F_q) ]
    //   + sum_edges sum_pools w c s sum_j m_j (N(q_j) - R_j) / (1 + s M),
    // R_j = (1 - u_j) N_a + u_j N_b, M = sum_j m_j over the edge's EE
    // samples. The fixed samples (vertices, face-rule points) have constant
    // coefficients, so each is evaluated and projected on its own stencil.
    // Each edge's EE term couples only that edge's samples and its two
    // endpoints, so it is one block. A corner gives an edge at most
    // c s M / (1 + s M) < c <= half its weight, so every coefficient stays
    // positive.
    template <int Order>
    Result<Order> evaluate(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params,
        PSDProjectionMethod projection = PSDProjectionMethod::NONE)
    {
        const auto cache = build_cache<Order>(collisions, mesh, X, params);
        const Eigen::VectorXd vertex_coef = vertex_coefficients(mesh, params);
        const NearFarBarrier<> nf(params.barrier, params.dbar_factor());
        const int nv = int(mesh.num_vertices());
        const auto& maps = collisions.maps();
        tbb::enumerable_thread_specific<Result<Order>> storage(
            [&] { return Result<Order>(X.size()); });

        // Vertex samples.
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, cache.vertices.size()),
            [&](const auto& range) {
                auto& out = storage.local();
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const auto& s = cache.vertices[i];
                    const double c = vertex_coef[s.vertex];
                    Eigen::VectorXd g;
                    Eigen::MatrixXd h;
                    if constexpr (Order >= 1) {
                        g = c * (s.near.grad + s.far.grad);
                    }
                    if constexpr (Order == 2) {
                        h = c * (s.near.hess + s.far.hess);
                    }
                    out.add(
                        s.dict->vertex_ids(), c * (s.near.value + s.far.value),
                        g, std::move(h), projection, nv);
                }
            });

        // Face-rule points other than the corners.
        std::vector<typename std::decay_t<
            decltype(maps.face_collisions)>::const_pointer>
            faces;
        faces.reserve(maps.face_collisions.size());
        for (const auto& entry : maps.face_collisions) {
            faces.push_back(&entry);
        }
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, faces.size()),
            [&](const auto& range) {
                auto& out = storage.local();
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const index_t f = faces[i]->first;
                    const auto& dicts = faces[i]->second;
                    const auto& rule =
                        *params.face_rule(mesh.face_areas()(f)).non_corner;
                    const double w = face_weight(mesh, params, f);
                    for (size_t qi = 0;
                         qi < std::min(rule.size(), dicts.size()); ++qi) {
                        const auto& dict = *dicts[qi];
                        if (dict.size() == 0) {
                            continue;
                        }
                        const auto& qp = rule[qi];
                        Eigen::RowVector3d q = Eigen::RowVector3d::Zero();
                        for (int k = 0; k < 3; ++k) {
                            q += qp.lambda[k] * X.row(mesh.faces()(f, k));
                        }
                        const VertexMatrixView<3> Xq(X, q);
                        const double c = w * qp.weight;
                        const auto [pn, pf] = PointPotentialHelper::
                            evaluate_potential_at_face_interior_point_with_cached_collisions_nearfar(
                                Xq, dict, params, nf);
                        Eigen::VectorXd g;
                        Eigen::MatrixXd h;
                        if constexpr (Order >= 1) {
                            const auto [gn, gf] = PointPotentialHelper::
                                evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
                                    Xq, dict, params, qp.lambda, nf);
                            g = c * (gn + gf);
                        }
                        if constexpr (Order == 2) {
                            const auto [hn, hf] = PointPotentialHelper::
                                evaluate_potential_hessian_at_face_interior_point_with_cached_collisions_nearfar(
                                    Xq, dict, params, qp.lambda,
                                    PSDProjectionMethod::NONE, nf);
                            h = c * (hn + hf);
                        }
                        out.add(
                            dict.vertex_ids(), c * (pn + pf), g, std::move(h),
                            projection, nv);
                    }
                }
            });

        // EE samples, one block per source edge.
        std::vector<index_t> active_edges;
        for (index_t e = 0; e < index_t(mesh.num_edges()); ++e) {
            if (cache.edge_offsets[e + 1] > cache.edge_offsets[e]) {
                active_edges.push_back(e);
            }
        }
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, active_edges.size()),
            [&](const auto& range) {
                auto& out = storage.local();
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const index_t e = active_edges[i];
                    const size_t begin = cache.edge_offsets[e];
                    const size_t end = cache.edge_offsets[e + 1];
                    bool any = false;
                    for (size_t k = begin; k < end; ++k) {
                        any = any || value(cache.edges[k].m) != 0;
                    }
                    if (!any) {
                        continue; // every m, and its derivatives, vanish
                    }
                    edge_pools(mesh, params, e, out.pools);

                    const std::array<index_t, 2> ends = {
                        { mesh.edges()(e, 0), mesh.edges()(e, 1) }
                    };
                    auto& ids = out.ids;
                    ids.assign(ends.begin(), ends.end());
                    for (const index_t v : ends) {
                        if (const auto* s = cache.vertex_at[v]) {
                            const auto& more = s->dict->vertex_ids();
                            ids.insert(ids.end(), more.begin(), more.end());
                        }
                    }
                    for (size_t k = begin; k < end; ++k) {
                        const auto& s = cache.edges[k];
                        const auto& more = s.dict->vertex_ids();
                        ids.insert(ids.end(), more.begin(), more.end());
                        ids.push_back(mesh.edges()(s.eb, 0));
                        ids.push_back(mesh.edges()(s.eb, 1));
                    }
                    std::sort(ids.begin(), ids.end());
                    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
                    const int ndof = 3 * int(ids.size());
                    const auto local_index = [&](index_t v) {
                        const auto it =
                            std::lower_bound(ids.begin(), ids.end(), v);
                        assert(it != ids.end() && *it == v);
                        return int(it - ids.begin());
                    };
                    // A dict-local term as a jet on the block's vertices.
                    const auto scattered = [&](double v, const auto& vids,
                                               const Eigen::VectorXd& g,
                                               const Eigen::MatrixXd& h) {
                        Jet<Order> jet(ndof, v);
                        if constexpr (Order >= 1) {
                            const int n = int(vids.size());
                            std::vector<int> at(n);
                            for (int a = 0; a < n; ++a) {
                                at[a] = local_index(vids[a]);
                                jet.grad.template segment<3>(3 * at[a]) +=
                                    g.template segment<3>(3 * a);
                            }
                            if constexpr (Order == 2) {
                                for (int a = 0; a < n; ++a) {
                                    for (int b = 0; b < n; ++b) {
                                        jet.hess.template block<3, 3>(
                                            3 * at[a], 3 * at[b]) +=
                                            h.template block<3, 3>(
                                                3 * a, 3 * b);
                                    }
                                }
                            }
                        }
                        return jet;
                    };
                    // A function of an EE pair's four vertices as a jet.
                    const auto geometry =
                        [&](const GeometryScalar<Order>& x,
                            const std::array<index_t, 4>& prim) {
                            Jet<Order> jet(ndof, value(x));
                            if constexpr (Order >= 1) {
                                for (int a = 0; a < 4; ++a) {
                                    const int ia = local_index(prim[a]);
                                    jet.grad.template segment<3>(3 * ia) +=
                                        x.grad.template segment<3>(3 * a);
                                    if constexpr (Order == 2) {
                                        for (int b = 0; b < 4; ++b) {
                                            jet.hess.template block<3, 3>(
                                                3 * ia,
                                                3 * local_index(prim[b])) +=
                                                x.Hess.template block<3, 3>(
                                                    3 * a, 3 * b);
                                        }
                                    }
                                }
                            }
                            return jet;
                        };

                    std::array<Jet<Order>, 2> near;
                    for (int k = 0; k < 2; ++k) {
                        near[k].set_zero(ndof);
                        if (const auto* s = cache.vertex_at[ends[k]]) {
                            near[k] = scattered(
                                s->near.value, s->dict->vertex_ids(),
                                s->near.grad, s->near.hess);
                        }
                    }
                    // A = sum_j m_j (N(q_j) - R_j), B = sum_j m_j.
                    Jet<Order> A(ndof), B(ndof);
                    for (size_t k = begin; k < end; ++k) {
                        const auto& s = cache.edges[k];
                        const std::array<index_t, 4> prim = { {
                            ends[0],
                            ends[1],
                            mesh.edges()(s.eb, 0),
                            mesh.edges()(s.eb, 1),
                        } };
                        Jet<Order> residual = scattered(
                            s.potential.value, s.dict->vertex_ids(),
                            s.potential.grad, s.potential.hess);
                        const Jet<Order> u = geometry(s.u, prim);
                        Jet<Order> one_minus_u = scale(u, -1);
                        one_minus_u.value += 1;
                        residual += scale(product(one_minus_u, near[0]), -1);
                        residual += scale(product(u, near[1]), -1);
                        const Jet<Order> m = geometry(s.m, prim);
                        A += product(m, residual);
                        B += m;
                    }
                    Jet<Order> term(ndof);
                    for (const Pool& pool : out.pools) {
                        Jet<Order> denominator = scale(B, pool.scale);
                        denominator.value += 1;
                        term += scale(quotient(A, denominator), pool.coef);
                    }
                    out.add(
                        ids, term.value, term.grad, std::move(term.hess),
                        projection, nv);
                }
            });

        Result<Order> out(X.size());
        for (auto& local : storage) {
            if constexpr (Order == 0) {
                out.energy += local.energy;
            }
            if constexpr (Order == 1) {
                out.gradient += local.gradient;
            }
            if constexpr (Order == 2) {
                out.hessian.insert(
                    out.hessian.end(), local.hessian.begin(),
                    local.hessian.end());
            }
        }
        return out;
    }
} // namespace

double local_transfer_value(
    const ESPCollisions& c,
    const CollisionMesh& m,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const ESPParameters& p)
{
    return c.empty() ? 0.0 : evaluate<0>(c, m, X, p).energy;
}
Eigen::VectorXd local_transfer_gradient(
    const ESPCollisions& c,
    const CollisionMesh& m,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const ESPParameters& p)
{
    return c.empty() ? Eigen::VectorXd::Zero(X.size()).eval()
                     : evaluate<1>(c, m, X, p).gradient;
}
Eigen::SparseMatrix<double> local_transfer_hessian(
    const ESPCollisions& c,
    const CollisionMesh& m,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const ESPParameters& p,
    PSDProjectionMethod psd)
{
    if (c.empty()) {
        return Eigen::SparseMatrix<double>(X.size(), X.size());
    }
    auto r = evaluate<2>(c, m, X, p, psd);
    Eigen::SparseMatrix<double> H(X.size(), X.size());
    H.setFromTriplets(r.hessian.begin(), r.hessian.end());
    return H;
}
} // namespace ipc::detail
