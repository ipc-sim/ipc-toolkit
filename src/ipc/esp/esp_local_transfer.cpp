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
            if constexpr (Order >= 1)
                grad.setZero(ndof);
            if constexpr (Order == 2)
                hess.setZero(ndof, ndof);
        }
        Jet& operator+=(const Jet& b)
        {
            value += b.value;
            if constexpr (Order >= 1)
                grad += b.grad;
            if constexpr (Order == 2)
                hess += b.hess;
            return *this;
        }
    };

    template <int O> Jet<O> scale(Jet<O> a, double s)
    {
        a.value *= s;
        if constexpr (O >= 1)
            a.grad *= s;
        if constexpr (O == 2)
            a.hess *= s;
        return a;
    }
    template <int O> Jet<O> product(const Jet<O>& a, const Jet<O>& b)
    {
        Jet<O> out(a.grad.size(), a.value * b.value);
        if constexpr (O >= 1)
            out.grad = b.value * a.grad + a.value * b.grad;
        if constexpr (O == 2)
            out.hess = b.value * a.hess + a.value * b.hess
                + a.grad * b.grad.transpose() + b.grad * a.grad.transpose();
        return out;
    }
    template <int O> Jet<O> quotient(const Jet<O>& a, const Jet<O>& b)
    {
        Jet<O> inv(b.grad.size(), 1.0 / b.value);
        if constexpr (O >= 1)
            inv.grad = -b.grad / (b.value * b.value);
        if constexpr (O == 2)
            inv.hess = -b.hess / (b.value * b.value)
                + (2.0 / (b.value * b.value * b.value)) * b.grad
                    * b.grad.transpose();
        return product(a, inv);
    }

    struct Element {
        std::array<index_t, 3> vertices {}, edges {};
        int num_vertices = 0, num_edges = 0;
        std::array<double, 3> vertex_weights {};
        index_t face = -1;
        const FaceQuadRule* non_corner = nullptr;
        double ee_scale = 1.0; // Element-local prototype area factor.
        double weight;
    };
    // PROTOTYPE: reference length for area-consistent EE weights (0 = off).
    // An EE sample takes at most the smaller endpoint weight c of the
    // element's mass, so the measure it can stand for is c times the
    // element's measure (area / 3 on faces with the vertex rule, length / 2
    // on rods).
    double area_ref(const ESPParameters& p)
    {
        const auto& cfg = ee_mollifier_config();
        return cfg.area_ref_dbar ? p.dbar : cfg.area_ref;
    }

    Element element(const CollisionMesh& mesh, const ESPParameters& p, size_t k)
    {
        Element out;
        if (k < mesh.num_faces()) {
            const auto rule = p.face_rule(mesh.face_areas()(k));
            const auto& cw = rule.corner_weights;
            out.non_corner = rule.non_corner;
            out.face = static_cast<index_t>(k);
            out.vertex_weights = { cw[0], cw[1], cw[2] };
            out.num_vertices = out.num_edges = 3;
            for (int i = 0; i < 3; ++i) {
                out.vertices[i] = mesh.faces()(k, i);
                out.edges[i] = mesh.faces_to_edges()(k, i);
            }
            out.weight = p.area_weights ? mesh.face_areas()(k) : 1.0;
            const double c = std::min({ cw[0], cw[1], cw[2] });
            out.ee_scale = ee_area_scale_factor(
                c * mesh.face_areas()(k), area_ref(p), /*face=*/true);
        } else if ((k -= mesh.num_faces()) < mesh.num_codim_edges()) {
            const auto e = mesh.codim_edges()[k];
            out.vertex_weights = { 0.5, 0.5, 0.0 };
            out.num_vertices = 2;
            out.num_edges = 1;
            out.vertices = { mesh.edges()(e, 0), mesh.edges()(e, 1) };
            out.edges = { e };
            out.weight = p.area_weights ? mesh.edge_area(e) : 1.0;
            out.ee_scale = ee_area_scale_factor(
                0.5 * mesh.edge_area(e), area_ref(p), /*face=*/false);
        } else {
            const auto v = mesh.codim_vertices()[k - mesh.num_codim_edges()];
            out.vertex_weights = { 1.0, 0.0, 0.0 };
            out.num_vertices = 1;
            out.vertices = { v };
            out.weight = p.area_weights ? mesh.vertex_area(v) : 1.0;
        }
        return out;
    }

    template <typename T> double value(const T& x)
    {
        if constexpr (std::is_same_v<T, double>)
            return x;
        else
            return x.val;
    }

    template <int Order>
    using GeometryScalar = std::conditional_t<
        Order == 0,
        double,
        std::conditional_t<Order == 1, ADGrad<12>, ADHessian<12>>>;

    template <int Order> struct VertexSample {
        const ESPCollisionDict<PointType::VERTEX>* dict;
        Jet<Order> near, far;
        explicit VertexSample(const ESPCollisionDict<PointType::VERTEX>* d)
            : dict(d)
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
    // derivative order. Only dict-local jets are cached: face reference and
    // denominator terms, cross-stencil coupling, and PSD remain face-local.
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
    EvaluationCache<Order> build_cache(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params)
    {
        using T = GeometryScalar<Order>;
        EvaluationCache<Order> cache;
        const auto& maps = collisions.maps();
        cache.vertex_at.resize(mesh.num_vertices(), nullptr);
        cache.vertices.reserve(maps.vertex_collisions.size());
        for (const auto& [v, dict] : maps.vertex_collisions) {
            cache.vertices.emplace_back(dict.get());
            cache.vertex_at[v] = &cache.vertices.back();
        }

        struct Pair {
            index_t ea, eb;
            const ESPCollisionDict<PointType::EDGE>* dict;
        };
        std::vector<Pair> pairs;
        pairs.reserve(maps.edge_edge_collisions.size());
        for (const auto& [edges, dict] : maps.edge_edge_collisions)
            pairs.push_back({ edges.first, edges.second, dict.get() });
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
                    auto& sample = cache.vertices[i];
                    const auto& dict = *sample.dict;
                    std::tie(
                        sample.near.value,
                        sample.far.value) = PointPotentialHelper::
                        evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                            X, dict, params, nf);
                    if constexpr (Order >= 1)
                        std::tie(
                            sample.near.grad,
                            sample.far.grad) = PointPotentialHelper::
                            evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                                X, dict, params, nf);
                    if constexpr (Order == 2)
                        std::tie(
                            sample.near.hess,
                            sample.far.hess) = PointPotentialHelper::
                            evaluate_potential_hessian_at_vertex_with_cached_collisions_nearfar(
                                X, dict, params, PSDProjectionMethod::NONE, nf);
                }
            });
        int power = 1;
        if (const auto* b = dynamic_cast<const InversePowerBarrier<>*>(
                params.barrier.get()))
            power = std::max(1, int(std::lround(b->power())) + 1);
        if (ee_mollifier_config().order_override > 0)
            power = ee_mollifier_config().order_override;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, cache.edges.size()),
            [&](const auto& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    auto& sample = cache.edges[i];
                    const std::array<index_t, 4> prim = {
                        mesh.edges()(sample.ea, 0), mesh.edges()(sample.ea, 1),
                        mesh.edges()(sample.eb, 0), mesh.edges()(sample.eb, 1)
                    };
                    Eigen::Matrix<T, 4, 3> positions;
                    if constexpr (Order == 0) {
                        for (int j = 0; j < 4; ++j)
                            positions.row(j) = X.row(prim[j]);
                    } else {
                        Eigen::Vector<double, 12> flat;
                        for (int j = 0; j < 4; ++j)
                            flat.segment<3>(3 * j) = X.row(prim[j]);
                        positions = slice_positions<T, 4, 3>(flat);
                    }
                    sample.u = closest_point_uv<T>(
                        positions.row(0), positions.row(1), positions.row(2),
                        positions.row(3));
                    const T d2 = line_line_sqr_distance<T>(
                        positions.row(0), positions.row(1), positions.row(2),
                        positions.row(3));
                    const auto types = edge_edge_mollifier_type(
                        X.row(prim[0]), X.row(prim[1]), X.row(prim[2]),
                        X.row(prim[3]), value(d2));
                    sample.m =
                        Math<T>::cubic_spline(sqrt(d2) / params.ee_support())
                        * ee_mollifier_config().ee_weight
                        * edge_edge_mollifier<T>(
                            positions.row(0).transpose(),
                            positions.row(1).transpose(),
                            positions.row(2).transpose(),
                            positions.row(3).transpose(), types, d2);
                    const T initial = sample.m;
                    for (int j = 1; j < power; ++j)
                        sample.m = sample.m * initial;
                    const Eigen::Vector3<T> q =
                        (positions.row(0)
                         + sample.u * (positions.row(1) - positions.row(0)))
                            .transpose();
                    Eigen::RowVector3d qvalue;
                    for (int j = 0; j < 3; ++j)
                        qvalue(j) = value(q(j));
                    const VertexMatrixView<3> extended(X, qvalue);
                    const auto& dict = *sample.dict;
                    sample.potential.value = PointPotentialHelper::
                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                            extended, dict, params, dict.ee_dtype(), nf);
                    if constexpr (Order >= 1)
                        sample.potential.grad = PointPotentialHelper::
                            evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                                T>(extended, dict, params, q, nf);
                    if constexpr (Order == 2)
                        sample.potential.hess = PointPotentialHelper::
                            evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions_near(
                                extended, dict, params, q, nf);
                }
            });
        return cache;
    }

    template <int O> struct ElementResult {
        std::vector<index_t> vertices;
        Jet<O> energy;
        ElementResult(std::vector<index_t> ids = {})
            : vertices(std::move(ids))
            , energy(3 * vertices.size())
        {
        }
    };

    // Private to a worker and this evaluation. Resize/zero every active jet;
    // only vector/matrix storage survives between source elements. The result
    // is scattered before the next element reuses it, including after PSD.
    template <int Order> struct ElementScratch {
        ElementResult<Order> result;
        std::vector<const EdgeSample<Order>*> pairs;
        std::vector<int> scatter_ids;
        std::array<Jet<Order>, 3> near_vertices;
        Jet<Order> far, numerator, denominator, potential;
    };

    template <int Order>
    ElementResult<Order>& evaluate_element(
        const Element& elem,
        const CollisionMesh& mesh,
        const EvaluationCache<Order>& cache,
        ElementScratch<Order>& scratch,
        const ESPCollisions& collisions,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params)
    {
        using T = GeometryScalar<Order>;
        auto& pairs = scratch.pairs;
        pairs.clear();
        auto& result = scratch.result;
        auto& ids = result.vertices;
        ids.assign(
            elem.vertices.begin(), elem.vertices.begin() + elem.num_vertices);
        for (int j = 0; j < elem.num_vertices; ++j) {
            const auto v = elem.vertices[j];
            if (const auto* sample = cache.vertex_at[v]) {
                const auto& more = sample->dict->vertex_ids();
                ids.insert(ids.end(), more.begin(), more.end());
            }
        }
        const auto& maps = collisions.maps();
        // The face rule's non-corner points (corner points are the vertex
        // samples above).
        const auto face_it = elem.face >= 0
            ? maps.face_collisions.find(elem.face)
            : maps.face_collisions.end();
        if (face_it != maps.face_collisions.end()) {
            for (const auto& dict : face_it->second) {
                const auto& more = dict->vertex_ids();
                ids.insert(ids.end(), more.begin(), more.end());
            }
        }
        for (int j = 0; j < elem.num_edges; ++j) {
            const auto ea = elem.edges[j];
            for (size_t k = cache.edge_offsets[ea];
                 k < cache.edge_offsets[ea + 1]; ++k) {
                const auto& sample = cache.edges[k];
                pairs.push_back(&sample);
                const auto& more = sample.dict->vertex_ids();
                ids.insert(ids.end(), more.begin(), more.end());
                for (const auto v : sample.dict->primary_vertex_ids())
                    if (v >= 0)
                        ids.push_back(v);
            }
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        const int ndof = 3 * result.vertices.size();
        const auto local_index = [&](index_t v) {
            const auto it = std::lower_bound(ids.begin(), ids.end(), v);
            assert(it != ids.end() && *it == v);
            return int(it - ids.begin());
        };

        auto scatter = [&](Jet<Order>& out, const auto& vids,
                           const Eigen::VectorXd& g, const Eigen::MatrixXd& h) {
            if constexpr (Order >= 1) {
                auto& indices = scratch.scatter_ids;
                indices.resize(vids.size());
                for (int i = 0; i < int(vids.size()); ++i) {
                    indices[i] = local_index(vids[i]);
                    out.grad.template segment<3>(3 * indices[i]) +=
                        g.segment<3>(3 * i);
                }
                if constexpr (Order == 2) {
                    for (int i = 0; i < int(vids.size()); ++i)
                        for (int j = 0; j < int(vids.size()); ++j)
                            out.hess.template block<3, 3>(
                                3 * indices[i], 3 * indices[j]) +=
                                h.block<3, 3>(3 * i, 3 * j);
                }
            }
        };
        auto geometry_jet = [&](const T& x,
                                const std::array<index_t, 4>& prim) {
            Jet<Order> out(ndof, value(x));
            if constexpr (Order >= 1) {
                std::array<int, 4> indices;
                for (int i = 0; i < 4; ++i)
                    indices[i] = local_index(prim[i]);
                for (int i = 0; i < 4; ++i)
                    out.grad.template segment<3>(3 * indices[i]) +=
                        x.grad.template segment<3>(3 * i);
                if constexpr (Order == 2) {
                    for (int i = 0; i < 4; ++i)
                        for (int j = 0; j < 4; ++j)
                            out.hess.template block<3, 3>(
                                3 * indices[i], 3 * indices[j]) +=
                                x.Hess.template block<3, 3>(3 * i, 3 * j);
                }
            }
            return out;
        };

        auto& base = result.energy;
        auto& numerator = scratch.numerator;
        auto& denominator = scratch.denominator;
        base.set_zero(ndof);
        numerator.set_zero(ndof);
        denominator.set_zero(ndof, 1.0);
        for (int j = 0; j < elem.num_vertices; ++j) {
            const auto v = elem.vertices[j];
            auto& near = scratch.near_vertices[j];
            auto& far = scratch.far;
            near.set_zero(ndof);
            far.set_zero(ndof);
            if (const auto* sample = cache.vertex_at[v]) {
                near.value = sample->near.value;
                far.value = sample->far.value;
                scatter(
                    near, sample->dict->vertex_ids(), sample->near.grad,
                    sample->near.hess);
                scatter(
                    far, sample->dict->vertex_ids(), sample->far.grad,
                    sample->far.hess);
            }
            base += scale(near, elem.vertex_weights[j]);
            base += scale(far, elem.vertex_weights[j]);
        }
        const NearFarBarrier<> nf(params.barrier, params.dbar_factor());
        if (face_it != maps.face_collisions.end()) {
            const auto& rule = *elem.non_corner;
            for (size_t qi = 0; qi < rule.size(); ++qi) {
                const auto& qp = rule[qi];
                const auto& dict = *face_it->second[qi];
                if (dict.size() == 0)
                    continue;
                Eigen::RowVector3d q_pos = Eigen::RowVector3d::Zero();
                for (int i = 0; i < 3; ++i)
                    q_pos += qp.lambda[i] * X.row(elem.vertices[i]);
                const VertexMatrixView<3> X_qp(X, q_pos);
                Jet<Order> near(ndof), far(ndof);
                std::tie(near.value, far.value) = PointPotentialHelper::
                    evaluate_potential_at_face_interior_point_with_cached_collisions_nearfar(
                        X_qp, dict, params, nf);
                if constexpr (Order >= 1) {
                    auto [gn, gf] = PointPotentialHelper::
                        evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
                            X_qp, dict, params, qp.lambda, nf);
                    Eigen::MatrixXd hn, hf;
                    if constexpr (Order == 2)
                        std::tie(hn, hf) = PointPotentialHelper::
                            evaluate_potential_hessian_at_face_interior_point_with_cached_collisions_nearfar(
                                X_qp, dict, params, qp.lambda,
                                PSDProjectionMethod::NONE, nf);
                    scatter(near, dict.vertex_ids(), gn, hn);
                    scatter(far, dict.vertex_ids(), gf, hf);
                }
                base += scale(near, qp.weight);
                base += scale(far, qp.weight);
            }
        }
        const auto near_vertex = [&](index_t v) -> const Jet<Order>& {
            const auto end = elem.vertices.begin() + elem.num_vertices;
            const auto it = std::find(elem.vertices.begin(), end, v);
            assert(it != end);
            return scratch.near_vertices[it - elem.vertices.begin()];
        };
        for (const auto* pair : pairs) {
            const std::array<index_t, 4> prim = { mesh.edges()(pair->ea, 0),
                                                  mesh.edges()(pair->ea, 1),
                                                  mesh.edges()(pair->eb, 0),
                                                  mesh.edges()(pair->eb, 1) };
            auto& potential = scratch.potential;
            potential.set_zero(ndof, pair->potential.value);
            scatter(
                potential, pair->dict->vertex_ids(), pair->potential.grad,
                pair->potential.hess);
            const auto uj = geometry_jet(pair->u, prim);
            Jet<Order> one_minus_u = scale(uj, -1);
            one_minus_u.value += 1;
            auto reference = product(one_minus_u, near_vertex(prim[0]));
            reference += product(uj, near_vertex(prim[1]));
            potential += scale(reference, -1);
            // Cache only edge geometry common to all incident elements.
            // Area scaling and endpoint quadrature weights remain
            // element-local.
            const auto mj = scale(geometry_jet(pair->m, prim), elem.ee_scale);
            const auto weight_at = [&](index_t v) {
                const auto end = elem.vertices.begin() + elem.num_vertices;
                const auto it = std::find(elem.vertices.begin(), end, v);
                assert(it != end);
                return elem.vertex_weights[it - elem.vertices.begin()];
            };
            const double c = std::min(weight_at(prim[0]), weight_at(prim[1]));
            numerator += scale(product(mj, potential), c);
            denominator += mj;
        }
        // E = w [sum_i w_i (P_near + P_far)(vertex i)
        //        + sum_q w_q (P_near + P_far)(non-corner point q)
        //        + sum_j c_j m_j (P(q_j) - P_ref,j) / (1 + sum_j m_j)].
        // Each vertex donates at most c M / (1 + M) < c <= w_i, so all
        // coefficients stay nonnegative, including with multiple EE samples.
        // With the vertex rule (w_i = c = 1/3) this is the original
        // w/3 [...] formula; rods use w_i = c = 1/2.
        base += quotient(numerator, denominator);
        const double weight = elem.weight;
        base.value *= weight;
        if constexpr (Order >= 1)
            base.grad *= weight;
        if constexpr (Order == 2)
            base.hess *= weight;
        return result;
    }

    template <int Order> struct Result {
        double energy = 0;
        Eigen::VectorXd gradient;
        std::vector<Eigen::Triplet<double>> hessian;
        ElementScratch<Order> scratch;
        explicit Result(int n)
        {
            if constexpr (Order == 1)
                gradient = Eigen::VectorXd::Zero(n);
        }
    };
    template <int Order>
    Result<Order> evaluate(
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params,
        PSDProjectionMethod projection = PSDProjectionMethod::NONE)
    {
        const auto cache = build_cache<Order>(collisions, mesh, X, params);
        tbb::enumerable_thread_specific<Result<Order>> storage(
            [&] { return Result<Order>(X.size()); });
        // Contact is usually localized on a much larger surface. An element
        // with no vertex, face-interior, or source-edge EE dict contributes
        // identically zero at every derivative order. Reject it before building
        // vertex unions, allocating jets, projecting a zero Hessian, or
        // emitting zero triplets. Rebuild these masks on every evaluation: the
        // collision set can change between line-search trials.
        const auto& maps = collisions.maps();
        std::vector<unsigned char> vertex_active(mesh.num_vertices(), 0);
        std::vector<unsigned char> edge_active(mesh.num_edges(), 0);
        std::vector<unsigned char> face_active(mesh.num_faces(), 0);
        for (const auto& [v, dict] : maps.vertex_collisions)
            vertex_active[v] = 1;
        for (const auto& [edges, dict] : maps.edge_edge_collisions)
            edge_active[edges.first] = 1;
        for (const auto& [face, dicts] : maps.face_collisions)
            face_active[face] = 1;

        const auto is_active = [&](size_t k) {
            if (k < mesh.num_faces()) {
                if (face_active[k])
                    return true;
                for (int j = 0; j < 3; ++j)
                    if (vertex_active[mesh.faces()(k, j)]
                        || edge_active[mesh.faces_to_edges()(k, j)])
                        return true;
                return false;
            }
            k -= mesh.num_faces();
            if (k < mesh.num_codim_edges()) {
                const auto e = mesh.codim_edges()[k];
                return bool(
                    edge_active[e] || vertex_active[mesh.edges()(e, 0)]
                    || vertex_active[mesh.edges()(e, 1)]);
            }
            return bool(
                vertex_active
                    [mesh.codim_vertices()[k - mesh.num_codim_edges()]]);
        };
        const size_t count = mesh.num_faces() + mesh.num_codim_edges()
            + mesh.num_codim_vertices();
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, count), [&](const auto& range) {
                auto& local = storage.local();
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    if (!is_active(i))
                        continue;
                    auto& r = evaluate_element<Order>(
                        element(mesh, params, i), mesh, cache, local.scratch,
                        collisions, X, params);
                    if constexpr (Order == 0)
                        local.energy += r.energy.value;
                    if constexpr (Order == 1)
                        local_gradient_to_global_gradient(
                            r.energy.grad, r.vertices, 3, local.gradient);
                    if constexpr (Order == 2) {
                        if (projection != PSDProjectionMethod::NONE)
                            r.energy.hess =
                                project_to_psd(r.energy.hess, projection);
                        local_hessian_to_global_triplets(
                            r.energy.hess, r.vertices, 3, local.hessian,
                            mesh.num_vertices());
                    }
                }
            });
        Result<Order> out(X.size());
        for (auto& local : storage) {
            if constexpr (Order == 0)
                out.energy += local.energy;
            if constexpr (Order == 1)
                out.gradient += local.gradient;
            if constexpr (Order == 2)
                out.hessian.insert(
                    out.hessian.end(), local.hessian.begin(),
                    local.hessian.end());
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
    if (c.empty())
        return Eigen::SparseMatrix<double>(X.size(), X.size());
    auto r = evaluate<2>(c, m, X, p, psd);
    Eigen::SparseMatrix<double> H(X.size(), X.size());
    H.setFromTriplets(r.hessian.begin(), r.hessian.end());
    return H;
}
} // namespace ipc::detail
