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
#include <map>
#include <type_traits>
#include <unordered_map>

namespace ipc::detail {
namespace {
    // A scalar's value and derivatives on one element's union of real vertices.
    // Small-stencil AD handles geometry; these jets combine existing pointwise
    // derivatives, retaining cross-stencil terms in the quotient-rule Hessian.
    template <int Order> struct Jet {
        double value;
        Eigen::VectorXd grad;
        Eigen::MatrixXd hess;

        Jet(int ndof, double v = 0) : value(v)
        {
            if constexpr (Order >= 1)
                grad = Eigen::VectorXd::Zero(ndof);
            if constexpr (Order == 2)
                hess = Eigen::MatrixXd::Zero(ndof, ndof);
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
        std::vector<index_t> vertices, edges;
        double weight;
        double ee_scale = 1.0; // PROTOTYPE: area-consistent EE weight factor
    };
    // PROTOTYPE: reference length for area-consistent EE weights (0 = off).
    // Transfer requires the vertex rule (unit weight per vertex), so the
    // measure per unit of weight is area / 3 on faces and length / 2 on rods.
    double area_ref(const ESPParameters& p)
    {
        const auto& cfg = ee_mollifier_config();
        return cfg.area_ref_dbar ? p.dbar : cfg.area_ref;
    }

    Element element(const CollisionMesh& mesh, const ESPParameters& p, size_t k)
    {
        Element out;
        if (k < mesh.num_faces()) {
            for (int i = 0; i < 3; ++i) {
                out.vertices.push_back(mesh.faces()(k, i));
                out.edges.push_back(mesh.faces_to_edges()(k, i));
            }
            out.weight = p.area_weights ? mesh.face_areas()(k) : 1.0;
            out.ee_scale = ee_area_scale_factor(
                mesh.face_areas()(k) / 3.0, area_ref(p), /*face=*/true);
        } else if ((k -= mesh.num_faces()) < mesh.num_codim_edges()) {
            const auto e = mesh.codim_edges()[k];
            out.vertices = { mesh.edges()(e, 0), mesh.edges()(e, 1) };
            out.edges = { e };
            out.weight = p.area_weights ? mesh.edge_area(e) : 1.0;
            out.ee_scale = ee_area_scale_factor(
                mesh.edge_area(e) / 2.0, area_ref(p), /*face=*/false);
        } else {
            const auto v = mesh.codim_vertices()[k - mesh.num_codim_edges()];
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

    template <int O> struct ElementResult {
        std::vector<index_t> vertices;
        Jet<O> energy;
        ElementResult(std::vector<index_t> ids)
            : vertices(std::move(ids))
            , energy(3 * vertices.size())
        {
        }
    };

    template <int Order>
    ElementResult<Order> evaluate_element(
        const Element& elem,
        const ESPCollisions& collisions,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const ESPParameters& params)
    {
        using T = std::conditional_t<
            Order == 0, double,
            std::conditional_t<Order == 1, ADGrad<12>, ADHessian<12>>>;
        using Dict = ESPCollisionDict<PointType::EDGE>;
        struct Pair {
            index_t ea, eb;
            const Dict* dict;
        };
        std::vector<Pair> pairs;
        std::vector<index_t> ids = elem.vertices;
        const auto& maps = collisions.maps();
        for (const auto v : elem.vertices) {
            const auto it = maps.vertex_collisions.find(v);
            if (it != maps.vertex_collisions.end()) {
                const auto& more = it->second->vertex_ids();
                ids.insert(ids.end(), more.begin(), more.end());
            }
        }
        for (const auto ea : elem.edges) {
            for (const auto eb : collisions.m_candidates.ee_set(mesh, ea)) {
                const auto it = maps.edge_edge_collisions.find({ ea, eb });
                if (it == maps.edge_edge_collisions.end())
                    continue;
                pairs.push_back({ ea, eb, it->second.get() });
                const auto& more = it->second->vertex_ids();
                ids.insert(ids.end(), more.begin(), more.end());
                for (const auto v : it->second->primary_vertex_ids())
                    if (v >= 0)
                        ids.push_back(v);
            }
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        ElementResult<Order> result(std::move(ids));
        const int ndof = 3 * result.vertices.size();
        std::unordered_map<index_t, int> local;
        for (int i = 0; i < int(result.vertices.size()); ++i)
            local[result.vertices[i]] = i;

        auto scatter = [&](Jet<Order>& out, const auto& vids,
                           const Eigen::VectorXd& g, const Eigen::MatrixXd& h) {
            if constexpr (Order >= 1) {
                for (int i = 0; i < int(vids.size()); ++i)
                    out.grad.template segment<3>(3 * local.at(vids[i])) +=
                        g.segment<3>(3 * i);
            }
            if constexpr (Order == 2) {
                for (int i = 0; i < int(vids.size()); ++i)
                    for (int j = 0; j < int(vids.size()); ++j)
                        out.hess.template block<3, 3>(
                            3 * local.at(vids[i]), 3 * local.at(vids[j])) +=
                            h.block<3, 3>(3 * i, 3 * j);
            }
        };
        auto geometry_jet = [&](const T& x,
                                const std::array<index_t, 4>& prim) {
            Jet<Order> out(ndof, value(x));
            if constexpr (Order >= 1) {
                for (int i = 0; i < 4; ++i)
                    out.grad.template segment<3>(3 * local.at(prim[i])) +=
                        x.grad.template segment<3>(3 * i);
            }
            if constexpr (Order == 2) {
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                        out.hess.template block<3, 3>(
                            3 * local.at(prim[i]), 3 * local.at(prim[j])) +=
                            x.Hess.template block<3, 3>(3 * i, 3 * j);
            }
            return out;
        };

        NearFarBarrier<> nf(params.barrier, params.dbar_factor());
        Jet<Order> base(ndof), numerator(ndof), denominator(ndof, 1.0);
        std::map<index_t, Jet<Order>> near_vertices;
        for (const auto v : elem.vertices) {
            Jet<Order> near(ndof), far(ndof);
            const auto it = maps.vertex_collisions.find(v);
            if (it != maps.vertex_collisions.end()) {
                const auto& dict = *it->second;
                std::tie(near.value, far.value) = PointPotentialHelper::
                    evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                        X, dict, params, nf);
                if constexpr (Order >= 1) {
                    auto [gn, gf] = PointPotentialHelper::
                        evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                            X, dict, params, nf);
                    Eigen::MatrixXd hn, hf;
                    if constexpr (Order == 2)
                        std::tie(hn, hf) = PointPotentialHelper::
                            evaluate_potential_hessian_at_vertex_with_cached_collisions_nearfar(
                                X, dict, params, PSDProjectionMethod::NONE, nf);
                    scatter(near, dict.vertex_ids(), gn, hn);
                    scatter(far, dict.vertex_ids(), gf, hf);
                }
            }
            base += near;
            base += far;
            near_vertices.emplace(v, std::move(near));
        }
        int power = 1;
        if (const auto* b = dynamic_cast<const InversePowerBarrier<>*>(
                params.barrier.get()))
            power = std::max(1, int(std::lround(b->power())) + 1);
        if (ee_mollifier_config().order_override > 0) // PROTOTYPE
            power = ee_mollifier_config().order_override;

        for (const auto& pair : pairs) {
            const std::array<index_t, 4> prim = { {
                mesh.edges()(pair.ea, 0),
                mesh.edges()(pair.ea, 1),
                mesh.edges()(pair.eb, 0),
                mesh.edges()(pair.eb, 1),
            } };
            Eigen::Matrix<T, 4, 3> positions;
            if constexpr (Order == 0) {
                for (int i = 0; i < 4; ++i)
                    positions.row(i) = X.row(prim[i]);
            } else {
                Eigen::Vector<double, 12> flat;
                for (int i = 0; i < 4; ++i)
                    flat.segment<3>(3 * i) = X.row(prim[i]);
                positions = slice_positions<T, 4, 3>(flat);
            }
            const T u = closest_point_uv<T>(
                positions.row(0), positions.row(1), positions.row(2),
                positions.row(3));
            const T d2 = line_line_sqr_distance<T>(
                positions.row(0), positions.row(1), positions.row(2),
                positions.row(3));
            const auto types = edge_edge_mollifier_type(
                X.row(prim[0]), X.row(prim[1]), X.row(prim[2]), X.row(prim[3]),
                value(d2));
            T m =
                Math<T>::cubic_spline(sqrt(d2) / params.ee_support()) * ee_mollifier_config().ee_weight /* PROTOTYPE */
                * edge_edge_mollifier<T>(
                    positions.row(0).transpose(), positions.row(1).transpose(),
                    positions.row(2).transpose(), positions.row(3).transpose(),
                    types, d2);
            const T initial = m;
            for (int k = 1; k < power; ++k)
                m = m * initial;
            m = m * elem.ee_scale; // PROTOTYPE
            const Eigen::Vector3<T> q =
                (positions.row(0) + u * (positions.row(1) - positions.row(0)))
                    .transpose();
            Eigen::RowVector3d qvalue;
            for (int i = 0; i < 3; ++i)
                qvalue(i) = value(q(i));
            const VertexMatrixView<3> extended(X, qvalue);
            const auto& dict = *pair.dict;
            Jet<Order> potential(
                ndof,
                PointPotentialHelper::
                    evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                        extended, dict, params, dict.ee_dtype(), nf));
            if constexpr (Order >= 1) {
                const auto g = PointPotentialHelper::
                    evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                        T>(extended, dict, params, q, nf);
                Eigen::MatrixXd h;
                if constexpr (Order == 2)
                    h = PointPotentialHelper::
                        evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions_near(
                            extended, dict, params, q, nf);
                scatter(potential, dict.vertex_ids(), g, h);
            }
            const auto uj = geometry_jet(u, prim);
            Jet<Order> one_minus_u = scale(uj, -1);
            one_minus_u.value += 1;
            auto reference = product(one_minus_u, near_vertices.at(prim[0]));
            reference += product(uj, near_vertices.at(prim[1]));
            potential += scale(reference, -1);
            const auto mj = geometry_jet(m, prim);
            numerator += product(mj, potential);
            denominator += mj;
        }
        // E = w/n [sum(P_near(vertex)+P_far(vertex))
        //           + sum(m_j*(P(q_j)-P_ref,j))/(1+sum(m_j))].
        // Total donated mass M/[n(1+M)] < 1/n, so all vertex
        // coefficients remain nonnegative, including with multiple EE samples.
        base += quotient(numerator, denominator);
        result.energy = scale(base, elem.weight / double(elem.vertices.size()));
        return result;
    }

    template <int Order> struct Result {
        double energy = 0;
        Eigen::VectorXd gradient;
        std::vector<Eigen::Triplet<double>> hessian;
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
        tbb::enumerable_thread_specific<Result<Order>> storage(
            [&] { return Result<Order>(X.size()); });
        const size_t count = mesh.num_faces() + mesh.num_codim_edges()
            + mesh.num_codim_vertices();
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, count), [&](const auto& range) {
                auto& local = storage.local();
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    auto r = evaluate_element<Order>(
                        element(mesh, params, i), collisions, mesh, X, params);
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
