#include "esp_potential.hpp"

#include "esp_local_transfer.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/esp/collisions/esp_quadrature.hpp>
#include <ipc/esp/collisions/vertex_matrix_view.hpp>
#include <ipc/esp/esp_collision_maps.hpp>
#include <ipc/esp/esp_distance.hpp>
#include <ipc/esp/quadrature_potential.hpp>
#include <ipc/gcp/distance/mollifier.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>

namespace ipc {

constexpr double FACE_QUADRATURE_WEIGHT_SCALE = 1.0;

namespace {
    template <class T> inline T pow_int(T x, int n)
    {
        T r = T(1);
        for (int i = 0; i < n; ++i) {
            r = r * x;
        }
        return r;
    }

    // An element we integrate the 3D potential over: a face, an edge in no
    // face (e.g., a rod), or a vertex in no edge (e.g., a particle). Each
    // element's potential is the weighted average of the pointwise potential
    // at its samples: its vertices, its face rule's other points (faces only),
    // and the closest points of its edges' edge-edge pairs. A face's vertices
    // carry the rule's corner weights; a rod's two vertices 1/2 each; a
    // particle's vertex 1. So every element's fixed samples weigh 1 in total.
    //
    // A rod or particle has no face, so without elements of its own it would
    // get exactly zero potential, however close it comes to anything.
    struct IntegrationElement {
        index_t face = -1; ///< Face id, or -1 for an edge or vertex element
        std::array<index_t, 3> edges {}; ///< Edges whose edge-edge pairs count
        int num_edges = 0;
        std::array<index_t, 3> vertices {}; ///< Vertex samples
        std::array<double, 3> vertex_weights {}; ///< Their weights
        int num_vertices = 0;
        /// The face rule's other points (faces only)
        const FaceQuadRule* non_corner = nullptr;
        double weight = 1; ///< Outer weight of the element's average
    };

    const FaceQuadRule NO_FACE_QUADRATURE;

    // We number the elements faces first, then edges in no face, then
    // vertices in no edge.
    size_t num_integration_elements(const CollisionMesh& mesh)
    {
        return mesh.num_faces() + mesh.num_codim_edges()
            + mesh.num_codim_vertices();
    }

    IntegrationElement integration_element(
        const CollisionMesh& mesh, const ESPParameters& params, size_t k)
    {
        // With area weights, an element weighs its rest measure: a face its
        // area, an edge in no face its edge area, and an isolated vertex its
        // vertex area (one by default), each times
        // params.element_weight_scale(). These have different units, so on a
        // mesh that mixes faces with rods or particles the barrier stiffness
        // has to be chosen with that in mind.
        IntegrationElement elem;
        if (k < mesh.num_faces()) {
            const auto f = static_cast<index_t>(k);
            elem.face = f;
            elem.num_edges = 3;
            for (int j = 0; j < 3; j++) {
                elem.edges[j] = mesh.faces_to_edges()(f, j);
            }
            // The rule's corner points are the face's vertex samples.
            const auto rule = params.face_rule(mesh.face_areas()(f));
            elem.num_vertices = 3;
            for (int j = 0; j < 3; j++) {
                elem.vertices[j] = mesh.faces()(f, j);
                elem.vertex_weights[j] = rule.corner_weights[j];
            }
            elem.non_corner = rule.non_corner;
            elem.weight = params.area_weights
                ? mesh.face_areas()(f) * params.element_weight_scale(3)
                : 1.;
            return elem;
        }
        k -= mesh.num_faces();

        if (k < mesh.num_codim_edges()) {
            const index_t e = mesh.codim_edges()[k];
            elem.num_edges = 1;
            elem.edges[0] = e;
            elem.num_vertices = 2;
            elem.vertices[0] = mesh.edges()(e, 0);
            elem.vertices[1] = mesh.edges()(e, 1);
            elem.vertex_weights = { { 0.5, 0.5, 0 } };
            elem.weight = params.area_weights
                ? mesh.edge_area(e) * params.element_weight_scale(2)
                : 1.;
            return elem;
        }
        k -= mesh.num_codim_edges();

        assert(k < mesh.num_codim_vertices());
        const index_t v = mesh.codim_vertices()[k];
        elem.num_vertices = 1;
        elem.vertices[0] = v;
        elem.vertex_weights = { { 1, 0, 0 } };
        elem.weight = params.area_weights ? mesh.vertex_area(v) : 1.;
        return elem;
    }

    // Make an EE sample stand for about one r x r patch (length r on rods)
    // instead of a fixed fraction of its element (params.ee_area_scaling).
    // An EE sample of weight m takes the share m / (1 + m) of the element's
    // average, since the element's fixed samples weigh 1 in total, so the
    // scale compares the element's rest measure with the patch. Rest
    // geometry only, so a per-element constant.
    double ee_area_scale(
        const CollisionMesh& mesh,
        const IntegrationElement& elem,
        const ESPParameters& params)
    {
        if (!params.ee_area_scaling) {
            return 1.0;
        }
        const double r = params.ee_area_reference();
        if (elem.face >= 0) {
            return ee_area_scale_factor(
                mesh.face_areas()(elem.face), r, /*face=*/true);
        }
        if (elem.num_edges == 1) {
            return ee_area_scale_factor(
                mesh.edge_area(elem.edges[0]), r, /*face=*/false);
        }
        return 1.0;
    }

    // ---- Lagged normalization (see ESPPotential::set_lagged_weights) ----

    // Edge-edge mollifier weight of the sample on edge e0 closest to e1,
    // without the element's area scale: (cubic spline * ee_weight *
    // mollifier)^power, as a scalar T (double or an AD type over the four
    // edge endpoints).
    template <typename T>
    T ee_sample_weight(
        const Eigen::Matrix<T, 4, 3>& positionsT,
        const Eigen::Ref<const Eigen::Matrix<double, 4, 3>>& positions,
        const ESPParameters& params)
    {
        using std::sqrt;
        const T dist_sqr = line_line_sqr_distance<T>(
            positionsT.row(0), positionsT.row(1), positionsT.row(2),
            positionsT.row(3));
        const T dist = sqrt(dist_sqr);
        double dist_sqr_val;
        if constexpr (std::is_same_v<T, double>) {
            dist_sqr_val = dist_sqr;
        } else {
            dist_sqr_val = dist_sqr.val;
        }
        const auto mtypes = edge_edge_mollifier_type(
            positions.row(0).transpose(), positions.row(1).transpose(),
            positions.row(2).transpose(), positions.row(3).transpose(),
            dist_sqr_val, params.ee_mollifier);
        T m = Math<T>::cubic_spline(dist / params.ee_support())
            * params.ee_weight;
        m *= edge_edge_mollifier<T>(
            positionsT.row(0).transpose(), positionsT.row(1).transpose(),
            positionsT.row(2).transpose(), positionsT.row(3).transpose(),
            mtypes, dist_sqr, params.ee_mollifier);
        return pow_int(m, params.ee_mollifier_exponent());
    }

    // Per-sample coefficients of the lagged energy
    //   E = sum_e w_e [ (sum_j m_j P_near(q_j) + sum_s w_s P_near(s)) / D_e
    //                 + (sum_s w_s P_far(s)) / W_e ],
    // where s runs over the element's vertex and face samples (weights w_s,
    // summing to W_e), j over its edge-edge samples, and D_e = W_e + sum_j
    // m_j is held fixed. Each sample then contributes coefficient * its own
    // potential, so it is evaluated once instead of once per element.
    struct LaggedCoefficients {
        Eigen::VectorXd vertex_near, vertex_far; // per vertex
        Eigen::VectorXd face_near, face_far;     // per face, times qp weight
        Eigen::VectorXd edge; // per edge: sum_e w_e * area scale / D_e
    };

    LaggedCoefficients lagged_coefficients(
        const CollisionMesh& mesh,
        const ESPParameters& params,
        const std::vector<double>& lagged_weights)
    {
        LaggedCoefficients c;
        c.vertex_near.setZero(mesh.num_vertices());
        c.vertex_far.setZero(mesh.num_vertices());
        c.face_near.setZero(mesh.num_faces());
        c.face_far.setZero(mesh.num_faces());
        c.edge.setZero(mesh.num_edges());
        for (size_t k = 0; k < num_integration_elements(mesh); k++) {
            const IntegrationElement elem =
                integration_element(mesh, params, k);
            double fixed_weight = 0; // W_e
            for (int lv = 0; lv < elem.num_vertices; lv++) {
                fixed_weight += elem.vertex_weights[lv];
            }
            if (elem.face >= 0) {
                for (const auto& qp : *elem.non_corner) {
                    fixed_weight += FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                }
            }
            const double near_scale = elem.weight / lagged_weights[k];
            const double far_scale = elem.weight / fixed_weight;
            for (int lv = 0; lv < elem.num_vertices; lv++) {
                const double vw = elem.vertex_weights[lv];
                c.vertex_near[elem.vertices[lv]] += vw * near_scale;
                c.vertex_far[elem.vertices[lv]] += vw * far_scale;
            }
            if (elem.face >= 0) {
                c.face_near[elem.face] += near_scale;
                c.face_far[elem.face] += far_scale;
            }
            const double ee_scale = ee_area_scale(mesh, elem, params);
            for (int le = 0; le < elem.num_edges; le++) {
                c.edge[elem.edges[le]] += near_scale * ee_scale;
            }
        }
        return c;
    }

    template <typename Map>
    std::vector<typename Map::const_pointer> map_entries(const Map& map)
    {
        std::vector<typename Map::const_pointer> entries;
        entries.reserve(map.size());
        for (const auto& entry : map) {
            entries.push_back(&entry);
        }
        return entries;
    }

    Eigen::Matrix<double, 4, 3> ee_positions(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> X,
        const std::pair<index_t, index_t>& pair)
    {
        Eigen::Matrix<double, 4, 3> p;
        p.row(0) = X.row(mesh.edges()(pair.first, 0));
        p.row(1) = X.row(mesh.edges()(pair.first, 1));
        p.row(2) = X.row(mesh.edges()(pair.second, 0));
        p.row(3) = X.row(mesh.edges()(pair.second, 1));
        return p;
    }
} // namespace

std::vector<double> ESPPotential::element_near_weights(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    // D_e = W_e + sum_j m_j at X: the element's fixed sample weights plus
    // its edge-edge samples' mollifier weights (with the area scale).
    const size_t n = num_integration_elements(mesh);
    std::vector<double> weights(n, 0.0);
    std::vector<std::vector<index_t>> edge_elements(mesh.num_edges());
    std::vector<double> elem_ee_scale(n, 1.0);
    for (size_t k = 0; k < n; k++) {
        const IntegrationElement elem = integration_element(mesh, m_params, k);
        for (int lv = 0; lv < elem.num_vertices; lv++) {
            weights[k] += elem.vertex_weights[lv];
        }
        if (elem.face >= 0) {
            for (const auto& qp : *elem.non_corner) {
                weights[k] += FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
            }
        }
        for (int le = 0; le < elem.num_edges; le++) {
            edge_elements[elem.edges[le]].push_back(static_cast<index_t>(k));
        }
        elem_ee_scale[k] = ee_area_scale(mesh, elem, m_params);
    }
    if (mesh.dim() == 3 && m_params.dbar_factor() > 0) {
        for (const auto& [pair, dict] :
             collisions.maps().edge_edge_collisions) {
            const Eigen::Matrix<double, 4, 3> p = ee_positions(mesh, X, pair);
            const double m = ee_sample_weight<double>(p, p, m_params);
            for (const index_t k : edge_elements[pair.first]) {
                weights[k] += m * elem_ee_scale[k];
            }
        }
    }
    return weights;
}

void ESPPotential::set_lagged_weights(std::vector<double> weights)
{
    for (const double w : weights) {
        if (!(w > 0) || !std::isfinite(w)) {
            throw std::invalid_argument(
                "Lagged ESP element weights must be positive and finite.");
        }
    }
    m_lagged_weights = std::move(weights);
}

bool ESPPotential::uses_lagged_weights(const CollisionMesh& mesh) const
{
    if (m_lagged_weights.empty() || mesh.dim() != 3 || !m_use_near_far
        || m_params.dbar_factor() <= 0) {
        return false;
    }
    if (m_lagged_weights.size() != num_integration_elements(mesh)) {
        throw std::invalid_argument(
            "Lagged ESP element weights do not match the mesh's elements.");
    }
    if (m_params.local_quadrature_transfer) {
        throw std::invalid_argument(
            "Lagged ESP normalization does not support local quadrature "
            "transfer.");
    }
    return true;
}

double ESPPotential::lagged_value(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    const LaggedCoefficients c =
        lagged_coefficients(mesh, m_params, m_lagged_weights);
    const NearFarBarrier<> nf(m_params.barrier, m_params.dbar_factor());
    const auto& maps = collisions.maps();
    tbb::combinable<double> total([] { return 0.0; });

    const auto vertices = map_entries(maps.vertex_collisions);
    tbb::parallel_for(size_t(0), vertices.size(), [&](size_t i) {
        const index_t v = vertices[i]->first;
        const auto [pn, pf] = PointPotentialHelper::
            evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                X, *vertices[i]->second, m_params, nf);
        total.local() += c.vertex_near[v] * pn + c.vertex_far[v] * pf;
    });

    const auto faces = map_entries(maps.face_collisions);
    tbb::parallel_for(size_t(0), faces.size(), [&](size_t i) {
        const index_t f = faces[i]->first;
        const auto& rule = *m_params.face_rule(mesh.face_areas()(f)).non_corner;
        const auto& dicts = faces[i]->second;
        for (size_t qi = 0; qi < std::min(rule.size(), dicts.size()); qi++) {
            const auto& qp = rule[qi];
            const Eigen::RowVector3d q =
                qp.lambda[0] * X.row(mesh.faces()(f, 0))
                + qp.lambda[1] * X.row(mesh.faces()(f, 1))
                + qp.lambda[2] * X.row(mesh.faces()(f, 2));
            const auto [pn, pf] = PointPotentialHelper::
                evaluate_potential_at_face_interior_point_with_cached_collisions_nearfar(
                    VertexMatrixView<3>(X, q), *dicts[qi], m_params, nf);
            const double s = FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
            total.local() += s * (c.face_near[f] * pn + c.face_far[f] * pf);
        }
    });

    const auto pairs = map_entries(maps.edge_edge_collisions);
    tbb::parallel_for(size_t(0), pairs.size(), [&](size_t i) {
        const auto& [pair, dict] = *pairs[i];
        const double coef = c.edge[pair.first];
        if (coef == 0) {
            return;
        }
        const Eigen::Matrix<double, 4, 3> p = ee_positions(mesh, X, pair);
        const double m = ee_sample_weight<double>(p, p, m_params);
        if (m == 0) {
            return;
        }
        const double u =
            closest_point_uv<double>(p.row(0), p.row(1), p.row(2), p.row(3));
        const Eigen::RowVector3d q = u * (p.row(1) - p.row(0)) + p.row(0);
        const double pn = PointPotentialHelper::
            evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                VertexMatrixView<3>(X, q), *dict, m_params, dict->ee_dtype(),
                nf);
        total.local() += coef * m * pn;
    });
    return total.combine(std::plus<double>());
}

Eigen::VectorXd ESPPotential::lagged_gradient(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    const LaggedCoefficients c =
        lagged_coefficients(mesh, m_params, m_lagged_weights);
    const NearFarBarrier<> nf(m_params.barrier, m_params.dbar_factor());
    const auto& maps = collisions.maps();
    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(X.size()));

    const auto vertices = map_entries(maps.vertex_collisions);
    tbb::parallel_for(size_t(0), vertices.size(), [&](size_t i) {
        const index_t v = vertices[i]->first;
        const auto& dict = *vertices[i]->second;
        const auto [gn, gf] = PointPotentialHelper::
            evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                X, dict, m_params, nf);
        storage.local()(dict.dofs()) +=
            c.vertex_near[v] * gn + c.vertex_far[v] * gf;
    });

    const auto faces = map_entries(maps.face_collisions);
    tbb::parallel_for(size_t(0), faces.size(), [&](size_t i) {
        const index_t f = faces[i]->first;
        const auto& rule = *m_params.face_rule(mesh.face_areas()(f)).non_corner;
        const auto& dicts = faces[i]->second;
        for (size_t qi = 0; qi < std::min(rule.size(), dicts.size()); qi++) {
            const auto& qp = rule[qi];
            const Eigen::RowVector3d q =
                qp.lambda[0] * X.row(mesh.faces()(f, 0))
                + qp.lambda[1] * X.row(mesh.faces()(f, 1))
                + qp.lambda[2] * X.row(mesh.faces()(f, 2));
            const auto [gn, gf] = PointPotentialHelper::
                evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
                    VertexMatrixView<3>(X, q), *dicts[qi], m_params, qp.lambda,
                    nf);
            const double s = FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
            storage.local()(dicts[qi]->dofs()) +=
                s * (c.face_near[f] * gn + c.face_far[f] * gf);
        }
    });

    const auto pairs = map_entries(maps.edge_edge_collisions);
    tbb::parallel_for(size_t(0), pairs.size(), [&](size_t i) {
        using T = ADGrad<12>;
        const auto& [pair, dict] = *pairs[i];
        const double coef = c.edge[pair.first];
        if (coef == 0) {
            return;
        }
        const Eigen::Matrix<double, 4, 3> p = ee_positions(mesh, X, pair);
        Eigen::Vector<double, 12> flat;
        for (int r = 0; r < 4; r++) {
            flat.segment<3>(3 * r) = p.row(r).transpose();
        }
        const Eigen::Matrix<T, 4, 3> pT = slice_positions<T, 4, 3>(flat);
        const T m = ee_sample_weight<T>(pT, p, m_params);
        if (m.val == 0) {
            return;
        }
        const T u =
            closest_point_uv<T>(pT.row(0), pT.row(1), pT.row(2), pT.row(3));
        const Eigen::RowVector3<T> qT = u * (pT.row(1) - pT.row(0)) + pT.row(0);
        const Eigen::RowVector3d q(qT(0).val, qT(1).val, qT(2).val);
        const VertexMatrixView<3> Xq(X, q);
        const double pn = PointPotentialHelper::
            evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                Xq, *dict, m_params, dict->ee_dtype(), nf);
        const Eigen::VectorXd gp = PointPotentialHelper::
            evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                T>(Xq, *dict, m_params, qT, nf);
        auto& grad = storage.local();
        grad(dict->dofs()) += (coef * m.val) * gp;
        grad(dict->primary_dofs()) += (coef * pn) * m.grad;
    });

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(X.size());
    for (const auto& g : storage) {
        grad += g;
    }
    return grad;
}

Eigen::SparseMatrix<double> ESPPotential::lagged_hessian(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const PSDProjectionMethod psd) const
{
    // Every sample's coefficient is nonnegative, so projecting each sample's
    // Hessian keeps the sum PSD; no per-element dense block is needed.
    const LaggedCoefficients c =
        lagged_coefficients(mesh, m_params, m_lagged_weights);
    const NearFarBarrier<> nf(m_params.barrier, m_params.dbar_factor());
    const auto& maps = collisions.maps();
    const int n_verts = X.rows();
    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;
    const auto project = [&](Eigen::MatrixXd& h) {
        if (psd != PSDProjectionMethod::NONE) {
            h = project_to_psd(h, psd);
        }
    };

    const auto vertices = map_entries(maps.vertex_collisions);
    tbb::parallel_for(size_t(0), vertices.size(), [&](size_t i) {
        const index_t v = vertices[i]->first;
        const auto& dict = *vertices[i]->second;
        auto [hn, hf] = PointPotentialHelper::
            evaluate_potential_hessian_at_vertex_with_cached_collisions_nearfar(
                X, dict, m_params, PSDProjectionMethod::NONE, nf);
        Eigen::MatrixXd h = c.vertex_near[v] * hn + c.vertex_far[v] * hf;
        project(h);
        local_hessian_to_global_triplets(
            h, dict.vertex_ids(), 3, storage.local(), n_verts);
    });

    const auto faces = map_entries(maps.face_collisions);
    tbb::parallel_for(size_t(0), faces.size(), [&](size_t i) {
        const index_t f = faces[i]->first;
        const auto& rule = *m_params.face_rule(mesh.face_areas()(f)).non_corner;
        const auto& dicts = faces[i]->second;
        for (size_t qi = 0; qi < std::min(rule.size(), dicts.size()); qi++) {
            const auto& qp = rule[qi];
            const Eigen::RowVector3d q =
                qp.lambda[0] * X.row(mesh.faces()(f, 0))
                + qp.lambda[1] * X.row(mesh.faces()(f, 1))
                + qp.lambda[2] * X.row(mesh.faces()(f, 2));
            auto [hn, hf] = PointPotentialHelper::
                evaluate_potential_hessian_at_face_interior_point_with_cached_collisions_nearfar(
                    VertexMatrixView<3>(X, q), *dicts[qi], m_params, qp.lambda,
                    PSDProjectionMethod::NONE, nf);
            const double s = FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
            Eigen::MatrixXd h = s * (c.face_near[f] * hn + c.face_far[f] * hf);
            project(h);
            local_hessian_to_global_triplets(
                h, dicts[qi]->vertex_ids(), 3, storage.local(), n_verts);
        }
    });

    const auto pairs = map_entries(maps.edge_edge_collisions);
    tbb::parallel_for(size_t(0), pairs.size(), [&](size_t i) {
        using T = ADHessian<12>;
        const auto& [pair, dict] = *pairs[i];
        const double coef = c.edge[pair.first];
        if (coef == 0) {
            return;
        }
        const Eigen::Matrix<double, 4, 3> p = ee_positions(mesh, X, pair);
        Eigen::Vector<double, 12> flat;
        for (int r = 0; r < 4; r++) {
            flat.segment<3>(3 * r) = p.row(r).transpose();
        }
        const Eigen::Matrix<T, 4, 3> pT = slice_positions<T, 4, 3>(flat);
        const T m = ee_sample_weight<T>(pT, p, m_params);
        if (m.val == 0) {
            return;
        }
        const T u =
            closest_point_uv<T>(pT.row(0), pT.row(1), pT.row(2), pT.row(3));
        const Eigen::RowVector3<T> qT = u * (pT.row(1) - pT.row(0)) + pT.row(0);
        const Eigen::RowVector3d q(qT(0).val, qT(1).val, qT(2).val);
        const VertexMatrixView<3> Xq(X, q);
        const double pn = PointPotentialHelper::
            evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                Xq, *dict, m_params, dict->ee_dtype(), nf);
        const Eigen::VectorXd gp = PointPotentialHelper::
            evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                T>(Xq, *dict, m_params, qT, nf);
        // H(m P) = m H(P) + P H(m) + sym(grad m (x) grad P), with m on the
        // four primary vertices.
        Eigen::MatrixXd h = m.val
            * PointPotentialHelper::
                evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions_near(
                                Xq, *dict, m_params, qT, nf);
        const auto& ids = dict->primary_local_ids();
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                h.block<3, 3>(ids[a] * 3, ids[b] * 3) +=
                    pn * m.Hess.block<3, 3>(a * 3, b * 3);
            }
            const Eigen::MatrixXd cross =
                m.grad.segment<3>(a * 3) * gp.transpose();
            h.middleRows(ids[a] * 3, 3) += cross;
            h.middleCols(ids[a] * 3, 3) += cross.transpose();
        }
        h *= coef;
        project(h);
        local_hessian_to_global_triplets(
            h, dict->vertex_ids(), 3, storage.local(), n_verts);
    });

    std::vector<Eigen::Triplet<double>> triplets;
    for (const auto& t : storage) {
        triplets.insert(triplets.end(), t.begin(), t.end());
    }
    Eigen::SparseMatrix<double> hess(X.size(), X.size());
    hess.setFromTriplets(triplets.begin(), triplets.end());
    return hess;
}

double ESPPotential::operator()(
    const ESPCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    assert(static_cast<size_t>(X.rows()) == mesh.num_vertices());
    collisions.check_quadrature(m_params);

    if (uses_lagged_weights(mesh)) {
        return lagged_value(collisions, mesh, X);
    }
    if (m_params.local_quadrature_transfer && mesh.dim() == 3
        && m_params.dbar_factor() > 0) {
        return detail::local_transfer_value(collisions, mesh, X, m_params);
    }

    if (collisions.empty()) {
        return 0;
    }

    double result = 0;

    if (mesh.dim() == 2) {
        tbb::enumerable_thread_specific<double> potential_storage(0.0);
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(m_params.quad_order);

        // Collect active edge ids into a flat vector for parallel indexing.
        std::vector<index_t> active_edges;
        active_edges.reserve(collisions.maps().edge_collisions_2d.size());
        for (const auto& [ei, _] : collisions.maps().edge_collisions_2d) {
            active_edges.push_back(ei);
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, active_edges.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                double& total = potential_storage.local();
                for (size_t k = r.begin(); k < r.end(); ++k) {
                    const index_t ei = active_edges[k];
                    const auto& qp_dicts =
                        collisions.maps().edge_collisions_2d.at(ei);
                    const double L = mesh.edge_area(ei);
                    const double w_edge = m_params.area_weights ? L : 1.;
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
                                     X_ext, dict, m_params);
                    }
                }
            });

        for (const double v : potential_storage) {
            result += v;
        }
    } else if (mesh.dim() == 3) {
        {
            tbb::enumerable_thread_specific<double> potential_storage(0.0);

            // Split near and far for every dbar_factor in (0, 2]. At 0 the
            // edge-edge terms are skipped, so the weights are constant.
            const double dbar_factor = m_params.dbar_factor();
            const bool use_nf = m_use_near_far && dbar_factor > 0;
            const bool skip_ee =
                (dbar_factor == 0); // Skip EE pairs when dbar_factor == 0

            std::unique_ptr<NearFarBarrier<>> nf_barrier;
            if (use_nf) {
                nf_barrier = std::make_unique<NearFarBarrier<>>(
                    m_params.barrier, dbar_factor);
            }

            auto loop_body = [&](const tbb::blocked_range<size_t>& r) {
                double& total = potential_storage.local();
                for (size_t k = r.begin(); k < r.end(); k++) {
                    const IntegrationElement elem =
                        integration_element(mesh, m_params, k);
                    const double w = elem.weight;

                    double total_w = 0, total_w_near = 0, total_w_far = 0;
                    double total_p = 0, total_p_near = 0, total_p_far = 0;

                    for (int le = 0; le < elem.num_edges; le++) {
                        const index_t edge_id = elem.edges[le];
                        const index_t ea = mesh.edges()(edge_id, 0);
                        const index_t eb = mesh.edges()(edge_id, 1);

                        for (index_t other_edge_id :
                             collisions.m_candidates.ee_set(mesh, edge_id)) {
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
                                    collisions.maps().edge_edge_collisions.find(
                                        std::make_pair(edge_id, other_edge_id));
                                iter
                                != collisions.maps()
                                       .edge_edge_collisions.end()) {

                                const auto dtype = iter->second->ee_dtype();
                                // The builder keeps only EA_EB pairs.
                                assert(dtype == EdgeEdgeDistanceType::EA_EB);

                                const double dist =
                                    sqrt(edge_edge_distance_parallel_safe(
                                        X.row(ea), X.row(eb), X.row(ec),
                                        X.row(ed), dtype));

                                const double uv = closest_point_uv<double>(
                                    X.row(ea), X.row(eb), X.row(ec), X.row(ed));

                                const Eigen::RowVector3d ee_closest_point =
                                    uv * (X.row(eb) - X.row(ea)) + X.row(ea);

                                const double dist_sqr =
                                    edge_edge_distance_parallel_safe(
                                        X.row(ea), X.row(eb), X.row(ec),
                                        X.row(ed), dtype);
                                const auto mtypes = edge_edge_mollifier_type(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), dist_sqr,
                                    m_params.ee_mollifier);

                                double mollifier =
                                    Math<double>::cubic_spline(
                                        dist / m_params.ee_support())
                                    * m_params.ee_weight;
                                mollifier *= edge_edge_mollifier<double>(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), mtypes, dist_sqr,
                                    m_params.ee_mollifier);
                                mollifier = pow_int(
                                    mollifier,
                                    m_params.ee_mollifier_exponent());
                                mollifier = mollifier
                                    * ee_area_scale(mesh, elem, m_params);

                                if (use_nf) {
                                    const double p_near = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                                            VertexMatrixView<3>(
                                                X, ee_closest_point),
                                            *(iter->second), m_params, dtype,
                                            *nf_barrier);
                                    total_w_near += mollifier;
                                    total_p_near += mollifier * p_near;
                                } else {
                                    const double P_val = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions(
                                            VertexMatrixView<3>(
                                                X, ee_closest_point),
                                            *(iter->second), m_params, dtype);
                                    total_w += mollifier;
                                    total_p += mollifier * P_val;
                                }
                            }
                        }
                    }

                    // Face-interior quadrature points from
                    // params.quad_rule() (faces only).
                    const FaceQuadRule& face_quad_rule =
                        elem.face >= 0 ? *elem.non_corner : NO_FACE_QUADRATURE;
                    {
                        auto iter =
                            collisions.maps().face_collisions.find(elem.face);
                        for (size_t qi = 0; qi < face_quad_rule.size(); qi++) {
                            const auto& qp = face_quad_rule[qi];
                            if (use_nf) {
                                total_w_near +=
                                    FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                                total_w_far +=
                                    FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                            } else {
                                total_w +=
                                    FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                            }
                            if (iter
                                != collisions.maps().face_collisions.end()) {
                                const Eigen::RowVector3d q_pos = qp.lambda[0]
                                        * X.row(mesh.faces()(elem.face, 0))
                                    + qp.lambda[1]
                                        * X.row(mesh.faces()(elem.face, 1))
                                    + qp.lambda[2]
                                        * X.row(mesh.faces()(elem.face, 2));
                                if (use_nf) {
                                    auto [fq_near, fq_far] = PointPotentialHelper::
                                        evaluate_potential_at_face_interior_point_with_cached_collisions_nearfar(
                                            VertexMatrixView<3>(X, q_pos),
                                            *iter->second[qi], m_params,
                                            *nf_barrier);
                                    total_p_near += FACE_QUADRATURE_WEIGHT_SCALE
                                        * qp.weight * fq_near;
                                    total_p_far += FACE_QUADRATURE_WEIGHT_SCALE
                                        * qp.weight * fq_far;
                                } else {
                                    const double fq_val = PointPotentialHelper::
                                        evaluate_potential_at_face_interior_point_with_cached_collisions(
                                            VertexMatrixView<3>(X, q_pos),
                                            *iter->second[qi], m_params);
                                    total_p += FACE_QUADRATURE_WEIGHT_SCALE
                                        * qp.weight * fq_val;
                                }
                            }
                        }
                    }

                    // Vertex samples (see integration_element())
                    for (int lv = 0; lv < elem.num_vertices; lv++) {
                        const index_t v = elem.vertices[lv];
                        const double vw = elem.vertex_weights[lv];
                        if (vw == 0) {
                            continue;
                        }
                        if (use_nf) {
                            total_w_near += vw;
                            total_w_far += vw;
                        } else {
                            total_w += vw;
                        }
                        if (auto iter =
                                collisions.maps().vertex_collisions.find(v);
                            iter != collisions.maps().vertex_collisions.end()) {
                            if (use_nf) {
                                auto [vt_near, vt_far] = PointPotentialHelper::
                                    evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                                        X, *(iter->second), m_params,
                                        *nf_barrier);
                                total_p_near += vw * vt_near;
                                total_p_far += vw * vt_far;
                            } else {
                                const double vt_val = PointPotentialHelper::
                                    evaluate_potential_at_vertex_with_cached_collisions(
                                        X, *(iter->second), m_params);
                                total_p += vw * vt_val;
                            }
                        }
                    }

                    if (use_nf) {
                        // The element's own vertices or quadrature points
                        // always add a positive weight to both totals.
                        assert(total_w_near >= total_w_far - 1e-14);
                        assert(total_w_far > 0);
                        total += w
                            * (total_p_near / total_w_near
                               + total_p_far / total_w_far);
                    } else if (m_use_near_far) {
                        assert(total_w > 0);
                        total += w * (total_p / total_w);
                    } else {
                        total += w * total_p;
                    }
                }
            };

            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, num_integration_elements(mesh)),
                loop_body);

            for (const auto& local_potential : potential_storage) {
                result += local_potential;
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
    assert(static_cast<size_t>(X.rows()) == mesh.num_vertices());
    collisions.check_quadrature(m_params);

    if (uses_lagged_weights(mesh)) {
        return lagged_gradient(collisions, mesh, X);
    }
    if (m_params.local_quadrature_transfer && mesh.dim() == 3
        && m_params.dbar_factor() > 0) {
        return detail::local_transfer_gradient(collisions, mesh, X, m_params);
    }

    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(X.size());
    }

    const int dim = X.cols();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(X.size()));

    if (mesh.dim() == 2) {
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(m_params.quad_order);

        std::vector<index_t> active_edges;
        active_edges.reserve(collisions.maps().edge_collisions_2d.size());
        for (const auto& [ei, _] : collisions.maps().edge_collisions_2d) {
            active_edges.push_back(ei);
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, active_edges.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                Eigen::VectorXd& global_grad = storage.local();

                for (size_t k = r.begin(); k < r.end(); ++k) {
                    const index_t ei = active_edges[k];
                    const auto& qp_dicts =
                        collisions.maps().edge_collisions_2d.at(ei);
                    const double L = mesh.edge_area(ei);
                    const double w_edge = m_params.area_weights ? L : 1.;
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
                                    X_ext, dict, m_params, lambda);

                        local_gradient_to_global_gradient(
                            local_grad, dict.vertex_ids(), dim, global_grad);
                    }
                }
            });
    } else if (mesh.dim() == 3) {
        {
            using T = ADGrad<12>;

            const double dbar_factor = m_params.dbar_factor();
            const bool use_nf_grad = m_use_near_far && dbar_factor > 0;
            const bool skip_ee_grad = (dbar_factor == 0);

            auto loop_body = [&](const tbb::blocked_range<size_t>& r) {
                Eigen::VectorXd& grad = storage.local();
                for (size_t k = r.begin(); k < r.end(); k++) {
                    const IntegrationElement elem =
                        integration_element(mesh, m_params, k);
                    const double w = elem.weight;
                    // Pass 1: collect all quadrature contributions for this
                    // face
                    struct EEGradEntry {
                        const ESPCollisionDict<PointType::EDGE>* dict;
                        double mol_val;
                        Eigen::Vector<double, 12> mol_grad;
                        double P;
                        Eigen::VectorXd grad_p;
                    };
                    struct ConstGradEntry {
                        const std::vector<index_t>* dofs;
                        Eigen::VectorXd grad_p_near,
                            grad_p_far;       // near/far or single (for non-nf)
                        double p_near, p_far; // near/far or single (for non-nf)
                    };
                    std::vector<EEGradEntry> ee_cache;
                    std::vector<ConstGradEntry> const_cache;
                    double total_w = 0;
                    double total_w_near = 0, total_p_near = 0;
                    double total_w_far = 0;

                    std::unique_ptr<NearFarBarrier<>> nf_barrier;
                    if (use_nf_grad) {
                        nf_barrier = std::make_unique<NearFarBarrier<>>(
                            m_params.barrier, m_params.dbar_factor());
                    }

                    for (int le = 0; le < elem.num_edges; le++) {
                        const index_t edge_id = elem.edges[le];
                        const index_t ea = mesh.edges()(edge_id, 0);
                        const index_t eb = mesh.edges()(edge_id, 1);

                        for (index_t other_edge_id :
                             collisions.m_candidates.ee_set(mesh, edge_id)) {
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
                                    collisions.maps().edge_edge_collisions.find(
                                        std::make_pair(edge_id, other_edge_id));
                                iter
                                != collisions.maps()
                                       .edge_edge_collisions.end()) {

                                const auto dtype = iter->second->ee_dtype();
                                // The builder keeps only EA_EB pairs.
                                assert(dtype == EdgeEdgeDistanceType::EA_EB);

                                Eigen::Vector<double, 12> positions;
                                positions << X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose();

                                Eigen::Matrix<T, 4, 3> positionsT =
                                    slice_positions<T, 4, 3>(positions);

                                const T dist = sqrt(
                                    line_line_sqr_distance<T>(
                                        positionsT.row(0), positionsT.row(1),
                                        positionsT.row(2), positionsT.row(3)));

                                const T uv = closest_point_uv<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3));

                                const Eigen::RowVector3<T> ee_closest_point_T =
                                    uv * (positionsT.row(1) - positionsT.row(0))
                                    + positionsT.row(0);
                                const Eigen::RowVector3<double>
                                    ee_closest_point(
                                        ee_closest_point_T(0).val,
                                        ee_closest_point_T(1).val,
                                        ee_closest_point_T(2).val);

                                const T dist_sqr = line_line_sqr_distance<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3));
                                const auto mtypes = edge_edge_mollifier_type(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), dist_sqr.val,
                                    m_params.ee_mollifier);

                                T mollifier = Math<T>::cubic_spline(
                                                  dist / m_params.ee_support())
                                    * m_params.ee_weight;
                                mollifier *= edge_edge_mollifier<T>(
                                    positionsT.row(0).transpose(),
                                    positionsT.row(1).transpose(),
                                    positionsT.row(2).transpose(),
                                    positionsT.row(3).transpose(), mtypes,
                                    dist_sqr, m_params.ee_mollifier);
                                mollifier = pow_int(
                                    mollifier,
                                    m_params.ee_mollifier_exponent());
                                mollifier = mollifier
                                    * ee_area_scale(mesh, elem, m_params);

                                const ESPCollisionDict<PointType::EDGE>& dict =
                                    *(iter->second);

                                VertexMatrixView<3> X_extended(
                                    X, ee_closest_point);
                                assert(X_extended.rows() == X.rows() + 1);
                                assert(
                                    X_extended.m_a == X.data()
                                    && "VertexMatrixView has made a deepcopy!");

                                double P;
                                if (use_nf_grad) {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                                            X_extended, dict, m_params, dtype,
                                            *nf_barrier);
                                } else {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions(
                                            X_extended, dict, m_params, dtype);
                                }
                                Eigen::VectorXd grad_p;
                                if (use_nf_grad) {
                                    grad_p = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                                            T>(
                                            X_extended, dict, m_params,
                                            ee_closest_point_T, *nf_barrier);
                                } else {
                                    grad_p = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions<
                                            T>(
                                            X_extended, dict, m_params,
                                            ee_closest_point_T);
                                }

                                ee_cache.push_back(
                                    { &dict, mollifier.val, mollifier.grad, P,
                                      grad_p });
                                total_w += mollifier.val;
                                if (use_nf_grad) {
                                    total_w_near += mollifier.val;
                                    total_p_near += mollifier.val * P;
                                }
                            }
                        }
                    }

                    // Face-interior quadrature points (faces only)
                    const FaceQuadRule& face_quad_rule =
                        elem.face >= 0 ? *elem.non_corner : NO_FACE_QUADRATURE;
                    {
                        auto iter =
                            collisions.maps().face_collisions.find(elem.face);
                        for (size_t qi = 0; qi < face_quad_rule.size(); qi++) {
                            const auto& qp = face_quad_rule[qi];
                            const double qp_weight_scale =
                                FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                            total_w += qp_weight_scale;
                            if (use_nf_grad) {
                                total_w_near += qp_weight_scale;
                                total_w_far += qp_weight_scale;
                            }
                            if (iter != collisions.maps().face_collisions.end()
                                && qi < iter->second.size()) {
                                const auto& dict = *iter->second[qi];
                                const Eigen::RowVector3d q_pos = qp.lambda[0]
                                        * X.row(mesh.faces()(elem.face, 0))
                                    + qp.lambda[1]
                                        * X.row(mesh.faces()(elem.face, 1))
                                    + qp.lambda[2]
                                        * X.row(mesh.faces()(elem.face, 2));
                                VertexMatrixView<3> X_qp(X, q_pos);

                                if (use_nf_grad) {
                                    auto [P_n, P_f] = PointPotentialHelper::
                                        evaluate_potential_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, m_params, *nf_barrier);
                                    auto [grad_n, grad_f] = PointPotentialHelper::
                                        evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, m_params, qp.lambda,
                                            *nf_barrier);
                                    const_cache.push_back(
                                        ConstGradEntry {
                                            &dict.dofs(),
                                            qp_weight_scale * grad_n,
                                            qp_weight_scale * grad_f,
                                            qp_weight_scale * P_n,
                                            qp_weight_scale * P_f });
                                    total_p_near += qp_weight_scale * P_n;
                                } else {
                                    const double P = PointPotentialHelper::
                                        evaluate_potential_at_face_interior_point_with_cached_collisions(
                                            X_qp, dict, m_params);
                                    const Eigen::VectorXd grad_p =
                                        PointPotentialHelper::
                                            evaluate_potential_gradient_at_face_interior_point_with_cached_collisions(
                                                X_qp, dict, m_params,
                                                qp.lambda);
                                    const_cache.push_back(
                                        ConstGradEntry {
                                            &dict.dofs(),
                                            qp_weight_scale * grad_p,
                                            Eigen::VectorXd::Zero(0),
                                            qp_weight_scale * P, 0 });
                                }
                            }
                        }
                    }

                    // Vertex samples (see integration_element())
                    for (int lv = 0; lv < elem.num_vertices; lv++) {
                        const index_t v = elem.vertices[lv];
                        const double vw = elem.vertex_weights[lv];
                        if (vw == 0) {
                            continue;
                        }
                        total_w += vw;
                        if (use_nf_grad) {
                            total_w_near += vw;
                            total_w_far += vw;
                        }
                        if (auto iter =
                                collisions.maps().vertex_collisions.find(v);
                            iter != collisions.maps().vertex_collisions.end()) {
                            if (use_nf_grad) {
                                auto [P_n, P_f] = PointPotentialHelper::
                                    evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                                        X, (*iter->second), m_params,
                                        *nf_barrier);
                                auto [grad_n, grad_f] = PointPotentialHelper::
                                    evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                                        X, (*iter->second), m_params,
                                        *nf_barrier);
                                const_cache.push_back(ConstGradEntry {
                                    &(*iter->second).dofs(), vw * grad_n,
                                    vw * grad_f, vw * P_n, vw * P_f });
                                total_p_near += vw * P_n;
                            } else {
                                const double P = PointPotentialHelper::
                                    evaluate_potential_at_vertex_with_cached_collisions(
                                        X, (*iter->second), m_params);
                                const Eigen::VectorXd grad_p = PointPotentialHelper::
                                    evaluate_potential_gradient_at_vertex_with_cached_collisions(
                                        X, (*iter->second), m_params);
                                const_cache.push_back(ConstGradEntry {
                                    &(*iter->second).dofs(), vw * grad_p,
                                    Eigen::VectorXd::Zero(0), vw * P, 0 });
                            }
                        }
                    }

                    // Pass 2: apply gradient
                    if (use_nf_grad) {
                        assert(total_w_near > 0);
                        assert(total_w_far > 0);
                        const double avg_P_near = total_p_near / total_w_near;
                        for (const auto& e : ee_cache) {
                            grad(e.dict->dofs()) +=
                                (w / total_w_near * e.mol_val) * e.grad_p;
                            grad(e.dict->primary_dofs()) +=
                                (w / total_w_near * (e.P - avg_P_near))
                                * e.mol_grad;
                        }
                        for (const auto& e : const_cache) {
                            grad(*e.dofs) += (w / total_w_near) * e.grad_p_near
                                + (w / total_w_far) * e.grad_p_far;
                        }
                    } else if (m_use_near_far) {
                        // Normalized but without NearFarBarrier splitting
                        // (only reached at dbar_factor = 0, which skips the
                        // edge-edge terms, so the weights are constant)
                        assert(total_w > 0);
                        assert(ee_cache.empty());
                        for (const auto& e : const_cache) {
                            grad(*e.dofs) += (w / total_w) * e.grad_p_near;
                        }
                    } else {
                        // Unnormalized
                        for (const auto& e : ee_cache) {
                            grad(e.dict->dofs()) += w * e.mol_val * e.grad_p;
                            grad(e.dict->primary_dofs()) +=
                                w * e.P * e.mol_grad;
                        }
                        for (const auto& e : const_cache) {
                            grad(*e.dofs) += w * e.grad_p_near;
                        }
                    }
                }
            };

            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, num_integration_elements(mesh)),
                loop_body);
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
    assert(static_cast<size_t>(X.rows()) == mesh.num_vertices());
    collisions.check_quadrature(m_params);

    if (uses_lagged_weights(mesh)) {
        return lagged_hessian(collisions, mesh, X, project_hessian_to_psd);
    }
    if (m_params.local_quadrature_transfer && mesh.dim() == 3
        && m_params.dbar_factor() > 0) {
        return detail::local_transfer_hessian(
            collisions, mesh, X, m_params, project_hessian_to_psd);
    }

    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(X.size(), X.size());
    }

    const int dim = X.cols();
    const int n_verts = X.rows();
    const int ndof = X.size();

    const int max_triplets_size = int(1e7);
    const int buffer_size = std::min(max_triplets_size, ndof);
    tbb::enumerable_thread_specific<LocalThreadMatStorage> storage(
        LocalThreadMatStorage(buffer_size, ndof, ndof));

    if (mesh.dim() == 2) {
        const GaussLobatto::Rule& rule =
            GaussLobatto::get_rule(m_params.quad_order);

        std::vector<index_t> active_edges;
        active_edges.reserve(collisions.maps().edge_collisions_2d.size());
        for (const auto& [ei, _] : collisions.maps().edge_collisions_2d) {
            active_edges.push_back(ei);
        }

        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, active_edges.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                auto& hess_triplets = storage.local();

                for (size_t k = r.begin(); k < r.end(); ++k) {
                    const index_t ei = active_edges[k];
                    const auto& qp_dicts =
                        collisions.maps().edge_collisions_2d.at(ei);
                    const double L = mesh.edge_area(ei);
                    const double w_edge = m_params.area_weights ? L : 1.;
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
                                    X_ext, dict, m_params, lambda,
                                    project_hessian_to_psd);

                        local_hessian_to_global_triplets(
                            local_hess, dict.vertex_ids(), dim,
                            *(hess_triplets.cache), n_verts);
                    }
                }
            });
    } else if (mesh.dim() == 3) {
        // With normalized weights (use_near_far) and 0 < dbar_factor <= 2, the
        // per-face hessian is assembled as
        //   Term A (sum of per-stencil H(p_i))
        // + Term B (negative weighted sum of H(mol_i))
        // + Term C (sign-indefinite cross terms ~ sym(G⊗∇Z)).
        // Terms B and C are indefinite, so projecting the per-stencil blocks
        // is not enough. To still guarantee a PSD per-face contribution (and
        // thus a PSD global hessian), defer all per-stencil projections and
        // project the assembled per-face block once, over the union of
        // involved DOFs. At dbar_factor = 0 the edge-edge terms are skipped,
        // so B = C = 0 and the per-stencil projections suffice; with
        // use_near_far = false the local-projection path is used.
        const double dbar_factor_hess = m_params.dbar_factor();
        const bool use_nf_hess = m_use_near_far && dbar_factor_hess > 0;
        const bool skip_ee_hess = (dbar_factor_hess == 0);
        const bool combined_psd_projection =
            use_nf_hess && project_hessian_to_psd != PSDProjectionMethod::NONE;
        const PSDProjectionMethod inner_psd_method = combined_psd_projection
            ? PSDProjectionMethod::NONE
            : project_hessian_to_psd;
        {
            using T = ADHessian<12>;

            auto loop_body = [&](const tbb::blocked_range<size_t>& r) {
                auto& hess_triplets = storage.local();
                for (size_t k = r.begin(); k < r.end(); k++) {
                    const IntegrationElement elem =
                        integration_element(mesh, m_params, k);
                    const double w = elem.weight;

                    // Pass 1: collect all quadrature contributions for this
                    // face
                    struct EEHessEntry {
                        const ESPCollisionDict<PointType::EDGE>* dict;
                        double mol_val;
                        Eigen::Vector<double, 12> mol_grad; // on primary_dofs
                        Eigen::Matrix<double, 12, 12>
                            mol_hess; // H(mol) on primary_dofs
                        double P;     // only near component for use_nf_hess
                        Eigen::VectorXd grad_p;     // indexed by dict->dofs()
                        Eigen::MatrixXd local_hess; // H(mol*P), PSD-projected
                    };
                    struct ConstHessEntry {
                        const std::vector<index_t>* vertex_ids;
                        const std::vector<index_t>* dofs;
                        double p_near, p_far;
                        Eigen::VectorXd grad_p_near,
                            grad_p_far; // indexed by dofs
                        Eigen::MatrixXd local_hess_near,
                            local_hess_far; // H(p_near), H(p_far)
                    };
                    std::vector<EEHessEntry> ee_cache;
                    std::vector<ConstHessEntry> const_cache;
                    double total_w = 0;
                    double total_w_near = 0, total_p_near = 0;
                    double total_w_far = 0;
                    [[maybe_unused]] double total_p_far = 0; // asserts only

                    // Construct NearFarBarrier if needed
                    std::unique_ptr<NearFarBarrier<>> nf_barrier;
                    if (use_nf_hess) {
                        nf_barrier = std::make_unique<NearFarBarrier<>>(
                            m_params.barrier, m_params.dbar_factor());
                    }

                    for (int le = 0; le < elem.num_edges; le++) {
                        const index_t edge_id = elem.edges[le];
                        const index_t ea = mesh.edges()(edge_id, 0);
                        const index_t eb = mesh.edges()(edge_id, 1);

                        for (index_t other_edge_id :
                             collisions.m_candidates.ee_set(mesh, edge_id)) {
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
                                    collisions.maps().edge_edge_collisions.find(
                                        std::make_pair(edge_id, other_edge_id));
                                iter
                                != collisions.maps()
                                       .edge_edge_collisions.end()) {

                                const auto dtype = iter->second->ee_dtype();
                                // The builder keeps only EA_EB pairs.
                                assert(dtype == EdgeEdgeDistanceType::EA_EB);

                                Eigen::Vector<double, 12> positions;
                                positions << X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose();

                                Eigen::Matrix<T, 4, 3> positionsT =
                                    slice_positions<T, 4, 3>(positions);

                                const T dist = sqrt(
                                    line_line_sqr_distance<T>(
                                        positionsT.row(0), positionsT.row(1),
                                        positionsT.row(2), positionsT.row(3)));

                                const T uv = closest_point_uv<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3));

                                const Eigen::RowVector3<T> ee_closest_point_T =
                                    uv * (positionsT.row(1) - positionsT.row(0))
                                    + positionsT.row(0);
                                const Eigen::RowVector3<double>
                                    ee_closest_point(
                                        ee_closest_point_T(0).val,
                                        ee_closest_point_T(1).val,
                                        ee_closest_point_T(2).val);

                                const T dist_sqr = line_line_sqr_distance<T>(
                                    positionsT.row(0), positionsT.row(1),
                                    positionsT.row(2), positionsT.row(3));
                                const auto mtypes = edge_edge_mollifier_type(
                                    X.row(ea).transpose(),
                                    X.row(eb).transpose(),
                                    X.row(ec).transpose(),
                                    X.row(ed).transpose(), dist_sqr.val,
                                    m_params.ee_mollifier);

                                T mollifier = Math<T>::cubic_spline(
                                                  dist / m_params.ee_support())
                                    * m_params.ee_weight;
                                mollifier *= edge_edge_mollifier<T>(
                                    positionsT.row(0).transpose(),
                                    positionsT.row(1).transpose(),
                                    positionsT.row(2).transpose(),
                                    positionsT.row(3).transpose(), mtypes,
                                    dist_sqr, m_params.ee_mollifier);
                                mollifier = pow_int(
                                    mollifier,
                                    m_params.ee_mollifier_exponent());
                                mollifier = mollifier
                                    * ee_area_scale(mesh, elem, m_params);

                                const ESPCollisionDict<PointType::EDGE>& dict =
                                    *(iter->second);

                                VertexMatrixView<3> X_extended(
                                    X, ee_closest_point);
                                assert(
                                    X_extended.m_a == X.data()
                                    && "VertexMatrixView has made a deepcopy!");

                                double P;
                                Eigen::VectorXd grad_p;
                                Eigen::MatrixXd base_hess;
                                if (use_nf_hess) {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions_near(
                                            X_extended, dict, m_params, dtype,
                                            *nf_barrier);
                                    grad_p = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions_near<
                                            T>(
                                            X_extended, dict, m_params,
                                            ee_closest_point_T, *nf_barrier);
                                    base_hess = PointPotentialHelper::
                                        evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions_near(
                                            X_extended, dict, m_params,
                                            ee_closest_point_T, *nf_barrier);
                                } else {
                                    P = PointPotentialHelper::
                                        evaluate_potential_at_edge_edge_closest_point_with_cached_collisions(
                                            X_extended, dict, m_params, dtype);
                                    grad_p = PointPotentialHelper::
                                        evaluate_potential_gradient_at_edge_edge_closest_point_with_cached_collisions<
                                            T>(
                                            X_extended, dict, m_params,
                                            ee_closest_point_T);
                                    base_hess = PointPotentialHelper::
                                        evaluate_potential_hessian_at_edge_edge_closest_point_with_cached_collisions(
                                            X_extended, dict, m_params,
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
                                        * grad_p.transpose();
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
                                    local_hess = project_to_psd(
                                        local_hess, project_hessian_to_psd);
                                }

                                ee_cache.push_back(
                                    { &dict, mollifier.val, mollifier.grad,
                                      mollifier.Hess, P, grad_p,
                                      std::move(local_hess) });
                                total_w += mollifier.val;
                                if (use_nf_hess) {
                                    total_w_near += mollifier.val;
                                    total_p_near += mollifier.val * P;
                                }
                            }
                        }
                    }

                    // Face-interior quadrature points (faces only)
                    const FaceQuadRule& face_quad_rule =
                        elem.face >= 0 ? *elem.non_corner : NO_FACE_QUADRATURE;
                    {
                        auto iter =
                            collisions.maps().face_collisions.find(elem.face);
                        for (size_t qi = 0; qi < face_quad_rule.size(); qi++) {
                            const auto& qp = face_quad_rule[qi];
                            total_w += FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                            if (use_nf_hess) {
                                total_w_near +=
                                    FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                                total_w_far +=
                                    FACE_QUADRATURE_WEIGHT_SCALE * qp.weight;
                            }
                            if (iter != collisions.maps().face_collisions.end()
                                && qi < iter->second.size()) {
                                const auto& dict = *iter->second[qi];
                                const Eigen::RowVector3d q_pos = qp.lambda[0]
                                        * X.row(mesh.faces()(elem.face, 0))
                                    + qp.lambda[1]
                                        * X.row(mesh.faces()(elem.face, 1))
                                    + qp.lambda[2]
                                        * X.row(mesh.faces()(elem.face, 2));
                                VertexMatrixView<3> X_qp(X, q_pos);
                                ConstHessEntry entry;
                                entry.vertex_ids = &dict.vertex_ids();
                                entry.dofs = &dict.dofs();

                                if (use_nf_hess) {
                                    auto [P_n, P_f] = PointPotentialHelper::
                                        evaluate_potential_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, m_params, *nf_barrier);
                                    auto [grad_n, grad_f] = PointPotentialHelper::
                                        evaluate_potential_gradient_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, m_params, qp.lambda,
                                            *nf_barrier);
                                    auto [hess_n, hess_f] = PointPotentialHelper::
                                        evaluate_potential_hessian_at_face_interior_point_with_cached_collisions_nearfar(
                                            X_qp, dict, m_params, qp.lambda,
                                            inner_psd_method, *nf_barrier);
                                    entry.p_near = FACE_QUADRATURE_WEIGHT_SCALE
                                        * qp.weight * P_n;
                                    entry.p_far = FACE_QUADRATURE_WEIGHT_SCALE
                                        * qp.weight * P_f;
                                    entry.grad_p_near =
                                        FACE_QUADRATURE_WEIGHT_SCALE * qp.weight
                                        * grad_n;
                                    entry.grad_p_far =
                                        FACE_QUADRATURE_WEIGHT_SCALE * qp.weight
                                        * grad_f;
                                    entry.local_hess_near =
                                        FACE_QUADRATURE_WEIGHT_SCALE * qp.weight
                                        * hess_n;
                                    entry.local_hess_far =
                                        FACE_QUADRATURE_WEIGHT_SCALE * qp.weight
                                        * hess_f;
                                    total_p_near += entry.p_near;
                                    total_p_far += entry.p_far;
                                } else {
                                    entry.p_near = FACE_QUADRATURE_WEIGHT_SCALE
                                        * qp.weight
                                        * PointPotentialHelper::
                                            evaluate_potential_at_face_interior_point_with_cached_collisions(
                                                       X_qp, dict, m_params);
                                    entry.grad_p_near =
                                        FACE_QUADRATURE_WEIGHT_SCALE * qp.weight
                                        * PointPotentialHelper::
                                            evaluate_potential_gradient_at_face_interior_point_with_cached_collisions(
                                                X_qp, dict, m_params,
                                                qp.lambda);
                                    entry.local_hess_near =
                                        FACE_QUADRATURE_WEIGHT_SCALE * qp.weight
                                        * PointPotentialHelper::
                                            evaluate_potential_hessian_at_face_interior_point_with_cached_collisions(
                                                X_qp, dict, m_params, qp.lambda,
                                                inner_psd_method);
                                    entry.p_far = 0;
                                    entry.grad_p_far = Eigen::VectorXd::Zero(0);
                                    entry.local_hess_far =
                                        Eigen::MatrixXd::Zero(0, 0);
                                }
                                const_cache.push_back(std::move(entry));
                            }
                        }
                    }

                    // Vertex samples (see integration_element())
                    for (int lv = 0; lv < elem.num_vertices; lv++) {
                        const index_t v = elem.vertices[lv];
                        const double vw = elem.vertex_weights[lv];
                        if (vw == 0) {
                            continue;
                        }
                        total_w += vw;
                        if (use_nf_hess) {
                            total_w_near += vw;
                            total_w_far += vw;
                        }
                        if (auto iter =
                                collisions.maps().vertex_collisions.find(v);
                            iter != collisions.maps().vertex_collisions.end()) {
                            const auto& dict = *iter->second;
                            ConstHessEntry entry;
                            entry.vertex_ids = &dict.vertex_ids();
                            entry.dofs = &dict.dofs();

                            if (use_nf_hess) {
                                auto [P_n, P_f] = PointPotentialHelper::
                                    evaluate_potential_at_vertex_with_cached_collisions_nearfar(
                                        X, dict, m_params, *nf_barrier);
                                auto [grad_n, grad_f] = PointPotentialHelper::
                                    evaluate_potential_gradient_at_vertex_with_cached_collisions_nearfar(
                                        X, dict, m_params, *nf_barrier);
                                auto [hess_n, hess_f] = PointPotentialHelper::
                                    evaluate_potential_hessian_at_vertex_with_cached_collisions_nearfar(
                                        X, dict, m_params, inner_psd_method,
                                        *nf_barrier);
                                entry.p_near = vw * P_n;
                                entry.p_far = vw * P_f;
                                entry.grad_p_near = vw * grad_n;
                                entry.grad_p_far = vw * grad_f;
                                entry.local_hess_near = vw * hess_n;
                                entry.local_hess_far = vw * hess_f;
                                total_p_near += entry.p_near;
                                total_p_far += entry.p_far;
                            } else {
                                entry.p_near = vw
                                    * PointPotentialHelper::
                                        evaluate_potential_at_vertex_with_cached_collisions(
                                                   X, dict, m_params);
                                entry.grad_p_near = vw
                                    * PointPotentialHelper::
                                        evaluate_potential_gradient_at_vertex_with_cached_collisions(
                                                        X, dict, m_params);
                                entry.local_hess_near = vw
                                    * PointPotentialHelper::
                                        evaluate_potential_hessian_at_vertex_with_cached_collisions(
                                                            X, dict, m_params,
                                                            inner_psd_method);
                                entry.p_far = 0;
                                entry.grad_p_far = Eigen::VectorXd::Zero(0);
                                entry.local_hess_far =
                                    Eigen::MatrixXd::Zero(0, 0);
                            }
                            const_cache.push_back(std::move(entry));
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

                            // Maps a global DOF index (vid * dim + d in
                            // RowMajor, d * n_verts + vid in ColMajor) to its
                            // local index inside H_face, which is always
                            // vertex-interleaved.
                            auto global_dof_to_local = [&](index_t gd) {
                                if constexpr (
                                    VERTEX_DERIVATIVE_LAYOUT
                                    == Eigen::RowMajor) {
                                    return vid_to_local.at(gd / dim) * dim
                                        + static_cast<int>(gd % dim);
                                } else {
                                    return vid_to_local.at(gd % n_verts) * dim
                                        + static_cast<int>(gd / n_verts);
                                }
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
                                        ek.dict->dofs(), ek.mol_val * ek.grad_p,
                                        prim_dofs_i, mol_grad_i, scale_C_near);
                                }
                                // EE-const near interactions
                                for (const auto& ej : const_cache) {
                                    add_sym_correction_dense(
                                        *ej.dofs, ej.grad_p_near, prim_dofs_i,
                                        mol_grad_i, scale_C_near);
                                }
                            }

                            H_face =
                                project_to_psd(H_face, project_hessian_to_psd);

                            local_hessian_to_global_triplets(
                                H_face, union_vids, dim, *(hess_triplets.cache),
                                n_verts);
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
                                local_hessian_to_global_triplets(
                                    (w / total_w_near) * e.local_hess,
                                    e.dict->vertex_ids(), dim,
                                    *(hess_triplets.cache), n_verts);
                            }
                            for (const auto& e : const_cache) {
                                local_hessian_to_global_triplets(
                                    (w / total_w_near) * e.local_hess_near,
                                    *e.vertex_ids, dim, *(hess_triplets.cache),
                                    n_verts);
                                if (total_w_far > 0) {
                                    local_hessian_to_global_triplets(
                                        (w / total_w_far) * e.local_hess_far,
                                        *e.vertex_ids, dim,
                                        *(hess_triplets.cache), n_verts);
                                }
                            }

                            // Term B: -(w*avg_P_near/total_w_near) * Σ_i
                            // H(mol_i) (Z_far has no V-dependence, so no far
                            // Term B.)
                            for (const auto& e : ee_cache) {
                                local_hessian_to_global_triplets(
                                    -(w * avg_P_near / total_w_near)
                                        * e.mol_hess,
                                    e.dict->primary_vertex_ids(), dim,
                                    *(hess_triplets.cache), n_verts);
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
                                        ek.dict->dofs(), ek.mol_val * ek.grad_p,
                                        prim_dofs_i, mol_grad_i, scale_C_near);
                                }
                                // EE-const near interactions
                                for (const auto& ej : const_cache) {
                                    add_sym_correction(
                                        *ej.dofs, ej.grad_p_near, prim_dofs_i,
                                        mol_grad_i, scale_C_near);
                                }
                            }
                        }
                    } else if (m_use_near_far) {
                        // Normalized but without NearFarBarrier splitting
                        // (only reached at dbar_factor = 0, which skips the
                        // edge-edge terms, so the weights are constant)
                        assert(total_w > 0);
                        assert(ee_cache.empty());
                        for (const auto& e : const_cache) {
                            local_hessian_to_global_triplets(
                                (w / total_w) * e.local_hess_near,
                                *e.vertex_ids, dim, *(hess_triplets.cache),
                                n_verts);
                        }
                    } else {
                        // Unnormalized: H(w*p_sum) = w * H(p_sum)
                        for (const auto& e : ee_cache) {
                            local_hessian_to_global_triplets(
                                w * e.local_hess, e.dict->vertex_ids(), dim,
                                *(hess_triplets.cache), n_verts);
                        }
                        for (const auto& e : const_cache) {
                            local_hessian_to_global_triplets(
                                w * e.local_hess_near, *e.vertex_ids, dim,
                                *(hess_triplets.cache), n_verts);
                        }
                    }
                }
            };

            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, num_integration_elements(mesh)),
                loop_body);
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
        return hess;
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
    return collision.weight * collision(positions, m_params);
}

Eigen::VectorXd ESPPotential::gradient(
    const ESPCollision& collision,
    Eigen::ConstRef<Eigen::VectorXd> positions) const
{
    return collision.weight * collision.gradient(positions, m_params);
}

Eigen::MatrixXd ESPPotential::hessian(
    const ESPCollision& collision,
    Eigen::ConstRef<Eigen::VectorXd> positions,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    Eigen::MatrixXd hess =
        collision.weight * collision.hessian(positions, m_params);
    return project_to_psd(hess, project_hessian_to_psd);
}

} // namespace ipc
