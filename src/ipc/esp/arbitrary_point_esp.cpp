#include "arbitrary_point_esp.hpp"

#include <ipc/distance/distance_type_exact.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/esp/collisions/esp_collision_template.hpp>
#include <ipc/esp/collisions/vertex_matrix_view.hpp>
#include <ipc/esp/inclusion_exclusion_weights.hpp>

#include <array>
#include <cassert>
#include <type_traits>
#include <vector>

namespace ipc {

namespace {

    // The collision set at one query point, in per-thread lists that outlive
    // the call: build_point_collisions() clears them and refills them, so a
    // query allocates nothing once each list has grown to its working size.
    // One list per collision type, in the order ESPCollisionDict
    // evaluates them (vertex-vertex, edge-vertex, face-vertex), each holding
    // the collisions by value -- they carry only fixed-size arrays of ids.
    template <int dim> struct PointCollisions {
        using VV = typename ESPCollisionDict<PointType::VERTEX, dim>::VVType;
        using EV = typename ESPCollisionDict<PointType::VERTEX, dim>::EVType;
        using FV = ESPCollisionTemplate<Face3P1, Vertex3>;

        std::vector<index_t> vertex_ids, edge_ids, face_ids; // broad phase
        std::vector<VV> vv;
        std::vector<EV> ev;
        std::vector<FV> fv; // unused in 2D

        template <typename F> void for_each(F&& f) const
        {
            for (const VV& c : vv) {
                f(c);
            }
            for (const EV& c : ev) {
                f(c);
            }
            for (const FV& c : fv) {
                f(c);
            }
        }
    };

    template <int dim> PointCollisions<dim>& point_collisions_scratch()
    {
        thread_local PointCollisions<dim> scratch;
        return scratch;
    }

    // Merges collisions that resolve to the same underlying feature (same
    // typed hash) by accumulating their weight, and drops the entry
    // entirely if the accumulated weight is exactly zero. This is the
    // symbolic-cancellation step: redundant +1/-1 contributions from
    // different primitives converging on the same feature (e.g. two
    // adjacent faces and their shared edge) cancel here, as integers,
    // before any barrier value is ever computed -- never as a
    // floating-point subtraction of two independently-rounded nearly-equal
    // barrier values (which is numerically unstable near dhat -> 0, where
    // the barrier and its derivatives blow up).
    //
    // The same merge the hash map of insert_pair() in quadrature_potential.cpp
    // performs, as a linear search over one type's list: the typed hash
    // starts with the type, so two collisions of different types never
    // merge, and a query point has only a handful of collisions.
    template <typename Collision>
    void insert_by_value(
        std::vector<Collision>& list,
        Collision&& collision,
        const double weight)
    {
        collision.weight = weight;
        const std::array<index_t, 3> hash = collision.get_typed_hash();
        for (auto it = list.begin(); it != list.end(); ++it) {
            if (it->get_typed_hash() == hash) {
                it->weight += collision.weight;
                if (it->weight == 0) {
                    list.erase(it);
                }
                return;
            }
        }
        list.push_back(std::move(collision));
    }

    // ESPCollisionsBuilder<3>::reduce_point_triangle_collision() and
    // reduce_point_edge_collision() (esp_collisions_builder.cpp),
    // emitting into the by-value lists instead of a new shared_ptr: the same
    // distance-type classification, the same dhat cut, the same collision
    // type for each case.
    void reduce_point_triangle(
        PointCollisions<3>& s,
        const index_t fi,
        const index_t vid,
        const double weight,
        const ESPParameters& params,
        const CollisionMesh& mesh,
        const VertexMatrixView<3>& vertices)
    {
        using VV = PointCollisions<3>::VV;
        using EV = PointCollisions<3>::EV;
        using FV = PointCollisions<3>::FV;

        const index_t t0 = mesh.faces()(fi, 0);
        const index_t t1 = mesh.faces()(fi, 1);
        const index_t t2 = mesh.faces()(fi, 2);

        const index_t e0 = mesh.faces_to_edges()(fi, 0);
        const index_t e1 = mesh.faces_to_edges()(fi, 1);
        const index_t e2 = mesh.faces_to_edges()(fi, 2);

        assert(vid != t0 && vid != t1 && vid != t2);

        const PointTriangleDistanceType dtype =
            point_triangle_distance_type_exact(
                vertices(vid), vertices(t0), vertices(t1), vertices(t2));

        const double dist_sqr = point_triangle_distance(
            vertices(vid), vertices(t0), vertices(t1), vertices(t2), dtype);
        if (dist_sqr >= params.dhat * params.dhat) {
            return;
        }

        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            insert_by_value(s.vv, VV(t0, vid, mesh), weight);
            break;
        case PointTriangleDistanceType::P_T1:
            insert_by_value(s.vv, VV(t1, vid, mesh), weight);
            break;
        case PointTriangleDistanceType::P_T2:
            insert_by_value(s.vv, VV(t2, vid, mesh), weight);
            break;
        case PointTriangleDistanceType::P_E0:
            insert_by_value(s.ev, EV(e0, vid, mesh), weight);
            break;
        case PointTriangleDistanceType::P_E1:
            insert_by_value(s.ev, EV(e1, vid, mesh), weight);
            break;
        case PointTriangleDistanceType::P_E2:
            insert_by_value(s.ev, EV(e2, vid, mesh), weight);
            break;
        default: // point_triangle_distance_type_exact never returns AUTO
            assert(dtype == PointTriangleDistanceType::P_T);
            insert_by_value(s.fv, FV(fi, vid, mesh), weight);
            break;
        }
    }

    void reduce_point_edge(
        PointCollisions<3>& s,
        const index_t ei,
        const index_t vid,
        const double weight,
        const ESPParameters& params,
        const CollisionMesh& mesh,
        const VertexMatrixView<3>& vertices)
    {
        using VV = PointCollisions<3>::VV;
        using EV = PointCollisions<3>::EV;

        const index_t t0 = mesh.edges()(ei, 0);
        const index_t t1 = mesh.edges()(ei, 1);

        const PointEdgeDistanceType dtype = point_edge_distance_type_exact(
            vertices(vid), vertices(t0), vertices(t1));

        const double dist_sqr = point_edge_distance(
            vertices(vid), vertices(t0), vertices(t1), dtype);
        if (dist_sqr >= params.dhat * params.dhat) {
            return;
        }

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            insert_by_value(s.vv, VV(t0, vid, mesh), weight);
            break;
        case PointEdgeDistanceType::P_E1:
            insert_by_value(s.vv, VV(t1, vid, mesh), weight);
            break;
        case PointEdgeDistanceType::P_E:
        default:
            assert(dtype == PointEdgeDistanceType::P_E);
            insert_by_value(s.ev, EV(ei, vid, mesh), weight);
            break;
        }
    }

    // 2D counterpart of ESPCollisionsBuilder<3>::
    // reduce_point_edge_collision (esp_collisions_builder.cpp), which
    // only exists on the <3> specialization: classify which sub-feature of
    // edge ei the query's closest point falls on and emit the
    // correspondingly-typed collision, so endpoint cases hash-merge with the
    // direct vertex collisions instead of double-counting a corner.
    //
    // Query vertex first in both templates, matching the convention of the
    // 2D edge-QP builder in quadrature_potential.cpp; Vertex2-Edge2P1's
    // evaluators assume that layout ([q, e0, e1]).
    void reduce_point_edge_2d(
        PointCollisions<2>& s,
        const index_t ei,
        const index_t vid,
        const double weight,
        const ESPParameters& params,
        const CollisionMesh& mesh,
        const VertexMatrixView<2>& vertices)
    {
        using VV = PointCollisions<2>::VV;
        using EV = PointCollisions<2>::EV;

        const index_t e0 = mesh.edges()(ei, 0);
        const index_t e1 = mesh.edges()(ei, 1);

        const PointEdgeDistanceType dtype = point_edge_distance_type_exact(
            vertices(vid), vertices(e0), vertices(e1));

        const double dist_sqr = point_edge_distance(
            vertices(vid), vertices(e0), vertices(e1), dtype);
        if (dist_sqr >= params.dhat * params.dhat) {
            return;
        }

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            insert_by_value(s.vv, VV(vid, e0, mesh), weight);
            break;
        case PointEdgeDistanceType::P_E1:
            insert_by_value(s.vv, VV(vid, e1, mesh), weight);
            break;
        case PointEdgeDistanceType::P_E:
            insert_by_value(s.ev, EV(vid, ei, mesh), weight);
            break;
        default:
            assert(false);
            break;
        }
    }

    // The collisions at query point q, in the calling thread's scratch.
    template <int dim>
    const PointCollisions<dim>& build_point_collisions(
        const LBVH& point_bvh,
        const CollisionMesh& mesh,
        const ESPParameters& params,
        const std::vector<int>& edge_weights,
        const std::vector<int>& vertex_weights,
        Eigen::ConstRef<Eigen::MatrixXd> V,
        Eigen::ConstRef<Eigen::RowVector<double, dim>> q)
    {
        PointCollisions<dim>& s = point_collisions_scratch<dim>();
        s.vv.clear();
        s.ev.clear();
        s.fv.clear();

        const index_t vid = static_cast<index_t>(V.rows()); // virtual vertex id
        const VertexMatrixView<dim> V_view(V, q);

        point_bvh.query_point(
            q, params.dhat, s.vertex_ids, s.edge_ids, s.face_ids);

        // Inclusion-exclusion: every element whose offset region can contain q
        // contributes a term with the element's own weight (see the
        // constructor) -- on a closed mesh, in 3D faces +1, edges -1, vertices
        // +1; in 2D edges +1, vertices -1. Each element is first *reduced* to
        // the sub-feature its closest point to q actually lies on (never just
        // "this face's interior" regardless of where the closest point falls),
        // so redundant terms converging on the same feature share a typed hash
        // and cancel as integers in insert_by_value() above. An element of
        // weight zero (a boundary edge or vertex) contributes nothing.
        if constexpr (dim == 3) {
            for (const index_t fi : s.face_ids) {
                reduce_point_triangle(s, fi, vid, +1, params, mesh, V_view);
            }
            for (const index_t ei : s.edge_ids) {
                if (edge_weights[ei] != 0) {
                    reduce_point_edge(
                        s, ei, vid, edge_weights[ei], params, mesh, V_view);
                }
            }
        } else {
            // In 2D edges are the top-dimensional elements, so they take the
            // weight faces take in 3D and there is no face loop (mesh.faces()
            // is empty and the face BVH is never built).
            for (const index_t ei : s.edge_ids) {
                reduce_point_edge_2d(
                    s, ei, vid, edge_weights[ei], params, mesh, V_view);
            }
        }

        // Vertices, with their own weights. These merge (and symbolically
        // cancel) with the vertex-typed collisions the reductions above emit
        // when both resolve to the same corner.
        //
        // The (query, mesh vertex) argument order is load-bearing in 2D and
        // only in 2D: get_typed_hash() is {type, primitive_a.id(),
        // primitive_b.id()}, and Vertex2-Vertex2 goes through the generic
        // constructor, which stores the ids as given -- so this has to match
        // what reduce_point_edge_2d emits or the two never merge.
        // Vertex3-Vertex3 has a specialized constructor that sorts its two ids
        // (esp_collision_template.cpp), so 3D merges either way.
        using VV = typename PointCollisions<dim>::VV;
        for (const index_t vi : s.vertex_ids) {
            if (vertex_weights[vi] == 0
                || (V.row(vi) - q).squaredNorm() >= params.dhat * params.dhat) {
                continue;
            }
            insert_by_value(s.vv, VV(vid, vi, mesh), vertex_weights[vi]);
        }
        return s;
    }

    // The collision's stencil positions, as ESPCollision::dof() gives
    // them, on the stack rather than in a new Eigen::VectorXd.
    template <int dim, typename Collision>
    VectorMax<double, ESPCollision::ELEMENT_SIZE>
    stencil_positions(const Collision& cc, const VertexMatrixView<dim>& V_view)
    {
        VectorMax<double, ESPCollision::ELEMENT_SIZE> x(
            cc.num_vertices() * dim);
        for (int i = 0; i < cc.num_vertices(); i++) {
            x.template segment<dim>(i * dim) = V_view(cc.vertex_id(i));
        }
        return x;
    }

    // Where the query vertex sits in the collision's stencil. Every collision
    // at a query point has the query vertex as exactly one of its vertices,
    // and only its block of a gradient or Hessian is ever returned, so that
    // block is all that is accumulated: the same per-collision terms, added
    // in the same collision order, as assembling the whole stencil and then
    // extracting the query's block.
    template <typename Collision>
    int query_slot(const Collision& cc, const index_t vid)
    {
        for (int i = 0; i < cc.num_vertices(); i++) {
            if (cc.vertex_id(i) == vid) {
                return i;
            }
        }
        assert(false);
        return -1;
    }

} // namespace

template <int dim>
ArbitraryPointESP<dim>::ArbitraryPointESP(
    const CollisionMesh& _mesh, ESPParameters _params)
    : m_mesh(_mesh)
    , m_params(std::move(_params))
{
    if (m_mesh.dim() != dim) {
        log_and_throw_error(
            "ArbitraryPointESP<{}> requires a {}D mesh (got {}D)!", dim, dim,
            m_mesh.dim());
    }

    // Inclusion-exclusion weights from the mesh connectivity: the alternating
    // signs on a closed mesh, zero on boundary edges and vertices, and one on
    // edges in no face and isolated vertices.
    compute_inclusion_exclusion_weights(
        m_mesh, m_edge_weights, m_vertex_weights);
}

template <int dim>
void ArbitraryPointESP<dim>::update(Eigen::ConstRef<Eigen::MatrixXd> V)
{
    // Zero inflation: LBVH::query_point applies the radius to the query box, so
    // one build serves queries at any radius.
    m_point_bvh.build(
        V, m_mesh.edges(), m_mesh.faces(), /*inflation_radius=*/0.0);
}

template <int dim>
double ArbitraryPointESP<dim>::operator()(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
{
    const PointCollisions<dim>& collisions = build_point_collisions<dim>(
        m_point_bvh, m_mesh, m_params, m_edge_weights, m_vertex_weights, V, q);
    const VertexMatrixView<dim> V_view(V, q);

    double value = 0.0;
    collisions.for_each([&](const auto& cc) {
        value += cc.weight * cc(stencil_positions<dim>(cc, V_view), m_params);
    });
    return value;
}

template <int dim>
auto ArbitraryPointESP<dim>::gradient(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
    -> Gradient
{
    const PointCollisions<dim>& collisions = build_point_collisions<dim>(
        m_point_bvh, m_mesh, m_params, m_edge_weights, m_vertex_weights, V, q);
    const VertexMatrixView<dim> V_view(V, q);
    const index_t vid = static_cast<index_t>(V.rows());

    Gradient grad = Gradient::Zero();
    collisions.for_each([&](const auto& cc) {
        const VectorMax<double, ESPCollision::ELEMENT_SIZE> g = cc.weight
            * cc.gradient(stencil_positions<dim>(cc, V_view), m_params);
        grad += g.template segment<dim>(dim * query_slot(cc, vid));
    });
    return grad;
}

template <int dim>
auto ArbitraryPointESP<dim>::hessian(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
    -> Hessian
{
    const PointCollisions<dim>& collisions = build_point_collisions<dim>(
        m_point_bvh, m_mesh, m_params, m_edge_weights, m_vertex_weights, V, q);
    const VertexMatrixView<dim> V_view(V, q);
    const index_t vid = static_cast<index_t>(V.rows());

    Hessian H = Hessian::Zero();
    collisions.for_each([&](const auto& cc) {
        const int j = dim * query_slot(cc, vid);
        // One Hessian temporary (ELEMENT_SIZE^2 doubles on the stack), read
        // through a block of the scaled expression rather than copied.
        H += (cc.hessian(stencil_positions<dim>(cc, V_view), m_params)
              * cc.weight)
                 .template block<dim, dim>(j, j);
    });
    return H;
}

template <int dim>
auto ArbitraryPointESP<dim>::evaluate(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
    -> std::tuple<double, Gradient, Hessian>
{
    const PointCollisions<dim>& collisions = build_point_collisions<dim>(
        m_point_bvh, m_mesh, m_params, m_edge_weights, m_vertex_weights, V, q);
    const VertexMatrixView<dim> V_view(V, q);
    const index_t vid = static_cast<index_t>(V.rows());

    double value = 0.0;
    Gradient grad = Gradient::Zero();
    Hessian H = Hessian::Zero();

    collisions.for_each([&](const auto& cc) {
        const VectorMax<double, ESPCollision::ELEMENT_SIZE> dof =
            stencil_positions<dim>(cc, V_view);
        const int j = dim * query_slot(cc, vid);

        value += cc.weight * cc(dof, m_params);

        const VectorMax<double, ESPCollision::ELEMENT_SIZE> g =
            cc.weight * cc.gradient(dof, m_params);
        grad += g.template segment<dim>(j);

        H += (cc.weight * cc.hessian(dof, m_params))
                 .template block<dim, dim>(j, j);
    });

    return { value, grad, H };
}

template class ArbitraryPointESP<2>;
template class ArbitraryPointESP<3>;

} // namespace ipc
