#include "arbitrary_point_potential.hpp"

#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/face_vertex.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/high_order_contact/collisions/high_order_collision_template.hpp>
#include <ipc/high_order_contact/collisions/vertex_matrix_view.hpp>
#include <ipc/high_order_contact/high_order_collisions_builder.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

#include <array>
#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

namespace ipc {

namespace {

    // Mirrors the local insert_pair() helper in quadrature_potential.cpp:
    // merges collisions that resolve to the same underlying feature (same
    // typed hash) by accumulating their weight, and drops the entry
    // entirely if the accumulated weight is exactly zero. This is the
    // symbolic-cancellation step: redundant +1/-1 contributions from
    // different primitives converging on the same feature (e.g. two
    // adjacent faces and their shared edge) cancel here, as integers,
    // before any barrier value is ever computed -- never as a
    // floating-point subtraction of two independently-rounded nearly-equal
    // barrier values (which is numerically unstable near dhat -> 0, where
    // the barrier and its derivatives blow up).
    template <typename KeyType, typename ValueType>
    void
    insert_pair(unordered_map<KeyType, ValueType>& map, ValueType&& collision)
    {
        if (auto iter = map.find(collision->get_typed_hash());
            iter != map.end()) {
            iter->second->weight += collision->weight;
            if (iter->second->weight == 0) {
                map.erase(iter);
            }
        } else {
            map[collision->get_typed_hash()] = std::move(collision);
        }
    }

    // 2D counterpart of HighOrderCollisionsBuilder<3>::
    // reduce_point_edge_collision (high_order_collisions_builder.cpp), which
    // only exists on the <3> specialization: classify which sub-feature of
    // edge ei the query's closest point falls on and emit the
    // correspondingly-typed collision, so endpoint cases hash-merge with the
    // direct vertex collisions instead of double-counting a corner.
    //
    // Query vertex first in both templates, matching the convention of the
    // 2D edge-QP builder in quadrature_potential.cpp; Vertex2-Edge2P1's
    // evaluators assume that layout ([q, e0, e1]).
    std::shared_ptr<HighOrderCollision> reduce_point_edge_collision_2d(
        const index_t ei,
        const index_t vid,
        const HighOrderContactParameters& params,
        const CollisionMesh& mesh,
        const VertexMatrixView<2>& vertices)
    {
        const index_t e0 = mesh.edges()(ei, 0);
        const index_t e1 = mesh.edges()(ei, 1);

        const PointEdgeDistanceType dtype =
            point_edge_distance_type(vertices(vid), vertices(e0), vertices(e1));

        const double dist_sqr = point_edge_distance(
            vertices(vid), vertices(e0), vertices(e1), dtype);
        if (dist_sqr >= params.dhat * params.dhat) {
            return nullptr;
        }

        switch (dtype) {
        case PointEdgeDistanceType::P_E0:
            return std::make_shared<
                HighOrderCollisionTemplate<Vertex2, Vertex2>>(vid, e0, mesh);
        case PointEdgeDistanceType::P_E1:
            return std::make_shared<
                HighOrderCollisionTemplate<Vertex2, Vertex2>>(vid, e1, mesh);
        case PointEdgeDistanceType::P_E:
            return std::make_shared<
                HighOrderCollisionTemplate<Vertex2, Edge2P1>>(vid, ei, mesh);
        default:
            assert(false);
            return nullptr;
        }
    }

} // namespace

template <int dim>
ArbitraryPointPotential<dim>::ArbitraryPointPotential(
    const CollisionMesh& _mesh, HighOrderContactParameters _params)
    : mesh(_mesh)
    , params(std::move(_params))
{
    if (mesh.dim() != dim) {
        log_and_throw_error(
            "ArbitraryPointPotential<{}> requires a {}D mesh (got {}D)!", dim,
            dim, mesh.dim());
    }
}

template <int dim>
void ArbitraryPointPotential<dim>::update(Eigen::ConstRef<Eigen::MatrixXd> V)
{
    point_bvh.update(V, mesh);
}

template <int dim>
std::unique_ptr<HighOrderCollisionDict<PointType::VERTEX, dim>>
ArbitraryPointPotential<dim>::build_collisions_at_point(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
{
    using VertexP = std::conditional_t<dim == 2, Vertex2, Vertex3>;

    const index_t vid = static_cast<index_t>(V.rows()); // virtual vertex id
    const VertexMatrixView<dim> V_view(V, q);

    std::vector<index_t> vertex_ids, edge_ids, face_ids;
    point_bvh.query_point(q, params.dhat, vertex_ids, edge_ids, face_ids);

    unordered_map<std::array<index_t, 3>, std::shared_ptr<HighOrderCollision>>
        pairs;

    // Inclusion-exclusion over codimension: every primitive whose offset
    // region can contain q contributes a term signed (-1)^(codim-1) -- in 3D
    // faces +1, edges -1, vertices +1; in 2D edges +1, vertices -1. Each
    // codim-1 primitive is first *reduced* to the sub-feature its closest
    // point to q actually lies on (never just "this face's interior"
    // regardless of where the closest point falls), so redundant terms
    // converging on the same feature share a typed hash and cancel as
    // integers in insert_pair() above.
    if constexpr (dim == 3) {
        for (const index_t fi : face_ids) {
            if (std::shared_ptr<HighOrderCollision> pair =
                    HighOrderCollisionsBuilder<3>::
                        reduce_point_triangle_collision(
                            FaceVertexCandidate(fi, vid), params, mesh,
                            V_view)) {
                // Weight stays at the class default (+1).
                insert_pair(pairs, std::move(pair));
            }
        }
        for (const index_t ei : edge_ids) {
            if (std::shared_ptr<HighOrderCollision> pair =
                    HighOrderCollisionsBuilder<3>::reduce_point_edge_collision(
                        EdgeVertexCandidate(ei, vid), params, mesh, V_view)) {
                pair->weight = -1;
                insert_pair(pairs, std::move(pair));
            }
        }
    } else {
        // In 2D edges are the codim-1 primitives, so they take the +1 the
        // faces take in 3D and there is no face loop (mesh.faces() is empty
        // and the face BVH is never built).
        for (const index_t ei : edge_ids) {
            if (std::shared_ptr<HighOrderCollision> pair =
                    reduce_point_edge_collision_2d(
                        ei, vid, params, mesh, V_view)) {
                insert_pair(pairs, std::move(pair));
            }
        }
    }

    // Vertices: the highest codimension, so +1 in 3D and -1 in 2D. These merge
    // (and symbolically cancel) with the vertex-typed collisions the
    // reductions above emit when both resolve to the same corner.
    //
    // The (query, mesh vertex) argument order is load-bearing in 2D and only
    // in 2D: get_typed_hash() is {type, primitive_a.id(), primitive_b.id()},
    // and Vertex2-Vertex2 goes through the generic constructor, which stores
    // the ids as given -- so this has to match what
    // reduce_point_edge_collision_2d emits or the two never merge.
    // Vertex3-Vertex3 has a specialized constructor that sorts its two ids
    // (high_order_collision_template.cpp), so 3D merges either way.
    for (const index_t vi : vertex_ids) {
        if ((V.row(vi) - q).squaredNorm() >= params.dhat * params.dhat) {
            continue;
        }
        std::shared_ptr<HighOrderCollision> pair =
            std::make_shared<HighOrderCollisionTemplate<VertexP, VertexP>>(
                vid, vi, mesh);
        if constexpr (dim == 2) {
            pair->weight = -1;
        }
        insert_pair(pairs, std::move(pair));
    }

    auto collisions =
        std::make_unique<HighOrderCollisionDict<PointType::VERTEX, dim>>();
    collisions->initialize(
        std::vector<index_t> { vid }, std::vector<index_t> { vid }, pairs);
    return collisions;
}

template <int dim>
double ArbitraryPointPotential<dim>::operator()(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
{
    const auto collisions = build_collisions_at_point(V, q);
    const VertexMatrixView<dim> V_view(V, q);

    double value = 0.0;
    for (int ci = 0; ci < collisions->size(); ci++) {
        const auto& cc = (*collisions)[ci];
        value += cc.weight * cc(cc.dof(V_view), params, /*adaptive=*/nullptr);
    }
    return value;
}

template <int dim>
auto ArbitraryPointPotential<dim>::gradient(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
    -> Gradient
{
    const auto collisions = build_collisions_at_point(V, q);
    const VertexMatrixView<dim> V_view(V, q);
    const index_t vid = static_cast<index_t>(V.rows());

    Eigen::VectorXd grad =
        Eigen::VectorXd::Zero(collisions->vertex_ids().size() * dim);
    for (int ci = 0; ci < collisions->size(); ci++) {
        const auto& cc = (*collisions)[ci];
        const Eigen::VectorXd g = cc.weight
            * cc.gradient(cc.dof(V_view), params, /*adaptive=*/nullptr);
        for (int j = 0; j < cc.num_vertices(); j++) {
            grad.template segment<dim>(
                dim * collisions->vertex_ids_inverse(cc.vertex_id(j))) +=
                g.template segment<dim>(dim * j);
        }
    }

    const index_t q_local = collisions->vertex_ids_inverse(vid);
    return grad.template segment<dim>(dim * q_local);
}

template <int dim>
auto ArbitraryPointPotential<dim>::hessian(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
    -> Hessian
{
    const auto collisions = build_collisions_at_point(V, q);
    const VertexMatrixView<dim> V_view(V, q);
    const index_t vid = static_cast<index_t>(V.rows());

    const int m = static_cast<int>(collisions->vertex_ids().size());
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m * dim, m * dim);
    for (int ci = 0; ci < collisions->size(); ci++) {
        const auto& cc = (*collisions)[ci];
        const Eigen::MatrixXd h =
            cc.hessian(cc.dof(V_view), params, /*adaptive=*/nullptr)
            * cc.weight;
        for (int i = 0; i < cc.num_vertices(); i++) {
            for (int j = 0; j < cc.num_vertices(); j++) {
                H.template block<dim, dim>(
                    dim * collisions->vertex_ids_inverse(cc.vertex_id(i)),
                    dim * collisions->vertex_ids_inverse(cc.vertex_id(j))) +=
                    h.template block<dim, dim>(dim * i, dim * j);
            }
        }
    }

    const index_t q_local = collisions->vertex_ids_inverse(vid);
    return H.template block<dim, dim>(dim * q_local, dim * q_local);
}

template <int dim>
auto ArbitraryPointPotential<dim>::evaluate(
    Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Point> q) const
    -> std::tuple<double, Gradient, Hessian>
{
    const auto collisions = build_collisions_at_point(V, q);
    const VertexMatrixView<dim> V_view(V, q);
    const index_t vid = static_cast<index_t>(V.rows());

    double value = 0.0;
    const int m = static_cast<int>(collisions->vertex_ids().size());
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(m * dim);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m * dim, m * dim);

    for (int ci = 0; ci < collisions->size(); ci++) {
        const auto& cc = (*collisions)[ci];
        const auto dof = cc.dof(V_view);

        value += cc.weight * cc(dof, params, /*adaptive=*/nullptr);

        const Eigen::VectorXd g =
            cc.weight * cc.gradient(dof, params, /*adaptive=*/nullptr);
        const Eigen::MatrixXd h =
            cc.weight * cc.hessian(dof, params, /*adaptive=*/nullptr);

        for (int i = 0; i < cc.num_vertices(); i++) {
            const index_t gi = collisions->vertex_ids_inverse(cc.vertex_id(i));
            grad.template segment<dim>(dim * gi) +=
                g.template segment<dim>(dim * i);
            for (int j = 0; j < cc.num_vertices(); j++) {
                const index_t gj =
                    collisions->vertex_ids_inverse(cc.vertex_id(j));
                H.template block<dim, dim>(dim * gi, dim * gj) +=
                    h.template block<dim, dim>(dim * i, dim * j);
            }
        }
    }

    const index_t q_local = collisions->vertex_ids_inverse(vid);
    return { value, grad.template segment<dim>(dim * q_local),
             H.template block<dim, dim>(dim * q_local, dim * q_local) };
}

template class ArbitraryPointPotential<2>;
template class ArbitraryPointPotential<3>;

} // namespace ipc
