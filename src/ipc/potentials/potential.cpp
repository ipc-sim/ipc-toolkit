#include "potential.hpp"

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/utils/hessian_assembler.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/meshfem_hessian_assembler.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <atomic>
#include <vector>

namespace ipc {

namespace {

    /// @brief Gather-based global gradient assembly.
    ///
    /// Sums per-stencil gradient contributions into the global gradient by
    /// gathering per *vertex* instead of scattering per stencil: a
    /// vertex→slot adjacency is built with a parallel counting sort, then
    /// each vertex sums its contributions independently. This avoids the
    /// dense per-thread accumulators of a scatter+reduce (whose zero+combine
    /// cost grows with ndof, not with the number of collisions).
    ///
    /// @param out_ndof Size of the output gradient.
    /// @param dim Spatial dimension.
    /// @param local_grads Flat buffer of slot contributions; slot s occupies
    ///     [dim·s, dim·(s+1)). Slots with an invalid vertex are never read.
    /// @param slot_vertex Global vertex of each slot (negative = unused).
    Eigen::VectorXd gather_global_gradient(
        const int out_ndof,
        const int dim,
        Eigen::ConstRef<Eigen::VectorXd> local_grads,
        const std::vector<index_t>& slot_vertex)
    {
        const size_t num_slots = slot_vertex.size();
        const size_t out_num_verts = size_t(out_ndof) / dim;
        assert(local_grads.size() == Eigen::Index(dim * num_slots));

        // --- Vertex → slot adjacency (parallel counting sort) -----------
        std::vector<std::atomic<std::ptrdiff_t>> cursor(out_num_verts);
        tbb::parallel_for(size_t(0), out_num_verts, [&](const size_t v) {
            cursor[v].store(0, std::memory_order_relaxed);
        });
        tbb::parallel_for(size_t(0), num_slots, [&](const size_t s) {
            if (slot_vertex[s] >= 0) {
                cursor[slot_vertex[s]].fetch_add(1, std::memory_order_relaxed);
            }
        });

        std::vector<std::ptrdiff_t> bucket_start(out_num_verts + 1);
        bucket_start[0] = 0;
        for (size_t v = 0; v < out_num_verts; v++) {
            bucket_start[v + 1] =
                bucket_start[v] + cursor[v].load(std::memory_order_relaxed);
            cursor[v].store(bucket_start[v], std::memory_order_relaxed);
        }

        std::vector<std::ptrdiff_t> bucket_slots(bucket_start.back());
        tbb::parallel_for(size_t(0), num_slots, [&](const size_t s) {
            if (slot_vertex[s] >= 0) {
                bucket_slots[cursor[slot_vertex[s]].fetch_add(
                    1, std::memory_order_relaxed)] = s;
            }
        });

        // --- Per-vertex gather -------------------------------------------
        Eigen::VectorXd grad(out_ndof);
        tbb::parallel_for(size_t(0), out_num_verts, [&](const size_t v) {
            VectorMax3d acc = VectorMax3d::Zero(dim);
            for (auto bi = bucket_start[v]; bi < bucket_start[v + 1]; bi++) {
                acc += local_grads.segment(dim * bucket_slots[bi], dim);
            }

            if constexpr (VERTEX_DERIVATIVE_LAYOUT == Eigen::RowMajor) {
                grad.segment(dim * v, dim) = acc;
            } else {
                for (int d = 0; d < dim; d++) {
                    grad[out_num_verts * d + v] = acc[d];
                }
            }
        });
        return grad;
    }

} // namespace

template <class TCollisions>
double Potential<TCollisions>::operator()(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X) const
{
    assert(X.rows() == mesh.num_vertices());
    IPC_TOOLKIT_PROFILE_BLOCK("Potential<T>::operator()");

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double partial_sum) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                // Quadrature weight is premultiplied by local potential
                partial_sum += (*this)(
                    collisions[i],
                    collisions[i].dof(X, mesh.edges(), mesh.faces()));
            }
            return partial_sum;
        },
        std::plus<double>());
}

template <class TCollisions>
Eigen::VectorXd Potential<TCollisions>::gradient(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const bool in_full_dof) const
{
    assert(X.rows() == mesh.num_vertices());
    IPC_TOOLKIT_PROFILE_BLOCK("Potential<T>::gradient()");

    // Assemble directly in full-mesh DOF when the DOF map is a pure selection
    // (remapping stencil vertex IDs is then equivalent to to_full_dof());
    // otherwise assemble in collision DOF and apply the map at the end.
    const bool fold_to_full = in_full_dof && mesh.is_selection_dof_map();
    const bool map_to_full = in_full_dof && !fold_to_full;

    const int out_ndof = fold_to_full ? mesh.full_ndof() : X.size();

    if (collisions.empty()) {
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(out_ndof);
        if (map_to_full) {
            return mesh.to_full_dof(grad);
        }
        return grad;
    }

    const int dim = X.cols();
    const size_t num_collisions = collisions.size();
    constexpr int STENCIL_SIZE = TCollision::STENCIL_SIZE;
    const size_t num_slots = STENCIL_SIZE * num_collisions;

    // Two assembly strategies with complementary scaling:
    //  - gather: per-vertex reduction through a slot adjacency. Cost scales
    //    with the number of contributions.
    //  - scatter+reduce: thread-local dense accumulators. Cost scales with
    //    out_ndof (zeroing + combining a gradient-sized vector per thread).
    // The crossover is at out_ndof ≈ num_slots on the benchmarked scenes.
    if (size_t(out_ndof) > num_slots) {
        // Evaluate the local gradients into a flat buffer of per-vertex
        // slots: vertex k of collision i is slot s = STENCIL_SIZE·i + k,
        // with its dim gradient entries at local_grads[dim·s].
        Eigen::VectorXd local_grads(dim * num_slots);
        std::vector<index_t> slot_vertex(num_slots);

        {
            IPC_TOOLKIT_PROFILE_BLOCK("Compute Local Gradients");
            tbb::parallel_for(size_t(0), num_collisions, [&](size_t i) {
                const TCollision& collision = collisions[i];

                const VectorMaxNd local_grad = this->gradient(
                    collision, collision.dof(X, mesh.edges(), mesh.faces()));
                local_grads.segment(STENCIL_SIZE * dim * i, local_grad.size()) =
                    local_grad;

                const auto ids =
                    collision.vertex_ids(mesh.edges(), mesh.faces());
                const int n_verts = local_grad.size() / dim;
                for (int k = 0; k < STENCIL_SIZE; k++) {
                    index_t id = k < n_verts ? ids[k] : index_t(-1);
                    if (fold_to_full && id >= 0) {
                        id = mesh.to_full_vertex_id(id);
                    }
                    slot_vertex[STENCIL_SIZE * i + k] = id;
                }
            });
        }

        IPC_TOOLKIT_PROFILE_BLOCK("Gather Local Gradients");
        Eigen::VectorXd grad =
            gather_global_gradient(out_ndof, dim, local_grads, slot_vertex);
        if (map_to_full) {
            return mesh.to_full_dof(grad);
        }
        return grad; // implicitly moved
    }

    tbb::combinable<Eigen::VectorXd> grad(Eigen::VectorXd::Zero(out_ndof));

    {
        IPC_TOOLKIT_PROFILE_BLOCK("Compute Local Gradients");
        tbb::parallel_for(size_t(0), num_collisions, [&](size_t i) {
            const TCollision& collision = collisions[i];

            const VectorMaxNd local_grad = this->gradient(
                collision, collision.dof(X, mesh.edges(), mesh.faces()));

            auto ids = collision.vertex_ids(mesh.edges(), mesh.faces());
            if (fold_to_full) {
                for (auto& id : ids) {
                    if (id >= 0) {
                        id = mesh.to_full_vertex_id(id);
                    }
                }
            }

            local_gradient_to_global_gradient(
                local_grad, ids, dim, grad.local());
        });
    }

    {
        IPC_TOOLKIT_PROFILE_BLOCK("Combine Local Gradients");
        Eigen::VectorXd combined_grad = grad.combine(
            [](const Eigen::VectorXd& a,
               const Eigen::VectorXd& b) -> Eigen::VectorXd { return a + b; });
        if (map_to_full) {
            return mesh.to_full_dof(combined_grad);
        }
        return combined_grad; // implicitly moved
    }
}

template <class TCollisions>
Eigen::SparseMatrix<double> Potential<TCollisions>::hessian(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const PSDProjectionMethod project_hessian_to_psd,
    const bool in_full_dof) const
{
    IPC_TOOLKIT_PROFILE_BLOCK("Potential<T>::hessian()");

    // Assemble directly in full-mesh DOF when the DOF map is a pure selection
    // (remapping stencil vertex IDs is then equivalent to to_full_dof());
    // otherwise assemble in collision DOF and apply the map at the end.
    const bool fold_to_full = in_full_dof && mesh.is_selection_dof_map();
    const bool map_to_full = in_full_dof && !fold_to_full;

    if (collisions.empty()) {
        // Short-circuit: building a sparsity pattern for an empty contact set
        // costs O(ndof) (or more) work to produce an all-zero matrix.
        assert(X.rows() == mesh.num_vertices());
        const int out_ndof = fold_to_full ? mesh.full_ndof() : X.size();
        Eigen::SparseMatrix<double> hess(out_ndof, out_ndof);
        if (map_to_full) {
            return mesh.to_full_dof(hess);
        }
        return hess;
    }

#ifdef IPC_TOOLKIT_WITH_MESHFEM_SPARSE
    MeshFEMHessianAssembler assembler;
    assemble_hessian(
        collisions, mesh, X, assembler, project_hessian_to_psd, fold_to_full);
    Eigen::SparseMatrix<double> hess = assembler.take_matrix();
#else
    TripletHessianAssembler assembler;
    assemble_hessian(
        collisions, mesh, X, assembler, project_hessian_to_psd, fold_to_full);
    Eigen::SparseMatrix<double> hess = assembler.get_matrix();
#endif

    if (map_to_full) {
        return mesh.to_full_dof(hess);
    }
    return hess; // implicitly moved
}

template <class TCollisions>
void Potential<TCollisions>::assemble_hessian(
    const TCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    HessianAssembler& assembler,
    const PSDProjectionMethod project_hessian_to_psd,
    const bool in_full_dof) const
{
    assert(X.rows() == mesh.num_vertices());
    IPC_TOOLKIT_PROFILE_BLOCK("Potential<T>::assemble_hessian()");

    // The HessianAssembler interface is sized for the universal collision
    // stencil (≤ 4 vertices ⇒ ≤ 12 DOF).
    static_assert(TCollision::STENCIL_SIZE == 4);
    static_assert(std::is_same_v<MatrixMaxNd, MatrixMax12d>);

    if (in_full_dof && !mesh.is_selection_dof_map()) {
        log_and_throw_error(
            "assemble_hessian: in_full_dof requires the mesh's DOF map to be "
            "a pure selection (see CollisionMesh::is_selection_dof_map); "
            "assemble in collision DOF and apply to_full_dof instead.");
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    const int dim = X.cols();
    const int ndof = in_full_dof ? mesh.full_ndof() : X.size(); // NOLINT

    // Stencil vertex IDs as assembled (remapped to full-mesh IDs if folding).
    const auto stencil_vertex_ids = [&](const size_t i) {
        auto ids = collisions[i].vertex_ids(edges, faces);
        if (in_full_dof) {
            for (auto& id : ids) {
                if (id >= 0) {
                    id = mesh.to_full_vertex_id(id);
                }
            }
        }
        return ids;
    };

    assembler.begin(ndof, dim, collisions.size(), stencil_vertex_ids);

    {
        IPC_TOOLKIT_PROFILE_BLOCK("compute and assemble local hessians");
        tbb::parallel_for(size_t(0), collisions.size(), [&](size_t i) {
            const TCollision& collision = collisions[i];

            MatrixMaxNd local_hess;
            {
                IPC_TOOLKIT_PROFILE_BLOCK("Compute Local Hessian");
                local_hess = this->hessian(
                    collision, collision.dof(X, edges, faces),
                    project_hessian_to_psd);
            }

            assembler.add_local_hessian(local_hess, stencil_vertex_ids(i));
        });
    }

    assembler.end();
}

template class Potential<NormalCollisions>;
template class Potential<TangentialCollisions>;

} // namespace ipc