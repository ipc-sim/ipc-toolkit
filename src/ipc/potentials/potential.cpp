#include "potential.hpp"

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/utils/gradient_assembler.hpp>
#include <ipc/utils/hessian_assembler.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/meshfem_hessian_assembler.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <array>

namespace ipc {

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

    // Local gradient of collision i (quadrature weight premultiplied).
    const auto local_gradient = [&](const size_t i) -> VectorMaxNd {
        const TCollision& collision = collisions[i];
        return this->gradient(
            collision, collision.dof(X, mesh.edges(), mesh.faces()));
    };

    // Stencil vertex ids of collision i (remapped to full ids if folding).
    using IDArray = std::array<index_t, TCollision::STENCIL_SIZE>;
    const auto stencil_ids = [&](const size_t i) -> IDArray {
        IDArray ids = collisions[i].vertex_ids(mesh.edges(), mesh.faces());
        if (fold_to_full) {
            for (auto& id : ids) {
                if (id >= 0) {
                    id = mesh.to_full_vertex_id(id);
                }
            }
        }
        return ids;
    };

    Eigen::VectorXd grad = assemble_gradient(
        out_ndof, dim, num_collisions, local_gradient, stencil_ids);

    if (map_to_full) {
        return mesh.to_full_dof(grad);
    }

    return grad;
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