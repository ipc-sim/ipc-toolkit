#include "rigid_contact_hessian.hpp"

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <array>

namespace ipc::rigid {

template <class TCollisions, class TPotential>
Eigen::SparseMatrix<double> per_contact_rigid_hessian(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses,
    const TCollisions& collisions,
    const TPotential& potential,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const double gn_scale,
    const double geom_scale,
    const PSDProjectionMethod project)
{
    const int dim = bodies.dim();
    const int rb_ndof = poses.empty() ? 0 : poses[0].ndof();
    const int n = int(poses.size()) * rb_ndof;
    const auto& E = bodies.edges();
    const auto& F = bodies.faces();

    if (collisions.size() == 0) {
        return Eigen::SparseMatrix<double>(n, n);
    }

    using Triplets = std::vector<Eigen::Triplet<double>>;
    tbb::enumerable_thread_specific<Triplets> storage;

    tbb::parallel_for(size_t(0), collisions.size(), [&](size_t ci) {
        Triplets& triplets = storage.local();
        const auto& c = collisions[ci];
        const int nv = c.num_vertices();

        // Vertex-space stencil derivatives; the interior stencil Hessian is
        // left UNPROJECTED (mirrors Rigid IPC's compute_potential_hessian with
        // project_hessian_to_psd=false) — only the completed rigid-DOF block is
        // projected below.
        const auto s_dof = c.dof(X, E, F);
        const auto grad_f = potential.gradient(c, s_dof);
        const auto hess_f =
            potential.hessian(c, s_dof, PSDProjectionMethod::NONE);

        const std::array<index_t, 4> vids = c.vertex_ids(E, F);

        // Map each stencil vertex to its body; collect the ≤2 unique bodies.
        std::array<index_t, 2> body_ids = { { -1, -1 } };
        std::array<int, 4> local_slot = { { 0, 0, 0, 0 } };
        int nb = 0;
        for (int i = 0; i < nv; ++i) {
            const index_t b = bodies.vertex_to_body(vids[i]);
            int slot = -1;
            for (int k = 0; k < nb; ++k) {
                if (body_ids[k] == b) {
                    slot = k;
                    break;
                }
            }
            if (slot < 0) {
                slot = nb;
                body_ids[nb++] = b;
            }
            local_slot[i] = slot;
        }

        // Jacobian of the stencil vertices w.r.t. the (≤2) bodies' rigid DOFs.
        // A stencil has at most 4 vertices × 3 coordinates = 12 rows, and the
        // ≤2 bodies have at most 2 × 6 rigid DOFs = 12 columns, so every local
        // matrix here fits in a stack-allocated MatrixMax12d.
        MatrixMax12d jac_Vi = MatrixMax12d::Zero(nv * dim, nb * rb_ndof);
        for (int i = 0; i < nv; ++i) {
            const index_t b = bodies.vertex_to_body(vids[i]);
            const index_t lv = vids[i] - bodies.body_vertex_start(b);
            jac_Vi.block(i * dim, local_slot[i] * rb_ndof, dim, rb_ndof) =
                poses[b].transform_vertex_jacobian(
                    bodies.body_rest_positions(b).row(lv).transpose());
        }

        // Gauss–Newton term (keeps the two-body coupling).
        MatrixMax12d local = gn_scale * (jac_Vi.transpose() * hess_f * jac_Vi);

        // Second-order geometric term on the diagonal body blocks only (a
        // vertex of body A has zero second derivative w.r.t. body B's DOFs).
        for (int i = 0; i < nv; ++i) {
            const auto gi = grad_f.segment(i * dim, dim);
            if (gi.isZero()) {
                continue;
            }
            const index_t b = bodies.vertex_to_body(vids[i]);
            const index_t lv = vids[i] - bodies.body_vertex_start(b);
            const MatrixMax<double, 3, 36> d2V =
                poses[b].transform_vertex_hessian(
                    bodies.body_rest_positions(b).row(lv).transpose());
            const int off = local_slot[i] * rb_ndof;
            for (int jj = 0; jj < rb_ndof; ++jj) {
                for (int ii = 0; ii < rb_ndof; ++ii) {
                    local(off + ii, off + jj) +=
                        geom_scale * gi.dot(d2V.col(ii * rb_ndof + jj));
                }
            }
        }

        // Project the completed local rigid-DOF block (mirrors
        // apply_chain_rule).
        if (project != PSDProjectionMethod::NONE) {
            local = project_to_psd(local, project);
        }

        // Scatter into the global rigid-DOF Hessian.
        for (int a = 0; a < nb; ++a) {
            for (int b2 = 0; b2 < nb; ++b2) {
                for (int r = 0; r < rb_ndof; ++r) {
                    for (int col = 0; col < rb_ndof; ++col) {
                        const double v =
                            local(a * rb_ndof + r, b2 * rb_ndof + col);
                        if (v != 0.0) {
                            triplets.emplace_back(
                                body_ids[a] * rb_ndof + r,
                                body_ids[b2] * rb_ndof + col, v);
                        }
                    }
                }
            }
        }
    });

    Triplets all;
    size_t total = 0;
    for (const auto& t : storage) {
        total += t.size();
    }
    all.reserve(total);
    for (const auto& t : storage) {
        all.insert(all.end(), t.begin(), t.end());
    }

    Eigen::SparseMatrix<double> H(n, n);
    H.setFromTriplets(all.begin(), all.end());
    return H;
}

// Explicit instantiations for the two contact potentials used by the demo.
template Eigen::SparseMatrix<double>
per_contact_rigid_hessian<NormalCollisions, BarrierPotential>(
    const RigidBodies&,
    const std::vector<Pose>&,
    const NormalCollisions&,
    const BarrierPotential&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const double,
    const double,
    const PSDProjectionMethod);

template Eigen::SparseMatrix<double>
per_contact_rigid_hessian<TangentialCollisions, FrictionPotential>(
    const RigidBodies&,
    const std::vector<Pose>&,
    const TangentialCollisions&,
    const FrictionPotential&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const double,
    const double,
    const PSDProjectionMethod);

} // namespace ipc::rigid
