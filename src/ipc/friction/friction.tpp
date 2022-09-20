#pragma once
#include "friction.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {

template <typename T>
T compute_friction_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const MatrixX<T>& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h)
{
    if (friction_constraint_set.empty()) {
        return T(0);
    }
    assert(epsv_times_h > 0);

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    // absolute linear dislacement of each point
    MatrixX<T> U = V1 - V0.cast<T>();

    tbb::enumerable_thread_specific<T> storage(0);

    const size_t num_vv = friction_constraint_set.vv_constraints.size();
    const size_t num_ev = friction_constraint_set.ev_constraints.size();
    const size_t num_ee = friction_constraint_set.ee_constraints.size();
    const size_t num_fv = friction_constraint_set.fv_constraints.size();

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), friction_constraint_set.size()),
        [&](tbb::blocked_range<size_t> r) {
            auto& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {

                size_t ci = i;
                if (ci < num_vv) {
                    local_potential +=
                        friction_constraint_set.vv_constraints[ci]
                            .compute_potential(U, E, F, epsv_times_h);
                } else if ((ci -= num_vv) < num_ev) {
                    local_potential +=
                        friction_constraint_set.ev_constraints[ci]
                            .compute_potential(U, E, F, epsv_times_h);
                } else if ((ci -= num_ev) < num_ee) {
                    local_potential +=
                        friction_constraint_set.ee_constraints[ci]
                            .compute_potential(U, E, F, epsv_times_h);
                } else {
                    ci -= num_ee;
                    assert(ci < num_fv);
                    local_potential +=
                        friction_constraint_set.fv_constraints[ci]
                            .compute_potential(U, E, F, epsv_times_h);
                }
            }
        });

    T potential(0);
    for (const auto& local_potential : storage) {
        potential += local_potential;
    }
    return potential;
}

} // namespace ipc
