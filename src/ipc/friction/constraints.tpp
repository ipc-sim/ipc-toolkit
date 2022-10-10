#pragma once
#include "constraints.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {

template <typename T>
T FrictionConstraints::compute_potential(
    const CollisionMesh& mesh,
    const MatrixX<T>& velocity,
    double epsv_times_h) const
{
    if (empty()) {
        return T(0);
    }
    assert(epsv_times_h > 0);

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    tbb::enumerable_thread_specific<T> storage(0);

    const size_t num_vv = vv_constraints.size();
    const size_t num_ev = ev_constraints.size();
    const size_t num_ee = ee_constraints.size();
    const size_t num_fv = fv_constraints.size();

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), size()),
        [&](tbb::blocked_range<size_t> r) {
            auto& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {

                size_t ci = i;
                if (ci < num_vv) {
                    local_potential += vv_constraints[ci].compute_potential(
                        velocity, E, F, epsv_times_h);
                } else if ((ci -= num_vv) < num_ev) {
                    local_potential += ev_constraints[ci].compute_potential(
                        velocity, E, F, epsv_times_h);
                } else if ((ci -= num_ev) < num_ee) {
                    local_potential += ee_constraints[ci].compute_potential(
                        velocity, E, F, epsv_times_h);
                } else {
                    ci -= num_ee;
                    assert(ci < num_fv);
                    local_potential += fv_constraints[ci].compute_potential(
                        velocity, E, F, epsv_times_h);
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
