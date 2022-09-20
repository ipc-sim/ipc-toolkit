#include "plane.hpp"

#include <ipc/distance/point_plane.hpp>
#include <ipc/ccd/point_static_plane.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {

void construct_point_plane_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& plane_origins,
    const Eigen::MatrixXd& plane_normals,
    const double dhat,
    std::vector<PlaneVertexConstraint>& pv_constraints,
    const double dmin,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    pv_constraints.clear();

    double dhat_squared = dhat * dhat;
    double dmin_squared = dmin * dmin;

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.

    size_t n_planes = plane_origins.rows();
    assert(plane_normals.rows() == n_planes);

    for (size_t vi = 0; vi < V.rows(); vi++) {
        for (size_t pi = 0; pi < n_planes; pi++) {
            if (!can_collide(vi, pi)) {
                continue;
            }

            const auto& plane_origin = plane_origins.row(pi);
            const auto& plane_normal = plane_normals.row(pi);

            double distance_sqr =
                point_plane_distance(V.row(vi), plane_origin, plane_normal);

            if (distance_sqr - dmin_squared < 2 * dmin * dhat + dhat_squared) {
                pv_constraints.emplace_back(plane_origin, plane_normal, vi);
                pv_constraints.back().minimum_distance = dmin;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

bool is_step_point_plane_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXd& plane_origins,
    const Eigen::MatrixXd& plane_normals,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    size_t n_planes = plane_origins.rows();
    assert(plane_normals.rows() == n_planes);
    assert(V0.rows() == V1.rows());

    for (size_t vi = 0; vi < V0.rows(); vi++) {
        for (size_t pi = 0; pi < n_planes; pi++) {
            if (!can_collide(vi, pi)) {
                continue;
            }

            const auto& plane_origin = plane_origins.row(pi);
            const auto& plane_normal = plane_normals.row(pi);

            double toi;
            bool is_collision = point_static_plane_ccd(
                V0.row(vi), V1.row(vi), plane_origin, plane_normal, toi);

            if (is_collision) {
                return false;
            }
        }
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////

double compute_point_plane_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXd& plane_origins,
    const Eigen::MatrixXd& plane_normals,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    size_t n_planes = plane_origins.rows();
    assert(plane_normals.rows() == n_planes);
    assert(V0.rows() == V1.rows());

    tbb::enumerable_thread_specific<double> storage(1);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, V0.rows()),
        [&](tbb::blocked_range<size_t> r) {
            double& earliest_toi = storage.local();

            for (size_t vi = r.begin(); vi < r.end(); vi++) {
                for (size_t pi = 0; pi < n_planes; pi++) {
                    if (!can_collide(vi, pi)) {
                        continue;
                    }

                    const auto& plane_origin = plane_origins.row(pi);
                    const auto& plane_normal = plane_normals.row(pi);

                    double toi;
                    bool are_colliding = point_static_plane_ccd(
                        V0.row(vi), V1.row(vi), plane_origin, plane_normal,
                        toi);

                    if (are_colliding) {
                        if (toi < earliest_toi) {
                            earliest_toi = toi;
                        }
                    }
                }
            }
        });

    double earliest_toi = 1;
    for (const auto& local_earliest_toi : storage) {
        earliest_toi = std::min(earliest_toi, local_earliest_toi);
    }
    assert(earliest_toi >= 0 && earliest_toi <= 1.0);
    return earliest_toi;
}

} // namespace ipc
