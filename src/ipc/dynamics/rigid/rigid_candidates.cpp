#include "rigid_candidates.hpp"

#include <ipc/config.hpp>
#include <ipc/ipc.hpp>
#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/dynamics/rigid/rigid_trajectory.hpp>
#include <ipc/io/write_candidates_obj.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <igl/remove_unreferenced.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <atomic>

namespace ipc::rigid {

void RigidCandidates::build(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses,
    const double inflation_radius,
    BroadPhase* broad_phase)
{
    build(bodies, poses, poses, inflation_radius, broad_phase);
}

void RigidCandidates::build(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double inflation_radius,
    BroadPhase* broad_phase)
{
    assert(poses_t0.size() == bodies.num_bodies());
    assert(poses_t1.size() == bodies.num_bodies());

    clear();

    // World-space vertex positions at the interval endpoints (for the
    // plane-vertex candidates below; body-body candidates work in relative
    // body space and never need world-space vertices).
    const Eigen::MatrixXd vertices_t0 = bodies.vertices(poses_t0);
    const Eigen::MatrixXd vertices_t1 = bodies.vertices(poses_t1);

    // Planes to vertices: the plane distance is affine along the chord, so it
    // suffices to check the endpoint distances inflated by the deviation bound.
    for (const auto& plane : bodies.planes) {
        for (index_t vi = 0; vi < bodies.num_vertices(); ++vi) {
            const index_t b = bodies.vertex_to_body(vi);
            const double deviation = rigid_chord_deviation(
                bodies.rest_positions().row(vi).norm(),
                (poses_t1[b].rotation - poses_t0[b].rotation).norm());
            const double d_t0 = plane.signedDistance(vertices_t0.row(vi));
            const double d_t1 = plane.signedDistance(vertices_t1.row(vi));
            if (d_t0 < inflation_radius + deviation
                || d_t1 < inflation_radius + deviation) {
                pv_candidates.emplace_back(plane, vi);
            }
        }
    }

    if (bodies.num_bodies() < 2) {
        return; // No body-body collisions possible
    }

    std::unique_ptr<BroadPhase> default_broad_phase;
    if (broad_phase == nullptr) {
        default_broad_phase = make_default_broad_phase();
        broad_phase = default_broad_phase.get();
    }

    const int dim = poses_t0[0].position.size();

    // 1. Broad phase between bodies

    // a. Build body AABBs (bounding sphere swept between the two positions)
    AABBs body_boxes(bodies.num_bodies());
    for (int i = 0; i < bodies.num_bodies(); i++) {
        const double r = bodies[i].bounding_radius() + inflation_radius;
        body_boxes[i] = AABB(
            AABB::from_point(poses_t0[i].position, r),
            AABB::from_point(poses_t1[i].position, r));
    }

    // b. Build the broad phase. The shared broad-phase object may carry a
    // *vertex*-level filter from a previous standard build (Candidates::build
    // sets can_vertices_collide to mesh.can_collide and leaves it set). Here
    // the "vertices" are bodies, so that filter must not be applied: feeding
    // body indices to a vertex filter silently drops body pairs (e.g., the
    // same-body vertex exclusion rejects low body indices because the
    // corresponding vertices belong to the same body). Save the caller's
    // filter, use an accept-all filter for the body-level pass, and restore
    // it afterwards.
    const CollisionFilter caller_filter = broad_phase->can_vertices_collide;
    // TODO: Implement a mechanism to allow the caller to specify a body-level
    //       filter (e.g., for jointed bodies).
    broad_phase->can_vertices_collide = CollisionFilter(); // accept all
    broad_phase->build(
        body_boxes, /*edges=*/Eigen::MatrixXi(), /*faces=*/Eigen::MatrixXi(),
        dim);

    // c. Detect body-body candidates: bodies are stored as "vertices," so we
    // can reuse vertex-vertex candidate detection.
    std::vector<VertexVertexCandidate> body_candidates;
    broad_phase->detect_vertex_vertex_candidates(body_candidates);
    broad_phase->can_vertices_collide = caller_filter;

    // 2. Broad phase between colliding bodies
    // NOTE: The loop over body pairs is serial because each two-tree query is
    // internally parallelized.
    for (const VertexVertexCandidate& body_candidate : body_candidates) {
        auto [body_i, body_j] = body_candidate; // Body indices

        // Ensure bi has more vertices than bj to take advantage of log(n)
        // traversal (bj's leaves are iterated against bi's tree).
        if (bodies.body_num_vertices(body_i)
            < bodies.body_num_vertices(body_j)) {
            std::swap(body_i, body_j);
        }

        // a. Build the boxes for body bj's swept trajectory inside the (rest)
        //    body space of body bi.
        const Pose body_j_to_i_t0 =
            poses_t0[body_i].inverse() * poses_t0[body_j];
        const Eigen::MatrixXd j_vertices_t0 =
            bodies.body_vertices(body_j, body_j_to_i_t0);

        const Pose body_j_to_i_t1 =
            poses_t1[body_i].inverse() * poses_t1[body_j];
        const Eigen::MatrixXd j_vertices_t1 =
            bodies.body_vertices(body_j, body_j_to_i_t1);

        // The relative trajectory is nonlinear, so the endpoint boxes are
        // inflated by a conservative bound on the deviation of the relative
        // curve from its chord (see relative_rigid_chord_deviation). Body i's
        // rest-frame BVH boxes are uninflated, so the inflation radius is
        // applied twice on this side (one-sided Minkowski growth).
        const double dtheta_i =
            (poses_t1[body_i].rotation - poses_t0[body_i].rotation).norm();
        const double dtheta_j =
            (poses_t1[body_j].rotation - poses_t0[body_j].rotation).norm();
        const VectorMax3d dq_t0 =
            poses_t0[body_j].position - poses_t0[body_i].position;
        const VectorMax3d dq_t1 =
            poses_t1[body_j].position - poses_t1[body_i].position;
        const double max_translation = std::max(dq_t0.norm(), dq_t1.norm());
        const double translation_change = (dq_t1 - dq_t0).norm();

        const auto j_rest_positions = bodies.body_rest_positions(body_j);

        AABBs body_j_vertex_boxes(bodies.body_num_vertices(body_j));
        for (size_t v = 0; v < body_j_vertex_boxes.size(); ++v) {
            const double deviation = relative_rigid_chord_deviation(
                j_rest_positions.row(v).norm(), dtheta_i, dtheta_j,
                max_translation, translation_change);
            body_j_vertex_boxes[v] = AABB::from_point(
                j_vertices_t0.row(v).transpose(),
                j_vertices_t1.row(v).transpose(),
                deviation + 2 * inflation_radius);
        }

        // b. Build a temporary BVH for body bj's swept boxes. Edge/face boxes
        //    are unions of the vertex boxes, which is conservative because
        //    every point of a moving edge/face is a convex combination of its
        //    vertices' trajectories.
        // TODO: We should skip this BVH build and directly pass over the boxes
        // to the two-tree traversal, but the current LBVH interface does not
        // support that.
        LBVH body_j_bvh;
        body_j_bvh.build(
            body_j_vertex_boxes, bodies.body_edges(body_j),
            bodies.body_faces(body_j), dim);

        // c. Detect candidates between body bi and bj using the static
        //    rest-frame BVH of body bi. Local per-body element ids are mapped
        //    to global collision mesh ids using the body start offsets.
        const LBVH& body_i_bvh = *bodies[body_i].bvh();

        const index_t vi_start = bodies.body_vertex_start(body_i);
        const index_t ei_start = bodies.body_edge_start(body_i);
        const index_t fi_start = bodies.body_face_start(body_i);
        const index_t vj_start = bodies.body_vertex_start(body_j);
        const index_t ej_start = bodies.body_edge_start(body_j);
        const index_t fj_start = bodies.body_face_start(body_j);

        // Vertex-level filters mirroring the single-tree traversal's
        // can_*_collide semantics, but with global ids so the collision mesh's
        // can_collide (e.g., a vertex-patch filter between jointed bodies) is
        // honored. Elements of distinct bodies never share endpoints, so the
        // shared-endpoint exclusions are unnecessary here.
        const Eigen::MatrixXi& E = bodies.edges();
        const Eigen::MatrixXi& F = bodies.faces();
        const auto can_vv = [&](const index_t va, const index_t vb) {
            return bodies.can_collide(va, vb);
        };
        const auto can_ev = [&](const index_t e, const index_t v) {
            return can_vv(E(e, 0), v) || can_vv(E(e, 1), v);
        };
        const auto can_ee = [&](const index_t ea, const index_t eb) {
            return can_ev(ea, E(eb, 0)) || can_ev(ea, E(eb, 1));
        };
        const auto can_fv = [&](const index_t f, const index_t v) {
            return can_vv(F(f, 0), v) || can_vv(F(f, 1), v)
                || can_vv(F(f, 2), v);
        };

        if (dim == 2) {
            // (eᵢ, vⱼ): receives (this_edge_id, other_vertex_id)
            std::vector<EdgeVertexCandidate> ev;
            body_i_bvh.detect_edge_vertex_candidates(
                body_j_bvh,
                [&](size_t ei, size_t vj) {
                    return can_ev(ei_start + ei, vj_start + vj);
                },
                ev);
            for (const EdgeVertexCandidate& c : ev) {
                ev_candidates.emplace_back(
                    ei_start + c.edge_id, vj_start + c.vertex_id);
            }

            // (eⱼ, vᵢ): receives (other_edge_id, this_vertex_id)
            ev.clear();
            body_i_bvh.detect_vertex_edge_candidates(
                body_j_bvh,
                [&](size_t ej, size_t vi) {
                    return can_ev(ej_start + ej, vi_start + vi);
                },
                ev);
            for (const EdgeVertexCandidate& c : ev) {
                ev_candidates.emplace_back(
                    ej_start + c.edge_id, vi_start + c.vertex_id);
            }
        } else {
            // (eᵢ, eⱼ): receives (this_edge_id, other_edge_id)
            std::vector<EdgeEdgeCandidate> ee;
            body_i_bvh.detect_edge_edge_candidates(
                body_j_bvh,
                [&](size_t ei, size_t ej) {
                    return can_ee(ei_start + ei, ej_start + ej);
                },
                ee);
            for (const EdgeEdgeCandidate& c : ee) {
                ee_candidates.emplace_back(
                    ei_start + c.edge0_id, ej_start + c.edge1_id);
            }

            // (fᵢ, vⱼ): receives (this_face_id, other_vertex_id)
            std::vector<FaceVertexCandidate> fv;
            body_i_bvh.detect_face_vertex_candidates(
                body_j_bvh,
                [&](size_t fi, size_t vj) {
                    return can_fv(fi_start + fi, vj_start + vj);
                },
                fv);
            for (const FaceVertexCandidate& c : fv) {
                fv_candidates.emplace_back(
                    fi_start + c.face_id, vj_start + c.vertex_id);
            }

            // (fⱼ, vᵢ): receives (other_face_id, this_vertex_id)
            fv.clear();
            body_i_bvh.detect_vertex_face_candidates(
                body_j_bvh,
                [&](size_t fj, size_t vi) {
                    return can_fv(fj_start + fj, vi_start + vi);
                },
                fv);
            for (const FaceVertexCandidate& c : fv) {
                fv_candidates.emplace_back(
                    fj_start + c.face_id, vi_start + c.vertex_id);
            }
        }

        // Codimensional candidates, mirroring the single-mesh logic of
        // Candidates::build: codim-vertex × codim-vertex pairs (narrow-phased
        // with point-point CCD) and, in 3D, codim-edge × codim-vertex pairs
        // (a codim edge against a non-codim vertex is already covered by the
        // ee/fv passes above, and a codim vertex against a face-adjacent edge
        // by fv). The passes are skipped entirely for surface-mesh pairs.
        const bool i_has_cv = bodies.body_num_codim_vertices(body_i) > 0;
        const bool j_has_cv = bodies.body_num_codim_vertices(body_j) > 0;

        if (i_has_cv && j_has_cv) {
            // (vᵢ, vⱼ): receives (this_vertex_id, other_vertex_id)
            std::vector<VertexVertexCandidate> vv;
            body_i_bvh.detect_vertex_vertex_candidates(
                body_j_bvh,
                [&](size_t vi, size_t vj) {
                    return bodies.is_codim_vertex(vi_start + vi)
                        && bodies.is_codim_vertex(vj_start + vj)
                        && can_vv(vi_start + vi, vj_start + vj);
                },
                vv);
            for (const VertexVertexCandidate& c : vv) {
                vv_candidates.emplace_back(
                    vi_start + c.vertex0_id, vj_start + c.vertex1_id);
            }
        }

        if (dim == 3) {
            if (bodies.body_num_codim_edges(body_i) > 0 && j_has_cv) {
                // (eᵢ, vⱼ): receives (this_edge_id, other_vertex_id)
                std::vector<EdgeVertexCandidate> ev;
                body_i_bvh.detect_edge_vertex_candidates(
                    body_j_bvh,
                    [&](size_t ei, size_t vj) {
                        return bodies.is_codim_edge(ei_start + ei)
                            && bodies.is_codim_vertex(vj_start + vj)
                            && can_ev(ei_start + ei, vj_start + vj);
                    },
                    ev);
                for (const EdgeVertexCandidate& c : ev) {
                    ev_candidates.emplace_back(
                        ei_start + c.edge_id, vj_start + c.vertex_id);
                }
            }

            if (i_has_cv && bodies.body_num_codim_edges(body_j) > 0) {
                // (eⱼ, vᵢ): receives (other_edge_id, this_vertex_id)
                std::vector<EdgeVertexCandidate> ev;
                body_i_bvh.detect_vertex_edge_candidates(
                    body_j_bvh,
                    [&](size_t ej, size_t vi) {
                        return bodies.is_codim_edge(ej_start + ej)
                            && bodies.is_codim_vertex(vi_start + vi)
                            && can_ev(ej_start + ej, vi_start + vi);
                    },
                    ev);
                for (const EdgeVertexCandidate& c : ev) {
                    ev_candidates.emplace_back(
                        ej_start + c.edge_id, vi_start + c.vertex_id);
                }
            }
        }
    }
}

// ============================================================================
// Narrow phase

namespace {
    /// @brief Nonlinear CCD between a rigidly-moving point and a static plane.
    ///
    /// Uses the conservative piecewise-linear CCD driver: within each
    /// linearized piece the position is affine in t, so the plane distance is
    /// affine and attains its minimum at an endpoint; the crossing time of the
    /// (deviation-inflated) minimum separation is computed exactly.
    bool rigid_point_plane_ccd(
        const Eigen::Hyperplane<double, 3>& plane,
        const RigidTrajectory& p,
        double& toi,
        const double min_distance,
        const double tmax,
        const double conservative_rescaling)
    {
        return NonlinearCCD::conservative_piecewise_linear_ccd(
            [&](const double t) {
                return std::abs(plane.signedDistance(p(t)));
            },
            [&](const double t0, const double t1) {
                return p.max_distance_from_linear(t0, t1);
            },
            [&](const double ti0, const double ti1, const double min_d,
                const bool no_zero_toi, double& piece_toi) {
                const double d0 = plane.signedDistance(p(ti0));
                const double d1 = plane.signedDistance(p(ti1));
                // Measure from the side the point starts on
                const double sign = d0 >= 0 ? 1.0 : -1.0;
                const double e0 = sign * d0;
                const double e1 = sign * d1;
                if (e1 >= min_d) {
                    return false; // Affine distance is minimized at an endpoint
                }
                piece_toi = e0 > e1
                    ? std::max((e0 - min_d) / (e0 - e1), 0.0)
                    : 0.0; // Let the driver subdivide degenerate cases
                return true;
            },
            toi, min_distance, tmax, conservative_rescaling);
    }
} // namespace

bool RigidCandidates::candidate_ccd(
    const size_t i,
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double min_distance,
    const double tmax,
    const NonlinearCCD& nonlinear_ccd,
    double& toi) const
{
    const auto trajectory = [&](const index_t vi) {
        const index_t b = bodies.vertex_to_body(vi);
        return RigidTrajectory(
            bodies.rest_positions().row(vi).transpose(), poses_t0[b],
            poses_t1[b]);
    };

    const Eigen::MatrixXi& E = bodies.edges();
    const Eigen::MatrixXi& F = bodies.faces();

    size_t j = i; // Match the ordering of Candidates::operator[]
    if (j < vv_candidates.size()) {
        const VertexVertexCandidate& c = vv_candidates[j];
        return nonlinear_ccd.point_point_ccd(
            trajectory(c.vertex0_id), trajectory(c.vertex1_id), toi,
            min_distance, tmax);
    }
    j -= vv_candidates.size();
    if (j < ev_candidates.size()) {
        const EdgeVertexCandidate& c = ev_candidates[j];
        return nonlinear_ccd.point_edge_ccd(
            trajectory(c.vertex_id), trajectory(E(c.edge_id, 0)),
            trajectory(E(c.edge_id, 1)), toi, min_distance, tmax);
    }
    j -= ev_candidates.size();
    if (j < ee_candidates.size()) {
        const EdgeEdgeCandidate& c = ee_candidates[j];
        return nonlinear_ccd.edge_edge_ccd(
            trajectory(E(c.edge0_id, 0)), trajectory(E(c.edge0_id, 1)),
            trajectory(E(c.edge1_id, 0)), trajectory(E(c.edge1_id, 1)), toi,
            min_distance, tmax);
    }
    j -= ee_candidates.size();
    if (j < fv_candidates.size()) {
        const FaceVertexCandidate& c = fv_candidates[j];
        return nonlinear_ccd.point_triangle_ccd(
            trajectory(c.vertex_id), trajectory(F(c.face_id, 0)),
            trajectory(F(c.face_id, 1)), trajectory(F(c.face_id, 2)), toi,
            min_distance, tmax);
    }
    j -= fv_candidates.size();
    assert(j < pv_candidates.size());
    const PlaneVertexCandidate& c = pv_candidates[j];
    return rigid_point_plane_ccd(
        c.plane, trajectory(c.vertex_id), toi, min_distance, tmax,
        nonlinear_ccd.conservative_rescaling);
}

bool RigidCandidates::is_step_collision_free(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double min_distance,
    const NonlinearCCD& nonlinear_ccd) const
{
    assert(poses_t0.size() == bodies.num_bodies());
    assert(poses_t1.size() == bodies.num_bodies());

    // Narrow phase
    for (size_t i = 0; i < size(); i++) {
        double toi;
        if (candidate_ccd(
                i, bodies, poses_t0, poses_t1, min_distance, /*tmax=*/1.0,
                nonlinear_ccd, toi)) {
            return false;
        }
    }

    return true;
}

double RigidCandidates::compute_collision_free_stepsize(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double min_distance,
    const NonlinearCCD& nonlinear_ccd) const
{
    assert(poses_t0.size() == bodies.num_bodies());
    assert(poses_t1.size() == bodies.num_bodies());

    if (empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    std::atomic<double> earliest_toi(1.0);

    tbb::parallel_for(size_t(0), size(), [&](size_t i) {
        const double tmax = earliest_toi.load(std::memory_order_relaxed);

        double toi = std::numeric_limits<double>::infinity(); // output
        const bool are_colliding = candidate_ccd(
            i, bodies, poses_t0, poses_t1, min_distance, tmax, nonlinear_ccd,
            toi);

        if (are_colliding) {
            // Update the earliest time of impact (TOI) atomically
            double prev = earliest_toi.load(std::memory_order_relaxed);
            while (toi < prev
                   && !earliest_toi.compare_exchange_weak(
                       prev, toi, std::memory_order_relaxed)) { }
        }
    });

    const double result = earliest_toi.load(std::memory_order_relaxed);
    assert(result >= 0 && result <= 1.0);
    return result;
}

double RigidCandidates::compute_noncandidate_conservative_stepsize(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double dhat) const
{
    assert(poses_t0.size() == bodies.num_bodies());
    assert(poses_t1.size() == bodies.num_bodies());

    // Bound each vertex's path length by ‖Δp‖ + ‖Δθ‖‖x̄‖ (see
    // rigid_chord_deviation for the rotational velocity bound).
    Eigen::MatrixXd displacement_bounds(bodies.num_vertices(), 1);
    for (index_t vi = 0; vi < bodies.num_vertices(); ++vi) {
        const index_t b = bodies.vertex_to_body(vi);
        displacement_bounds(vi, 0) =
            (poses_t1[b].position - poses_t0[b].position).norm()
            + (poses_t1[b].rotation - poses_t0[b].rotation).norm()
                * bodies.rest_positions().row(vi).norm();
    }

    return Candidates::compute_noncandidate_conservative_stepsize(
        bodies, displacement_bounds, dhat);
}

double RigidCandidates::compute_cfl_stepsize(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double dhat,
    const double min_distance,
    BroadPhase* broad_phase,
    const NonlinearCCD& nonlinear_ccd) const
{
    const double alpha_C = this->compute_collision_free_stepsize(
        bodies, poses_t0, poses_t1, min_distance, nonlinear_ccd);

    const double alpha_F = this->compute_noncandidate_conservative_stepsize(
        bodies, poses_t0, poses_t1, dhat);

    // If alpha_F < 0.5 * alpha_C, then we should do full CCD.
    if (alpha_F < 0.5 * alpha_C) {
        RigidCandidates full_candidates;
        full_candidates.build(
            bodies, poses_t0, poses_t1, /*inflation_radius=*/0, broad_phase);
        return full_candidates.compute_collision_free_stepsize(
            bodies, poses_t0, poses_t1, min_distance, nonlinear_ccd);
    }
    return std::min(alpha_C, alpha_F);
}

} // namespace ipc::rigid
