#include "candidates.hpp"

#include <ipc/config.hpp>
#include <ipc/ipc.hpp>
#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/io/write_candidates_obj.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <igl/remove_unreferenced.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <fstream>
#include <shared_mutex>

namespace ipc {

void RigidCandidates::build(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double inflation_radius,
    const std::shared_ptr<BroadPhase> broad_phase)
{
    assert(poses_t0.size() == bodies.num_bodies());
    assert(poses_t1.size() == bodies.num_bodies());
    assert(broad_phase != nullptr);

    clear();

    const int dim = poses_t0[0].position.size();

    // 1. Broad phase between bodies
    std::vector<std::pair<index_t, index_t>> body_pairs;
    {
        Eigen::MatrixXd V(2 * bodies.num_bodies(), dim);
        Eigen::MatrixXi E(bodies.num_bodies(), 2);
        // Eigen::VectorXi group_ids(2 * num_bodies());
        double max_radius = 0;
        for (int i = 0; i < num_bodies(); i++) {
            V.row(2 * i + 0) = poses_t0[i].position;
            V.row(2 * i + 1) = poses_t1[i].position;
            E(i, 0) = 2 * i + 0;
            E(i, 1) = 2 * i + 1;
            max_radius = std::max(m_rbs[i].r_max, max_radius);
            group_ids[2 * i + 1] = group_ids[2 * i] = m_rbs[i].group_id;
        }
    }
}

bool RigidCandidates::is_step_collision_free(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    const double min_distance,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    // Narrow phase
    for (size_t i = 0; i < size(); i++) {
        const CollisionStencil& candidate = (*this)[i];

        double toi;
        bool is_collision = candidate.ccd(
            candidate.dof(vertices_t0, mesh.edges(), mesh.faces()),
            candidate.dof(vertices_t1, mesh.edges(), mesh.faces()), //
            toi, min_distance, /*tmax=*/1.0, narrow_phase_ccd);

        if (is_collision) {
            return false;
        }
    }

    return true;
}

double RigidCandidates::compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    const double min_distance,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    if (empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    double earliest_toi = 1;
    std::shared_mutex earliest_toi_mutex;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, size()),
        [&](tbb::blocked_range<size_t> r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                // Use the mutex to read as well in case writing double takes
                // more than one clock cycle.
                double tmax;
                {
                    std::shared_lock lock(earliest_toi_mutex);
                    tmax = earliest_toi;
                }

                const CollisionStencil& candidate = (*this)[i];

                double toi = std::numeric_limits<double>::infinity(); // output
                const bool are_colliding = candidate.ccd(
                    candidate.dof(vertices_t0, mesh.edges(), mesh.faces()),
                    candidate.dof(vertices_t1, mesh.edges(), mesh.faces()), //
                    toi, min_distance, tmax, narrow_phase_ccd);

                if (are_colliding) {
                    std::unique_lock lock(earliest_toi_mutex);
                    if (toi < earliest_toi) {
                        earliest_toi = toi;
                    }
                }
            }
        });

    assert(earliest_toi >= 0 && earliest_toi <= 1.0);
    return earliest_toi;
}

double RigidCandidates::compute_noncandidate_conservative_stepsize(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> displacements,
    const double dhat) const
{
    assert(displacements.rows() == mesh.num_vertices());

    if (empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    std::vector<bool> is_vertex_a_candidates(mesh.num_vertices(), false);
    for (size_t i = 0; i < size(); i++) {
        for (const index_t vid : (*this)[i].vertex_ids(E, F)) {
            if (vid < 0) {
                break;
            }
            is_vertex_a_candidates[vid] = true;
        }
    }

    double max_displacement = 0;
    for (size_t i = 0; i < displacements.rows(); i++) {
        if (!is_vertex_a_candidates[i]) {
            continue;
        }
        max_displacement =
            std::max(max_displacement, displacements.row(i).norm());
    }

    return 0.5 * dhat / max_displacement;
}

double RigidCandidates::compute_cfl_stepsize(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    const double dhat,
    const double min_distance,
    const std::shared_ptr<BroadPhase> broad_phase,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    const double alpha_C = this->compute_collision_free_stepsize(
        mesh, vertices_t0, vertices_t1, min_distance, narrow_phase_ccd);

    const double alpha_F = this->compute_noncandidate_conservative_stepsize(
        mesh, vertices_t1 - vertices_t0, dhat);

    // If alpha_F < 0.5 * alpha_C, then we should do full CCD.
    if (alpha_F < 0.5 * alpha_C) {
        return ipc::compute_collision_free_stepsize(
            mesh, vertices_t0, vertices_t1, min_distance, broad_phase,
            narrow_phase_ccd);
    }
    return std::min(alpha_C, alpha_F);
}

size_t RigidCandidates::size() const
{
    return vv_candidates.size() + ev_candidates.size() + ee_candidates.size()
        + fv_candidates.size();
}

bool RigidCandidates::empty() const
{
    return vv_candidates.empty() && ev_candidates.empty()
        && ee_candidates.empty() && fv_candidates.empty();
}

void RigidCandidates::clear()
{
    vv_candidates.clear();
    ev_candidates.clear();
    ee_candidates.clear();
    fv_candidates.clear();
}

CollisionStencil& RigidCandidates::operator[](size_t i)
{
    if (i < vv_candidates.size()) {
        return vv_candidates[i];
    }
    i -= vv_candidates.size();
    if (i < ev_candidates.size()) {
        return ev_candidates[i];
    }
    i -= ev_candidates.size();
    if (i < ee_candidates.size()) {
        return ee_candidates[i];
    }
    i -= ee_candidates.size();
    if (i < fv_candidates.size()) {
        return fv_candidates[i];
    }
    throw std::out_of_range("Candidate index is out of range!");
}

const CollisionStencil& RigidCandidates::operator[](size_t i) const
{
    if (i < vv_candidates.size()) {
        return vv_candidates[i];
    }
    i -= vv_candidates.size();
    if (i < ev_candidates.size()) {
        return ev_candidates[i];
    }
    i -= ev_candidates.size();
    if (i < ee_candidates.size()) {
        return ee_candidates[i];
    }
    i -= ee_candidates.size();
    if (i < fv_candidates.size()) {
        return fv_candidates[i];
    }
    throw std::out_of_range("Candidate index is out of range!");
}

bool RigidCandidates::write_obj(
    const std::string& filename,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces) const
{
    std::ofstream obj(filename, std::ios::out);
    if (!obj.is_open()) {
        return false;
    }
    int v_offset = 0;
    ipc::write_candidates_obj(
        obj, vertices, edges, faces, vv_candidates, v_offset);
    v_offset += vv_candidates.size() * 2;
    ipc::write_candidates_obj(
        obj, vertices, edges, faces, ev_candidates, v_offset);
    v_offset += ev_candidates.size() * 3;
    ipc::write_candidates_obj(
        obj, vertices, edges, faces, ee_candidates, v_offset);
    v_offset += ee_candidates.size() * 4;
    ipc::write_candidates_obj(
        obj, vertices, faces, faces, fv_candidates, v_offset);
    return true;
}

} // namespace ipc
