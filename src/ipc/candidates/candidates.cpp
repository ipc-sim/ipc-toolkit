#include "candidates.hpp"

#include <ipc/ipc.hpp>
#include <ipc/utils/save_obj.hpp>

#include <ipc/config.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <shared_mutex>

#include <fstream>

namespace ipc {

void Candidates::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double inflation_radius,
    const BroadPhaseMethod broad_phase_method)
{
    const int dim = vertices.cols();

    clear();

    std::unique_ptr<BroadPhase> broad_phase =
        BroadPhase::make_broad_phase(broad_phase_method);
    broad_phase->can_vertices_collide = mesh.can_collide;
    broad_phase->build(vertices, mesh.edges(), mesh.faces(), inflation_radius);
    broad_phase->detect_collision_candidates(dim, *this);
    broad_phase->clear();
}

void Candidates::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double inflation_radius,
    const BroadPhaseMethod broad_phase_method)
{
    const int dim = vertices_t0.cols();

    clear();

    std::unique_ptr<BroadPhase> broad_phase =
        BroadPhase::make_broad_phase(broad_phase_method);
    broad_phase->can_vertices_collide = mesh.can_collide;
    broad_phase->build(
        vertices_t0, vertices_t1, mesh.edges(), mesh.faces(), inflation_radius);
    broad_phase->detect_collision_candidates(dim, *this);
    broad_phase->clear();
}

bool Candidates::is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double min_distance,
    const double tolerance,
    const long max_iterations) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    // Narrow phase
    for (size_t i = 0; i < size(); i++) {
        double toi;
        bool is_collision = (*this)[i].ccd(
            vertices_t0, vertices_t1, mesh.edges(), mesh.faces(), toi,
            min_distance,
            /*tmax=*/1.0, tolerance, max_iterations);

        if (is_collision) {
            return false;
        }
    }

    return true;
}

double Candidates::compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double min_distance,
    const double tolerance,
    const long max_iterations) const
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

                double toi = std::numeric_limits<double>::infinity(); // output
                const bool are_colliding = (*this)[i].ccd(
                    vertices_t0, vertices_t1, mesh.edges(), mesh.faces(), toi,
                    min_distance, tmax, tolerance, max_iterations);

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

double Candidates::compute_noncandidate_conservative_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& displacements,
    const double dhat) const
{
    assert(displacements.rows() == mesh.num_vertices());

    if (empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    const auto& E = mesh.edges();
    const auto& F = mesh.faces();

    std::vector<bool> is_vertex_a_candidates(mesh.num_vertices(), false);
    for (size_t i = 0; i < size(); i++) {
        for (const long vid : (*this)[i].vertex_ids(E, F)) {
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

double Candidates::compute_cfl_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double dhat,
    const BroadPhaseMethod broad_phase_method,
    const double min_distance,
    const double tolerance,
    const long max_iterations) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    const double alpha_C = this->compute_collision_free_stepsize(
        mesh, vertices_t0, vertices_t1, min_distance, tolerance,
        max_iterations);

    const double alpha_F = this->compute_noncandidate_conservative_stepsize(
        mesh, vertices_t1 - vertices_t0, dhat);

    // If alpha_F < 0.5 * alpha_C, then we should do full CCD.
    if (alpha_F < 0.5 * alpha_C) {
        return ipc::compute_collision_free_stepsize(
            mesh, vertices_t0, vertices_t1, broad_phase_method, min_distance,
            tolerance, max_iterations);
    }
    return std::min(alpha_C, alpha_F);
}

size_t Candidates::size() const
{
    return ev_candidates.size() + ee_candidates.size() + fv_candidates.size();
}

bool Candidates::empty() const
{
    return ev_candidates.empty() && ee_candidates.empty()
        && fv_candidates.empty();
}

void Candidates::clear()
{
    ev_candidates.clear();
    ee_candidates.clear();
    fv_candidates.clear();
}

ContinuousCollisionCandidate& Candidates::operator[](size_t idx)
{
    if (idx < ev_candidates.size()) {
        return ev_candidates[idx];
    }
    idx -= ev_candidates.size();
    if (idx < ee_candidates.size()) {
        return ee_candidates[idx];
    }
    idx -= ee_candidates.size();
    if (idx < fv_candidates.size()) {
        return fv_candidates[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

const ContinuousCollisionCandidate& Candidates::operator[](size_t idx) const
{
    if (idx < ev_candidates.size()) {
        return ev_candidates[idx];
    }
    idx -= ev_candidates.size();
    if (idx < ee_candidates.size()) {
        return ee_candidates[idx];
    }
    idx -= ee_candidates.size();
    if (idx < fv_candidates.size()) {
        return fv_candidates[idx];
    }
    throw std::out_of_range("Constraint index is out of range!");
}

bool Candidates::save_obj(
    const std::string& filename,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    std::ofstream obj(filename, std::ios::out);
    if (!obj.is_open()) {
        return false;
    }
    int v_offset = 0;
    ipc::save_obj(obj, vertices, edges, faces, ev_candidates, v_offset);
    v_offset += ev_candidates.size() * 3;
    ipc::save_obj(obj, vertices, edges, faces, ee_candidates, v_offset);
    v_offset += ee_candidates.size() * 4;
    ipc::save_obj(obj, vertices, faces, faces, fv_candidates, v_offset);
    return true;
}

} // namespace ipc
