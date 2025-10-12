#include "smooth_collisions.hpp"

#include "smooth_collisions_builder.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/maybe_parallel_for.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

void SmoothCollisions::compute_adaptive_dhat(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices, // set to zero for rest pose
    const ParameterType param,
    const std::shared_ptr<BroadPhase> broad_phase)
{
    assert(vertices.rows() == mesh.num_vertices());

    const double dhat = param.dhat;
    double inflation_radius = dhat / 2;

    // Candidates candidates;
    candidates.build(mesh, vertices, inflation_radius, broad_phase);
    this->build(
        candidates, mesh, vertices, param,
        false /*disable adaptive dhat to compute true pairs*/);

    vert_adaptive_dhat.setConstant(mesh.num_vertices(), dhat);
    edge_adaptive_dhat.setConstant(mesh.num_edges(), dhat);
    if (mesh.dim() == 3)
        face_adaptive_dhat.setConstant(mesh.num_faces(), dhat);
    else
        face_adaptive_dhat.resize(0);

    auto assign_min = [](double& a, const double& b) -> void {
        a = std::min(a, b);
    };

    for (const auto& cc : collisions) {
        const double dist =
            param.adaptive_dhat_ratio() * sqrt(cc->compute_distance(vertices));
        switch (cc->type()) {
        case CollisionType::EDGE_EDGE:
            assign_min(edge_adaptive_dhat((*cc)[0]), dist);
            assign_min(edge_adaptive_dhat((*cc)[1]), dist);
            break;
        case CollisionType::EDGE_VERTEX:
            assign_min(edge_adaptive_dhat((*cc)[0]), dist);
            assign_min(vert_adaptive_dhat((*cc)[1]), dist);
            break;
        case CollisionType::FACE_VERTEX:
            assign_min(face_adaptive_dhat((*cc)[0]), dist);
            assign_min(vert_adaptive_dhat((*cc)[1]), dist);
            break;
        case CollisionType::VERTEX_VERTEX:
            assign_min(vert_adaptive_dhat((*cc)[0]), dist);
            assign_min(vert_adaptive_dhat((*cc)[1]), dist);
            break;
        default:
            throw std::runtime_error("Invalid collision type!");
        }
    }

    // face adaptive dhat should be minimum of all its adjacent vertices and
    // edges
    if (mesh.dim() == 3)
        for (int f = 0; f < mesh.num_faces(); f++) {
            for (int lv = 0; lv < 3; lv++) {
                face_adaptive_dhat(f) = std::min(
                    face_adaptive_dhat(f),
                    vert_adaptive_dhat(mesh.faces()(f, lv)));
                face_adaptive_dhat(f) = std::min(
                    face_adaptive_dhat(f),
                    edge_adaptive_dhat(mesh.faces_to_edges()(f, lv)));
            }
        }

    // edge adaptive dhat should be minimum of all its adjacent vertices
    for (int e = 0; e < mesh.num_edges(); e++) {
        for (int lv = 0; lv < 2; lv++) {
            edge_adaptive_dhat(e) = std::min(
                edge_adaptive_dhat(e), vert_adaptive_dhat(mesh.edges()(e, lv)));
        }
    }

    logger().debug(
        "Adaptive dhat: vert dhat min {:.2e}, max {:.2e}",
        vert_adaptive_dhat.minCoeff(), vert_adaptive_dhat.maxCoeff());
    logger().debug(
        "Adaptive dhat: edge dhat min {:.2e}, max {:.2e}",
        edge_adaptive_dhat.minCoeff(), edge_adaptive_dhat.maxCoeff());
    if (mesh.dim() == 3)
        logger().debug(
            "Adaptive dhat: face dhat min {:.2e}, max {:.2e}",
            face_adaptive_dhat.minCoeff(), face_adaptive_dhat.maxCoeff());
}

void SmoothCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ParameterType param,
    const bool use_adaptive_dhat,
    const std::shared_ptr<BroadPhase> broad_phase)
{
    assert(vertices.rows() == mesh.num_vertices());

    double inflation_radius = param.dhat / 2;

    // Candidates candidates;
    candidates.build(mesh, vertices, inflation_radius, broad_phase);
    // std::cout << "Candidate Memory " << getCurrentRSS() / (1024.*1024) <<
    // "MB\n";
    this->build(candidates, mesh, vertices, param, use_adaptive_dhat);
}

void SmoothCollisions::build(
    const Candidates& candidates_,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ParameterType param,
    const bool use_adaptive_dhat)
{
    assert(vertices.rows() == mesh.num_vertices());

    clear();

    const double dhat = param.dhat;
    if (!use_adaptive_dhat) {
        vert_adaptive_dhat.resize(1);
        vert_adaptive_dhat(0) = dhat;
        edge_adaptive_dhat.resize(1);
        edge_adaptive_dhat(0) = dhat;
        if (mesh.dim() == 3) {
            face_adaptive_dhat.resize(1);
            face_adaptive_dhat(0) = dhat;
        } else
            face_adaptive_dhat.resize(0);
    }

    auto vert_dhat = [&](const long& v_id) {
        return this->get_vert_dhat(v_id);
    };
    auto edge_dhat = [&](const long& e_id) {
        return this->get_edge_dhat(e_id);
    };
    auto face_dhat = [&](const long& f_id) {
        return this->get_face_dhat(f_id);
    };

    if (mesh.dim() == 2) {
        auto storage = create_thread_storage<SmoothCollisionsBuilder<2>>(
            SmoothCollisionsBuilder<2>());
        maybe_parallel_for(
            candidates_.ev_candidates.size(),
            [&](int start, int end, int thread_id) {
                SmoothCollisionsBuilder<2>& local_storage =
                    get_local_thread_storage(storage, thread_id);
                local_storage.add_edge_vertex_collisions(
                    mesh, vertices, candidates_.ev_candidates, param, vert_dhat,
                    edge_dhat, start, end);
            });
        SmoothCollisionsBuilder<2>::merge(storage, *this);
    } else {
        auto storage = create_thread_storage<SmoothCollisionsBuilder<3>>(
            SmoothCollisionsBuilder<3>());
        maybe_parallel_for(
            candidates_.ee_candidates.size(),
            [&](int start, int end, int thread_id) {
                SmoothCollisionsBuilder<3>& local_storage =
                    get_local_thread_storage(storage, thread_id);
                local_storage.add_edge_edge_collisions(
                    mesh, vertices, candidates_.ee_candidates, param, vert_dhat,
                    edge_dhat, start, end);
            });

        maybe_parallel_for(
            candidates_.fv_candidates.size(),
            [&](int start, int end, int thread_id) {
                SmoothCollisionsBuilder<3>& local_storage =
                    get_local_thread_storage(storage, thread_id);
                local_storage.add_face_vertex_collisions(
                    mesh, vertices, candidates_.fv_candidates, param, vert_dhat,
                    edge_dhat, face_dhat, start, end);
            });
        SmoothCollisionsBuilder<3>::merge(storage, *this);
    }
    candidates = candidates_;
}

// ============================================================================
size_t SmoothCollisions::size() const { return collisions.size(); }
bool SmoothCollisions::empty() const { return collisions.empty(); }
void SmoothCollisions::clear() { collisions.clear(); }

SmoothCollision& SmoothCollisions::operator[](size_t i)
{
    if (i < collisions.size()) {
        return *collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

const SmoothCollision& SmoothCollisions::operator[](size_t i) const
{
    if (i < collisions.size()) {
        return *collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

std::string SmoothCollisions::to_string(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const ParameterType& params) const
{
    std::stringstream ss;
    for (const auto& cc : collisions) {
        ss << "\n";
        {
            ss << fmt::format(
                "[{}]: ({} {}) dist {} potential {} grad {}", cc->name(),
                (*cc)[0], (*cc)[1], cc->compute_distance(vertices),
                (*cc)(cc->dof(vertices), params),
                (*cc).gradient(cc->dof(vertices), params).norm());
        }
    }
    return ss.str();
}

// NOTE: Actually distance squared
double SmoothCollisions::compute_minimum_distance(
    const CollisionMesh& mesh, Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (candidates.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, candidates.size()),
        [&](tbb::blocked_range<size_t> r) {
            double& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist = candidates[i].compute_distance(
                    candidates[i].dof(vertices, edges, faces));

                if (dist < local_min_dist) {
                    local_min_dist = dist;
                }
            }
        });

    return storage.combine([](double a, double b) { return std::min(a, b); });
}

double SmoothCollisions::compute_active_minimum_distance(
    const CollisionMesh& mesh, Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (collisions.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, collisions.size()),
        [&](tbb::blocked_range<size_t> r) {
            double& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist = collisions[i]->compute_distance(vertices);

                if (collisions[i]->is_active() && dist < local_min_dist) {
                    local_min_dist = dist;
                }
            }
        });

    return storage.combine([](double a, double b) { return std::min(a, b); });
}

} // namespace ipc
