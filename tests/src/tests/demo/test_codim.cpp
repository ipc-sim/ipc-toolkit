// Codimensional (point-cloud / segment-net) contact tests for the rigid
// body-pair broad phase.

#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/ccd/nonlinear_ccd.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/rigid_candidates.hpp>
#include <simulator.hpp>

#include <set>

using namespace ipc;
using namespace ipc::rigid;

namespace {

/// A 3×3 grid of isolated points in the z = 0 plane (3D) with a small tetra
/// offset so the inertia is nonsingular.
void point_grid(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F)
{
    V.resize(10, 3);
    int r = 0;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            V.row(r++) = Eigen::RowVector3d(0.5 * i, 0, 0.5 * j);
        }
    }
    V.row(r++) = Eigen::RowVector3d(0, 0.3, 0); // break planarity
    E.resize(0, 2);
    F.resize(0, 3);
}

/// A zig-zag segment net (edges only, no faces).
void segment_net(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F)
{
    V.resize(5, 3);
    V << -1, 0, 0, -0.5, 0.2, 0, 0, 0, 0, 0.5, 0.2, 0, 1, 0, 0;
    E.resize(4, 2);
    E << 0, 1, 1, 2, 2, 3, 3, 4;
    F.resize(0, 3);
}

} // namespace

TEST_CASE(
    "Codim broad phase completeness", "[codim][rigid][candidates]")
{
    Eigen::MatrixXd V_cloud;
    Eigen::MatrixXi E_cloud, F_cloud;
    point_grid(V_cloud, E_cloud, F_cloud);

    Eigen::MatrixXd V_net;
    Eigen::MatrixXi E_net, F_net;
    segment_net(V_net, E_net, F_net);

    // Cloud falls straight down through the net.
    std::vector<Pose> poses_t0(2, Pose::Identity(3));
    poses_t0[0].position = Eigen::Vector3d(0, 1.0, 0); // cloud
    poses_t0[1].position = Eigen::Vector3d::Zero();    // net

    auto bodies = RigidBodies::build_from_meshes(
        { V_cloud, V_net }, { E_cloud, E_net }, { F_cloud, F_net },
        { 1000.0, 1000.0 }, poses_t0);

    REQUIRE(bodies->body_num_codim_vertices(0) == 10);
    REQUIRE(bodies->body_num_codim_edges(1) == 4);

    std::vector<Pose> poses_t1 = poses_t0;
    poses_t1[0].position.y() -= 2.0; // sweep through the net

    const double inflation = 1e-3;
    RigidCandidates candidates;
    candidates.build(*bodies, poses_t0, poses_t1, inflation);

    // Sampled ground truth: every (codim edge, cloud vertex) pair that comes
    // within the inflation radius at some sampled time must be a candidate.
    std::set<std::pair<index_t, index_t>> close_ev;
    const int n_samples = 200;
    for (int s = 0; s <= n_samples; ++s) {
        const double t = s / double(n_samples);
        std::vector<Pose> poses(2);
        for (int b = 0; b < 2; ++b) {
            poses[b] = Pose(
                (1 - t) * poses_t0[b].position + t * poses_t1[b].position,
                poses_t0[b].rotation);
        }
        const Eigen::MatrixXd V = bodies->vertices(poses);
        for (index_t e = 0; e < bodies->edges().rows(); ++e) {
            for (index_t v = 0; v < 10; ++v) { // cloud vertices are 0..9
                const double d_sq = point_edge_distance(
                    V.row(v), V.row(bodies->edges()(e, 0)),
                    V.row(bodies->edges()(e, 1)));
                if (d_sq <= inflation * inflation) {
                    close_ev.emplace(e, v);
                }
            }
        }
    }
    REQUIRE(!close_ev.empty());

    std::set<std::pair<index_t, index_t>> found_ev;
    for (const auto& c : candidates.ev_candidates) {
        found_ev.emplace(c.edge_id, c.vertex_id);
    }
    for (const auto& pair : close_ev) {
        CAPTURE(pair.first, pair.second);
        CHECK(found_ev.count(pair) == 1);
    }
}

TEST_CASE(
    "Codim vertex-vertex candidates and CCD", "[codim][rigid][candidates]")
{
    Eigen::MatrixXd V_cloud;
    Eigen::MatrixXi E_cloud, F_cloud;
    point_grid(V_cloud, E_cloud, F_cloud);

    // Two identical clouds, one falling straight onto the other: aligned
    // points must generate vertex-vertex candidates and point-point CCD must
    // stop the sweep.
    std::vector<Pose> poses_t0(2, Pose::Identity(3));
    poses_t0[0].position = Eigen::Vector3d(0, 1.0, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V_cloud, V_cloud }, { E_cloud, E_cloud }, { F_cloud, F_cloud },
        { 1000.0, 1000.0 }, poses_t0);

    std::vector<Pose> poses_t1 = poses_t0;
    poses_t1[0].position.y() -= 2.0; // full pass-through if unimpeded

    RigidCandidates candidates;
    candidates.build(*bodies, poses_t0, poses_t1, /*inflation_radius=*/1e-3);
    CHECK(!candidates.vv_candidates.empty());

    NonlinearCCD ccd;
    const double alpha = candidates.compute_collision_free_stepsize(
        *bodies, poses_t0, poses_t1, /*min_distance=*/0, ccd);
    CHECK(alpha < 1.0); // the clouds may not pass through each other

    // The clamped step keeps all point pairs separated.
    std::vector<Pose> poses_alpha(2);
    for (int b = 0; b < 2; ++b) {
        poses_alpha[b] = Pose(
            (1 - alpha) * poses_t0[b].position
                + alpha * poses_t1[b].position,
            poses_t0[b].rotation);
    }
    const Eigen::MatrixXd V = bodies->vertices(poses_alpha);
    for (int v = 0; v < 10; ++v) {
        // Strictly separated (squared distance > 0); the conservative CCD
        // rescaling leaves only a small gap.
        CHECK(point_point_distance(V.row(v), V.row(10 + v)) > 0);
    }

    // The discrete (distance) overload also produces the vv candidates.
    RigidCandidates discrete;
    discrete.build(*bodies, poses_alpha, /*inflation_radius=*/1e-2);
    CHECK(!discrete.vv_candidates.empty());
}

TEST_CASE("Codim point cloud drop", "[codim][simulator]")
{
    using ipc::demo::Simulator;

    Eigen::MatrixXd V_cloud;
    Eigen::MatrixXi E_cloud, F_cloud;
    point_grid(V_cloud, E_cloud, F_cloud);

    // A point cloud dropped onto an identical static cloud: point-point
    // contacts must stop it.
    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[0].position = Eigen::Vector3d(0, 0.2, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V_cloud, V_cloud }, { E_cloud, E_cloud }, { F_cloud, F_cloud },
        { 1000.0, 1000.0 }, initial_poses);
    (*bodies)[1].set_type(RigidBody::Type::STATIC);

    Simulator sim(bodies, initial_poses, /*dt=*/0.01);
    REQUIRE(sim.run(/*t_end=*/0.5));

    // The falling cloud rests on (strictly above) the static one.
    const Eigen::MatrixXd V = bodies->vertices(sim.rigid_pose_history().back());
    double min_d_sq = std::numeric_limits<double>::infinity();
    for (int v = 0; v < 10; ++v) {
        min_d_sq = std::min(
            min_d_sq, point_point_distance(V.row(v), V.row(10 + v)));
    }
    CHECK(min_d_sq > 0);
    CHECK(sim.poses()[0].position(1) > 0); // did not fall through
}

TEST_CASE("2D codim point drop", "[2d][codim][simulator]")
{
    using ipc::demo::Simulator;

    // Two isolated 2D points (one falling, one static below).
    Eigen::MatrixXd V(3, 2);
    V << -0.5, 0, 0.5, 0, 0, 0.4; // three points (nonsingular inertia)
    const Eigen::MatrixXi E(0, 2), F(0, 3);

    std::vector<Pose> initial_poses(2, Pose::Identity(2));
    initial_poses[0].position = Eigen::Vector2d(0, 0.5);

    auto bodies = RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);
    (*bodies)[1].set_type(RigidBody::Type::STATIC);

    REQUIRE(bodies->body_num_codim_vertices(0) == 3);

    Simulator sim(bodies, initial_poses, /*dt=*/0.01);
    REQUIRE(sim.run(/*t_end=*/0.4));

    // Aligned point pairs collide (vertex-vertex): the falling triple rests
    // above the static one.
    CHECK(sim.poses()[0].position(1) > 0);
}
