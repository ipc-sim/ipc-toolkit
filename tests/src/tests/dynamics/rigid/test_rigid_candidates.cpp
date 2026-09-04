#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/collision_filter.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/rigid_candidates.hpp>
#include <ipc/dynamics/rigid/rigid_trajectory.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <igl/PI.h>
#include <igl/edges.h>
#include <igl/upsample.h>

#include <set>

using namespace ipc;
using namespace ipc::rigid;

namespace {

/// A long thin rod along x (half-length 2) with distinct y/z extents so the
/// principal axes are uniquely determined.
void make_rod(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F)
{
    tests::load_mesh("cube.ply", V, E, F);
    V.col(0) *= 4.0;
    V.col(1) *= 0.2;
    V.col(2) *= 0.1;
}

std::vector<Pose> lerp_poses(
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double t)
{
    std::vector<Pose> poses(poses_t0.size());
    for (size_t i = 0; i < poses.size(); i++) {
        poses[i] = Pose(
            (1 - t) * poses_t0[i].position + t * poses_t1[i].position,
            (1 - t) * poses_t0[i].rotation + t * poses_t1[i].rotation);
    }
    return poses;
}

/// Compose a world-frame rotation on top of a pose's rotation.
/// NOTE: build_from_meshes folds each body's principal-axes rotation R₀ into
/// the initial poses, so world-frame rotations must be composed with — not
/// assigned over — the post-build rotation.
VectorMax3d
compose_rotation(const Eigen::Matrix3d& world_rotation, const Pose& pose)
{
    return rotation_matrix_to_vector(world_rotation * pose.rotation_matrix());
}

/// Collect the inter-body element pairs that come within `radius` at any
/// sampled time. These are pairs the broad phase must report.
void sampled_close_pairs(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double radius,
    std::set<std::pair<index_t, index_t>>& ee_pairs,
    std::set<std::pair<index_t, index_t>>& fv_pairs,
    const int n_samples = 21)
{
    const Eigen::MatrixXi& E = bodies.edges();
    const Eigen::MatrixXi& F = bodies.faces();

    const auto edge_body = [&](index_t e) {
        return bodies.vertex_to_body(E(e, 0));
    };
    const auto face_body = [&](index_t f) {
        return bodies.vertex_to_body(F(f, 0));
    };

    for (int s = 0; s <= n_samples; s++) {
        const double t = s / double(n_samples);
        const Eigen::MatrixXd V =
            bodies.vertices(lerp_poses(poses_t0, poses_t1, t));

        for (index_t ea = 0; ea < E.rows(); ea++) {
            for (index_t eb = ea + 1; eb < E.rows(); eb++) {
                if (edge_body(ea) == edge_body(eb)) {
                    continue;
                }
                const double d = edge_edge_distance(
                    V.row(E(ea, 0)), V.row(E(ea, 1)), V.row(E(eb, 0)),
                    V.row(E(eb, 1)));
                if (d <= radius * radius) {
                    ee_pairs.emplace(std::min(ea, eb), std::max(ea, eb));
                }
            }
        }

        for (index_t f = 0; f < F.rows(); f++) {
            for (index_t v = 0; v < V.rows(); v++) {
                if (face_body(f) == bodies.vertex_to_body(v)) {
                    continue;
                }
                const double d = point_triangle_distance(
                    V.row(v), V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)));
                if (d <= radius * radius) {
                    fv_pairs.emplace(f, v);
                }
            }
        }
    }
}

/// Sampled first time at which any inter-body element pair comes within
/// `threshold`. A small positive threshold with fine sampling is required
/// because exactly-touching instants fall between samples, and once a vertex
/// penetrates the other body's interior all surface distances are positive
/// again.
double sampled_first_contact(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses_t0,
    const std::vector<Pose>& poses_t1,
    const double threshold = 1e-3,
    const int n_samples = 20'000)
{
    const Eigen::MatrixXi& E = bodies.edges();
    const Eigen::MatrixXi& F = bodies.faces();
    for (int s = 0; s <= n_samples; s++) {
        const double t = s / double(n_samples);
        const Eigen::MatrixXd V =
            bodies.vertices(lerp_poses(poses_t0, poses_t1, t));
        for (index_t f = 0; f < F.rows(); f++) {
            for (index_t v = 0; v < V.rows(); v++) {
                if (bodies.vertex_to_body(F(f, 0))
                    == bodies.vertex_to_body(v)) {
                    continue;
                }
                const double d = point_triangle_distance(
                    V.row(v), V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)));
                if (d <= threshold * threshold) {
                    return t;
                }
            }
        }
        for (index_t ea = 0; ea < E.rows(); ea++) {
            for (index_t eb = ea + 1; eb < E.rows(); eb++) {
                if (bodies.vertex_to_body(E(ea, 0))
                    == bodies.vertex_to_body(E(eb, 0))) {
                    continue;
                }
                const double d = edge_edge_distance(
                    V.row(E(ea, 0)), V.row(E(ea, 1)), V.row(E(eb, 0)),
                    V.row(E(eb, 1)));
                if (d <= threshold * threshold) {
                    return t;
                }
            }
        }
    }
    return 1.0;
}

/// Build the spinning-rod-and-cube scene: a rod spinning half a revolution in
/// place with a small static cube placed on the rod tip's mid-sweep position.
/// The endpoint positions of the rod hug its initial axis, so linearized
/// vertex trajectories never approach the cube — only the mid-step sweep does.
std::shared_ptr<RigidBodies> make_spinning_rod_scene(
    std::vector<Pose>& poses_t0, std::vector<Pose>& poses_t1)
{
    Eigen::MatrixXd V_rod, V_cube;
    Eigen::MatrixXi E_rod, F_rod, E_cube, F_cube;
    make_rod(V_rod, E_rod, F_rod);
    tests::load_mesh("cube.ply", V_cube, E_cube, F_cube);

    // Stage 1: build the rod alone to learn its post-build pose (R₀ folded
    // in) and locate the tip's mid-sweep position.
    std::vector<Pose> rod_poses(1, Pose::Identity(3));
    auto rod_only = RigidBodies::build_from_meshes(
        { V_rod }, { E_rod }, { F_rod }, { 1000.0 }, rod_poses);

    Pose rod_t1 = rod_poses[0];
    rod_t1.rotation = compose_rotation(
        Eigen::AngleAxisd(igl::PI, Eigen::Vector3d::UnitZ()).toRotationMatrix(),
        rod_poses[0]);

    // Tip = rest vertex with the largest radius
    index_t tip = 0;
    rod_only->body_rest_positions(0).rowwise().norm().maxCoeff(&tip);
    const VectorMax3d tip_mid = RigidTrajectory(
        rod_only->body_rest_positions(0).row(tip).transpose(), rod_poses[0],
        rod_t1)(0.5);

    // Stage 2: build the combined scene with the cube at the mid-sweep point,
    // offset slightly outward so the endpoints are collision free.
    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[1].position = 1.2 * tip_mid;

    auto bodies = RigidBodies::build_from_meshes(
        { V_rod, V_cube }, { E_rod, E_cube }, { F_rod, F_cube },
        { 1000.0, 1000.0 }, initial_poses);

    poses_t0 = initial_poses;
    poses_t1 = initial_poses;
    poses_t1[0] = rod_t1; // The rod's build is identical in both stages

    return bodies;
}

} // namespace

TEST_CASE("Rigid broad phase completeness", "[rigid][candidates]")
{
    // Two bodies with different vertex counts (exercises the body swap logic):
    // a cube and a once-upsampled cube.
    Eigen::MatrixXd V_cube;
    Eigen::MatrixXi E_cube, F_cube;
    tests::load_mesh("cube.ply", V_cube, E_cube, F_cube);

    Eigen::MatrixXd V_fine;
    Eigen::MatrixXi F_fine, E_fine;
    igl::upsample(V_cube, F_cube, V_fine, F_fine, 1);
    igl::edges(F_fine, E_fine);

    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        { V_fine, V_cube }, { E_fine, E_cube }, { F_fine, F_cube },
        { 1000.0, 1000.0 }, initial_poses);

    std::vector<Pose> poses_t0 = initial_poses;
    std::vector<Pose> poses_t1 = initial_poses;

    SECTION("translation")
    {
        poses_t0[1].position = Eigen::Vector3d(2.5, 0, 0);
        poses_t1[1].position = Eigen::Vector3d(1.05, 0, 0);
    }
    SECTION("rotation and translation")
    {
        poses_t0[1].position = Eigen::Vector3d(1.5, 1.5, 0);
        poses_t1[1].position = Eigen::Vector3d(1.05, 0.25, 0);
        poses_t1[1].rotation = compose_rotation(
            Eigen::AngleAxisd(igl::PI / 2, Eigen::Vector3d::UnitZ())
                .toRotationMatrix(),
            poses_t0[1]);
        poses_t1[0].rotation = compose_rotation(
            Eigen::AngleAxisd(igl::PI / 4, Eigen::Vector3d::UnitY())
                .toRotationMatrix(),
            poses_t0[0]);
    }

    const double inflation_radius = 0.1;

    RigidCandidates candidates;
    candidates.build(*bodies, poses_t0, poses_t1, inflation_radius);

    std::set<std::pair<index_t, index_t>> ee_pairs, fv_pairs;
    sampled_close_pairs(
        *bodies, poses_t0, poses_t1, inflation_radius, ee_pairs, fv_pairs);

    REQUIRE((ee_pairs.size() > 0 || fv_pairs.size() > 0));

    std::set<std::pair<index_t, index_t>> found_ee, found_fv;
    for (const EdgeEdgeCandidate& c : candidates.ee_candidates) {
        found_ee.emplace(
            std::min(c.edge0_id, c.edge1_id), std::max(c.edge0_id, c.edge1_id));
    }
    for (const FaceVertexCandidate& c : candidates.fv_candidates) {
        found_fv.emplace(c.face_id, c.vertex_id);
    }

    for (const auto& p : ee_pairs) {
        CHECK(found_ee.count(p) == 1);
    }
    for (const auto& p : fv_pairs) {
        CHECK(found_fv.count(p) == 1);
    }
}

TEST_CASE(
    "Rigid broad phase catches rotational tunneling", "[rigid][candidates]")
{
    std::vector<Pose> poses_t0, poses_t1;
    auto bodies = make_spinning_rod_scene(poses_t0, poses_t1);

    const double inflation_radius = 0.01;

    // Sanity check the scene: the bodies are far apart at both endpoints but
    // touch mid-sweep.
    std::set<std::pair<index_t, index_t>> ee_pairs, fv_pairs;
    sampled_close_pairs(
        *bodies, poses_t0, poses_t1, inflation_radius, ee_pairs, fv_pairs,
        /*n_samples=*/1); // only t=0 and t=1
    REQUIRE(ee_pairs.empty());
    REQUIRE(fv_pairs.empty());
    REQUIRE(sampled_first_contact(*bodies, poses_t0, poses_t1) < 1.0);

    RigidCandidates candidates;
    candidates.build(*bodies, poses_t0, poses_t1, inflation_radius);
    CHECK(candidates.size() > 0);

    // Completeness: all sampled close pairs are found. The rod tip sweeps
    // quickly, so it is only near the cube for a short time window — sample
    // finely enough to catch it.
    sampled_close_pairs(
        *bodies, poses_t0, poses_t1, inflation_radius, ee_pairs, fv_pairs,
        /*n_samples=*/1001);
    REQUIRE((ee_pairs.size() > 0 || fv_pairs.size() > 0));

    std::set<std::pair<index_t, index_t>> found_ee, found_fv;
    for (const EdgeEdgeCandidate& c : candidates.ee_candidates) {
        found_ee.emplace(
            std::min(c.edge0_id, c.edge1_id), std::max(c.edge0_id, c.edge1_id));
    }
    for (const FaceVertexCandidate& c : candidates.fv_candidates) {
        found_fv.emplace(c.face_id, c.vertex_id);
    }
    for (const auto& p : ee_pairs) {
        CHECK(found_ee.count(p) == 1);
    }
    for (const auto& p : fv_pairs) {
        CHECK(found_fv.count(p) == 1);
    }
}

TEST_CASE("Rigid narrow phase catches tunneling", "[rigid][candidates][ccd]")
{
    std::vector<Pose> poses_t0, poses_t1;
    auto bodies = make_spinning_rod_scene(poses_t0, poses_t1);

    RigidCandidates candidates;
    candidates.build(*bodies, poses_t0, poses_t1, /*inflation_radius=*/0.01);
    REQUIRE(candidates.size() > 0);

    const double first_contact =
        sampled_first_contact(*bodies, poses_t0, poses_t1);
    REQUIRE(first_contact < 1.0);

    // Linearized CCD misses the collision entirely: the endpoint positions of
    // the rod hug its initial axis, so lerped vertex trajectories never
    // approach the cube.
    const double linear_alpha =
        candidates.Candidates::compute_collision_free_stepsize(
            *bodies, bodies->vertices(poses_t0), bodies->vertices(poses_t1));
    CHECK(linear_alpha == 1.0);

    // Rigid CCD stops before contact and near the true time of impact.
    // (conservative_rescaling controls how early the CCD conservatively
    // declares impact; 0.9 matches the toolkit's nonlinear CCD tests.)
    NonlinearCCD ccd;
    ccd.conservative_rescaling = 0.9;
    const double alpha = candidates.compute_collision_free_stepsize(
        *bodies, poses_t0, poses_t1, /*min_distance=*/0.0, ccd);
    CHECK(alpha < 1.0);
    CHECK(alpha <= first_contact);
    CHECK(alpha == Catch::Approx(first_contact).margin(0.1));
    CHECK(!candidates.is_step_collision_free(
        *bodies, poses_t0, poses_t1, /*min_distance=*/0.0, ccd));

    // The step to the collision-free fraction must itself be collision free.
    const std::vector<Pose> poses_alpha = lerp_poses(poses_t0, poses_t1, alpha);
    CHECK(candidates.is_step_collision_free(
        *bodies, poses_t0, poses_alpha, /*min_distance=*/0.0, ccd));
}

TEST_CASE("Rigid plane-vertex candidates and CCD", "[rigid][candidates][ccd]")
{
    // A rod spinning half a revolution above a ground plane: the tip sweeps
    // below the plane mid-rotation even though the endpoint poses are clear.
    Eigen::MatrixXd V_rod;
    Eigen::MatrixXi E_rod, F_rod;
    make_rod(V_rod, E_rod, F_rod);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    initial_poses[0].position = Eigen::Vector3d(0, 1.0, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V_rod }, { E_rod }, { F_rod }, { 1000.0 }, initial_poses);
    bodies->planes.emplace_back(Eigen::Vector3d(0, 1, 0), 0); // ground y = 0

    std::vector<Pose> poses_t0 = initial_poses;
    std::vector<Pose> poses_t1 = initial_poses;
    poses_t1[0].rotation = compose_rotation(
        Eigen::AngleAxisd(igl::PI, Eigen::Vector3d::UnitZ()).toRotationMatrix(),
        poses_t0[0]);

    // Sampled ground truth: first time any vertex dips below the plane
    double first_contact = 1.0;
    for (int s = 0; s <= 10'000; s++) {
        const double t = s / 10'000.0;
        const Eigen::MatrixXd V =
            bodies->vertices(lerp_poses(poses_t0, poses_t1, t));
        if ((V.col(1).array() <= 0.0).any()) {
            first_contact = t;
            break;
        }
    }
    REQUIRE(first_contact < 1.0);

    // Endpoints are clear of the plane
    REQUIRE((bodies->vertices(poses_t0).col(1).array() > 0.0).all());
    REQUIRE((bodies->vertices(poses_t1).col(1).array() > 0.0).all());

    RigidCandidates candidates;
    candidates.build(*bodies, poses_t0, poses_t1, /*inflation_radius=*/0.01);
    REQUIRE(candidates.size() > 0);

    NonlinearCCD ccd;
    ccd.conservative_rescaling = 0.9;
    const double alpha = candidates.compute_collision_free_stepsize(
        *bodies, poses_t0, poses_t1, /*min_distance=*/0.0, ccd);
    CHECK(alpha < 1.0);
    CHECK(alpha <= first_contact);
    CHECK(alpha == Catch::Approx(first_contact).margin(0.1));

    const std::vector<Pose> poses_alpha = lerp_poses(poses_t0, poses_t1, alpha);
    CHECK(candidates.is_step_collision_free(
        *bodies, poses_t0, poses_alpha, /*min_distance=*/0.0, ccd));
    CHECK((bodies->vertices(poses_alpha).col(1).array() > 0.0).all());
}

TEST_CASE(
    "Rigid broad phase with a reused broad-phase object", "[rigid][candidates]")
{
    // The simulator reuses a single broad-phase object for both the standard
    // vertex-level candidate builds (barrier potential) and the rigid CCD
    // build. Candidates::build sets broad_phase->can_vertices_collide to the
    // mesh's *vertex*-level filter and leaves it set. The rigid body-level
    // prefilter stores bodies as broad-phase "vertices", so applying the stale
    // vertex filter to *body* indices silently drops body pairs (low body
    // indices map to same-body vertex pairs, which the default filter
    // rejects). Regression test: this exact scenario caused a bunny to tunnel
    // through a bowl while every fresh-broad-phase test passed.
    Eigen::MatrixXd V_cube;
    Eigen::MatrixXi E_cube, F_cube;
    tests::load_mesh("cube.ply", V_cube, E_cube, F_cube);

    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(1.2, 0, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V_cube, V_cube }, { E_cube, E_cube }, { F_cube, F_cube },
        { 1000.0, 1000.0 }, initial_poses);

    std::vector<Pose> poses_t0 = initial_poses;
    std::vector<Pose> poses_t1 = initial_poses;
    poses_t1[1].position = Eigen::Vector3d(1.02, 0, 0);

    const double inflation_radius = 0.05;

    LBVH broad_phase;

    // Poison the filter exactly the way the simulator does: a standard
    // vertex-level build leaves mesh.can_collide set on the broad phase.
    Candidates linear_candidates;
    linear_candidates.build(
        *bodies, bodies->vertices(poses_t0), bodies->vertices(poses_t1),
        inflation_radius, &broad_phase);
    // The stale filter rejects vertices 0 and 1 (both in body 0), which are
    // the ids the body-level prefilter would feed it for the body pair (0, 1).
    REQUIRE(!broad_phase.can_vertices_collide(0, 1));

    RigidCandidates candidates;
    candidates.build(
        *bodies, poses_t0, poses_t1, inflation_radius, &broad_phase);
    CHECK(candidates.size() > 0); // With the bug: 0 body-body candidates

    // The caller's vertex-level filter must be restored afterwards.
    CHECK(!broad_phase.can_vertices_collide(0, 1));            // same body
    CHECK(broad_phase.can_vertices_collide(0, V_cube.rows())); // cross body
}

TEST_CASE(
    "Rigid broad phase honors vertex-level collision filter",
    "[rigid][candidates]")
{
    // Two overlapping cubes whose collisions are disabled by a vertex-patch
    // filter (the jointed-bodies use case): the rigid broad phase must apply
    // the collision mesh's can_collide in its two-tree queries, otherwise CCD
    // reports a false impact at t=0 and blocks all motion.
    Eigen::MatrixXd V_cube;
    Eigen::MatrixXi E_cube, F_cube;
    tests::load_mesh("cube.ply", V_cube, E_cube, F_cube);

    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(0.75, 0, 0); // Overlapping

    auto bodies = RigidBodies::build_from_meshes(
        { V_cube, V_cube }, { E_cube, E_cube }, { F_cube, F_cube },
        { 1000.0, 1000.0 }, initial_poses);

    std::vector<Pose> poses_t0 = initial_poses;
    std::vector<Pose> poses_t1 = initial_poses;
    poses_t1[1].position = Eigen::Vector3d(0.7, 0, 0);

    RigidCandidates candidates;
    candidates.build(*bodies, poses_t0, poses_t1, /*inflation_radius=*/0.01);
    CHECK(candidates.size() > 0); // Sanity: unfiltered pairs are found

    // Same patch for every vertex => every pair is blocked.
    bodies->can_collide = make_vertex_patches_filter(
        Eigen::VectorXi::Zero(bodies->num_vertices()));

    candidates.build(*bodies, poses_t0, poses_t1, /*inflation_radius=*/0.01);
    CHECK(candidates.empty());
}

TEST_CASE(
    "Rigid broad phase superset of linearized broad phase",
    "[rigid][candidates]")
{
    // For small rotations, the rigid broad phase must find at least the
    // candidates found by the linearized (vertex-trajectory) broad phase.
    Eigen::MatrixXd V_cube;
    Eigen::MatrixXi E_cube, F_cube;
    tests::load_mesh("cube.ply", V_cube, E_cube, F_cube);

    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(1.2, 0, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V_cube, V_cube }, { E_cube, E_cube }, { F_cube, F_cube },
        { 1000.0, 1000.0 }, initial_poses);

    std::vector<Pose> poses_t0 = initial_poses;
    std::vector<Pose> poses_t1 = initial_poses;
    poses_t1[1].position = Eigen::Vector3d(1.05, 0.05, 0);
    poses_t1[1].rotation = compose_rotation(
        Eigen::AngleAxisd(0.01, Eigen::Vector3d::UnitZ()).toRotationMatrix(),
        poses_t0[1]);

    const double inflation_radius = 0.05;

    RigidCandidates rigid_candidates;
    rigid_candidates.build(*bodies, poses_t0, poses_t1, inflation_radius);

    Candidates linear_candidates;
    linear_candidates.build(
        *bodies, bodies->vertices(poses_t0), bodies->vertices(poses_t1),
        inflation_radius);

    std::set<std::pair<index_t, index_t>> found_ee, found_fv;
    for (const EdgeEdgeCandidate& c : rigid_candidates.ee_candidates) {
        found_ee.emplace(
            std::min(c.edge0_id, c.edge1_id), std::max(c.edge0_id, c.edge1_id));
    }
    for (const FaceVertexCandidate& c : rigid_candidates.fv_candidates) {
        found_fv.emplace(c.face_id, c.vertex_id);
    }

    const Eigen::MatrixXi& E = bodies->edges();
    const Eigen::MatrixXi& F = bodies->faces();

    int num_interbody = 0;
    for (const EdgeEdgeCandidate& c : linear_candidates.ee_candidates) {
        if (bodies->vertex_to_body(E(c.edge0_id, 0))
            == bodies->vertex_to_body(E(c.edge1_id, 0))) {
            continue; // Rigid candidates exclude intra-body pairs
        }
        num_interbody++;
        CHECK(
            found_ee.count(
                { std::min(c.edge0_id, c.edge1_id),
                  std::max(c.edge0_id, c.edge1_id) })
            == 1);
    }
    for (const FaceVertexCandidate& c : linear_candidates.fv_candidates) {
        if (bodies->vertex_to_body(F(c.face_id, 0))
            == bodies->vertex_to_body(c.vertex_id)) {
            continue;
        }
        num_interbody++;
        CHECK(found_fv.count({ c.face_id, c.vertex_id }) == 1);
    }
    CHECK(num_interbody > 0);
}
