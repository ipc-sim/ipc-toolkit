// Correctness tests for Potential<T>::gradient's assembly.
//
// Which of the two assembly strategies runs is a function of the scene's
// shape, so the fixture below covers both by sizing scenes on either side of
// the crossover: scenes with more DOF than stencil slots take the serial
// scatter, the rest take the per-thread accumulators. That makes the coverage
// implicit, hence "Gradient assembly test scenes cover both branches", which
// fails if a future edit to the scene sizes strands one strategy untested.

#include "assembly_scene.hpp"
#include "synthetic_contact_scene.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/utils/local_to_global.hpp>

#include <tbb/global_control.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <thread>
#include <vector>

using namespace ipc;
using namespace ipc::tests;

namespace {

/// @brief Whether a scene's collision-DOF assembly takes the two-pass path.
bool is_two_pass(const AssemblyScene& scene)
{
    return assembly_is_two_pass(scene.num_collisions(), scene.ndof());
}

/// @brief Collision-DOF gradient assembled by a single-threaded scatter in
/// collision order. The ground truth every assembly is compared against.
Eigen::VectorXd serial_reference_gradient(const AssemblyScene& scene)
{
    const Eigen::MatrixXi& edges = scene.mesh().edges();
    const Eigen::MatrixXi& faces = scene.mesh().faces();
    const Eigen::MatrixXd& X = scene.vertices();
    const int dim = X.cols();

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(X.size());
    for (size_t i = 0; i < scene.collisions().size(); i++) {
        const NormalCollision& collision = scene.collisions()[i];
        const VectorMax12d local_grad = scene.potential().gradient(
            collision, collision.dof(X, edges, faces));
        local_gradient_to_global_gradient(
            local_grad, collision.vertex_ids(edges, faces), dim, grad);
    }
    return grad;
}

/// @brief The scenes shared by the assembly tests, built once.
///
/// Spans: a real mesh (reusing the assembly benchmark fixture), both assembly
/// branches, both 3D collision types, 2D, and a displacement-mapped mesh
/// (forces the map_to_full path).
const std::vector<AssemblyScene>& assembly_test_scenes()
{
    static const std::vector<AssemblyScene> scenes = [] {
        std::vector<AssemblyScene> s;

        // Real mesh scene; skip silently if the mesh is unavailable.
        std::optional<AssemblyScene> real =
            build_assembly_scene(assembly_scene_specs().front());
        if (real.has_value()) {
            s.push_back(std::move(*real));
        }

        SyntheticContactSpec spec;
        spec.seed = 42;

        // 3D face-vertex, thread-local branch (num_slots = 8000 > ndof).
        spec.num_collisions = 2000;
        spec.target_ndof = 1500;
        spec.dim = 3;
        spec.type = SyntheticContactType::FACE_VERTEX;
        s.push_back(build_synthetic_contact_scene(spec));

        // 3D face-vertex, two-pass branch (ndof = 30000 > num_slots = 800).
        spec.num_collisions = 200;
        spec.target_ndof = 30000;
        s.push_back(build_synthetic_contact_scene(spec));

        // 3D edge-edge (mollified local gradients).
        spec.num_collisions = 1000;
        spec.target_ndof = 6000;
        spec.type = SyntheticContactType::EDGE_EDGE;
        s.push_back(build_synthetic_contact_scene(spec));

        // 2D edge-vertex.
        spec.num_collisions = 1000;
        spec.target_ndof = 2000;
        spec.dim = 2;
        spec.type = SyntheticContactType::EDGE_VERTEX;
        s.push_back(build_synthetic_contact_scene(spec));

        // 3D face-vertex with a displacement map: in_full_dof goes through
        // map_to_full instead of folding.
        spec.num_collisions = 500;
        spec.target_ndof = 1500;
        spec.dim = 3;
        spec.type = SyntheticContactType::FACE_VERTEX;
        spec.with_displacement_map = true;
        s.push_back(build_synthetic_contact_scene(spec));

        return s;
    }();
    return scenes;
}

} // namespace

TEST_CASE(
    "Gradient assembly test scenes cover both branches",
    "[potential][assembly][gradient]")
{
    // Guards the fixture above: without this, a change to the scene sizes
    // could silently leave one assembly branch untested.
    const auto& scenes = assembly_test_scenes();
    CHECK(std::any_of(scenes.begin(), scenes.end(), [](const AssemblyScene& s) {
        return is_two_pass(s);
    }));
    CHECK(std::any_of(scenes.begin(), scenes.end(), [](const AssemblyScene& s) {
        return !is_two_pass(s);
    }));
}

TEST_CASE(
    "Gradient assembly matches the serial reference",
    "[potential][assembly][gradient]")
{
    const bool in_full_dof = GENERATE(false, true);

    for (const AssemblyScene& scene : assembly_test_scenes()) {
        CAPTURE(
            in_full_dof, scene.label(), scene.num_collisions(), scene.ndof(),
            is_two_pass(scene));

        Eigen::VectorXd ref = serial_reference_gradient(scene);
        if (in_full_dof) {
            ref = scene.mesh().to_full_dof(ref);
        }

        const Eigen::VectorXd grad = scene.potential().gradient(
            scene.collisions(), scene.mesh(), scene.vertices(), in_full_dof);

        REQUIRE(
            grad.size()
            == Eigen::Index(in_full_dof ? scene.full_ndof() : scene.ndof()));

        const double norm_scale = std::max(ref.norm(), 1.0);
        const double max_scale = std::max(ref.cwiseAbs().maxCoeff(), 1.0);
        CAPTURE(
            (grad - ref).norm() / norm_scale,
            (grad - ref).cwiseAbs().maxCoeff() / max_scale);
        CHECK((grad - ref).norm() <= 1e-14 * norm_scale);
        CHECK((grad - ref).cwiseAbs().maxCoeff() <= 1e-13 * max_scale);
    }
}

TEST_CASE(
    "Two-pass gradient assembly is bitwise deterministic",
    "[potential][assembly][gradient]")
{
    // ndof (30k) exceeds num_slots (20k), so this takes the two-pass branch,
    // whose fixed scatter order makes it reproducible bit for bit.
    SyntheticContactSpec spec;
    spec.num_collisions = 5000;
    spec.target_ndof = 30000;
    spec.dim = 3;
    spec.type = SyntheticContactType::FACE_VERTEX;
    spec.seed = 7;
    const AssemblyScene scene = build_synthetic_contact_scene(spec);
    REQUIRE(is_two_pass(scene));

    // The fixed scatter order must also match the serial reference exactly.
    const Eigen::VectorXd ref = serial_reference_gradient(scene);

    const std::size_t hw_threads =
        std::max(2U, std::thread::hardware_concurrency());
    for (const std::size_t num_threads : { std::size_t(2), hw_threads }) {
        CAPTURE(num_threads);
        tbb::global_control control(
            tbb::global_control::max_allowed_parallelism, num_threads);
        for (int rep = 0; rep < 5; rep++) {
            CAPTURE(rep);
            const Eigen::VectorXd grad = scene.potential().gradient(
                scene.collisions(), scene.mesh(), scene.vertices());
            REQUIRE(grad.size() == ref.size());
            REQUIRE(
                std::memcmp(
                    grad.data(), ref.data(), sizeof(double) * grad.size())
                == 0);
        }
    }
}
