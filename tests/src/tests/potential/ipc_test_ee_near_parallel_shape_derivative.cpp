// Regression test for the edge-edge parallel mollifier shape-derivative path.
//
// Reproduces a single mollified EE pair (near-parallel edges, distance near
// dhat) and FD-checks the analytic shape_derivative against a centered
// difference of `gradient` at fixed collision set — the same probe IPC's own
// shape_derivative test uses (test_barrier_potential.cpp:377-391).
//
// Prior to the fix in edge_edge_mollifier.cpp, the assertion failed by
// roughly two orders of magnitude on these geometries: the analytic
// per-pair shape_derivative was inflated when (a) the EE pair is mollified
// (cross_squared_norm < eps_x), AND (b) the closest distance is near dhat.
// The root cause was `edge_edge_mollifier_gradient_wrt_x` returning
// ∂m/∂eps_x · ∇_position m instead of ∂m/∂eps_x · ∇_rest eps_x. The fix
// uses edge_edge_mollifier_threshold_gradient as the second factor, matching
// what edge_edge_mollifier_gradient_jacobian_wrt_x already does.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <Eigen/Core>

using namespace ipc;

namespace {

// Build a 2-edge / 4-vertex CollisionMesh: two unit-length parallel-ish edges
// in 3D, separated in y by `y_sep`, with a small z-twist `z_twist` on the
// second edge to keep the cross product slightly nonzero (so the mollifier is
// triggered but s = |A x B|^2 != 0).
//
//   ea0 = (0, 0,        0)        ea1 = (1, 0,        0)
//   eb0 = (0, y_sep,    0)        eb1 = (1, y_sep,    z_twist)
//
// With z_twist > 0 small, |A x B|^2 = z_twist^2 > 0 but << eps_x = 1e-3
// (since both edges have length ~1), so the mollifier branch fires. The
// closest-point distance between the edges stays near y_sep.
struct EEScene {
    Eigen::MatrixXd rest;
    Eigen::MatrixXd displaced;
    Eigen::MatrixXi edges;
    Eigen::MatrixXi faces;
};

EEScene make_near_parallel_ee_scene(double y_sep, double z_twist)
{
    EEScene s;
    s.rest.resize(4, 3);
    s.rest << 0.0, 0.0,     0.0,
              1.0, 0.0,     0.0,
              0.0, y_sep,   0.0,
              1.0, y_sep,   z_twist;
    // displaced = rest: the shape-derivative bug is independent of any
    // additional displacement, so we test the simplest configuration.
    s.displaced = s.rest;
    s.edges.resize(2, 2);
    s.edges << 0, 1,
               2, 3;
    s.faces.resize(0, 3);
    return s;
}

} // namespace

TEST_CASE(
    "EE shape_derivative on near-parallel mollified pair",
    "[potential][barrier_potential][shape_derivative][mollifier]")
{
    const double dhat = 1e-3;

    // y_sep just below dhat (active contact), z_twist makes |A x B|^2 small
    // but nonzero so the mollifier branch is entered cleanly.
    const double y_sep   = GENERATE(9.99e-4, 9.9e-4);
    const double z_twist = GENERATE(1e-3);

    EEScene scene = make_near_parallel_ee_scene(y_sep, z_twist);
    CollisionMesh mesh(scene.rest, scene.edges, scene.faces);
    if (!mesh.are_area_jacobians_initialized())
        mesh.init_area_jacobians();

    NormalCollisions collisions;
    collisions.set_enable_shape_derivatives(true);
    collisions.build(mesh, scene.displaced, dhat);
    REQUIRE(collisions.size() == 1);
    REQUIRE(collisions.is_edge_edge(0));
    REQUIRE(collisions[0].is_mollified());

    // kappa = 1e8 keeps entries well above the precision floor of the
    // centered-FD probe.
    BarrierPotential bp(dhat, /*stiffness=*/1e8, /*use_physical_barrier=*/false);

    // -------- Analytic 12x12 (collision-mesh DOF order) --------
    Eigen::SparseMatrix<double> JF_an_sparse =
        bp.shape_derivative(collisions, mesh, scene.displaced);
    Eigen::MatrixXd JF_an(JF_an_sparse);

    // -------- FD 12x12 column-by-column on rest DOFs, material
    // displacement held constant, collision set held fixed. Build in
    // ROW-MAJOR DOF order (v * dim + d) to match the layout returned by
    // shape_derivative. --------
    const int n_verts = static_cast<int>(scene.rest.rows());
    const int dim     = static_cast<int>(scene.rest.cols());
    const int ndof    = n_verts * dim;
    const double eps  = 1e-7;
    Eigen::MatrixXd JF_fd(ndof, ndof);
    JF_fd.setZero();
    for (int v = 0; v < n_verts; ++v) {
        for (int d = 0; d < dim; ++d) {
            Eigen::MatrixXd rest_plus  = scene.rest;
            Eigen::MatrixXd rest_minus = scene.rest;
            rest_plus(v, d)  += eps;
            rest_minus(v, d) -= eps;
            // Hold material displacement (= displaced - rest) constant.
            const Eigen::MatrixXd disp_plus =
                rest_plus  + (scene.displaced - scene.rest);
            const Eigen::MatrixXd disp_minus =
                rest_minus + (scene.displaced - scene.rest);
            CollisionMesh mesh_plus(rest_plus, scene.edges, scene.faces);
            if (mesh.are_area_jacobians_initialized())
                mesh_plus.init_area_jacobians();
            CollisionMesh mesh_minus(rest_minus, scene.edges, scene.faces);
            if (mesh.are_area_jacobians_initialized())
                mesh_minus.init_area_jacobians();
            const Eigen::VectorXd F_plus =
                bp.gradient(collisions, mesh_plus,  disp_plus);
            const Eigen::VectorXd F_minus =
                bp.gradient(collisions, mesh_minus, disp_minus);
            JF_fd.col(v * dim + d) = (F_plus - F_minus) / (2.0 * eps);
        }
    }

    const double an   = JF_an.norm();
    const double fn   = JF_fd.norm();
    const double diff = (JF_an - JF_fd).norm();
    const double rel  = diff / std::max(std::max(an, fn), 1e-30);

    UNSCOPED_INFO(
        "y_sep=" << y_sep << " z_twist=" << z_twist
        << "  ||analytic||=" << an
        << "  ||fd||=" << fn
        << "  ||diff||=" << diff
        << "  rel=" << rel);

    CHECK(rel < 5e-2);
}
