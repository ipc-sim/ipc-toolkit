#include <catch2/catch.hpp>

#include <igl/edges.h>
#include <igl/pathinfo.h>
#include <igl/read_triangle_mesh.h>

#include <ipc.hpp>

using namespace ipc;

TEST_CASE("Dummy test for IPC compilation", "[ipc]")
{
    double dhat_squared;
    std::string mesh_name;

    SECTION("cube")
    {
        dhat_squared = 2.0;
        mesh_name = "cube.obj";
    }
    SECTION("bunny")
    {
        dhat_squared = 1e-4;
        mesh_name = "bunny.obj";
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    std::string dirname, basename, extension, filename;
    igl::pathinfo(
        std::string(__FILE__), dirname, basename, extension, filename);
#ifdef win32
    igl::read_triangle_mesh(dirname + "\\meshes\\" + mesh_name, V, F);
#else
    igl::read_triangle_mesh(dirname + "/meshes/" + mesh_name, V, F);
#endif
    REQUIRE(V.size());
    REQUIRE(F.size());
    igl::edges(F, E);
    REQUIRE(E.size());

    ccd::Candidates constraint_set;
    ipc::construct_constraint_set(V, E, F, dhat_squared, constraint_set);
    CAPTURE(mesh_name, dhat_squared);
    CHECK(constraint_set.ee_candidates.size() > 0);
    CHECK(constraint_set.fv_candidates.size() > 0);

    double b = ipc::compute_barrier_potential(
        V, V, E, F, constraint_set, dhat_squared);
    Eigen::VectorXd grad_b = ipc::compute_barrier_potential_gradient(
        V, V, E, F, constraint_set, dhat_squared);
    Eigen::MatrixXd hess_b = ipc::compute_barrier_potential_hessian(
        V, V, E, F, constraint_set, dhat_squared);
}
