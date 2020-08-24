#include <catch2/catch.hpp>

#if defined(WIN32)
#include <filesystem>
#endif

#include <igl/dirname.h>
#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>

#include <finitediff.hpp>

#include <ipc.hpp>

using namespace ipc;

// Flatten the matrix rowwise
Eigen::VectorXd flatten(const Eigen::MatrixXd& X)
{
    Eigen::MatrixXd XT = X.transpose();
    return Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(XT.data(), XT.size()));
}

/// Unflatten rowwise
Eigen::MatrixXd unflatten(const Eigen::VectorXd& x, int dim)
{
    assert(x.size() % dim == 0);
    Eigen::MatrixXd unflat_x(x.size() / dim, dim);
    for (int i = 0; i < x.size(); i++) {
        unflat_x(i / dim, i % dim) = x(i);
    }
    return unflat_x;
}

bool load_mesh(
    const std::string& mesh_name,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
#if defined(WIN32)
    std::string mesh_path =
        (std::filesystem::path(__FILE__).parent_path() / "meshes" / mesh_name)
            .string();
#else
    std::string mesh_path =
        igl::dirname(std::string(__FILE__)) + "/meshes/" + mesh_name;
#endif
    bool success = igl::read_triangle_mesh(mesh_path, V, F);
    if (F.size()) {
        igl::edges(F, E);
    }
    return success && V.size() && F.size() && E.size();
}

TEST_CASE("Flatten and unflatten", "[utils]")
{
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(1000, 3);
    Eigen::MatrixXd R = unflatten(flatten(X), X.cols());
    CHECK(X == R);
}

TEST_CASE("Dummy test for IPC compilation", "[ipc]")
{
    double dhat_squared = -1;
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
    bool success = load_mesh(mesh_name, V, E, F);
    REQUIRE(success);

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

TEST_CASE("Test IPC full gradient", "[ipc][grad]")
{
    double dhat_squared = -1;
    std::string mesh_name;

    SECTION("cube")
    {
        dhat_squared = 2.0;
        mesh_name = "cube.obj";
    }
    // SECTION("bunny")
    // {
    //     dhat_squared = 1e-4;
    //     mesh_name = "bunny.obj";
    // }

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    bool success = load_mesh(mesh_name, V, E, F);
    REQUIRE(success);

    ccd::Candidates constraint_set;
    ipc::construct_constraint_set(V, E, F, dhat_squared, constraint_set);
    CAPTURE(mesh_name, dhat_squared);
    CHECK(constraint_set.ee_candidates.size() > 0);
    CHECK(constraint_set.fv_candidates.size() > 0);

    Eigen::VectorXd grad_b = ipc::compute_barrier_potential_gradient(
        V, V, E, F, constraint_set, dhat_squared);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return ipc::compute_barrier_potential(
            V, unflatten(x, V.cols()), E, F, constraint_set, dhat_squared);
    };
    Eigen::VectorXd fgrad_b;
    fd::finite_gradient(flatten(V), f, fgrad_b);

    REQUIRE(grad_b.squaredNorm() > 1e-8);
    CHECK(fd::compare_gradient(grad_b, fgrad_b));
}

// TEST_CASE("Test IPC full hessian", "[ipc][hess]")
// {
//     double dhat_squared = -1;
//     std::string mesh_name;
//
//     // SECTION("cube")
//     // {
//     dhat_squared = 2.0;
//     mesh_name = "cube.obj";
//     // }
//     // SECTION("bunny")
//     // {
//     //     dhat_squared = 1e-4;
//     //     mesh_name = "bunny.obj";
//     // }
//
//     Eigen::MatrixXd V;
//     Eigen::MatrixXi E, F;
//     bool success = load_mesh(mesh_name, V, E, F);
//     REQUIRE(success);
//
//     ccd::Candidates constraint_set;
//     ipc::construct_constraint_set(V, E, F, dhat_squared, constraint_set);
//     CAPTURE(mesh_name, dhat_squared);
//     CHECK(constraint_set.ee_candidates.size() > 0);
//     CHECK(constraint_set.fv_candidates.size() > 0);
//
//     size_t num_fv_candidates = constraint_set.fv_candidates.size();
//     int saved_id = GENERATE_COPY(range(size_t(0), num_fv_candidates));
//     auto saved_candidate = constraint_set.fv_candidates[saved_id];
//     constraint_set.ee_candidates.clear();
//     constraint_set.fv_candidates.clear();
//     constraint_set.fv_candidates.push_back(saved_candidate);
//
//     Eigen::MatrixXd hess_b = ipc::compute_barrier_potential_hessian(
//         /*V_rest=*/V, V, E, F, constraint_set, dhat_squared);
//
//     // Compute the gradient using finite differences
//     auto f = [&](const Eigen::VectorXd& x) {
//         return ipc::compute_barrier_potential_gradient(
//             V, unflatten(x, V.cols()), E, F, constraint_set, dhat_squared);
//     };
//     Eigen::MatrixXd fhess_b;
//     fd::finite_jacobian(flatten(V), f, fhess_b);
//
//     REQUIRE(hess_b.squaredNorm() > 1e-8);
//     if (!fd::compare_hessian(hess_b, fhess_b, 1e-3)) {
//         std::cout << (hess_b - fhess_b).lpNorm<Eigen::Infinity>() <<
//         std::endl; std::cout
//             <<
//             "hess_b-------------------------------------------------------"
//                "-------------------------------------------------------------"
//                "---------"
//             << std::endl;
//         std::cout << hess_b << std::endl;
//         std::cout
//             <<
//             "fhess_b------------------------------------------------------"
//                "-------------------------------------------------------------"
//                "---------"
//             << std::endl;
//         std::cout << fhess_b << std::endl;
//     }
//     CHECK(fd::compare_hessian(hess_b, fhess_b, 1e-3));
// }
