#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <ipc/smooth_contact/smooth_contact_potential.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>
#include <igl/readCSV.h>
#include <ipc/ipc.hpp>

#include <ipc/utils/getRSS.h>

using namespace ipc;

// TEST_CASE(
//     "Smooth barrier potential convergence 2D",
//     "[smooth_potential]")
// {
//     const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID;
//     const bool adaptive_dhat = false;

//     double dhat = -1;
//     std::string mesh_name = "";
//     SECTION("cusp")
//     {
//         mesh_name = ("/home/zizhou/polyfem/result/smooth-ipc/examples/cusp-2d/out_surf_contact.vtu");
//         dhat = 0.2;
//     }

//     SECTION("convex")
//     {
//         mesh_name = ("/home/zizhou/polyfem/result/smooth-ipc/examples/visualize-contact/2D/circle-bar/out_surf_contact.vtu");
//         dhat = 0.5;
//     }

//     double min_dist_ratio = 0.8;
//     Eigen::MatrixXd vertices;
//     Eigen::MatrixXi edges, faces;
//     std::cout << mesh_name << "\n";
//     bool success = igl::readCSV(mesh_name + "-v.csv", vertices);
//     success = success && igl::readCSV(mesh_name + "-e.csv", edges);
//     CAPTURE(mesh_name);
//     REQUIRE(success);

//     CollisionMesh mesh;
//     ParameterType param(dhat*dhat, 0.8, 0.0, 0.1, 0.0, 1);
//     param.set_adaptive_dhat_ratio(min_dist_ratio);
//     SmoothCollisions<2> collisions;
//     mesh = CollisionMesh(vertices, edges, faces);
//     collisions.compute_adaptive_dhat(mesh, vertices, param, method);
//     collisions.build(mesh, vertices, param, adaptive_dhat, method);
//     CAPTURE(dhat, method, adaptive_dhat);
//     CHECK(collisions.size() > 0);
//     std::cout << "smooth collision candidate size " << collisions.size() << "\n";

//     CHECK(!has_intersections(mesh, vertices));

//     SmoothContactPotential<SmoothCollisions<2>> potential(param);
//     std::cout << "energy: " << potential(collisions, mesh, vertices) << "\n";

//     // -------------------------------------------------------------------------
//     // Gradient
//     // -------------------------------------------------------------------------

//     // const Eigen::VectorXd grad_b =
//     //     potential.gradient(collisions, mesh, vertices);
//     // std::cout << "grad: " << grad_b.norm() << "\n";

//     // -------------------------------------------------------------------------
//     // Hessian
//     // -------------------------------------------------------------------------

//     // Eigen::MatrixXd hess_b =
//     //     potential.hessian(collisions, mesh, vertices);
//     // std::cout << "hess: " << hess_b.norm() << "\n";
// }

// TEST_CASE(
//     "Smooth barrier potential convergence 3D",
//     "[smooth_potential]")
// {
//     const BroadPhaseMethod method = BroadPhaseMethod::BVH;
//     const bool adaptive_dhat = false;

//     double dhat = -1;
//     std::string mesh_name = "";
//     SECTION("debug")
//     {
//         mesh_name = ("2-cubes.obj");
//         dhat = 0.25;
//     }

//     double min_dist_ratio = 0.8;
//     Eigen::MatrixXd vertices;
//     Eigen::MatrixXi edges, faces;
//     bool success = tests::load_mesh(mesh_name, vertices, edges, faces);
//     CAPTURE(mesh_name);
//     REQUIRE(success);

//     CollisionMesh mesh;
//     ParameterType param(dhat, 0.8, 0, 0.1, 0, 2);
//     param.set_adaptive_dhat_ratio(min_dist_ratio);
//     SmoothCollisions<3> collisions;
//     mesh = CollisionMesh(vertices, edges, faces);
//     std::cout << "Mesh Memory " << getCurrentRSS() / (1024.*1024) << "MB\n";

//     // ignore self-contact
//     {
//         std::vector<int> vids(vertices.rows(), 0);
//         for (int i = 0; i < vertices.rows(); i++)
//             if (vertices(i, 0) < -1.1)
//                 vids[i] = 1;
//         mesh.can_collide = [vids](size_t vi, size_t vj) {
//             return vids[vi] != vids[vj];
//         };
//     }

//     // collisions.compute_adaptive_dhat(mesh, vertices, param, method);
//     collisions.build(mesh, vertices, param, adaptive_dhat, method);
//     std::cout << "Collision set Memory " << getCurrentRSS() / (1024.*1024) << "MB" << std::endl;
//     CAPTURE(dhat, method, adaptive_dhat);
//     CHECK(collisions.size() > 0);
//     std::cout << "smooth collision candidate size " << collisions.size() << std::endl;

//     CHECK(!has_intersections(mesh, vertices));

//     SmoothContactPotential<SmoothCollisions<3>> potential(param);
//     std::cout << std::setprecision(12) << "energy: " << potential(collisions, mesh, vertices) << std::endl;

//     // -------------------------------------------------------------------------
//     // Gradient
//     // -------------------------------------------------------------------------

//     // const Eigen::VectorXd grad_b =
//     //     potential.gradient(collisions, mesh, vertices);
//     // std::cout << "grad: " << grad_b.norm() << "\n";
// }
