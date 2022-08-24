#include <catch2/catch.hpp>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/utils/intersection.hpp>

#include "test_utils.hpp"

#include <tbb/parallel_sort.h>

#include <igl/predicates/segment_segment_intersect.h>
#include <igl/edges.h>

#include <Eigen/Geometry>

using namespace ipc;

Eigen::MatrixXi remove_faces_with_degenerate_edges(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    Eigen::MatrixXi F_new(F.rows(), F.cols());
    int num_faces = 0;
    for (int i = 0; i < F.rows(); ++i) {
        if (V.row(F(i, 0)) != V.row(F(i, 1)) && V.row(F(i, 0)) != V.row(F(i, 2))
            && V.row(F(i, 1)) != V.row(F(i, 2))) {
            F_new.row(num_faces++) = F.row(i);
        }
    }
    F_new.conservativeResize(num_faces, F.cols());
    return F_new;
}

bool combine_meshes(
    const std::string& mesh1_name,
    const std::string& mesh2_name,
    const Eigen::Matrix3d& R1,
    const Eigen::Matrix3d& R2,
    int dim,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& E,
    Eigen::MatrixXi& F)
{
    Eigen::MatrixXd V1, V2;
    Eigen::MatrixXi E1, E2, F1, F2;

    bool success =
        load_mesh(mesh1_name, V1, E1, F1) && load_mesh(mesh2_name, V2, E2, F2);
    if (!success) {
        return false;
    }

    V1 = V1 * R1.transpose(); // (RVᵀ)ᵀ = VRᵀ
    V2 = V2 * R2.transpose();

    REQUIRE(dim <= V1.cols());
    REQUIRE(dim <= V2.cols());
    V = Eigen::MatrixXd(V1.rows() + V2.rows(), dim);
    V.topRows(V1.rows()) = V1.leftCols(dim);
    V.bottomRows(V2.rows()) = V2.leftCols(dim);

    REQUIRE(E1.cols() == E2.cols());
    E = Eigen::MatrixXi(E1.rows() + E2.rows(), E1.cols());
    E.topRows(E1.rows()) = E1;
    E.bottomRows(E2.rows()) = E2;
    E.bottomRows(E2.rows()).array() += V1.rows();

    REQUIRE(F1.cols() == F2.cols());
    F = Eigen::MatrixXi(F1.rows() + F2.rows(), F1.cols());
    F.topRows(F1.rows()) = F1;
    F.bottomRows(F2.rows()) = F2;
    F.bottomRows(F2.rows()).array() += V1.rows();

    F = remove_faces_with_degenerate_edges(V, F);
    if (F.rows() != 0) {
        igl::edges(F, E);
    }

    return true;
}

// TEST_CASE("Test HashGrid EF vs brute-force", "[intersection][brute_force]")
// {
// #ifdef NDEBUG
//     std::string mesh1_name = GENERATE("cube.obj", "bunny.obj");
//     std::string mesh2_name = GENERATE("cube.obj", "bunny.obj");
// #else
//     std::string mesh1_name = "cube.obj";
//     std::string mesh2_name = "cube.obj";
// #endif
//     int dim = GENERATE(2, 3);
//     bool perform_aabb_check = true;

//     Eigen::Matrix3d R1 = GENERATE(take(4, RotationGenerator::create()));
//     Eigen::Matrix3d R2 = GENERATE(take(4, RotationGenerator::create()));

//     Eigen::MatrixXd V;
//     Eigen::MatrixXi E, F;
//     bool success = combine_meshes(mesh1_name, mesh2_name, R1, R2, dim, V, E,
//     F); REQUIRE(success);

//     HashGrid hash_grid;
//     hash_grid.resize(V, E);
//     hash_grid.addEdges(V, E);
//     if (V.cols() == 3) {
//         // These are not needed for 2D
//         hash_grid.addFaces(V, F);
//     }

//     if (V.cols() == 2) { // Need to check segment-segment intersections in 2D
//         std::vector<EdgeEdgeCandidate> hg_ee_candidates;
//         hash_grid.getEdgeEdgePairs(E, hg_ee_candidates);

//         std::vector<EdgeEdgeCandidate> bf_ee_candidates;
//         detect_edge_edge_collision_candidates_brute_force(
//             V, E, bf_ee_candidates, perform_aabb_check);

//         REQUIRE(hg_ee_candidates.size() <= bf_ee_candidates.size());

//         tbb::parallel_sort(hg_ee_candidates.begin(), hg_ee_candidates.end());
//         tbb::parallel_sort(bf_ee_candidates.begin(), bf_ee_candidates.end());

//         // narrow-phase using igl
//         igl::predicates::exactinit();
//         int hg_ci = 0;
//         for (int bf_ci = 0; bf_ci < bf_ee_candidates.size(); bf_ci++) {
//             if (hg_ee_candidates.size() <= hg_ci
//                 || bf_ee_candidates[bf_ci] != hg_ee_candidates[hg_ci]) {
//                 long eai = bf_ee_candidates[bf_ci].edge0_index;
//                 long ebi = bf_ee_candidates[bf_ci].edge1_index;
//                 bool intersects = igl::predicates::segment_segment_intersect(
//                     V.row(E(eai, 0)).head<2>(), V.row(E(eai, 1)).head<2>(),
//                     V.row(E(ebi, 0)).head<2>(), V.row(E(ebi, 1)).head<2>());

//                 CHECK(!intersects); // Check for FN

//             } else {
//                 hg_ci++;
//             }
//         }
//         CHECK(hg_ci >= hg_ee_candidates.size());
//     } else { // Need to check segment-triangle intersections in 3D
//         assert(V.cols() == 3);

//         std::vector<EdgeFaceCandidate> hg_ef_candidates;
//         hash_grid.getEdgeFacePairs(E, F, hg_ef_candidates);

//         std::vector<EdgeFaceCandidate> bf_ef_candidates;
//         detect_edge_face_collision_candidates_brute_force(
//             V, E, F, bf_ef_candidates, perform_aabb_check);

//         REQUIRE(hg_ef_candidates.size() <= bf_ef_candidates.size());

//         tbb::parallel_sort(hg_ef_candidates.begin(), hg_ef_candidates.end());
//         tbb::parallel_sort(bf_ef_candidates.begin(), bf_ef_candidates.end());

//         // narrow-phase using igl
//         int hg_ci = 0;
//         for (int bf_ci = 0; bf_ci < bf_ef_candidates.size(); bf_ci++) {
//             if (hg_ef_candidates.size() <= hg_ci
//                 || bf_ef_candidates[bf_ci] != hg_ef_candidates[hg_ci]) {
//                 long ei = bf_ef_candidates[bf_ci].edge_index;
//                 long fi = bf_ef_candidates[bf_ci].face_index;
//                 bool intersects = is_edge_intersecting_triangle(
//                     V.row(E(ei, 0)), V.row(E(ei, 1)), V.row(F(fi, 0)),
//                     V.row(F(fi, 1)), V.row(F(fi, 2)));

//                 CHECK(!intersects); // Check for FN
//             } else {
//                 hg_ci++;
//             }
//         }
//         CHECK(hg_ci >= hg_ef_candidates.size());
//     }

//     CAPTURE(mesh1_name, mesh2_name, R1, R2);
// }

TEST_CASE("Test has_intersections()", "[intersection][thisone]")
{
    std::string mesh1_name = GENERATE("cube.obj", "bunny.obj");
    std::string mesh2_name = GENERATE("cube.obj", "bunny.obj");
    int dim = GENERATE(2, 3);

#ifdef NDEBUG
    Eigen::Matrix3d R1 = GENERATE(take(4, RotationGenerator::create()));
    Eigen::Matrix3d R2 = GENERATE(take(4, RotationGenerator::create()));
#else
    Eigen::Matrix3d R1 = GENERATE(take(2, RotationGenerator::create()));
    Eigen::Matrix3d R2 = GENERATE(take(2, RotationGenerator::create()));
#endif

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    bool success = combine_meshes(mesh1_name, mesh2_name, R1, R2, dim, V, E, F);
    REQUIRE(success);

    CAPTURE(mesh1_name, mesh2_name, R1, R2);
    CHECK(has_intersections(CollisionMesh(V, E, F), V));
}
