#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <tbb/parallel_sort.h>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/ccd/ccd.hpp>

using namespace ipc;

TEST_CASE("Compare SpatialHash against brute force", "[spatial_hash]")
{
    using namespace ipc;

    Eigen::MatrixXd V0, U;
    Eigen::MatrixXi E, F;
    Eigen::VectorXi group_ids;

    SECTION("Simple")
    {
        V0.resize(4, 3);
        V0.row(0) << -1, -1, 0;
        V0.row(1) << 1, -1, 0;
        V0.row(2) << 0, 1, 1;
        V0.row(3) << 0, 1, -1;

        E.resize(2, 2);
        E.row(0) << 0, 1;
        E.row(1) << 2, 3;

        SECTION("Without group ids") {}
        // SECTION("With group ids")
        // {
        //     group_ids.resize(4);
        //     group_ids << 0, 0, 1, 1;
        // }

        F.resize(0, 3);

        U = Eigen::MatrixXd::Zero(V0.rows(), V0.cols());
        U.col(1).head(2).setConstant(2);
        U.col(1).tail(2).setConstant(-2);
    }
    SECTION("Complex")
    {
        // #ifdef NDEBUG
        //         std::string filename =
        //             GENERATE(std::string("cube.obj"),
        //             std::string("bunny.obj"));
        // #else
        std::string filename = "cube.obj";
        // #endif
        std::string mesh_path = std::string(TEST_DATA_DIR) + filename;
        bool success = igl::read_triangle_mesh(mesh_path, V0, F);
        REQUIRE(success);
        igl::edges(F, E);

        U = Eigen::MatrixXd::Zero(V0.rows(), V0.cols());
        U.col(1).setOnes();
    }

    SpatialHash sh;
    Candidates sh_candidates, bf_candidates;

    double inflation_radius = 1e-2; // GENERATE(0, 1e-4, 1e-3, 1e-2, 1e-1);

    auto can_collide = [&group_ids](size_t vi, size_t vj) {
        return group_ids.size() == 0 || group_ids(vi) != group_ids(vj);
    };

    for (int i = 0; i < 2; i++) {
        Eigen::MatrixXd V1 = V0 + U;
        sh.build(V0, V1, E, F, inflation_radius);

        sh_candidates.clear();
        // TODO: use can_collide
        sh.queryMeshForCandidates(
            V0, V1, E, F, sh_candidates,
            /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true);

        bf_candidates.clear();
        detect_collision_candidates_brute_force(
            V0, V1, E, F, bf_candidates,
            /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true,
            /*perform_aabb_check=*/false, inflation_radius,
            /*ignore_codimensional_vertices=*/false, can_collide);

        CHECK(
            sh_candidates.ev_candidates.size()
            <= bf_candidates.ev_candidates.size());
        CHECK(
            sh_candidates.ee_candidates.size()
            <= bf_candidates.ee_candidates.size());
        CHECK(
            sh_candidates.fv_candidates.size()
            <= bf_candidates.fv_candidates.size());

        tbb::parallel_sort(
            sh_candidates.ev_candidates.begin(),
            sh_candidates.ev_candidates.end());
        tbb::parallel_sort(
            bf_candidates.ev_candidates.begin(),
            bf_candidates.ev_candidates.end());
        int sh_ci = 0;
        for (int bf_ci = 0; bf_ci < bf_candidates.ev_candidates.size();
             bf_ci++) {
            if (sh_candidates.ev_candidates.size() <= sh_ci
                || bf_candidates.ev_candidates[bf_ci]
                    != sh_candidates.ev_candidates[sh_ci]) {

                long ei = bf_candidates.ev_candidates[bf_ci].edge_index;
                long vi = bf_candidates.ev_candidates[bf_ci].vertex_index;
                double toi;
                bool hit = point_edge_ccd(
                    V0.row(vi), V0.row(E(ei, 0)), V0.row(E(ei, 1)), // t = 0
                    V1.row(vi), V1.row(E(ei, 0)), V1.row(E(ei, 1)), // t = 1
                    toi, 1.0);
                CHECK(!hit); // Check for FN

            } else {
                sh_ci++;
            }
        }
        CHECK(sh_ci >= sh_candidates.ev_candidates.size());

        tbb::parallel_sort(
            sh_candidates.ee_candidates.begin(),
            sh_candidates.ee_candidates.end());
        tbb::parallel_sort(
            bf_candidates.ee_candidates.begin(),
            bf_candidates.ee_candidates.end());
        sh_ci = 0;
        for (int bf_ci = 0; bf_ci < bf_candidates.ee_candidates.size();
             bf_ci++) {
            if (sh_candidates.ee_candidates.size() <= sh_ci
                || bf_candidates.ee_candidates[bf_ci]
                    != sh_candidates.ee_candidates[sh_ci]) {

                long eai = bf_candidates.ee_candidates[bf_ci].edge0_index;
                long ebi = bf_candidates.ee_candidates[bf_ci].edge1_index;
                double toi;
                bool hit = edge_edge_ccd(
                    V0.row(E(eai, 0)), V0.row(E(eai, 1)), // t = 0
                    V0.row(E(ebi, 0)), V0.row(E(ebi, 1)), // t = 0
                    V1.row(E(eai, 0)), V1.row(E(eai, 1)), // t = 1
                    V1.row(E(ebi, 0)), V1.row(E(ebi, 1)), // t = 1
                    toi,
                    /*tmax=*/1.0,
                    /*tolerance=*/1e-6,
                    /*max_iterations=*/1e7,
                    /*conservative_rescaling=*/1.0);
                CHECK(!hit); // Check for FN

            } else {
                sh_ci++;
            }
        }
        CHECK(sh_ci >= sh_candidates.ee_candidates.size());

        tbb::parallel_sort(
            sh_candidates.fv_candidates.begin(),
            sh_candidates.fv_candidates.end());
        tbb::parallel_sort(
            bf_candidates.fv_candidates.begin(),
            bf_candidates.fv_candidates.end());
        sh_ci = 0;
        for (int bf_ci = 0; bf_ci < bf_candidates.fv_candidates.size();
             bf_ci++) {
            if (sh_candidates.fv_candidates.size() <= sh_ci
                || bf_candidates.fv_candidates[bf_ci]
                    != sh_candidates.fv_candidates[sh_ci]) {

                long fi = bf_candidates.fv_candidates[bf_ci].face_index;
                long vi = bf_candidates.fv_candidates[bf_ci].vertex_index;
                double toi;
                bool hit = point_triangle_ccd(
                    V0.row(vi),                                           //
                    V0.row(F(fi, 0)), V0.row(F(fi, 1)), V0.row(F(fi, 2)), //
                    V1.row(vi),                                           //
                    V1.row(F(fi, 0)), V1.row(F(fi, 1)), V1.row(F(fi, 2)), //
                    toi,
                    /*tmax=*/1.0,
                    /*tolerance=*/1e-6,
                    /*max_iterations=*/1e7,
                    /*conservative_rescaling=*/1.0);
                CHECK(!hit); // Check for FN

            } else {
                sh_ci++;
            }
        }
        CHECK(sh_ci >= sh_candidates.fv_candidates.size());

        U.setRandom();
        U *= 3;
    }
}
