#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <tbb/parallel_sort.h>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/ccd/ccd.hpp>

using namespace ipc;

// TEST_CASE("AABB initilization", "[hash_grid][AABB]")
// {
//     int dim = GENERATE(2, 3);
//     CAPTURE(dim);
//     AABB aabb;
//     ArrayMax3d actual_center(dim);
//     SECTION("Empty AABB")
//     {
//         aabb = AABB(ArrayMax3d::Zero(dim), ArrayMax3d::Zero(dim));
//         actual_center.setZero();
//     }
//     SECTION("Box centered at zero")
//     {
//         ArrayMax3d min =
//             ArrayMax3d::Random(dim).array() - 1; // in range [-2, 0]
//         ArrayMax3d max = -min;
//         aabb = AABB(min, max);
//         actual_center.setZero();
//     }
//     SECTION("Box not centered at zero")
//     {
//         ArrayMax3d min(dim), max(dim);
//         if (dim == 2) {
//             min << 5.1, 3.14;
//             max << 10.4, 7.89;
//             actual_center << 7.75, 5.515;
//         } else {
//             min << 5.1, 3.14, 7.94;
//             max << 10.4, 7.89, 10.89;
//             actual_center << 7.75, 5.515, 9.415;
//         }
//         aabb = AABB(min, max);
//     }
//     ArrayMax3d center_diff = aabb.getCenter() - actual_center;
//     CHECK(center_diff.matrix().norm() == Approx(0.0).margin(1e-12));
// }

TEST_CASE("AABB overlapping", "[has_grid][AABB]")
{
    AABB a, b;
    bool are_overlapping = false;
    SECTION("a to the right of b")
    {
        a = AABB(Eigen::Array2d(-1, 0), Eigen::Array2d(0, 1));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Array2d(-0.5, 0), Eigen::Array2d(0.5, 1));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Array2d(0.5, 0), Eigen::Array2d(1.5, 1));
            are_overlapping = false;
        }
    }
    SECTION("b to the right of a")
    {
        b = AABB(Eigen::Array2d(-1, 0), Eigen::Array2d(0, 1));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Array2d(-0.5, 0), Eigen::Array2d(0.5, 1));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Array2d(0.5, 0), Eigen::Array2d(1.5, 1));
            are_overlapping = false;
        }
    }
    SECTION("a above b")
    {
        a = AABB(Eigen::Array2d(0, -1), Eigen::Array2d(1, 0));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Array2d(0, -0.5), Eigen::Array2d(1, 0.5));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Array2d(0, 0.5), Eigen::Array2d(1, 1.5));
            are_overlapping = false;
        }
    }
    SECTION("a above b")
    {
        b = AABB(Eigen::Array2d(0, -1), Eigen::Array2d(1, 0));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Array2d(0, -0.5), Eigen::Array2d(1, 0.5));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Array2d(0, 0.5), Eigen::Array2d(1, 1.5));
            are_overlapping = false;
        }
    }
    CHECK(AABB::are_overlapping(a, b) == are_overlapping);
}

TEST_CASE("Compare HashGrid against brute force", "[hash_grid]")
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
        SECTION("With group ids")
        {
            group_ids.resize(4);
            group_ids << 0, 0, 1, 1;
        }

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

    HashGrid hashgrid;
    Candidates hg_candidates, bf_candidates;

    double inflation_radius = 1e-2; // GENERATE(0.0, 1e-4, 1e-3, 1e-2, 1e-1);

    auto can_collide = [&group_ids](size_t vi, size_t vj) {
        return group_ids.size() == 0 || group_ids(vi) != group_ids(vj);
    };

    for (int i = 0; i < 2; i++) {
        Eigen::MatrixXd V1 = V0 + U;
        hashgrid.resize(V0, V1, E, inflation_radius);
        hashgrid.addVertices(V0, V1, inflation_radius);
        hashgrid.addEdges(V0, V1, E, inflation_radius);
        hashgrid.addFaces(V0, V1, F, inflation_radius);

        hg_candidates.clear();
        hashgrid.getVertexEdgePairs(
            E, hg_candidates.ev_candidates, can_collide);
        hashgrid.getEdgeEdgePairs(E, hg_candidates.ee_candidates, can_collide);
        hashgrid.getFaceVertexPairs(
            F, hg_candidates.fv_candidates, can_collide);

        bf_candidates.clear();
        detect_collision_candidates_brute_force(
            V0, V1, E, F, bf_candidates,
            /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true,
            /*perform_aabb_check=*/false, inflation_radius,
            /*ignore_codimensional_vertices=*/false, can_collide);

        CHECK(
            hg_candidates.ev_candidates.size()
            <= bf_candidates.ev_candidates.size());
        CHECK(
            hg_candidates.ee_candidates.size()
            <= bf_candidates.ee_candidates.size());
        CHECK(
            hg_candidates.fv_candidates.size()
            <= bf_candidates.fv_candidates.size());

        tbb::parallel_sort(
            hg_candidates.ev_candidates.begin(),
            hg_candidates.ev_candidates.end());
        tbb::parallel_sort(
            bf_candidates.ev_candidates.begin(),
            bf_candidates.ev_candidates.end());
        int hg_ci = 0;
        for (int bf_ci = 0; bf_ci < bf_candidates.ev_candidates.size();
             bf_ci++) {
            if (hg_candidates.ev_candidates.size() <= hg_ci
                || bf_candidates.ev_candidates[bf_ci]
                    != hg_candidates.ev_candidates[hg_ci]) {
                long ei = bf_candidates.ev_candidates[bf_ci].edge_index;
                long vi = bf_candidates.ev_candidates[bf_ci].vertex_index;
                double toi;
                bool hit = point_edge_ccd(
                    V0.row(vi), V0.row(E(ei, 0)), V0.row(E(ei, 1)), // t = 0
                    V1.row(vi), V1.row(E(ei, 0)), V1.row(E(ei, 1)), // t = 1
                    toi, 1.0);
                CHECK(!hit); // Check for FN

            } else {
                hg_ci++;
            }
        }
        CHECK(hg_ci >= hg_candidates.ev_candidates.size());

        tbb::parallel_sort(
            hg_candidates.ee_candidates.begin(),
            hg_candidates.ee_candidates.end());
        tbb::parallel_sort(
            bf_candidates.ee_candidates.begin(),
            bf_candidates.ee_candidates.end());
        hg_ci = 0;
        for (int bf_ci = 0; bf_ci < bf_candidates.ee_candidates.size();
             bf_ci++) {
            if (hg_candidates.ee_candidates.size() <= hg_ci
                || bf_candidates.ee_candidates[bf_ci]
                    != hg_candidates.ee_candidates[hg_ci]) {
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
                hg_ci++;
            }
        }
        CHECK(hg_ci >= hg_candidates.ee_candidates.size());

        tbb::parallel_sort(
            hg_candidates.fv_candidates.begin(),
            hg_candidates.fv_candidates.end());
        tbb::parallel_sort(
            bf_candidates.fv_candidates.begin(),
            bf_candidates.fv_candidates.end());
        hg_ci = 0;
        for (int bf_ci = 0; bf_ci < bf_candidates.fv_candidates.size();
             bf_ci++) {
            if (hg_candidates.fv_candidates.size() <= hg_ci
                || bf_candidates.fv_candidates[bf_ci]
                    != hg_candidates.fv_candidates[hg_ci]) {
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
                hg_ci++;
            }
        }
        CHECK(hg_ci >= hg_candidates.fv_candidates.size());

        U.setRandom();
        U *= 3;
    }
}
