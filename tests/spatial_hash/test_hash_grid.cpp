#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/spatial_hash/hash_grid.hpp>
#include <ipc/spatial_hash/brute_force.hpp>
#include <ipc/ccd/ccd.hpp>

using namespace ipc;

TEST_CASE("AABB initilization", "[hash_grid][AABB]")
{
    int dim = GENERATE(2, 3);
    CAPTURE(dim);
    AABB aabb;
    Eigen::VectorXd actual_center(dim);
    SECTION("Empty AABB")
    {
        aabb = AABB(Eigen::VectorXd::Zero(dim), Eigen::VectorXd::Zero(dim));
        actual_center = Eigen::VectorXd::Zero(dim);
    }
    SECTION("Box centered at zero")
    {
        Eigen::VectorXd min =
            Eigen::VectorXd::Random(dim).array() - 1; // in range [-2, 0]
        Eigen::VectorXd max = -min;
        aabb = AABB(min, max);
        actual_center = Eigen::VectorXd::Zero(dim);
    }
    SECTION("Box not centered at zero")
    {
        Eigen::VectorXd min(dim), max(dim);
        if (dim == 2) {
            min << 5.1, 3.14;
            max << 10.4, 7.89;
            actual_center << 7.75, 5.515;
        } else {
            min << 5.1, 3.14, 7.94;
            max << 10.4, 7.89, 10.89;
            actual_center << 7.75, 5.515, 9.415;
        }
        aabb = AABB(min, max);
    }
    Eigen::VectorXd center_diff = aabb.getCenter() - actual_center;
    CHECK(center_diff.norm() == Approx(0.0).margin(1e-12));
}

TEST_CASE("AABB overlapping", "[has_grid][AABB]")
{
    AABB a, b;
    bool are_overlapping = false;
    SECTION("a to the right of b")
    {
        a = AABB(Eigen::Vector2d(-1, 0), Eigen::Vector2d(0, 1));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 1));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Vector2d(0.5, 0), Eigen::Vector2d(1.5, 1));
            are_overlapping = false;
        }
    }
    SECTION("b to the right of a")
    {
        b = AABB(Eigen::Vector2d(-1, 0), Eigen::Vector2d(0, 1));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 1));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Vector2d(0.5, 0), Eigen::Vector2d(1.5, 1));
            are_overlapping = false;
        }
    }
    SECTION("a above b")
    {
        a = AABB(Eigen::Vector2d(0, -1), Eigen::Vector2d(1, 0));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Vector2d(0, -0.5), Eigen::Vector2d(1, 0.5));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Vector2d(0, 0.5), Eigen::Vector2d(1, 1.5));
            are_overlapping = false;
        }
    }
    SECTION("a above b")
    {
        b = AABB(Eigen::Vector2d(0, -1), Eigen::Vector2d(1, 0));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Vector2d(0, -0.5), Eigen::Vector2d(1, 0.5));
            are_overlapping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Vector2d(0, 0.5), Eigen::Vector2d(1, 1.5));
            are_overlapping = false;
        }
    }
    CHECK(AABB::are_overlapping(a, b) == are_overlapping);
}

TEST_CASE("Vertex-Vertex Spatial Hash", "[ccd][has_grid]")
{
    Eigen::MatrixXd V_t0(4, 2);
    V_t0.row(0) << 1.11111, 0.5;  // edge 0 vertex 0
    V_t0.row(1) << 1.11111, 0.75; // edge 0 vertex 1
    V_t0.row(2) << 1, 0.5;        // edge 1 vertex 0
    V_t0.row(3) << 1, 0.75;       // edge 1 vertex 1

    Eigen::MatrixXd V_t1 = V_t0;
    V_t1.row(0) << 0.888889, 0.5;  // edge 0 vertex 0
    V_t1.row(1) << 0.888889, 0.75; // edge 0 vertex 1

    Eigen::MatrixXi E(2, 2);
    E.row(0) << 1, 0;
    E.row(1) << 2, 3;

    bool ignore_internal_vertices = GENERATE(false, true);

    bool is_valid_step = ipc::is_step_collision_free(
        V_t0, V_t1, E, /*F=*/Eigen::MatrixXi(), ignore_internal_vertices);

    CAPTURE(ignore_internal_vertices);
    CHECK(!is_valid_step);
}

TEST_CASE("Entire 2D Mesh", "[ccd][has_grid]")
{
    Eigen::MatrixXd V_t0;
    igl::readCSV(std::string(TEST_DATA_DIR) + "V_t0.txt", V_t0);

    Eigen::MatrixXd V_t1;
    igl::readCSV(std::string(TEST_DATA_DIR) + "V_t1.txt", V_t1);

    Eigen::MatrixXd E_double;
    igl::readCSV(std::string(TEST_DATA_DIR) + "E.txt", E_double);
    Eigen::MatrixXi E = E_double.cast<int>();

    bool ignore_internal_vertices = GENERATE(false, true);

    bool is_valid_step;
    SECTION("2D")
    {
        is_valid_step = ipc::is_step_collision_free(
            V_t0.leftCols(2), V_t1.rightCols(2), E, /*F=*/Eigen::MatrixXi(),
            ignore_internal_vertices);
    }
    SECTION("3D")
    {
        is_valid_step = ipc::is_step_collision_free(
            V_t0, V_t1, E, /*F=*/Eigen::MatrixXi(), ignore_internal_vertices);
    }

    CAPTURE(ignore_internal_vertices);
    CHECK(!is_valid_step);
}

TEST_CASE(
    "Test construct_constraint_set() with codimensional points",
    "[construct_constraint_set][has_grid]")
{
    double dhat = 0.1;
    // double dhat = 0.00173123;
    Eigen::MatrixXd V_rest, V;
    igl::readDMAT(
        std::string(TEST_DATA_DIR) + "codim-points/V_rest.dmat", V_rest);
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/V.dmat", V);
    Eigen::MatrixXi E, F;
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/E.dmat", E);
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/F.dmat", F);

    Constraints constraint_set;
    construct_constraint_set(
        V_rest, V, E, F, dhat, constraint_set,
        /*ignore_internal_vertices=*/false);

    CHECK(constraint_set.size() != 0);
}

TEST_CASE("Compare HashGrid against brute force", "[thisone][hash_grid]")
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
#ifdef NDEBUG
        std::string filename =
            GENERATE(std::string("cube.obj"), std::string("bunny.obj"));
#else
        std::string filename = "cube.obj";
#endif
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

    for (int i = 0; i < 2; i++) {
        Eigen::MatrixXd V1 = V0 + U;
        hashgrid.resize(V0, V1, E, inflation_radius);
        hashgrid.addVertices(V0, V1, inflation_radius);
        hashgrid.addEdges(V0, V1, E, inflation_radius);
        hashgrid.addFaces(V0, V1, F, inflation_radius);

        hg_candidates.clear();
        hashgrid.getVertexEdgePairs(E, group_ids, hg_candidates.ev_candidates);
        hashgrid.getEdgeEdgePairs(E, group_ids, hg_candidates.ee_candidates);
        hashgrid.getFaceVertexPairs(F, group_ids, hg_candidates.fv_candidates);

        bf_candidates.clear();
        detect_collision_candidates_brute_force(
            V0, V1, E, F, bf_candidates,
            /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true,
            /*perform_aabb_check=*/true, inflation_radius, group_ids);

        CHECK(
            hg_candidates.ev_candidates.size()
            <= bf_candidates.ev_candidates.size());
        CHECK(
            hg_candidates.ee_candidates.size()
            <= bf_candidates.ee_candidates.size());
        CHECK(
            hg_candidates.fv_candidates.size()
            <= bf_candidates.fv_candidates.size());

        std::sort(
            hg_candidates.ev_candidates.begin(),
            hg_candidates.ev_candidates.end());
        std::sort(
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

        std::sort(
            hg_candidates.ee_candidates.begin(),
            hg_candidates.ee_candidates.end());
        std::sort(
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

        std::sort(
            hg_candidates.fv_candidates.begin(),
            hg_candidates.fv_candidates.end());
        std::sort(
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
