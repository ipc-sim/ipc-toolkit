#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>

#include <ipc.hpp>
#include <spatial_hash/hash_grid.hpp>

using namespace ccd;

TEST_CASE("AABB initilization", "[spatial_hash][AABB]")
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

TEST_CASE("AABB overlapping", "[spatial_hash][AABB]")
{
    AABB a, b;
    bool are_overlaping = false;
    SECTION("a to the right of b")
    {
        a = AABB(Eigen::Vector2d(-1, 0), Eigen::Vector2d(0, 1));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 1));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Vector2d(0.5, 0), Eigen::Vector2d(1.5, 1));
            are_overlaping = false;
        }
    }
    SECTION("b to the right of a")
    {
        b = AABB(Eigen::Vector2d(-1, 0), Eigen::Vector2d(0, 1));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Vector2d(-0.5, 0), Eigen::Vector2d(0.5, 1));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Vector2d(0.5, 0), Eigen::Vector2d(1.5, 1));
            are_overlaping = false;
        }
    }
    SECTION("a above b")
    {
        a = AABB(Eigen::Vector2d(0, -1), Eigen::Vector2d(1, 0));
        SECTION("overlapping")
        {
            b = AABB(Eigen::Vector2d(0, -0.5), Eigen::Vector2d(1, 0.5));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            b = AABB(Eigen::Vector2d(0, 0.5), Eigen::Vector2d(1, 1.5));
            are_overlaping = false;
        }
    }
    SECTION("a above b")
    {
        b = AABB(Eigen::Vector2d(0, -1), Eigen::Vector2d(1, 0));
        SECTION("overlapping")
        {
            a = AABB(Eigen::Vector2d(0, -0.5), Eigen::Vector2d(1, 0.5));
            are_overlaping = true;
        }
        SECTION("not overlapping")
        {
            a = AABB(Eigen::Vector2d(0, 0.5), Eigen::Vector2d(1, 1.5));
            are_overlaping = false;
        }
    }
    CHECK(AABB::are_overlaping(a, b) == are_overlaping);
}

TEST_CASE("Vertex-Vertex Spatial Hash", "[ccd][spatial_hash]")
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

TEST_CASE("Entire 2D Mesh", "[ccd][spatial_hash]")
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

// TEST_CASE("Benchmark different spatiahashes", "[!benchmark][spatial_hash]")
// {
//     Eigen::MatrixXd vertices;
//     Eigen::MatrixXd displacements;
//     Eigen::MatrixXi edges;
//     Eigen::MatrixXi faces;
//     Eigen::VectorXi group_ids;
//
//     SECTION("Simple")
//     {
//         vertices.resize(4, 3);
//         vertices.row(0) << -1, -1, 0;
//         vertices.row(1) << 1, -1, 0;
//         vertices.row(2) << 0, 1, 1;
//         vertices.row(3) << 0, 1, -1;
//
//         edges.resize(2, 2);
//         edges.row(0) << 0, 1;
//         edges.row(1) << 2, 3;
//
//         SECTION("Without group ids") {}
//         SECTION("With group ids")
//         {
//             group_ids.resize(4);
//             group_ids << 0, 0, 1, 1;
//         }
//
//         faces.resize(0, 3);
//
//         displacements = Eigen::MatrixXd::Zero(vertices.rows(),
//         vertices.cols()); displacements.col(1).head(2).setConstant(2);
//         displacements.col(1).tail(2).setConstant(-2);
//     }
//     SECTION("Complex")
//     {
//         std::string fname = GENERATE(std::string("cube.obj")
//                                      /*, std::string("bunny-lowpoly.obj")*/);
//
//         boost::filesystem::path mesh_path = boost::filesystem::path(__FILE__)
//                                                 .parent_path()
//                                                 .parent_path()
//                                                 .parent_path()
//             / "meshes" / fname;
//         igl::read_triangle_mesh(mesh_path.string(), vertices, faces);
//         igl::edges(faces, edges);
//
//         displacements = Eigen::MatrixXd::Zero(vertices.rows(),
//         vertices.cols()); displacements.col(1).setOnes();
//     }
//
//     for (int i = 0; i < 10; i++) {
//         ConcurrentImpacts brute_force_impacts;
//         detect_collisions(
//             vertices, vertices + displacements, edges, faces, group_ids,
//             CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX,
//             brute_force_impacts, DetectionMethod::BRUTE_FORCE);
//         REQUIRE(brute_force_impacts.ev_impacts.size() == 0);
//
//         ConcurrentImpacts hash_impacts;
//         detect_collisions(
//             vertices, vertices + displacements, edges, faces, group_ids,
//             CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX,
//             hash_impacts, DetectionMethod::HASH_GRID);
//         REQUIRE(hash_impacts.ev_impacts.size() == 0);
//
//         REQUIRE(
//             brute_force_impacts.ee_impacts.size()
//             == hash_impacts.ee_impacts.size());
//         std::sort(
//             brute_force_impacts.ee_impacts.begin(),
//             brute_force_impacts.ee_impacts.end(),
//             compare_impacts_by_time<EdgeEdgeImpact>);
//         std::sort(
//             hash_impacts.ee_impacts.begin(), hash_impacts.ee_impacts.end(),
//             compare_impacts_by_time<EdgeEdgeImpact>);
//         bool is_equa=
//             brute_force_impacts.ee_impacts == hash_impacts.ee_impacts;
//         if (!is_equa&& brute_force_impacts.size() > 0) {
//             spdlog::error(
//                 "bf_impacts.ee[0]={{time={:g}, e0i={:d}, α₀={:g}, e0i={:d}, "
//                 "α₁={:g}}}",
//                 brute_force_impacts.ee_impacts[0].time,
//                 brute_force_impacts.ee_impacts[0].impacted_edge_index,
//                 brute_force_impacts.ee_impacts[0].impacted_alpha,
//                 brute_force_impacts.ee_impacts[0].impacting_edge_index,
//                 brute_force_impacts.ee_impacts[0].impacting_alpha);
//             spdlog::error(
//                 "hash_impacts[0]={{time={:g}, e0i={:d}, α₀={:g}, e0i={:d}, "
//                 "α₁={:g}}}",
//                 hash_impacts.ee_impacts[0].time,
//                 hash_impacts.ee_impacts[0].impacted_edge_index,
//                 hash_impacts.ee_impacts[0].impacted_alpha,
//                 hash_impacts.ee_impacts[0].impacting_edge_index,
//                 hash_impacts.ee_impacts[0].impacting_alpha);
//         }
//         CHECK(is_equal);
//
//         REQUIRE(
//             brute_force_impacts.fv_impacts.size()
//             == hash_impacts.fv_impacts.size());
//         std::sort(
//             brute_force_impacts.fv_impacts.begin(),
//             brute_force_impacts.fv_impacts.end(),
//             compare_impacts_by_time<FaceVertexImpact>);
//         std::sort(
//             hash_impacts.fv_impacts.begin(), hash_impacts.fv_impacts.end(),
//             compare_impacts_by_time<FaceVertexImpact>);
//         CHECK(brute_force_impacts.fv_impacts == hash_impacts.fv_impacts);
//
//         displacements.setRandom();
//         displacements *= 3;
//     }
//
//     BENCHMARK("IPC") {}
//     BENCHMARK("new") {}
// }
