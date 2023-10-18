#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/broad_phase/voxel_size_heuristic.hpp>

TEST_CASE("Voxel size heuristic", "[broad_phase][voxel_size]")
{
    Eigen::MatrixXd vertices_t0(4, 3);
    vertices_t0 << 0.0, 0.0, 0.0, //
        1.0, 0.0, 0.0,            //
        0.0, 1.0, 0.0,            //
        0.0, 0.0, 1.0;
    Eigen::MatrixXd vertices_t1 = vertices_t0;
    vertices_t1.row(3) += Eigen::RowVector3d(0.1, -0.1, -0.1);

    Eigen::MatrixXi edges(6, 2);
    edges << 0, 1, //
        0, 2,      //
        0, 3,      //
        1, 2,      //
        1, 3,      //
        2, 3;

    // Mean edge length
    {
        double std_deviation;
        const double mean = ipc::mean_edge_length(
            vertices_t0, vertices_t1, edges, std_deviation);

        CHECK(mean == Catch::Approx(1.189116069));
        CHECK(std_deviation == Catch::Approx(0.2085736678));
    }

    // Mean displacement length
    {
        double std_deviation;
        const double mean = ipc::mean_displacement_length(
            vertices_t1 - vertices_t0, std_deviation);

        CHECK(mean == Catch::Approx(0.04330127019));
        CHECK(std_deviation == Catch::Approx(0.075));
    }

    // Median edge length
    {
        const double median =
            ipc::median_edge_length(vertices_t0, vertices_t1, edges);

        CHECK(median == Catch::Approx(1.138357267));
    }

    // Median displacement length
    {
        const double median =
            ipc::median_displacement_length(vertices_t1 - vertices_t0);

        CHECK(median == Catch::Approx(0.0));
    }

    // Max edge length
    {
        const double max =
            ipc::max_edge_length(vertices_t0, vertices_t1, edges);

        CHECK(max == Catch::Approx(1.424780685));
    }

    // Max displacement length
    {
        const double max =
            ipc::max_displacement_length(vertices_t1 - vertices_t0);

        CHECK(max == Catch::Approx(0.1732050808));
    }
}