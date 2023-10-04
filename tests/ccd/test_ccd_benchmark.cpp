#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_TEST_CCD_BENCHMARK
#include <catch2/catch_all.hpp>

#include <ipc/ccd/ccd.hpp>
#include <ipc/ccd/additive_ccd.hpp>
#include <fmt/format.h>

#include <igl/Timer.h>
#include <ccd_io/read_ccd_queries.hpp>

#include <iostream>

using namespace ipc;

void run_benchmark(
    const std::function<bool(
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        double&)>& edge_edge_ccd,
    const std::function<bool(
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        const Eigen::Vector3d&,
        double&)>& point_triangle_ccd)
{
    logger().set_level(spdlog::level::err);

    namespace fs = std::filesystem;
    using Matrix8x3 = Eigen::Matrix<double, 8, 3, Eigen::RowMajor>;

    std::vector<fs::path> csv_dirs;
    bool is_data_split = true;
    // SECTION("Original Benchmark")
    // {
    //     is_data_split = false;

    //     const fs::path root_path(
    //         // CCD_IO_SAMPLE_QUERIES_DIR
    //         "/Users/zachary/Development/research/collision-detection/dataset");

    //     const std::vector<std::string> folders { {
    //         "chain",
    //         "cow-heads",
    //         "golf-ball",
    //         "mat-twist",
    //         "erleben-sliding-spike",
    //         "erleben-spike-wedge",
    //         "erleben-sliding-wedge",
    //         "erleben-wedge-crack",
    //         "erleben-spike-crack",
    //         "erleben-wedges",
    //         "erleben-cube-cliff-edges",
    //         "erleben-spike-hole",
    //         "erleben-cube-internal-edges",
    //         "erleben-spikes",
    //         "unit-tests",
    //     } };

    //     const std::array<std::string, 2> subfolders = { {
    //         "edge-edge",
    //         "vertex-face",
    //     } };

    //     for (const std::string& folder : folders) {
    //         for (const std::string& subfolder : subfolders) {
    //             csv_dirs.push_back(root_path / folder / subfolder);
    //         }
    //     }
    // }
    SECTION("New Benchmark")
    {
        bool is_data_split = true;

        const fs::path root_path(
            "/Users/zachary/Documents/Papers/Time of Impact Dataset for Continuous Collision Detection and a Scalable Conservative Algorithm/data/scalable-ccd-data");

        const std::vector<std::string> folders { {
            "armadillo-rollers",
            "cloth-ball",
            "cloth-funnel",
            "n-body-simulation",
            "rod-twist",
        } };

        for (const std::string& folder : folders) {
            csv_dirs.push_back(root_path / folder / "queries");
        }
    }

    int i = 0;
    double total_time = 0;
    int num_fp = 0, num_fn = 0;
    igl::Timer timer;

    for (const auto& csv_dir : csv_dirs) {
        for (const auto& csv : fs::directory_iterator(csv_dir)) {
            std::vector<ccd_io::CCDQuery> queries;
            if (is_data_split) {
                queries = ccd_io::read_ccd_queries(
                    csv.path().string(),
                    csv.path().parent_path().parent_path() / "mma_bool"
                        / (csv.path().stem().string() + "_mma_bool.json"));
            } else {
                queries = ccd_io::read_ccd_queries(csv.path().string());
            }

            bool is_edge_edge = is_data_split
                ? csv.path().stem().string().find("ee") != std::string::npos
                : csv.path().parent_path().filename().string() == "edge-edge";

            for (int qi = 0; qi < queries.size(); qi++) {
                INFO(fmt::format("{}(Q{})", csv.path().string(), qi));
                std::cout << "\r                      "
                          << "\rTesting query: " << ++i << std::flush;

                Eigen::Map<const Matrix8x3> V(&queries[qi].vertices[0][0]);
                const bool expected_result = queries[qi].ground_truth;

                bool result;
                double toi;
                timer.start();
                if (is_edge_edge) {
                    result = edge_edge_ccd(
                        V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
                        V.row(5), V.row(6), V.row(7), toi);
                } else {
                    result = point_triangle_ccd(
                        V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
                        V.row(5), V.row(6), V.row(7), toi);
                }
                timer.stop();
                total_time += timer.getElapsedTimeInMicroSec();

                if (result != expected_result) {
                    if (result) {
                        num_fp++;
                    } else {
                        num_fn++;
                    }
                }

                if (result < expected_result) {
                    fmt::print("\n");
                }
                CHECK(result >= expected_result); // false positive is ok
            }
        }
    }
    fmt::print("\n");
    fmt::print("Total time: {}s\n", total_time / 1e6);
    fmt::print("Average time: {}μs\n", total_time / i);
    fmt::print("False positives: {}\n", num_fp);
    fmt::print("False negatives: {}\n", num_fn);
}

#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
TEST_CASE(
    "Run CCD Benchmark on TI CCD",
    "[ccd][benchmark][3D][point-triangle][edge-edge]")
{
    fmt::print("Tight Inclusion CCD:\n");
#else
TEST_CASE(
    "Run CCD Benchmark on FP CCD",
    "[ccd][benchmark][3D][point-triangle][edge-edge][!mayfail]")
{
    fmt::print("Floating-Point CCD:\n");
#endif

    run_benchmark(
        [](const Eigen::Vector3d& ea0_t0, const Eigen::Vector3d& ea1_t0,
           const Eigen::Vector3d& eb0_t0, const Eigen::Vector3d& eb1_t0,
           const Eigen::Vector3d& ea0_t1, const Eigen::Vector3d& ea1_t1,
           const Eigen::Vector3d& eb0_t1, const Eigen::Vector3d& eb1_t1,
           double& toi) -> bool {
            return edge_edge_ccd(
                ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
                toi);
        },
        [](const Eigen::Vector3d& p_t0, const Eigen::Vector3d& t0_t0,
           const Eigen::Vector3d& t1_t0, const Eigen::Vector3d& t2_t0,
           const Eigen::Vector3d& p_t1, const Eigen::Vector3d& t0_t1,
           const Eigen::Vector3d& t1_t1, const Eigen::Vector3d& t2_t1,
           double& toi) -> bool {
            return point_triangle_ccd(
                p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);
        });
}

TEST_CASE(
    "Run CCD Benchmark on ACCD",
    "[ccd][benchmark][3D][point-triangle][edge-edge][!mayfail]")
{
    fmt::print("Additive CCD:\n");

    run_benchmark(
        [](const Eigen::Vector3d& ea0_t0, const Eigen::Vector3d& ea1_t0,
           const Eigen::Vector3d& eb0_t0, const Eigen::Vector3d& eb1_t0,
           const Eigen::Vector3d& ea0_t1, const Eigen::Vector3d& ea1_t1,
           const Eigen::Vector3d& eb0_t1, const Eigen::Vector3d& eb1_t1,
           double& toi) -> bool {
            return additive_ccd::edge_edge_ccd(
                ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
                toi);
        },

        [](const Eigen::Vector3d& p_t0, const Eigen::Vector3d& t0_t0,
           const Eigen::Vector3d& t1_t0, const Eigen::Vector3d& t2_t0,
           const Eigen::Vector3d& p_t1, const Eigen::Vector3d& t0_t1,
           const Eigen::Vector3d& t1_t1, const Eigen::Vector3d& t2_t1,
           double& toi) -> bool {
            return additive_ccd::point_triangle_ccd(
                p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);
        });
}

// TEST_CASE("Failing Benchmark Cases", "[ccd]")
// {
//     using Matrix8x3 = Eigen::Matrix<double, 8, 3, Eigen::RowMajor>;

//     const static std::vector<std::string> paths;
//     const static std::vector<std::vector<int>> qids;

//     for (int i = 0; i < paths.size(); i++) {
//         const std::vector<ccd_io::CCDQuery> queries =
//         ccd_io::read_ccd_queries(
//             std::string(CCD_IO_SAMPLE_QUERIES_DIR) + paths[i]);
//         for (const auto& qi : qids[i]) {
//             Eigen::Map<const Matrix8x3> V(&queries[qi].vertices[0][0]);
//             const bool expected_result = queries[qi].ground_truth;

//             bool result;
//             double toi;
//             // if (subfolder == "edge-edge") {
//             //     result = edge_edge_ccd(
//             //         V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
//             //         V.row(5), V.row(6), V.row(7), toi);
//             // } else {
//             result = additive_ccd::point_triangle_ccd(
//                 V.row(0), V.row(1), V.row(2), V.row(3), V.row(4), V.row(5),
//                 V.row(6), V.row(7), toi);
//             // }
//             CHECK(result >= expected_result); // false positive is ok
//         }
//     }
// }
#endif