#include <catch2/catch_test_macros.hpp>

#include <ipc/broad_phase/cuda/lbvh.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/utils/cuda/is_available.hpp>

TEST_CASE(
    "CUDA availability predicts whether the GPU broad phase runs",
    "[utils][cuda]")
{
#ifdef IPC_TOOLKIT_WITH_CUDA
    // Two parallel triangles 0.05 apart, so an inflation radius of 0.1 makes
    // every vertex of one overlap the face of the other.
    Eigen::MatrixXd vertices(6, 3);
    vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0, //
        0, 0, 0.05, 1, 0, 0.05, 0, 1, 0.05;
    Eigen::MatrixXi faces(2, 3);
    faces << 0, 1, 2, 3, 4, 5;
    Eigen::MatrixXi edges(6, 2);
    edges << 0, 1, 1, 2, 2, 0, 3, 4, 4, 5, 5, 3;

    ipc::cuda::LBVH lbvh;
    ipc::Candidates candidates;
    if (ipc::cuda::is_available()) {
        lbvh.build(vertices, edges, faces, /*inflation_radius=*/0.1);
        lbvh.detect_collision_candidates(candidates);
        CHECK(!candidates.fv_candidates.empty());
    } else {
        // Without a device, the first CUDA call fails and we report it as an
        // exception rather than crashing.
        CHECK_THROWS(lbvh.build(vertices, edges, faces, 0.1));
    }
#else
    CHECK_FALSE(ipc::cuda::is_available());
#endif
}
