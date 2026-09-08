// Exercises the __host__ __device__ distance functions from CUDA kernels.
//
// Compilation proves the framework works: nvcc must instantiate every
// converted distance function (including the distance-type classification and
// the signed variants) as device code. On a machine with a CUDA device the
// kernels also run and their outputs are checked against the host
// implementations — host/device parity is the contract. Without a device
// (e.g. a GPU-less container) the runtime checks are skipped but the device
// code has still been compiled.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/distance/signed/line_line.hpp>
#include <ipc/distance/signed/point_line.hpp>
#include <ipc/distance/signed/point_plane.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cuda_runtime.h>

#include <vector>

namespace {

#define REQUIRE_CUDA(expr) REQUIRE((expr) == cudaSuccess)

/// @brief Uniform kernel signature: unpack inputs from `in`, write results
/// flattened (row-major) to `out`.
using GpuTestKernel = void (*)(const double* in, double* out);

/// @brief Copy a fixed-size Eigen object into `out` row-major; returns the
/// number of doubles written.
template <typename Derived>
__device__ int
write_out(const Eigen::MatrixBase<Derived>& m, double* out, int offset)
{
    for (int r = 0; r < m.rows(); ++r) {
        for (int c = 0; c < m.cols(); ++c) {
            out[offset++] = m(r, c);
        }
    }
    return offset;
}

/// @brief Append a fixed-size Eigen object to a host vector row-major.
template <typename Derived>
void write_expected(
    const Eigen::MatrixBase<Derived>& m, std::vector<double>& expected)
{
    for (int r = 0; r < m.rows(); ++r) {
        for (int c = 0; c < m.cols(); ++c) {
            expected.push_back(m(r, c));
        }
    }
}

/// @brief Launch `kernel` with the given inputs; returns the device outputs.
std::vector<double> run_gpu_kernel(
    GpuTestKernel kernel, const std::vector<double>& in, const size_t n_out)
{
    double *d_in = nullptr, *d_out = nullptr;
    REQUIRE_CUDA(cudaMalloc(&d_in, in.size() * sizeof(double)));
    REQUIRE_CUDA(cudaMalloc(&d_out, n_out * sizeof(double)));
    REQUIRE_CUDA(cudaMemcpy(
        d_in, in.data(), in.size() * sizeof(double), cudaMemcpyHostToDevice));

    kernel<<<1, 1>>>(d_in, d_out);
    REQUIRE_CUDA(cudaGetLastError());
    REQUIRE_CUDA(cudaDeviceSynchronize());

    std::vector<double> out(n_out);
    REQUIRE_CUDA(cudaMemcpy(
        out.data(), d_out, n_out * sizeof(double), cudaMemcpyDeviceToHost));
    cudaFree(d_in);
    cudaFree(d_out);
    return out;
}

/// @brief SKIP the test if no CUDA device is present (kernels have still been
/// compiled as device code, which is the primary verification).
void skip_if_no_cuda_device()
{
    int device_count = 0;
    const cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }
}

void check_gpu_matches_host(
    GpuTestKernel kernel,
    const std::vector<double>& in,
    const std::vector<double>& expected)
{
    const std::vector<double> out = run_gpu_kernel(kernel, in, expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        CAPTURE(i);
        CHECK(out[i] == Catch::Approx(expected[i]).margin(1e-12));
    }
}

// ---------------------------------------------------------------------------

// Shared fixtures (kept in one place so kernels and host references agree).
// clang-format off
const std::vector<double> PP_IN  = { 0, 0, 0,  1, 2, 2 };
const std::vector<double> PL_IN  = { 0.1, 1, 0.2,  -1, 0, 0,  1, 0, 0.1 };
const std::vector<double> LL_IN  = { -1, 0, 0,  1, 0, 0,  0, 1, -1,  0.1, 1, 1 };
const std::vector<double> PLANE_IN = { 0, 1, 0,  0, 0, 0,  0, 1, 0,
                                       -1, 0, 1,  1, 0, 1,  0, 0, -1 };
const std::vector<double> PE_IN  = { 0.5, 1, 0,  -1, 0, 0,  1, 0, 0 };
const std::vector<double> EE_IN  = { -1, 0, 0,  1, 0, 0,  -0.5, 1, -1,  0.5, 1, 1 };
const std::vector<double> PT_IN  = { 0, 1, 0,  -1, 0, 1,  1, 0, 1,  0, 0, -1 };
const std::vector<double> SGN_IN = { 0.1, 1,  -1, 0,  1, 0 }; // 2D point-line
// clang-format on

__global__ void point_point_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const ipc::VectorMax3d p0 = Eigen::Vector3d(in[0], in[1], in[2]);
    const ipc::VectorMax3d p1 = Eigen::Vector3d(in[3], in[4], in[5]);

    int offset = 0;
    out[offset++] = ipc::point_point_distance(p0, p1);
    offset = write_out(ipc::point_point_distance_gradient(p0, p1), out, offset);
    offset = write_out(ipc::point_point_distance_hessian(p0, p1), out, offset);
}

__global__ void point_line_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const ipc::VectorMax3d p = Eigen::Vector3d(in[0], in[1], in[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(in[3], in[4], in[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(in[6], in[7], in[8]);

    int offset = 0;
    out[offset++] = ipc::point_line_distance(p, e0, e1);
    offset =
        write_out(ipc::point_line_distance_gradient(p, e0, e1), out, offset);
    offset =
        write_out(ipc::point_line_distance_hessian(p, e0, e1), out, offset);
}

__global__ void line_line_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d ea0(in[0], in[1], in[2]);
    const Eigen::Vector3d ea1(in[3], in[4], in[5]);
    const Eigen::Vector3d eb0(in[6], in[7], in[8]);
    const Eigen::Vector3d eb1(in[9], in[10], in[11]);

    int offset = 0;
    out[offset++] = ipc::line_line_distance(ea0, ea1, eb0, eb1);
    offset = write_out(
        ipc::line_line_distance_gradient(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(
        ipc::line_line_distance_hessian(ea0, ea1, eb0, eb1), out, offset);
}

__global__ void point_plane_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d p(in[0], in[1], in[2]);
    const Eigen::Vector3d origin(in[3], in[4], in[5]);
    const Eigen::Vector3d normal(in[6], in[7], in[8]);
    const Eigen::Vector3d t0(in[9], in[10], in[11]);
    const Eigen::Vector3d t1(in[12], in[13], in[14]);
    const Eigen::Vector3d t2(in[15], in[16], in[17]);

    int offset = 0;
    out[offset++] = ipc::point_plane_distance(p, origin, normal);
    offset = write_out(
        ipc::point_plane_distance_gradient(p, origin, normal), out, offset);
    offset = write_out(
        ipc::point_plane_distance_hessian(p, origin, normal), out, offset);
    out[offset++] = ipc::point_plane_distance(p, t0, t1, t2);
    offset = write_out(
        ipc::point_plane_distance_gradient(p, t0, t1, t2), out, offset);
    offset = write_out(
        ipc::point_plane_distance_hessian(p, t0, t1, t2), out, offset);
}

__global__ void point_edge_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const ipc::VectorMax3d p = Eigen::Vector3d(in[0], in[1], in[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(in[3], in[4], in[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(in[6], in[7], in[8]);

    int offset = 0;
    out[offset++] =
        static_cast<double>(ipc::point_edge_distance_type(p, e0, e1));
    out[offset++] = ipc::point_edge_distance(p, e0, e1);
    offset =
        write_out(ipc::point_edge_distance_gradient(p, e0, e1), out, offset);
    offset =
        write_out(ipc::point_edge_distance_hessian(p, e0, e1), out, offset);
}

__global__ void edge_edge_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d ea0(in[0], in[1], in[2]);
    const Eigen::Vector3d ea1(in[3], in[4], in[5]);
    const Eigen::Vector3d eb0(in[6], in[7], in[8]);
    const Eigen::Vector3d eb1(in[9], in[10], in[11]);

    int offset = 0;
    out[offset++] =
        static_cast<double>(ipc::edge_edge_distance_type(ea0, ea1, eb0, eb1));
    out[offset++] = static_cast<double>(
        ipc::edge_edge_parallel_distance_type(ea0, ea1, eb0, eb1));
    out[offset++] = ipc::edge_edge_distance(ea0, ea1, eb0, eb1);
    offset = write_out(
        ipc::edge_edge_distance_gradient(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(
        ipc::edge_edge_distance_hessian(ea0, ea1, eb0, eb1), out, offset);
}

__global__ void point_triangle_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d p(in[0], in[1], in[2]);
    const Eigen::Vector3d t0(in[3], in[4], in[5]);
    const Eigen::Vector3d t1(in[6], in[7], in[8]);
    const Eigen::Vector3d t2(in[9], in[10], in[11]);

    int offset = 0;
    out[offset++] =
        static_cast<double>(ipc::point_triangle_distance_type(p, t0, t1, t2));
    out[offset++] = ipc::point_triangle_distance(p, t0, t1, t2);
    offset = write_out(
        ipc::point_triangle_distance_gradient(p, t0, t1, t2), out, offset);
    offset = write_out(
        ipc::point_triangle_distance_hessian(p, t0, t1, t2), out, offset);
}

__global__ void mollifier_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d ea0(in[0], in[1], in[2]);
    const Eigen::Vector3d ea1(in[3], in[4], in[5]);
    const Eigen::Vector3d eb0(in[6], in[7], in[8]);
    const Eigen::Vector3d eb1(in[9], in[10], in[11]);
    const double eps_x = in[12];

    int offset = 0;
    const double c = ipc::edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
    out[offset++] = c;
    offset = write_out(
        ipc::edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1), out,
        offset);
    offset = write_out(
        ipc::edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1), out,
        offset);

    const double x = 0.5 * eps_x; // scalar mollifier in the active range
    out[offset++] = ipc::edge_edge_mollifier(x, eps_x);
    out[offset++] = ipc::edge_edge_mollifier_gradient(x, eps_x);
    out[offset++] = ipc::edge_edge_mollifier_hessian(x, eps_x);
    out[offset++] = ipc::edge_edge_mollifier_derivative_wrt_eps_x(x, eps_x);
    out[offset++] =
        ipc::edge_edge_mollifier_gradient_derivative_wrt_eps_x(x, eps_x);

    out[offset++] = ipc::edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    offset = write_out(
        ipc::edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x), out,
        offset);
    offset = write_out(
        ipc::edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x), out,
        offset);
    offset = write_out(
        ipc::edge_edge_mollifier_gradient_wrt_x(
            ea0, ea1, eb0, eb1, ea0, ea1, eb0, eb1),
        out, offset);
    offset = write_out(
        ipc::edge_edge_mollifier_gradient_jacobian_wrt_x(
            ea0, ea1, eb0, eb1, ea0, ea1, eb0, eb1),
        out, offset);

    out[offset++] = ipc::edge_edge_mollifier_threshold(ea0, ea1, eb0, eb1);
    offset = write_out(
        ipc::edge_edge_mollifier_threshold_gradient(ea0, ea1, eb0, eb1), out,
        offset);
}

__global__ void signed_distance_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    // 2D point-line (inputs 0-5)
    const Eigen::Vector2d p2(in[0], in[1]);
    const Eigen::Vector2d f0(in[2], in[3]);
    const Eigen::Vector2d f1(in[4], in[5]);
    // 3D line-line + point-plane (inputs 6-17 and 18-29)
    const Eigen::Vector3d ea0(in[6], in[7], in[8]);
    const Eigen::Vector3d ea1(in[9], in[10], in[11]);
    const Eigen::Vector3d eb0(in[12], in[13], in[14]);
    const Eigen::Vector3d eb1(in[15], in[16], in[17]);
    const Eigen::Vector3d p(in[18], in[19], in[20]);
    const Eigen::Vector3d t0(in[21], in[22], in[23]);
    const Eigen::Vector3d t1(in[24], in[25], in[26]);
    const Eigen::Vector3d t2(in[27], in[28], in[29]);

    int offset = 0;
    out[offset++] = ipc::point_line_signed_distance(p2, f0, f1);
    offset = write_out(
        ipc::point_line_signed_distance_gradient(p2, f0, f1), out, offset);
    offset = write_out(
        ipc::point_line_signed_distance_hessian(p2, f0, f1), out, offset);

    out[offset++] = ipc::line_line_signed_distance(ea0, ea1, eb0, eb1);
    offset = write_out(
        ipc::line_line_signed_distance_gradient(ea0, ea1, eb0, eb1), out,
        offset);
    offset = write_out(
        ipc::line_line_signed_distance_hessian(ea0, ea1, eb0, eb1), out,
        offset);

    out[offset++] = ipc::point_plane_signed_distance(p, t0, t1, t2);
    offset = write_out(
        ipc::point_plane_signed_distance_gradient(p, t0, t1, t2), out, offset);
    offset = write_out(
        ipc::point_plane_signed_distance_hessian(p, t0, t1, t2), out, offset);
}

} // namespace

TEST_CASE("GPU point-point distance", "[distance][point_point][gpu]")
{
    skip_if_no_cuda_device();

    const ipc::VectorMax3d p0 = Eigen::Vector3d(0, 0, 0);
    const ipc::VectorMax3d p1 = Eigen::Vector3d(1, 2, 2);

    std::vector<double> expected;
    expected.push_back(ipc::point_point_distance(p0, p1));
    write_expected(ipc::point_point_distance_gradient(p0, p1), expected);
    write_expected(ipc::point_point_distance_hessian(p0, p1), expected);
    REQUIRE(expected[0] == Catch::Approx(9.0)); // sanity-check the fixture

    check_gpu_matches_host(point_point_kernel, PP_IN, expected);
}

TEST_CASE("GPU point-line distance", "[distance][point_line][gpu]")
{
    skip_if_no_cuda_device();

    const ipc::VectorMax3d p = Eigen::Vector3d(PL_IN[0], PL_IN[1], PL_IN[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(PL_IN[3], PL_IN[4], PL_IN[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(PL_IN[6], PL_IN[7], PL_IN[8]);

    std::vector<double> expected;
    expected.push_back(ipc::point_line_distance(p, e0, e1));
    write_expected(ipc::point_line_distance_gradient(p, e0, e1), expected);
    write_expected(ipc::point_line_distance_hessian(p, e0, e1), expected);

    check_gpu_matches_host(point_line_kernel, PL_IN, expected);
}

TEST_CASE("GPU line-line distance", "[distance][line_line][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d ea0(LL_IN[0], LL_IN[1], LL_IN[2]);
    const Eigen::Vector3d ea1(LL_IN[3], LL_IN[4], LL_IN[5]);
    const Eigen::Vector3d eb0(LL_IN[6], LL_IN[7], LL_IN[8]);
    const Eigen::Vector3d eb1(LL_IN[9], LL_IN[10], LL_IN[11]);

    std::vector<double> expected;
    expected.push_back(ipc::line_line_distance(ea0, ea1, eb0, eb1));
    write_expected(
        ipc::line_line_distance_gradient(ea0, ea1, eb0, eb1), expected);
    write_expected(
        ipc::line_line_distance_hessian(ea0, ea1, eb0, eb1), expected);

    check_gpu_matches_host(line_line_kernel, LL_IN, expected);
}

TEST_CASE("GPU point-plane distance", "[distance][point_plane][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d p(PLANE_IN[0], PLANE_IN[1], PLANE_IN[2]);
    const Eigen::Vector3d origin(PLANE_IN[3], PLANE_IN[4], PLANE_IN[5]);
    const Eigen::Vector3d normal(PLANE_IN[6], PLANE_IN[7], PLANE_IN[8]);
    const Eigen::Vector3d t0(PLANE_IN[9], PLANE_IN[10], PLANE_IN[11]);
    const Eigen::Vector3d t1(PLANE_IN[12], PLANE_IN[13], PLANE_IN[14]);
    const Eigen::Vector3d t2(PLANE_IN[15], PLANE_IN[16], PLANE_IN[17]);

    std::vector<double> expected;
    expected.push_back(ipc::point_plane_distance(p, origin, normal));
    write_expected(
        ipc::point_plane_distance_gradient(p, origin, normal), expected);
    write_expected(
        ipc::point_plane_distance_hessian(p, origin, normal), expected);
    expected.push_back(ipc::point_plane_distance(p, t0, t1, t2));
    write_expected(ipc::point_plane_distance_gradient(p, t0, t1, t2), expected);
    write_expected(ipc::point_plane_distance_hessian(p, t0, t1, t2), expected);

    check_gpu_matches_host(point_plane_kernel, PLANE_IN, expected);
}

TEST_CASE("GPU point-edge distance", "[distance][point_edge][gpu]")
{
    skip_if_no_cuda_device();

    const ipc::VectorMax3d p = Eigen::Vector3d(PE_IN[0], PE_IN[1], PE_IN[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(PE_IN[3], PE_IN[4], PE_IN[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(PE_IN[6], PE_IN[7], PE_IN[8]);

    std::vector<double> expected;
    expected.push_back(
        static_cast<double>(ipc::point_edge_distance_type(p, e0, e1)));
    expected.push_back(ipc::point_edge_distance(p, e0, e1));
    write_expected(ipc::point_edge_distance_gradient(p, e0, e1), expected);
    write_expected(ipc::point_edge_distance_hessian(p, e0, e1), expected);

    check_gpu_matches_host(point_edge_kernel, PE_IN, expected);
}

TEST_CASE("GPU edge-edge distance", "[distance][edge_edge][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d ea0(EE_IN[0], EE_IN[1], EE_IN[2]);
    const Eigen::Vector3d ea1(EE_IN[3], EE_IN[4], EE_IN[5]);
    const Eigen::Vector3d eb0(EE_IN[6], EE_IN[7], EE_IN[8]);
    const Eigen::Vector3d eb1(EE_IN[9], EE_IN[10], EE_IN[11]);

    std::vector<double> expected;
    expected.push_back(
        static_cast<double>(ipc::edge_edge_distance_type(ea0, ea1, eb0, eb1)));
    expected.push_back(
        static_cast<double>(
            ipc::edge_edge_parallel_distance_type(ea0, ea1, eb0, eb1)));
    expected.push_back(ipc::edge_edge_distance(ea0, ea1, eb0, eb1));
    write_expected(
        ipc::edge_edge_distance_gradient(ea0, ea1, eb0, eb1), expected);
    write_expected(
        ipc::edge_edge_distance_hessian(ea0, ea1, eb0, eb1), expected);

    check_gpu_matches_host(edge_edge_kernel, EE_IN, expected);
}

TEST_CASE("GPU point-triangle distance", "[distance][point_triangle][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d p(PT_IN[0], PT_IN[1], PT_IN[2]);
    const Eigen::Vector3d t0(PT_IN[3], PT_IN[4], PT_IN[5]);
    const Eigen::Vector3d t1(PT_IN[6], PT_IN[7], PT_IN[8]);
    const Eigen::Vector3d t2(PT_IN[9], PT_IN[10], PT_IN[11]);

    std::vector<double> expected;
    expected.push_back(
        static_cast<double>(ipc::point_triangle_distance_type(p, t0, t1, t2)));
    expected.push_back(ipc::point_triangle_distance(p, t0, t1, t2));
    write_expected(
        ipc::point_triangle_distance_gradient(p, t0, t1, t2), expected);
    write_expected(
        ipc::point_triangle_distance_hessian(p, t0, t1, t2), expected);

    check_gpu_matches_host(point_triangle_kernel, PT_IN, expected);
}

TEST_CASE("GPU edge-edge mollifier", "[distance][edge_edge_mollifier][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d ea0(EE_IN[0], EE_IN[1], EE_IN[2]);
    const Eigen::Vector3d ea1(EE_IN[3], EE_IN[4], EE_IN[5]);
    const Eigen::Vector3d eb0(EE_IN[6], EE_IN[7], EE_IN[8]);
    const Eigen::Vector3d eb1(EE_IN[9], EE_IN[10], EE_IN[11]);
    const double eps_x = 1e-3;

    std::vector<double> in = EE_IN;
    in.push_back(eps_x);

    std::vector<double> expected;
    expected.push_back(ipc::edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1));
    write_expected(
        ipc::edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1),
        expected);
    write_expected(
        ipc::edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1), expected);

    const double x = 0.5 * eps_x;
    expected.push_back(ipc::edge_edge_mollifier(x, eps_x));
    expected.push_back(ipc::edge_edge_mollifier_gradient(x, eps_x));
    expected.push_back(ipc::edge_edge_mollifier_hessian(x, eps_x));
    expected.push_back(ipc::edge_edge_mollifier_derivative_wrt_eps_x(x, eps_x));
    expected.push_back(
        ipc::edge_edge_mollifier_gradient_derivative_wrt_eps_x(x, eps_x));

    expected.push_back(ipc::edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x));
    write_expected(
        ipc::edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x), expected);
    write_expected(
        ipc::edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x), expected);
    write_expected(
        ipc::edge_edge_mollifier_gradient_wrt_x(
            ea0, ea1, eb0, eb1, ea0, ea1, eb0, eb1),
        expected);
    write_expected(
        ipc::edge_edge_mollifier_gradient_jacobian_wrt_x(
            ea0, ea1, eb0, eb1, ea0, ea1, eb0, eb1),
        expected);

    expected.push_back(ipc::edge_edge_mollifier_threshold(ea0, ea1, eb0, eb1));
    write_expected(
        ipc::edge_edge_mollifier_threshold_gradient(ea0, ea1, eb0, eb1),
        expected);

    check_gpu_matches_host(mollifier_kernel, in, expected);
}

TEST_CASE("GPU signed distances", "[distance][signed][gpu]")
{
    skip_if_no_cuda_device();

    std::vector<double> in = SGN_IN;                 // 2D point-line
    in.insert(in.end(), LL_IN.begin(), LL_IN.end()); // 3D line-line
    in.insert(in.end(), PT_IN.begin(), PT_IN.end()); // 3D point-plane

    const Eigen::Vector2d p2(SGN_IN[0], SGN_IN[1]);
    const Eigen::Vector2d f0(SGN_IN[2], SGN_IN[3]);
    const Eigen::Vector2d f1(SGN_IN[4], SGN_IN[5]);
    const Eigen::Vector3d ea0(LL_IN[0], LL_IN[1], LL_IN[2]);
    const Eigen::Vector3d ea1(LL_IN[3], LL_IN[4], LL_IN[5]);
    const Eigen::Vector3d eb0(LL_IN[6], LL_IN[7], LL_IN[8]);
    const Eigen::Vector3d eb1(LL_IN[9], LL_IN[10], LL_IN[11]);
    const Eigen::Vector3d p(PT_IN[0], PT_IN[1], PT_IN[2]);
    const Eigen::Vector3d t0(PT_IN[3], PT_IN[4], PT_IN[5]);
    const Eigen::Vector3d t1(PT_IN[6], PT_IN[7], PT_IN[8]);
    const Eigen::Vector3d t2(PT_IN[9], PT_IN[10], PT_IN[11]);

    std::vector<double> expected;
    expected.push_back(ipc::point_line_signed_distance(p2, f0, f1));
    write_expected(
        ipc::point_line_signed_distance_gradient(p2, f0, f1), expected);
    write_expected(
        ipc::point_line_signed_distance_hessian(p2, f0, f1), expected);

    expected.push_back(ipc::line_line_signed_distance(ea0, ea1, eb0, eb1));
    write_expected(
        ipc::line_line_signed_distance_gradient(ea0, ea1, eb0, eb1), expected);
    write_expected(
        ipc::line_line_signed_distance_hessian(ea0, ea1, eb0, eb1), expected);

    expected.push_back(ipc::point_plane_signed_distance(p, t0, t1, t2));
    write_expected(
        ipc::point_plane_signed_distance_gradient(p, t0, t1, t2), expected);
    write_expected(
        ipc::point_plane_signed_distance_hessian(p, t0, t1, t2), expected);

    check_gpu_matches_host(signed_distance_kernel, in, expected);
}

#endif // IPC_TOOLKIT_WITH_CUDA
