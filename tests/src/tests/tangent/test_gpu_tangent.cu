// Exercises the __host__ __device__ tangent functions from CUDA kernels.
//
// Compilation proves the framework works: nvcc must instantiate every
// converted tangent function (tangent bases, closest points — including the
// std::array<Matrix12d, 2> hessians — and relative velocities) as device
// code. On a machine with a CUDA device the kernels also run and their
// outputs are checked against the host implementations — host/device parity
// is the contract. Without a device (e.g. a GPU-less container) the runtime
// checks are skipped but the device code has still been compiled.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/tangent/closest_point.hpp>
#include <ipc/tangent/relative_velocity.hpp>
#include <ipc/tangent/tangent_basis.hpp>

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
const std::vector<double> PE_IN = { 0.5, 1, 0.2,  -1, 0, 0,  1, 0, 0.1 };
const std::vector<double> EE_IN = { -1, 0, 0,  1, 0, 0,  -0.5, 1, -1,  0.5, 1, 1 };
const std::vector<double> PT_IN = { 0, 1, 0.1,  -1, 0, 1,  1, 0, 1,  0, 0, -1 };
// Velocities for the relative-velocity kernels (4 x 3D) + coords (2).
const std::vector<double> RV_IN = { 0.1, -0.2, 0.3,  0.4, 0.5, -0.6,
                                    -0.7, 0.8, 0.9,  1.0, -1.1, 1.2,
                                    0.25, 0.5 };
// clang-format on

__global__ void tangent_basis_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const ipc::VectorMax3d p = Eigen::Vector3d(in[0], in[1], in[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(in[3], in[4], in[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(in[6], in[7], in[8]);
    const Eigen::Vector3d ea0(in[9], in[10], in[11]);
    const Eigen::Vector3d ea1(in[12], in[13], in[14]);
    const Eigen::Vector3d eb0(in[15], in[16], in[17]);
    const Eigen::Vector3d eb1(in[18], in[19], in[20]);
    const Eigen::Vector3d q(in[21], in[22], in[23]);
    const Eigen::Vector3d t0(in[24], in[25], in[26]);
    const Eigen::Vector3d t1(in[27], in[28], in[29]);
    const Eigen::Vector3d t2(in[30], in[31], in[32]);

    int offset = 0;
    offset = write_out(ipc::point_point_tangent_basis(p, e0), out, offset);
    offset =
        write_out(ipc::point_point_tangent_basis_jacobian(p, e0), out, offset);
    offset = write_out(ipc::point_edge_tangent_basis(p, e0, e1), out, offset);
    offset = write_out(
        ipc::point_edge_tangent_basis_jacobian(p, e0, e1), out, offset);
    offset = write_out(
        ipc::edge_edge_tangent_basis(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(
        ipc::edge_edge_tangent_basis_jacobian(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(
        ipc::point_triangle_tangent_basis(q, t0, t1, t2), out, offset);
    offset = write_out(
        ipc::point_triangle_tangent_basis_jacobian(q, t0, t1, t2), out, offset);
}

__global__ void closest_point_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const ipc::VectorMax3d p = Eigen::Vector3d(in[0], in[1], in[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(in[3], in[4], in[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(in[6], in[7], in[8]);
    const Eigen::Vector3d ea0(in[9], in[10], in[11]);
    const Eigen::Vector3d ea1(in[12], in[13], in[14]);
    const Eigen::Vector3d eb0(in[15], in[16], in[17]);
    const Eigen::Vector3d eb1(in[18], in[19], in[20]);
    const Eigen::Vector3d q(in[21], in[22], in[23]);
    const Eigen::Vector3d t0(in[24], in[25], in[26]);
    const Eigen::Vector3d t1(in[27], in[28], in[29]);
    const Eigen::Vector3d t2(in[30], in[31], in[32]);

    int offset = 0;
    out[offset++] = ipc::point_edge_closest_point(p, e0, e1);
    offset = write_out(
        ipc::point_edge_closest_point_jacobian(p, e0, e1), out, offset);
    offset = write_out(
        ipc::point_edge_closest_point_hessian(p, e0, e1), out, offset);

    offset = write_out(
        ipc::edge_edge_closest_point(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(
        ipc::edge_edge_closest_point_jacobian(ea0, ea1, eb0, eb1), out, offset);
    const std::array<ipc::Matrix12d, 2> ee_hess =
        ipc::edge_edge_closest_point_hessian(ea0, ea1, eb0, eb1);
    offset = write_out(ee_hess[0], out, offset);
    offset = write_out(ee_hess[1], out, offset);

    offset = write_out(
        ipc::point_triangle_closest_point(q, t0, t1, t2), out, offset);
    offset = write_out(
        ipc::point_triangle_closest_point_jacobian(q, t0, t1, t2), out, offset);
    const std::array<ipc::Matrix12d, 2> pt_hess =
        ipc::point_triangle_closest_point_hessian(q, t0, t1, t2);
    offset = write_out(pt_hess[0], out, offset);
    offset = write_out(pt_hess[1], out, offset);
}

__global__ void relative_velocity_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d dp(in[0], in[1], in[2]);
    const Eigen::Vector3d dv0(in[3], in[4], in[5]);
    const Eigen::Vector3d dv1(in[6], in[7], in[8]);
    const Eigen::Vector3d dv2(in[9], in[10], in[11]);
    const Eigen::Vector2d coords(in[12], in[13]);
    const double alpha = in[12];

    int offset = 0;
    offset =
        write_out(ipc::point_point_relative_velocity(dp, dv0), out, offset);
    offset =
        write_out(ipc::point_point_relative_velocity_jacobian(3), out, offset);
    offset =
        write_out(ipc::point_point_relative_velocity_dx_dbeta(3), out, offset);

    offset = write_out(
        ipc::point_edge_relative_velocity(dp, dv0, dv1, alpha), out, offset);
    offset = write_out(
        ipc::point_edge_relative_velocity_jacobian(3, alpha), out, offset);
    offset = write_out(
        ipc::point_edge_relative_velocity_dx_dbeta(3, alpha), out, offset);

    offset = write_out(
        ipc::edge_edge_relative_velocity(dp, dv0, dv1, dv2, coords), out,
        offset);
    offset = write_out(
        ipc::edge_edge_relative_velocity_jacobian(coords), out, offset);
    offset = write_out(
        ipc::edge_edge_relative_velocity_dx_dbeta(coords), out, offset);

    offset = write_out(
        ipc::point_triangle_relative_velocity(dp, dv0, dv1, dv2, coords), out,
        offset);
    offset = write_out(
        ipc::point_triangle_relative_velocity_jacobian(coords), out, offset);
    offset = write_out(
        ipc::point_triangle_relative_velocity_dx_dbeta(coords), out, offset);
}

std::vector<double> tangent_fixture()
{
    std::vector<double> in = PE_IN;
    in.insert(in.end(), EE_IN.begin(), EE_IN.end());
    in.insert(in.end(), PT_IN.begin(), PT_IN.end());
    return in;
}

} // namespace

TEST_CASE("GPU tangent basis", "[tangent][tangent_basis][gpu]")
{
    skip_if_no_cuda_device();

    const ipc::VectorMax3d p = Eigen::Vector3d(PE_IN[0], PE_IN[1], PE_IN[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(PE_IN[3], PE_IN[4], PE_IN[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(PE_IN[6], PE_IN[7], PE_IN[8]);
    const Eigen::Vector3d ea0(EE_IN[0], EE_IN[1], EE_IN[2]);
    const Eigen::Vector3d ea1(EE_IN[3], EE_IN[4], EE_IN[5]);
    const Eigen::Vector3d eb0(EE_IN[6], EE_IN[7], EE_IN[8]);
    const Eigen::Vector3d eb1(EE_IN[9], EE_IN[10], EE_IN[11]);
    const Eigen::Vector3d q(PT_IN[0], PT_IN[1], PT_IN[2]);
    const Eigen::Vector3d t0(PT_IN[3], PT_IN[4], PT_IN[5]);
    const Eigen::Vector3d t1(PT_IN[6], PT_IN[7], PT_IN[8]);
    const Eigen::Vector3d t2(PT_IN[9], PT_IN[10], PT_IN[11]);

    std::vector<double> expected;
    write_expected(ipc::point_point_tangent_basis(p, e0), expected);
    write_expected(ipc::point_point_tangent_basis_jacobian(p, e0), expected);
    write_expected(ipc::point_edge_tangent_basis(p, e0, e1), expected);
    write_expected(ipc::point_edge_tangent_basis_jacobian(p, e0, e1), expected);
    write_expected(ipc::edge_edge_tangent_basis(ea0, ea1, eb0, eb1), expected);
    write_expected(
        ipc::edge_edge_tangent_basis_jacobian(ea0, ea1, eb0, eb1), expected);
    write_expected(ipc::point_triangle_tangent_basis(q, t0, t1, t2), expected);
    write_expected(
        ipc::point_triangle_tangent_basis_jacobian(q, t0, t1, t2), expected);

    check_gpu_matches_host(tangent_basis_kernel, tangent_fixture(), expected);
}

TEST_CASE("GPU closest point", "[tangent][closest_point][gpu]")
{
    skip_if_no_cuda_device();

    const ipc::VectorMax3d p = Eigen::Vector3d(PE_IN[0], PE_IN[1], PE_IN[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(PE_IN[3], PE_IN[4], PE_IN[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(PE_IN[6], PE_IN[7], PE_IN[8]);
    const Eigen::Vector3d ea0(EE_IN[0], EE_IN[1], EE_IN[2]);
    const Eigen::Vector3d ea1(EE_IN[3], EE_IN[4], EE_IN[5]);
    const Eigen::Vector3d eb0(EE_IN[6], EE_IN[7], EE_IN[8]);
    const Eigen::Vector3d eb1(EE_IN[9], EE_IN[10], EE_IN[11]);
    const Eigen::Vector3d q(PT_IN[0], PT_IN[1], PT_IN[2]);
    const Eigen::Vector3d t0(PT_IN[3], PT_IN[4], PT_IN[5]);
    const Eigen::Vector3d t1(PT_IN[6], PT_IN[7], PT_IN[8]);
    const Eigen::Vector3d t2(PT_IN[9], PT_IN[10], PT_IN[11]);

    std::vector<double> expected;
    expected.push_back(ipc::point_edge_closest_point(p, e0, e1));
    write_expected(ipc::point_edge_closest_point_jacobian(p, e0, e1), expected);
    write_expected(ipc::point_edge_closest_point_hessian(p, e0, e1), expected);
    write_expected(ipc::edge_edge_closest_point(ea0, ea1, eb0, eb1), expected);
    write_expected(
        ipc::edge_edge_closest_point_jacobian(ea0, ea1, eb0, eb1), expected);
    const std::array<ipc::Matrix12d, 2> ee_hess =
        ipc::edge_edge_closest_point_hessian(ea0, ea1, eb0, eb1);
    write_expected(ee_hess[0], expected);
    write_expected(ee_hess[1], expected);
    write_expected(ipc::point_triangle_closest_point(q, t0, t1, t2), expected);
    write_expected(
        ipc::point_triangle_closest_point_jacobian(q, t0, t1, t2), expected);
    const std::array<ipc::Matrix12d, 2> pt_hess =
        ipc::point_triangle_closest_point_hessian(q, t0, t1, t2);
    write_expected(pt_hess[0], expected);
    write_expected(pt_hess[1], expected);

    check_gpu_matches_host(closest_point_kernel, tangent_fixture(), expected);
}

TEST_CASE("GPU relative velocity", "[tangent][relative_velocity][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d dp(RV_IN[0], RV_IN[1], RV_IN[2]);
    const Eigen::Vector3d dv0(RV_IN[3], RV_IN[4], RV_IN[5]);
    const Eigen::Vector3d dv1(RV_IN[6], RV_IN[7], RV_IN[8]);
    const Eigen::Vector3d dv2(RV_IN[9], RV_IN[10], RV_IN[11]);
    const Eigen::Vector2d coords(RV_IN[12], RV_IN[13]);
    const double alpha = RV_IN[12];

    std::vector<double> expected;
    write_expected(ipc::point_point_relative_velocity(dp, dv0), expected);
    write_expected(ipc::point_point_relative_velocity_jacobian(3), expected);
    write_expected(ipc::point_point_relative_velocity_dx_dbeta(3), expected);
    write_expected(
        ipc::point_edge_relative_velocity(dp, dv0, dv1, alpha), expected);
    write_expected(
        ipc::point_edge_relative_velocity_jacobian(3, alpha), expected);
    write_expected(
        ipc::point_edge_relative_velocity_dx_dbeta(3, alpha), expected);
    write_expected(
        ipc::edge_edge_relative_velocity(dp, dv0, dv1, dv2, coords), expected);
    write_expected(ipc::edge_edge_relative_velocity_jacobian(coords), expected);
    write_expected(ipc::edge_edge_relative_velocity_dx_dbeta(coords), expected);
    write_expected(
        ipc::point_triangle_relative_velocity(dp, dv0, dv1, dv2, coords),
        expected);
    write_expected(
        ipc::point_triangle_relative_velocity_jacobian(coords), expected);
    write_expected(
        ipc::point_triangle_relative_velocity_dx_dbeta(coords), expected);

    check_gpu_matches_host(relative_velocity_kernel, RV_IN, expected);
}

#endif // IPC_TOOLKIT_WITH_CUDA
