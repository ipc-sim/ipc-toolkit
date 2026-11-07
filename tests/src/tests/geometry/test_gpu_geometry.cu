// Exercises the __host__ __device__ geometry functions from CUDA kernels.
//
// Compilation proves the framework works: nvcc must instantiate every
// converted geometry function (and its Eigen/std::tuple/std::array machinery)
// as device code. On a machine with a CUDA device the kernels also run and
// their outputs are checked against the host implementations — host/device
// parity is the contract. Without a device (e.g. a GPU-less container) the
// runtime checks are skipped but the device code has still been compiled.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/geometry/angle.hpp>
#include <ipc/geometry/area.hpp>
#include <ipc/geometry/normal.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <tests/gpu_utils.hpp>

#include <cuda_runtime.h>

#include <vector>

using namespace ipc::tests;

namespace {

// ---------------------------------------------------------------------------

__global__ void area_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d e0(in[0], in[1], in[2]);
    const Eigen::Vector3d e1(in[3], in[4], in[5]);
    const Eigen::Vector3d t0(in[6], in[7], in[8]);
    const Eigen::Vector3d t1(in[9], in[10], in[11]);
    const Eigen::Vector3d t2(in[12], in[13], in[14]);

    int offset = 0;
    out[offset++] = ipc::edge_length(e0, e1);
    offset = write_out(ipc::edge_length_gradient(e0, e1), out, offset);
    out[offset++] = ipc::triangle_area(t0, t1, t2);
    offset = write_out(ipc::triangle_area_gradient(t0, t1, t2), out, offset);
}

__global__ void angle_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d x0(in[0], in[1], in[2]);
    const Eigen::Vector3d x1(in[3], in[4], in[5]);
    const Eigen::Vector3d x2(in[6], in[7], in[8]);
    const Eigen::Vector3d x3(in[9], in[10], in[11]);

    int offset = 0;
    out[offset++] = ipc::dihedral_angle(x0, x1, x2, x3);
    offset =
        write_out(ipc::dihedral_angle_gradient(x0, x1, x2, x3), out, offset);
    offset =
        write_out(ipc::dihedral_angle_hessian(x0, x1, x2, x3), out, offset);
}

__global__ void normalization_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const ipc::VectorMax3d x = Eigen::Vector3d(in[0], in[1], in[2]);

    int offset = 0;
    const auto [xhat, J] = ipc::normalization_and_jacobian(x);
    offset = write_out(xhat, out, offset);
    offset = write_out(J, out, offset);

    const auto [xhat2, J2, H] = ipc::normalization_and_jacobian_and_hessian(x);
    offset = write_out(xhat2, out, offset);
    offset = write_out(J2, out, offset);
    for (const auto& Hi : H) {
        offset = write_out(Hi, out, offset);
    }

    offset = write_out(ipc::cross_product_matrix(x), out, offset);
    offset =
        write_out(ipc::cross_product_matrix_jacobian<double>(), out, offset);
}

__global__ void point_line_normal_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const ipc::VectorMax3d p = Eigen::Vector3d(in[0], in[1], in[2]);
    const ipc::VectorMax3d e0 = Eigen::Vector3d(in[3], in[4], in[5]);
    const ipc::VectorMax3d e1 = Eigen::Vector3d(in[6], in[7], in[8]);

    int offset = 0;
    offset =
        write_out(ipc::point_line_unnormalized_normal(p, e0, e1), out, offset);
    offset = write_out(ipc::point_line_normal(p, e0, e1), out, offset);
    offset = write_out(
        ipc::point_line_unnormalized_normal_jacobian(p, e0, e1), out, offset);
    offset = write_out(ipc::point_line_normal_jacobian(p, e0, e1), out, offset);
    offset = write_out(
        ipc::point_line_unnormalized_normal_hessian(p, e0, e1), out, offset);
    offset = write_out(ipc::point_line_normal_hessian(p, e0, e1), out, offset);
}

__global__ void triangle_normal_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d a(in[0], in[1], in[2]);
    const Eigen::Vector3d b(in[3], in[4], in[5]);
    const Eigen::Vector3d c(in[6], in[7], in[8]);

    int offset = 0;
    offset = write_out(ipc::triangle_unnormalized_normal(a, b, c), out, offset);
    offset = write_out(ipc::triangle_normal(a, b, c), out, offset);
    offset = write_out(
        ipc::triangle_unnormalized_normal_jacobian(a, b, c), out, offset);
    offset = write_out(ipc::triangle_normal_jacobian(a, b, c), out, offset);
    offset = write_out(
        ipc::triangle_unnormalized_normal_hessian(a, b, c), out, offset);
    offset = write_out(ipc::triangle_normal_hessian(a, b, c), out, offset);
}

__global__ void line_line_normal_kernel(const double* in, double* out)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) {
        return;
    }
    const Eigen::Vector3d ea0(in[0], in[1], in[2]);
    const Eigen::Vector3d ea1(in[3], in[4], in[5]);
    const Eigen::Vector3d eb0(in[6], in[7], in[8]);
    const Eigen::Vector3d eb1(in[9], in[10], in[11]);

    int offset = 0;
    offset = write_out(
        ipc::line_line_unnormalized_normal(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(ipc::line_line_normal(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(
        ipc::line_line_unnormalized_normal_jacobian(ea0, ea1, eb0, eb1), out,
        offset);
    offset = write_out(
        ipc::line_line_normal_jacobian(ea0, ea1, eb0, eb1), out, offset);
    offset = write_out(
        ipc::line_line_unnormalized_normal_hessian(ea0, ea1, eb0, eb1), out,
        offset);
    offset = write_out(
        ipc::line_line_normal_hessian(ea0, ea1, eb0, eb1), out, offset);
}

} // namespace

TEST_CASE("GPU edge length and triangle area", "[geometry][area][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d e0(0, 0, 0), e1(1, 2, 2);
    const Eigen::Vector3d t0(-1, 0, 1), t1(1, 0, 1), t2(0, 0, -1);

    std::vector<double> expected;
    expected.push_back(ipc::edge_length(e0, e1));
    write_expected(ipc::edge_length_gradient(e0, e1), expected);
    expected.push_back(ipc::triangle_area(t0, t1, t2));
    write_expected(ipc::triangle_area_gradient(t0, t1, t2), expected);

    check_gpu_matches_host(
        area_kernel, { 0, 0, 0, 1, 2, 2, -1, 0, 1, 1, 0, 1, 0, 0, -1 },
        expected);
}

TEST_CASE("GPU dihedral angle", "[geometry][angle][gpu]")
{
    skip_if_no_cuda_device();

    // Two triangles sharing the edge (x0, x1), folded at 90°.
    const Eigen::Vector3d x0(0, -1, 0), x1(0, 1, 0);
    const Eigen::Vector3d x2(0.5, 0, 0.5), x3(-0.5, 0, 0.5);

    std::vector<double> expected;
    expected.push_back(ipc::dihedral_angle(x0, x1, x2, x3));
    write_expected(ipc::dihedral_angle_gradient(x0, x1, x2, x3), expected);
    write_expected(ipc::dihedral_angle_hessian(x0, x1, x2, x3), expected);

    check_gpu_matches_host(
        angle_kernel, { 0, -1, 0, 0, 1, 0, 0.5, 0, 0.5, -0.5, 0, 0.5 },
        expected);
}

TEST_CASE("GPU normalization and cross product", "[geometry][normal][gpu]")
{
    skip_if_no_cuda_device();

    const ipc::VectorMax3d x = Eigen::Vector3d(1, 2, 2);

    std::vector<double> expected;
    const auto [xhat, J] = ipc::normalization_and_jacobian(x);
    write_expected(xhat, expected);
    write_expected(J, expected);
    const auto [xhat2, J2, H] = ipc::normalization_and_jacobian_and_hessian(x);
    write_expected(xhat2, expected);
    write_expected(J2, expected);
    for (const auto& Hi : H) {
        write_expected(Hi, expected);
    }
    write_expected(
        ipc::cross_product_matrix(Eigen::Vector3d(1, 2, 2)), expected);
    write_expected(ipc::cross_product_matrix_jacobian<double>(), expected);

    check_gpu_matches_host(normalization_kernel, { 1, 2, 2 }, expected);
}

TEST_CASE("GPU point-line normal", "[geometry][normal][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d p(0.1, 1, 0.2), e0(-1, 0, 0), e1(1, 0, 0.1);

    std::vector<double> expected;
    write_expected(ipc::point_line_unnormalized_normal(p, e0, e1), expected);
    write_expected(ipc::point_line_normal(p, e0, e1), expected);
    write_expected(
        ipc::point_line_unnormalized_normal_jacobian(p, e0, e1), expected);
    write_expected(ipc::point_line_normal_jacobian(p, e0, e1), expected);
    write_expected(
        ipc::point_line_unnormalized_normal_hessian(p, e0, e1), expected);
    write_expected(ipc::point_line_normal_hessian(p, e0, e1), expected);

    check_gpu_matches_host(
        point_line_normal_kernel, { 0.1, 1, 0.2, -1, 0, 0, 1, 0, 0.1 },
        expected);
}

TEST_CASE("GPU triangle normal", "[geometry][normal][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d a(-1, 0, 1), b(1, 0, 1), c(0, 0.1, -1);

    std::vector<double> expected;
    write_expected(ipc::triangle_unnormalized_normal(a, b, c), expected);
    write_expected(ipc::triangle_normal(a, b, c), expected);
    write_expected(
        ipc::triangle_unnormalized_normal_jacobian(a, b, c), expected);
    write_expected(ipc::triangle_normal_jacobian(a, b, c), expected);
    write_expected(
        ipc::triangle_unnormalized_normal_hessian(a, b, c), expected);
    write_expected(ipc::triangle_normal_hessian(a, b, c), expected);

    check_gpu_matches_host(
        triangle_normal_kernel, { -1, 0, 1, 1, 0, 1, 0, 0.1, -1 }, expected);
}

TEST_CASE("GPU line-line normal", "[geometry][normal][gpu]")
{
    skip_if_no_cuda_device();

    const Eigen::Vector3d ea0(-1, 0, 0), ea1(1, 0, 0);
    const Eigen::Vector3d eb0(0, 1, -1), eb1(0.1, 1, 1);

    std::vector<double> expected;
    write_expected(
        ipc::line_line_unnormalized_normal(ea0, ea1, eb0, eb1), expected);
    write_expected(ipc::line_line_normal(ea0, ea1, eb0, eb1), expected);
    write_expected(
        ipc::line_line_unnormalized_normal_jacobian(ea0, ea1, eb0, eb1),
        expected);
    write_expected(
        ipc::line_line_normal_jacobian(ea0, ea1, eb0, eb1), expected);
    write_expected(
        ipc::line_line_unnormalized_normal_hessian(ea0, ea1, eb0, eb1),
        expected);
    write_expected(ipc::line_line_normal_hessian(ea0, ea1, eb0, eb1), expected);

    check_gpu_matches_host(
        line_line_normal_kernel, { -1, 0, 0, 1, 0, 0, 0, 1, -1, 0.1, 1, 1 },
        expected);
}

#endif // IPC_TOOLKIT_WITH_CUDA
