// Exercises ipc::cuda::BarrierPotential against the CPU BarrierPotential.
//
// Compilation proves the kernels instantiate as device code. On a machine
// with a CUDA device the potential, gradient, and hessian are also computed
// on the GPU and checked against the CPU implementation — host/device parity
// is the contract. Without a device the runtime checks are skipped but the
// device code has still been compiled.

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/cuda/barrier_potential.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>

#include <cuda_runtime.h>

#include <cmath>
#include <string>

namespace {

bool has_cuda_device()
{
    int device_count = 0;
    const cudaError_t err = cudaGetDeviceCount(&device_count);
    return err == cudaSuccess && device_count > 0;
}

} // namespace

TEST_CASE(
    "GPU barrier potential parity", "[potential][barrier_potential][cuda][gpu]")
{
    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

    const bool use_area_weighting = GENERATE(true, false);
    const ipc::NormalCollisions::CollisionSetType collision_set_type = GENERATE(
        ipc::NormalCollisions::CollisionSetType::IPC,
        ipc::NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX);
    const bool use_physical_barrier = GENERATE(true, false);

    double dhat = -1;
    std::string mesh_name;
    SECTION("cube")
    {
        dhat = sqrt(2.0);
        mesh_name = "cube.ply";
    }
#ifdef NDEBUG
    SECTION("two cubes close")
    {
        dhat = 1e-1;
        mesh_name = "two-cubes-close.ply";
    }
#endif

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(ipc::tests::load_mesh(mesh_name, vertices, edges, faces));
    CAPTURE(
        mesh_name, use_area_weighting, collision_set_type,
        use_physical_barrier);

    const ipc::CollisionMesh mesh(vertices, edges, faces);

    ipc::NormalCollisions collisions;
    collisions.set_use_area_weighting(use_area_weighting);
    collisions.set_collision_set_type(collision_set_type);
    collisions.build(mesh, vertices, dhat);
    REQUIRE(!collisions.empty());

    const double kappa = 1.0;
    const ipc::BarrierPotential cpu_potential(
        dhat, kappa, use_physical_barrier);
    const ipc::cuda::BarrierPotential gpu_potential(
        dhat, kappa, use_physical_barrier);

    // ------------------------------------------------------------------
    // Drop-in API

    const double cpu_value = cpu_potential(collisions, mesh, vertices);
    const double gpu_value = gpu_potential(collisions, mesh, vertices);
    CHECK(gpu_value == Catch::Approx(cpu_value).epsilon(1e-12));

    const Eigen::VectorXd cpu_grad =
        cpu_potential.gradient(collisions, mesh, vertices);
    const Eigen::VectorXd gpu_grad =
        gpu_potential.gradient(collisions, mesh, vertices);
    REQUIRE(gpu_grad.size() == cpu_grad.size());
    CHECK((gpu_grad - cpu_grad).norm() <= 1e-10 * (1 + cpu_grad.norm()));

    for (const ipc::PSDProjectionMethod method :
         { ipc::PSDProjectionMethod::NONE, ipc::PSDProjectionMethod::CLAMP,
           ipc::PSDProjectionMethod::ABS }) {
        CAPTURE(static_cast<int>(method));
        const Eigen::SparseMatrix<double> cpu_hess =
            cpu_potential.hessian(collisions, mesh, vertices, method);
        const Eigen::SparseMatrix<double> gpu_hess =
            gpu_potential.hessian(collisions, mesh, vertices, method);
        REQUIRE(gpu_hess.rows() == cpu_hess.rows());
        REQUIRE(gpu_hess.cols() == cpu_hess.cols());
        // NONE: the GPU assembly sums the same device-computed blocks as the
        // CPU (just in a different order) -> tight. CLAMP/ABS: the GPU projects
        // with cuSOLVER's batched Jacobi while the CPU uses Eigen's QR
        // eigensolver, so the projected blocks differ by conditioning-amplified
        // rounding (worst near the near-zero eigenvalue clusters of contact
        // hessians) -> looser tolerance.
        const double tol =
            (method == ipc::PSDProjectionMethod::NONE) ? 1e-10 : 1e-6;
        CHECK((gpu_hess - cpu_hess).norm() <= tol * (1 + cpu_hess.norm()));
    }

    // ------------------------------------------------------------------
    // Fast-path API (prebuilt device objects, reused across evaluations)

    const ipc::cuda::DeviceNormalCollisions device_collisions(collisions, mesh);
    ipc::cuda::DeviceVertices device_vertices(vertices);

    CHECK(device_collisions.size() == collisions.size());

    CHECK(
        gpu_potential(device_collisions, device_vertices)
        == Catch::Approx(cpu_value).epsilon(1e-12));

    // Re-upload (exercises DeviceVertices::update) and re-evaluate.
    device_vertices.update(vertices);

    const Eigen::VectorXd gpu_grad_fast =
        gpu_potential.gradient(device_collisions, device_vertices);
    CHECK((gpu_grad_fast - cpu_grad).norm() <= 1e-10 * (1 + cpu_grad.norm()));

    const Eigen::SparseMatrix<double> cpu_hess = cpu_potential.hessian(
        collisions, mesh, vertices, ipc::PSDProjectionMethod::CLAMP);
    const Eigen::SparseMatrix<double> gpu_hess_fast = gpu_potential.hessian(
        device_collisions, device_vertices, ipc::PSDProjectionMethod::CLAMP);
    CHECK((gpu_hess_fast - cpu_hess).norm() <= 1e-6 * (1 + cpu_hess.norm()));
}

TEST_CASE(
    "GPU barrier potential unsupported inputs",
    "[potential][barrier_potential][cuda][gpu]")
{
    // These throw host-side regardless of whether a device is present.
    CHECK_THROWS(
        ipc::cuda::BarrierPotential(
            /*dhat=*/1e-3, /*stiffness=*/1.0, /*use_physical_barrier=*/false,
            ipc::BarrierType::CUBIC));

    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

    // 2D vertices are rejected.
    // Eigen::MatrixXd vertices_2d(3, 2);
    // vertices_2d << 0, 0, 1, 0, 0, 1;
    // CHECK_THROWS(ipc::cuda::DeviceVertices(vertices_2d));
}

TEST_CASE(
    "Benchmark GPU barrier potential",
    "[!benchmark][potential][barrier_potential][cuda][gpu]")
{
    if (!has_cuda_device()) {
        SKIP("No CUDA device available; kernels compiled but not executed.");
    }

    // A frame of the cloth-on-ball simulation scene (46.6k vertices, 92.2k
    // faces) with the cloth in resting contact with the ball.
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges, faces;
    REQUIRE(ipc::tests::load_mesh("cloth_ball93.ply", vertices, edges, faces));

    // Standard IPC activation distance: 10⁻³ of the bounding-box diagonal.
    const double dhat = 1e-3 * ipc::world_bbox_diagonal_length(vertices);

    const ipc::CollisionMesh mesh(vertices, edges, faces);

    ipc::NormalCollisions collisions;
    collisions.build(mesh, vertices, dhat);
    REQUIRE(!collisions.empty());
    CAPTURE(dhat, collisions.size());

    const double kappa = 1.0;
    const ipc::BarrierPotential cpu_potential(dhat, kappa);
    const ipc::cuda::BarrierPotential gpu_potential(dhat, kappa);

    // -- CPU reference -------------------------------------------------

    BENCHMARK("CPU potential")
    {
        return cpu_potential(collisions, mesh, vertices);
    };
    BENCHMARK("CPU gradient")
    {
        return cpu_potential.gradient(collisions, mesh, vertices);
    };
    BENCHMARK("CPU hessian")
    {
        return cpu_potential.hessian(collisions, mesh, vertices);
    };
    BENCHMARK("CPU hessian with PSD projection")
    {
        return cpu_potential.hessian(
            collisions, mesh, vertices, ipc::PSDProjectionMethod::CLAMP);
    };

    // -- GPU drop-in API (re-uploads the collisions and vertices) -------

    BENCHMARK("GPU potential (incl. upload)")
    {
        return gpu_potential(collisions, mesh, vertices);
    };

    // -- GPU fast-path API (device-resident inputs) ----------------------

    const ipc::cuda::DeviceNormalCollisions device_collisions(collisions, mesh);
    ipc::cuda::DeviceVertices device_vertices(vertices);

    BENCHMARK("Build DeviceNormalCollisions")
    {
        return ipc::cuda::DeviceNormalCollisions(collisions, mesh);
    };
    BENCHMARK("DeviceVertices::update") { device_vertices.update(vertices); };
    BENCHMARK("GPU potential")
    {
        return gpu_potential(device_collisions, device_vertices);
    };
    BENCHMARK("GPU gradient")
    {
        return gpu_potential.gradient(device_collisions, device_vertices);
    };
    BENCHMARK("GPU hessian")
    {
        return gpu_potential.hessian(device_collisions, device_vertices);
    };
    BENCHMARK("GPU hessian with PSD projection")
    {
        return gpu_potential.hessian(
            device_collisions, device_vertices,
            ipc::PSDProjectionMethod::CLAMP);
    };
}

#endif // IPC_TOOLKIT_WITH_CUDA
