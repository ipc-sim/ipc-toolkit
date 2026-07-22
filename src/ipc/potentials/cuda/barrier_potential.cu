#include "barrier_potential.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/barrier/barrier.hpp>
#include <ipc/collisions/normal/cuda/device_normal_collisions_impl.cuh>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/potentials/cuda/hessian_assembly.cuh>
#include <ipc/potentials/cuda/psd_projection.cuh>
#include <ipc/utils/cuda/device_utils.cuh>
#include <ipc/utils/logger.hpp>

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>

#include <array>
#include <limits>
#include <vector>

namespace ipc::cuda {

namespace {

    // Padded block layout consumed by the batched eigensolver and assembler:
    // each collision's ndof×ndof block lives in the top-left of a zero-padded
    // 12×12 column-major slot.
    constexpr int PADDED_DIM = BatchedPSDProjection::DIM;
    constexpr int PADDED_BLOCK = BatchedPSDProjection::BLOCK_SIZE;

    // ========================================================================
    // Device barrier dispatch (mirrors BarrierPotential's scalar hooks in
    // src/ipc/potentials/barrier_potential.cpp)

    /// Constant parameters shared by every collision in a launch.
    struct BarrierParams {
        double dhat;
        double stiffness;
        bool use_physical_barrier;
        BarrierType barrier_type;
    };

    __device__ double
    device_barrier(const BarrierType type, const double d, const double dhat)
    {
        switch (type) {
        case BarrierType::CLAMPED_LOG:
            return barrier(d, dhat);
        case BarrierType::NORMALIZED_CLAMPED_LOG:
            return barrier(d / dhat, 1.0);
        default:
            assert(false && "Unsupported barrier type on the GPU!");
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    __device__ double device_barrier_first_derivative(
        const BarrierType type, const double d, const double dhat)
    {
        switch (type) {
        case BarrierType::CLAMPED_LOG:
            return barrier_first_derivative(d, dhat);
        case BarrierType::NORMALIZED_CLAMPED_LOG:
            return barrier_first_derivative(d / dhat, 1.0) / dhat;
        default:
            assert(false && "Unsupported barrier type on the GPU!");
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    __device__ double device_barrier_second_derivative(
        const BarrierType type, const double d, const double dhat)
    {
        switch (type) {
        case BarrierType::CLAMPED_LOG:
            return barrier_second_derivative(d, dhat);
        case BarrierType::NORMALIZED_CLAMPED_LOG:
            return barrier_second_derivative(d / dhat, 1.0) / (dhat * dhat);
        default:
            assert(false && "Unsupported barrier type on the GPU!");
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    __device__ double
    device_barrier_units(const BarrierType type, const double dhat)
    {
        switch (type) {
        case BarrierType::CLAMPED_LOG:
            return dhat * dhat;
        case BarrierType::NORMALIZED_CLAMPED_LOG:
            return 1.0; // The normalized barrier is dimensionless.
        default:
            assert(false && "Unsupported barrier type on the GPU!");
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    /// f(d²): mirrors BarrierPotential::operator()(d², dmin).
    __device__ double potential_value(
        const double d_sqr, const double dmin, const BarrierParams& p)
    {
        const double dhat_eff = (2 * dmin + p.dhat) * p.dhat;
        double b =
            device_barrier(p.barrier_type, d_sqr - dmin * dmin, dhat_eff);
        if (p.use_physical_barrier) {
            b *= p.dhat / device_barrier_units(p.barrier_type, dhat_eff);
        }
        return p.stiffness * b;
    }

    /// f'(d²): mirrors BarrierPotential::gradient(d², dmin).
    __device__ double potential_gradient(
        const double d_sqr, const double dmin, const BarrierParams& p)
    {
        const double dhat_eff = (2 * dmin + p.dhat) * p.dhat;
        double db = device_barrier_first_derivative(
            p.barrier_type, d_sqr - dmin * dmin, dhat_eff);
        if (p.use_physical_barrier) {
            db *= p.dhat / device_barrier_units(p.barrier_type, dhat_eff);
        }
        return p.stiffness * db;
    }

    /// f″(d²): mirrors BarrierPotential::hessian(d², dmin).
    __device__ double potential_hessian(
        const double d_sqr, const double dmin, const BarrierParams& p)
    {
        const double dhat_eff = (2 * dmin + p.dhat) * p.dhat;
        double d2b = device_barrier_second_derivative(
            p.barrier_type, d_sqr - dmin * dmin, dhat_eff);
        if (p.use_physical_barrier) {
            d2b *= p.dhat / device_barrier_units(p.barrier_type, dhat_eff);
        }
        return p.stiffness * d2b;
    }

    // ========================================================================
    // Per-collision-type stencil evaluation

    enum class CollisionType : uint8_t { VV, EV, EE, FV };

    template <CollisionType CT> constexpr int num_stencil_vertices()
    {
        if constexpr (CT == CollisionType::VV) {
            return 2;
        } else if constexpr (CT == CollisionType::EV) {
            return 3;
        } else {
            return 4;
        }
    }

    /// Raw device pointers of a DeviceNormalCollisions::Impl.
    struct CollisionsView {
        const index_t* vertex_id_0;
        const index_t* vertex_id_1;
        const index_t* vertex_id_2;
        const index_t* vertex_id_3;
        const double* weight;
        const double* dmin;
        const double* ee_eps_x;  ///< Indexed by (i - ee_begin).
        const uint8_t* ee_dtype; ///< Indexed by (i - ee_begin).
        size_t ee_begin;
    };

    template <CollisionType CT>
    __device__ void gather_stencil(
        const CollisionsView& view,
        const double* positions,
        const size_t i,
        index_t ids[4],
        VectorMax12d& x)
    {
        constexpr int NV = num_stencil_vertices<CT>();
        ids[0] = view.vertex_id_0[i];
        ids[1] = view.vertex_id_1[i];
        ids[2] = view.vertex_id_2[i];
        ids[3] = view.vertex_id_3[i];
        x.resize(3 * NV);
        for (int k = 0; k < NV; k++) {
            assert(ids[k] >= 0);
            for (int d = 0; d < 3; d++) {
                x[3 * k + d] = positions[3 * ids[k] + d];
            }
        }
    }

    /// d(x)²: mirrors the candidates' compute_distance() (3D).
    template <CollisionType CT>
    __device__ double
    compute_distance(const VectorMax12d& x, const EdgeEdgeDistanceType ee_dtype)
    {
        if constexpr (CT == CollisionType::VV) {
            return point_point_distance(x.head<3>(), x.tail<3>());
        } else if constexpr (CT == CollisionType::EV) {
            return point_edge_distance(
                x.head<3>(), x.segment<3>(3), x.tail<3>(),
                PointEdgeDistanceType::P_E);
        } else if constexpr (CT == CollisionType::EE) {
            return edge_edge_distance(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                ee_dtype);
        } else {
            return point_triangle_distance(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                PointTriangleDistanceType::P_T);
        }
    }

    /// ∇d(x)²: mirrors the candidates' compute_distance_gradient() (3D).
    template <CollisionType CT>
    __device__ VectorMax12d compute_distance_gradient(
        const VectorMax12d& x, const EdgeEdgeDistanceType ee_dtype)
    {
        if constexpr (CT == CollisionType::VV) {
            return point_point_distance_gradient(x.head<3>(), x.tail<3>());
        } else if constexpr (CT == CollisionType::EV) {
            return point_edge_distance_gradient(
                x.head<3>(), x.segment<3>(3), x.tail<3>(),
                PointEdgeDistanceType::P_E);
        } else if constexpr (CT == CollisionType::EE) {
            return edge_edge_distance_gradient(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                ee_dtype);
        } else {
            return point_triangle_distance_gradient(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                PointTriangleDistanceType::P_T);
        }
    }

    /// ∇²d(x)²: mirrors the candidates' compute_distance_hessian() (3D).
    template <CollisionType CT>
    __device__ MatrixMax12d compute_distance_hessian(
        const VectorMax12d& x, const EdgeEdgeDistanceType ee_dtype)
    {
        if constexpr (CT == CollisionType::VV) {
            return point_point_distance_hessian(x.head<3>(), x.tail<3>());
        } else if constexpr (CT == CollisionType::EV) {
            return point_edge_distance_hessian(
                x.head<3>(), x.segment<3>(3), x.tail<3>(),
                PointEdgeDistanceType::P_E);
        } else if constexpr (CT == CollisionType::EE) {
            return edge_edge_distance_hessian(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                ee_dtype);
        } else {
            return point_triangle_distance_hessian(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                PointTriangleDistanceType::P_T);
        }
    }

    template <CollisionType CT>
    __device__ EdgeEdgeDistanceType
    collision_ee_dtype(const CollisionsView& view, const size_t i)
    {
        if constexpr (CT == CollisionType::EE) {
            return static_cast<EdgeEdgeDistanceType>(
                view.ee_dtype[i - view.ee_begin]);
        } else {
            return EdgeEdgeDistanceType::AUTO; // unused
        }
    }

    // ========================================================================
    // Kernels (mirror NormalPotential::operator()/gradient()/hessian() in
    // src/ipc/potentials/normal_potential.cpp; only edge-edge collisions are
    // mollified)

    template <CollisionType CT>
    __global__ void value_kernel(
        const CollisionsView view,
        const double* positions,
        const BarrierParams params,
        const size_t begin,
        const size_t count,
        double* values)
    {
        const size_t k =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (k >= count) {
            return;
        }
        const size_t i = begin + k;

        index_t ids[4];
        VectorMax12d x;
        gather_stencil<CT>(view, positions, i, ids, x);

        const EdgeEdgeDistanceType ee_dtype = collision_ee_dtype<CT>(view, i);
        const double d = compute_distance<CT>(x, ee_dtype);

        // w * m(x) * f(d(x))
        double m = 1.0;
        if constexpr (CT == CollisionType::EE) {
            m = edge_edge_mollifier(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                view.ee_eps_x[i - view.ee_begin]);
        }
        values[i] =
            view.weight[i] * m * potential_value(d, view.dmin[i], params);
    }

    template <CollisionType CT>
    __global__ void gradient_kernel(
        const CollisionsView view,
        const double* positions,
        const index_t n_total_vertices,
        const BarrierParams params,
        const size_t begin,
        const size_t count,
        double* grad)
    {
        const size_t k =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (k >= count) {
            return;
        }
        const size_t i = begin + k;

        index_t ids[4];
        VectorMax12d x;
        gather_stencil<CT>(view, positions, i, ids, x);

        const EdgeEdgeDistanceType ee_dtype = collision_ee_dtype<CT>(view, i);
        const double weight = view.weight[i];
        const double dmin = view.dmin[i];

        VectorMax12d local_grad;
        if constexpr (CT == CollisionType::EE) {
            // Mollified path: ∇[m(x) f(d(x))] = f ∇m + m f' ∇d
            const double eps_x = view.ee_eps_x[i - view.ee_begin];
            const double m = edge_edge_mollifier(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                eps_x);
            // When m = 0 the mollifier also satisfies ∇m = 0, so the full
            // gradient is zero without needing to evaluate f'(d).
            if (m <= 0) {
                return;
            }

            const double d = compute_distance<CT>(x, ee_dtype);
            const VectorMax12d grad_d =
                compute_distance_gradient<CT>(x, ee_dtype);
            const double f = potential_value(d, dmin, params);
            const double grad_f = potential_gradient(d, dmin, params);
            const Vector12d grad_m = edge_edge_mollifier_gradient(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                eps_x);

            local_grad = (weight * f) * grad_m + (weight * m * grad_f) * grad_d;
        } else {
            // ∇[f(d(x))] = f'(d(x)) ∇d(x)
            const double d = compute_distance<CT>(x, ee_dtype);
            const VectorMax12d grad_d =
                compute_distance_gradient<CT>(x, ee_dtype);
            const double grad_f = potential_gradient(d, dmin, params);

            local_grad = (weight * grad_f) * grad_d;
        }

        constexpr int NV = num_stencil_vertices<CT>();
        for (int v = 0; v < NV; v++) {
            for (int d = 0; d < 3; d++) {
                atomicAdd(
                    &grad[global_dof_index(ids[v], d, n_total_vertices)],
                    local_grad[3 * v + d]);
            }
        }
    }

    template <CollisionType CT>
    __global__ void hessian_kernel(
        const CollisionsView view,
        const double* positions,
        const BarrierParams params,
        const size_t begin,
        const size_t count,
        double* blocks) // n × PADDED_BLOCK doubles, zero-initialized
    {
        const size_t k =
            static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
        if (k >= count) {
            return;
        }
        const size_t i = begin + k;

        index_t ids[4];
        VectorMax12d x;
        gather_stencil<CT>(view, positions, i, ids, x);

        const EdgeEdgeDistanceType ee_dtype = collision_ee_dtype<CT>(view, i);
        const double weight = view.weight[i];
        const double dmin = view.dmin[i];

        constexpr int NDOF = 3 * num_stencil_vertices<CT>();
        Eigen::Matrix<double, NDOF, NDOF> hess;

        const double d = compute_distance<CT>(x, ee_dtype);

        if constexpr (CT == CollisionType::EE) {
            const double eps_x = view.ee_eps_x[i - view.ee_begin];
            const double m = edge_edge_mollifier(
                x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                eps_x);
            if (m <= 0) {
                // ∇²[m(x) f(d(x))] = f ∇²m (∇m = 0 when m = 0)
                const double f = potential_value(d, dmin, params);
                const Matrix12d hess_m = edge_edge_mollifier_hessian(
                    x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                    eps_x);
                hess = (weight * f) * hess_m;
            } else {
                const VectorMax12d grad_d =
                    compute_distance_gradient<CT>(x, ee_dtype);
                const MatrixMax12d hess_d =
                    compute_distance_hessian<CT>(x, ee_dtype);
                const double f = potential_value(d, dmin, params);
                const double grad_f = potential_gradient(d, dmin, params);
                const double hess_f = potential_hessian(d, dmin, params);
                const Vector12d grad_m = edge_edge_mollifier_gradient(
                    x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                    eps_x);
                const Matrix12d hess_m = edge_edge_mollifier_hessian(
                    x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>(),
                    eps_x);

                // ∇f(d(x)) ∇m(x)ᵀ
                const Matrix12d grad_f_grad_m =
                    (weight * grad_f) * grad_d * grad_m.transpose();

                // ∇²[m(x) f(d(x))] = ∇²m f + ∇f ∇mᵀ + ∇m ∇fᵀ + m ∇²f
                hess = (weight * f) * hess_m + grad_f_grad_m
                    + grad_f_grad_m.transpose()
                    + (weight * m * hess_f) * grad_d * grad_d.transpose()
                    + (weight * m * grad_f) * hess_d;
            }
        } else {
            const VectorMax12d grad_d =
                compute_distance_gradient<CT>(x, ee_dtype);
            const MatrixMax12d hess_d =
                compute_distance_hessian<CT>(x, ee_dtype);
            const double grad_f = potential_gradient(d, dmin, params);
            const double hess_f = potential_hessian(d, dmin, params);

            // ∇²[f(d(x))] = f″ ∇d ∇dᵀ + f' ∇²d
            hess = (weight * hess_f) * grad_d * grad_d.transpose()
                + (weight * grad_f) * hess_d;
        }

        // Write the weighted (unprojected) block into the top-left NDOF×NDOF of
        // a zero-padded 12×12 column-major slot (leading dim PADDED_DIM). The
        // batched eigensolver projects it and the GPU assembler scatters it.
        double* block = blocks + i * size_t(PADDED_BLOCK);
        for (int c = 0; c < NDOF; c++) {
            for (int r = 0; r < NDOF; r++) {
                block[c * PADDED_DIM + r] = hess(r, c);
            }
        }
    }

    // ========================================================================
    // Host-side launch helpers

    CollisionsView make_view(const DeviceNormalCollisions::Impl& impl)
    {
        CollisionsView view;
        view.vertex_id_0 = thrust::raw_pointer_cast(impl.vertex_id_0.data());
        view.vertex_id_1 = thrust::raw_pointer_cast(impl.vertex_id_1.data());
        view.vertex_id_2 = thrust::raw_pointer_cast(impl.vertex_id_2.data());
        view.vertex_id_3 = thrust::raw_pointer_cast(impl.vertex_id_3.data());
        view.weight = thrust::raw_pointer_cast(impl.weight.data());
        view.dmin = thrust::raw_pointer_cast(impl.dmin.data());
        view.ee_eps_x = thrust::raw_pointer_cast(impl.ee_eps_x.data());
        view.ee_dtype = thrust::raw_pointer_cast(impl.ee_dtype.data());
        view.ee_begin = impl.ee_begin();
        return view;
    }

    /// Per-type subranges of the flattened collision order (vv | ev | ee |
    /// fv), used to launch one kernel instantiation per collision type.
    struct TypeRange {
        CollisionType type;
        size_t begin;
        size_t end;
    };

    std::array<TypeRange, 4>
    type_ranges(const DeviceNormalCollisions::Impl& impl)
    {
        return { {
            { CollisionType::VV, impl.vv_begin(), impl.ev_begin() },
            { CollisionType::EV, impl.ev_begin(), impl.ee_begin() },
            { CollisionType::EE, impl.ee_begin(), impl.fv_begin() },
            { CollisionType::FV, impl.fv_begin(), impl.size() },
        } };
    }

} // namespace

// ============================================================================
// BarrierPotential

BarrierPotential::BarrierPotential(
    const double dhat,
    const double stiffness,
    const bool use_physical_barrier,
    const BarrierType barrier_type)
    : m_dhat(dhat)
    , m_stiffness(stiffness)
    , m_use_physical_barrier(use_physical_barrier)
{
    assert(dhat > 0);
    assert(stiffness > 0);
    set_barrier_type(barrier_type);
}

void BarrierPotential::set_barrier_type(const BarrierType barrier_type)
{
    if (barrier_type != BarrierType::CLAMPED_LOG
        && barrier_type != BarrierType::NORMALIZED_CLAMPED_LOG) {
        log_and_throw_error(
            "ipc::cuda::BarrierPotential only supports the (normalized) "
            "clamped log barrier!");
    }
    m_barrier_type = barrier_type;
}

// ----------------------------------------------------------------------------
// Drop-in API

double BarrierPotential::operator()(
    const NormalCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    if (collisions.empty()) {
        return 0.0;
    }
    return (*this)(
        DeviceNormalCollisions(collisions, mesh), DeviceVertices(vertices));
}

Eigen::VectorXd BarrierPotential::gradient(
    const NormalCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(vertices.size());
    }
    return gradient(
        DeviceNormalCollisions(collisions, mesh), DeviceVertices(vertices));
}

Eigen::SparseMatrix<double> BarrierPotential::hessian(
    const NormalCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(vertices.size(), vertices.size());
    }
    return hessian(
        DeviceNormalCollisions(collisions, mesh), DeviceVertices(vertices),
        project_hessian_to_psd);
}

// ----------------------------------------------------------------------------
// Fast-path API

double BarrierPotential::operator()(
    const DeviceNormalCollisions& collisions,
    const DeviceVertices& vertices) const
{
    if (collisions.empty()) {
        return 0.0;
    }

    const DeviceNormalCollisions::Impl& impl = collisions.impl();
    const CollisionsView view = make_view(impl);
    const double* positions =
        thrust::raw_pointer_cast(vertices.impl().positions.data());
    const BarrierParams params { dhat(), stiffness(), use_physical_barrier(),
                                 barrier_type() };

    thrust::device_vector<double> values(impl.size());
    double* d_values = thrust::raw_pointer_cast(values.data());

    for (const TypeRange& range : type_ranges(impl)) {
        const size_t count = range.end - range.begin;
        if (count == 0) {
            continue;
        }
        switch (range.type) {
        case CollisionType::VV:
            value_kernel<CollisionType::VV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_values);
            break;
        case CollisionType::EV:
            value_kernel<CollisionType::EV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_values);
            break;
        case CollisionType::EE:
            value_kernel<CollisionType::EE>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_values);
            break;
        case CollisionType::FV:
            value_kernel<CollisionType::FV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_values);
            break;
        }
    }
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    return thrust::reduce(values.begin(), values.end(), 0.0);
}

Eigen::VectorXd BarrierPotential::gradient(
    const DeviceNormalCollisions& collisions,
    const DeviceVertices& vertices) const
{
    const index_t n_total_vertices = vertices.num_vertices();
    const index_t ndof = 3 * n_total_vertices;

    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(ndof);
    }

    const DeviceNormalCollisions::Impl& impl = collisions.impl();
    const CollisionsView view = make_view(impl);
    const double* positions =
        thrust::raw_pointer_cast(vertices.impl().positions.data());
    const BarrierParams params { dhat(), stiffness(), use_physical_barrier(),
                                 barrier_type() };

    thrust::device_vector<double> grad(ndof, 0.0);
    double* d_grad = thrust::raw_pointer_cast(grad.data());

    for (const TypeRange& range : type_ranges(impl)) {
        const size_t count = range.end - range.begin;
        if (count == 0) {
            continue;
        }
        switch (range.type) {
        case CollisionType::VV:
            gradient_kernel<CollisionType::VV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, n_total_vertices, params, range.begin,
                    count, d_grad);
            break;
        case CollisionType::EV:
            gradient_kernel<CollisionType::EV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, n_total_vertices, params, range.begin,
                    count, d_grad);
            break;
        case CollisionType::EE:
            gradient_kernel<CollisionType::EE>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, n_total_vertices, params, range.begin,
                    count, d_grad);
            break;
        case CollisionType::FV:
            gradient_kernel<CollisionType::FV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, n_total_vertices, params, range.begin,
                    count, d_grad);
            break;
        }
    }
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    Eigen::VectorXd out(ndof);
    thrust::copy(grad.begin(), grad.end(), out.data());
    return out;
}

Eigen::SparseMatrix<double> BarrierPotential::hessian(
    const DeviceNormalCollisions& collisions,
    const DeviceVertices& vertices,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const index_t n_total_vertices = vertices.num_vertices();
    const index_t ndof = 3 * n_total_vertices;

    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(ndof, ndof);
    }

    const DeviceNormalCollisions::Impl& impl = collisions.impl();
    const CollisionsView view = make_view(impl);
    const double* positions =
        thrust::raw_pointer_cast(vertices.impl().positions.data());
    const BarrierParams params { dhat(), stiffness(), use_physical_barrier(),
                                 barrier_type() };

    const size_t n = impl.size();

    // Per-collision weighted blocks, each in the top-left of a zero-padded
    // 12×12 column-major slot (zero-init so the padding is exactly zero).
    thrust::device_vector<double> d_blocks(n * size_t(PADDED_BLOCK), 0.0);
    double* d_blocks_ptr = thrust::raw_pointer_cast(d_blocks.data());

    for (const TypeRange& range : type_ranges(impl)) {
        const size_t count = range.end - range.begin;
        if (count == 0) {
            continue;
        }
        switch (range.type) {
        case CollisionType::VV:
            hessian_kernel<CollisionType::VV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_blocks_ptr);
            break;
        case CollisionType::EV:
            hessian_kernel<CollisionType::EV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_blocks_ptr);
            break;
        case CollisionType::EE:
            hessian_kernel<CollisionType::EE>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_blocks_ptr);
            break;
        case CollisionType::FV:
            hessian_kernel<CollisionType::FV>
                <<<kernel_grid_size(count), KERNEL_BLOCK_SIZE>>>(
                    view, positions, params, range.begin, count, d_blocks_ptr);
            break;
        }
    }
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    // Project each block onto the PSD cone on device with a one-thread-per-
    // collision Eigen SelfAdjointEigenSolver. NONE is a no-op.
    if (project_hessian_to_psd != PSDProjectionMethod::NONE) {
        BatchedPSDProjection projector;
        projector.project(
            d_blocks_ptr, static_cast<int>(n), project_hessian_to_psd);
    }

    // Assemble the sparse hessian on device (Thrust sort + segmented reduce
    // into CSC arrays), then wrap as an Eigen sparse matrix.
    return assemble_sparse_hessian(
        d_blocks_ptr, view.vertex_id_0, view.vertex_id_1, view.vertex_id_2,
        view.vertex_id_3, n, n_total_vertices);
}

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
