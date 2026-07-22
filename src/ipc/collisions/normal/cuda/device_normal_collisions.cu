#include "device_normal_collisions.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include "device_normal_collisions_impl.cuh"

#include <ipc/utils/cuda/device_utils.cuh>
#include <ipc/utils/logger.hpp>

#include <thrust/copy.h>

namespace ipc::cuda {

// ============================================================================
// DeviceVertices

DeviceVertices::DeviceVertices(Eigen::ConstRef<Eigen::MatrixXd> vertices)
    : m_impl(std::make_unique<Impl>())
{
    update(vertices);
}

DeviceVertices::~DeviceVertices() = default;
DeviceVertices::DeviceVertices(DeviceVertices&&) noexcept = default;
DeviceVertices& DeviceVertices::operator=(DeviceVertices&&) noexcept = default;

const DeviceVertices::Impl& DeviceVertices::impl() const { return *m_impl; }

void DeviceVertices::update(Eigen::ConstRef<Eigen::MatrixXd> vertices)
{
    assert(vertices.cols() == 3);
    assert(m_num_vertices == 0 || vertices.rows() == m_num_vertices);
    m_num_vertices = static_cast<index_t>(vertices.rows());
    // Repack into a row-major buffer so vertex i occupies positions[3i..3i+2].
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> row_major =
        vertices;
    m_impl->positions.resize(row_major.size());
    thrust::copy(
        row_major.data(), row_major.data() + row_major.size(),
        m_impl->positions.begin());
}

// ============================================================================
// DeviceNormalCollisions

namespace {

    /// Gather the per-collision fields of one typed collision vector into the
    /// host staging buffers.
    template <typename Collision>
    void gather_collisions(
        const std::vector<Collision>& collisions,
        const CollisionMesh& mesh,
        std::vector<std::array<index_t, 4>>& vertex_ids,
        std::vector<double>& weight,
        std::vector<double>& dmin)
    {
        for (const Collision& collision : collisions) {
            vertex_ids.push_back(
                collision.vertex_ids(mesh.edges(), mesh.faces()));
            weight.push_back(collision.weight);
            dmin.push_back(collision.dmin);
        }
    }

} // namespace

DeviceNormalCollisions::DeviceNormalCollisions(
    const NormalCollisions& collisions, const CollisionMesh& mesh)
    : m_impl(std::make_unique<Impl>())
{
    assert(collisions.pv_collisions.empty());
    assert(mesh.dim() == 3);

    Impl& impl = *m_impl;
    impl.n_vv = collisions.vv_collisions.size();
    impl.n_ev = collisions.ev_collisions.size();
    impl.n_ee = collisions.ee_collisions.size();
    impl.n_fv = collisions.fv_collisions.size();
    const size_t n = impl.size();

    // Stage on the host in the flattened order vv | ev | ee | fv.
    std::vector<std::array<index_t, 4>> vertex_ids;
    vertex_ids.reserve(n);
    std::vector<double> weight, dmin;
    weight.reserve(n);
    dmin.reserve(n);

    gather_collisions(collisions.vv_collisions, mesh, vertex_ids, weight, dmin);
    gather_collisions(collisions.ev_collisions, mesh, vertex_ids, weight, dmin);
    gather_collisions(collisions.ee_collisions, mesh, vertex_ids, weight, dmin);
    gather_collisions(collisions.fv_collisions, mesh, vertex_ids, weight, dmin);

    std::vector<double> ee_eps_x;
    std::vector<uint8_t> ee_dtype;
    ee_eps_x.reserve(impl.n_ee);
    ee_dtype.reserve(impl.n_ee);
    for (const EdgeEdgeNormalCollision& ee : collisions.ee_collisions) {
        ee_eps_x.push_back(ee.eps_x);
        ee_dtype.push_back(static_cast<uint8_t>(ee.dtype));
    }

    // Transpose the vertex ids into SoA staging buffers.
    std::array<std::vector<index_t>, 4> ids_soa;
    for (auto& ids : ids_soa) {
        ids.resize(n);
    }
    for (size_t i = 0; i < n; i++) {
        for (int k = 0; k < 4; k++) {
            ids_soa[k][i] = vertex_ids[i][k];
        }
    }

    // Upload to the device.
    impl.vertex_id_0 = ids_soa[0];
    impl.vertex_id_1 = ids_soa[1];
    impl.vertex_id_2 = ids_soa[2];
    impl.vertex_id_3 = ids_soa[3];
    impl.weight = weight;
    impl.dmin = dmin;
    impl.ee_eps_x = ee_eps_x;
    impl.ee_dtype = ee_dtype;
}

DeviceNormalCollisions::~DeviceNormalCollisions() = default;
DeviceNormalCollisions::DeviceNormalCollisions(
    DeviceNormalCollisions&&) noexcept = default;
DeviceNormalCollisions&
DeviceNormalCollisions::operator=(DeviceNormalCollisions&&) noexcept = default;

const DeviceNormalCollisions::Impl& DeviceNormalCollisions::impl() const
{
    return *m_impl;
}

size_t DeviceNormalCollisions::size() const { return m_impl->size(); }

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
