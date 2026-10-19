#pragma once
#include "esp_collision.hpp"
#include "esp_primitives.hpp"

#include <ipc/barrier/barrier.hpp>

#include <array>

namespace ipc {

/// @brief Templated class for various types of contact pairs
template <typename PrimitiveA, typename PrimitiveB>
class ESPCollisionTemplate : public ESPCollision {
public:
    using Super = ESPCollision;
    static constexpr int N_CORE_POINTS =
        PrimitiveA::N_CORE_POINTS + PrimitiveB::N_CORE_POINTS;
    static constexpr int DIM = PrimitiveA::DIM;
    static constexpr int N_CORE_DOFS_A = PrimitiveA::N_CORE_POINTS * DIM;
    static constexpr int N_CORE_DOFS_B = PrimitiveB::N_CORE_POINTS * DIM;
    static constexpr int N_CORE_DOFS = N_CORE_POINTS * DIM;
    static constexpr int ELEMENT_SIZE = Super::ELEMENT_SIZE;
    static_assert(
        PrimitiveA::N_DOFS + PrimitiveB::N_DOFS <= ELEMENT_SIZE,
        "ESPCollision::ELEMENT_SIZE is too small for this stencil");

    ESPCollisionTemplate(
        index_t primitive0, index_t primitive1, const CollisionMesh& mesh);

    virtual ~ESPCollisionTemplate() = default;

    std::string name() const override;

    int n_dofs() const override
    {
        return m_primitive_a.n_dofs() + m_primitive_b.n_dofs();
    }
    ESPCollisionType type() const override;

    std::pair<index_t, index_t> get_hash() const override
    {
        return std::make_pair(m_primitive_a.id(), m_primitive_b.id());
    }

    std::array<index_t, 3> get_typed_hash() const override
    {
        return { { static_cast<index_t>(type()), m_primitive_a.id(),
                   m_primitive_b.id() } };
    }

    index_t operator[](int idx) const override
    {
        if (idx == 0) {
            return m_primitive_a.id();
        } else if (idx == 1) {
            return m_primitive_b.id();
        } else {
            throw std::runtime_error("Invalid index in ESP collision!");
        }
    }

    int num_vertices() const override
    {
        return m_primitive_a.n_vertices() + m_primitive_b.n_vertices();
    }

    index_t vertex_id(index_t i) const override;

    size_t n_vertices_a() const override { return m_primitive_a.n_vertices(); }
    size_t n_vertices_b() const override { return m_primitive_b.n_vertices(); }

    double operator()(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params) const override;

    VectorMax<double, ELEMENT_SIZE> gradient(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params) const override;

    MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> hessian(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params) const override;

    double
    compute_distance(Eigen::ConstRef<Eigen::MatrixXd> vertices) const override;

    // The near/far functions are only declared here. Their generic stubs and
    // the 3D specializations are both defined in esp_collision_template.cpp:
    // a generic definition in the class body can be instantiated (MSVC does)
    // before the specializations are seen, putting the zero stub in the
    // vtable.

    std::pair<double, double> operator_nearfar(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const NearFarBarrier* nf_barrier) const override;

    std::pair<VectorMax<double, ELEMENT_SIZE>, VectorMax<double, ELEMENT_SIZE>>
    gradient_nearfar(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const NearFarBarrier* nf_barrier) const override;

    std::pair<
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>,
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>>
    hessian_nearfar(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const NearFarBarrier* nf_barrier) const override;

private:
    PrimitiveA m_primitive_a;
    PrimitiveB m_primitive_b;
};

// Keep old name as alias for backward compatibility within this codebase
template <typename PrimitiveA, typename PrimitiveB>
using ESPCollision3DTemplate = ESPCollisionTemplate<PrimitiveA, PrimitiveB>;

// 2D alias (for use with 2D primitives)
template <typename PrimitiveA, typename PrimitiveB>
using ESPCollision2DTemplate = ESPCollisionTemplate<PrimitiveA, PrimitiveB>;

} // namespace ipc
