#pragma once
#include "esp_collision.hpp"
#include "esp_primitives.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/utils/logger.hpp>

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
            log_and_throw_error("Invalid index in ESP collision!");
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
        const NearFarBarrier<>* nf_barrier) const override;

    std::pair<VectorMax<double, ELEMENT_SIZE>, VectorMax<double, ELEMENT_SIZE>>
    gradient_nearfar(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const NearFarBarrier<>* nf_barrier) const override;

    std::pair<
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>,
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>>
    hessian_nearfar(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const NearFarBarrier<>* nf_barrier) const override;

private:
    PrimitiveA m_primitive_a;
    PrimitiveB m_primitive_b;
};

// We define the members below as explicit specializations in
// esp_collision_template.cpp, and every translation unit that uses one must
// see its declaration first; otherwise the program is ill-formed (no
// diagnostic required). In particular, the vertex-vertex constructor sorts its
// two ids, which ArbitraryPointESP relies on.

template <>
ESPCollisionTemplate<Vertex3, Vertex3>::ESPCollisionTemplate(
    index_t primitive0, index_t primitive1, const CollisionMesh& mesh);

#define IPC_DECLARE_ESP_COLLISION_SPECIALIZATIONS(A, B)                        \
    template <> std::string ESPCollisionTemplate<A, B>::name() const;          \
    template <> ESPCollisionType ESPCollisionTemplate<A, B>::type() const;     \
    template <>                                                                \
    double ESPCollisionTemplate<A, B>::compute_distance(                       \
        Eigen::ConstRef<Eigen::MatrixXd> vertices) const;                      \
    template <>                                                                \
    double ESPCollisionTemplate<A, B>::operator()(                             \
        Eigen::ConstRef<VectorMax<double, ESPCollision::ELEMENT_SIZE>>         \
            positions,                                                         \
        const ESPParameters& params) const;                                    \
    template <>                                                                \
    VectorMax<double, ESPCollision::ELEMENT_SIZE>                              \
    ESPCollisionTemplate<A, B>::gradient(                                      \
        Eigen::ConstRef<VectorMax<double, ESPCollision::ELEMENT_SIZE>>         \
            positions,                                                         \
        const ESPParameters& params) const;                                    \
    template <>                                                                \
    MatrixMax<double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>  \
    ESPCollisionTemplate<A, B>::hessian(                                       \
        Eigen::ConstRef<VectorMax<double, ESPCollision::ELEMENT_SIZE>>         \
            positions,                                                         \
        const ESPParameters& params) const

#define IPC_DECLARE_ESP_COLLISION_NEARFAR_SPECIALIZATIONS(A, B)                \
    template <>                                                                \
    std::pair<double, double> ESPCollisionTemplate<A, B>::operator_nearfar(    \
        Eigen::ConstRef<VectorMax<double, ESPCollision::ELEMENT_SIZE>>         \
            positions,                                                         \
        const ESPParameters& params, const NearFarBarrier<>* nf_barrier)       \
        const;                                                                 \
    template <>                                                                \
    std::pair<                                                                 \
        VectorMax<double, ESPCollision::ELEMENT_SIZE>,                         \
        VectorMax<double, ESPCollision::ELEMENT_SIZE>>                         \
    ESPCollisionTemplate<A, B>::gradient_nearfar(                              \
        Eigen::ConstRef<VectorMax<double, ESPCollision::ELEMENT_SIZE>>         \
            positions,                                                         \
        const ESPParameters& params, const NearFarBarrier<>* nf_barrier)       \
        const;                                                                 \
    template <>                                                                \
    std::pair<                                                                 \
        MatrixMax<                                                             \
            double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>,   \
        MatrixMax<                                                             \
            double, ESPCollision::ELEMENT_SIZE, ESPCollision::ELEMENT_SIZE>>   \
    ESPCollisionTemplate<A, B>::hessian_nearfar(                               \
        Eigen::ConstRef<VectorMax<double, ESPCollision::ELEMENT_SIZE>>         \
            positions,                                                         \
        const ESPParameters& params, const NearFarBarrier<>* nf_barrier) const

IPC_DECLARE_ESP_COLLISION_SPECIALIZATIONS(Vertex3, Vertex3);
IPC_DECLARE_ESP_COLLISION_SPECIALIZATIONS(Edge3P1, Vertex3);
IPC_DECLARE_ESP_COLLISION_SPECIALIZATIONS(Face3P1, Vertex3);
IPC_DECLARE_ESP_COLLISION_SPECIALIZATIONS(Vertex2, Vertex2);
IPC_DECLARE_ESP_COLLISION_SPECIALIZATIONS(Vertex2, Edge2P1);

IPC_DECLARE_ESP_COLLISION_NEARFAR_SPECIALIZATIONS(Vertex3, Vertex3);
IPC_DECLARE_ESP_COLLISION_NEARFAR_SPECIALIZATIONS(Edge3P1, Vertex3);
IPC_DECLARE_ESP_COLLISION_NEARFAR_SPECIALIZATIONS(Face3P1, Vertex3);

#undef IPC_DECLARE_ESP_COLLISION_SPECIALIZATIONS
#undef IPC_DECLARE_ESP_COLLISION_NEARFAR_SPECIALIZATIONS

// Keep old name as alias for backward compatibility within this codebase
template <typename PrimitiveA, typename PrimitiveB>
using ESPCollision3DTemplate = ESPCollisionTemplate<PrimitiveA, PrimitiveB>;

// 2D alias (for use with 2D primitives)
template <typename PrimitiveA, typename PrimitiveB>
using ESPCollision2DTemplate = ESPCollisionTemplate<PrimitiveA, PrimitiveB>;

} // namespace ipc
