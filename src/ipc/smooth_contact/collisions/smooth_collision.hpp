#pragma once

#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/smooth_contact/distance/primitive_distance.hpp>

namespace ipc {

enum class CollisionType : uint8_t {
    EDGE_VERTEX,
    VERTEX_VERTEX,
    FACE_VERTEX,
    EDGE_EDGE,
};

/// @brief Contact pair class for Geometric Contact Potential.
/// @note Unlike NormalCollision, SmoothCollision has to be reconstructed whenever vertices change position
class SmoothCollision {
public:
    static constexpr int ELEMENT_SIZE = 3 * MAX_VERT_3D;

    SmoothCollision(
        const index_t _primitive0,
        const index_t _primitive1,
        const double _dhat,
        const CollisionMesh& mesh)
        : primitive0(_primitive0)
        , primitive1(_primitive1)
        , m_dhat(_dhat)
    {
    }

    virtual ~SmoothCollision() = default;

    /// @brief Check if this contact pair is active (depending on both orientation and distance)
    bool is_active() const { return m_is_active; }

    /// @brief dhat value for this contact pair
    double dhat() const { return m_dhat; }

    /// @brief Name of the contact pair type
    virtual std::string name() const = 0;

    /// @brief Number of vertices involved times the dimension
    virtual int n_dofs() const = 0;

    /// @brief Contact pair type
    virtual CollisionType type() const = 0;

    /// @brief Get the number of vertices in the collision stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the vertex IDs of the collision stencil.
    /// @return The vertex IDs of the collision stencil. Size is always 4, but elements i > num_vertices() are -1.
    std::vector<index_t> vertex_ids() const { return m_vertex_ids; }

    /// @brief Get the vertex attributes of the collision stencil.
    /// @param vertices Vertex attributes
    /// @return The vertex positions of the collision stencil. Size is always 4, but elements i > num_vertices() are NaN.
    Eigen::MatrixXd vertices(Eigen::ConstRef<Eigen::MatrixXd> vertices) const
    {
        const int DIM = vertices.cols();
        Eigen::MatrixXd stencil_vertices(vertex_ids().size(), DIM);
        for (int i = 0; i < vertex_ids().size(); i++) {
            stencil_vertices.row(i) = vertices.row(vertex_ids()[i]);
        }

        return stencil_vertices;
    }

    /// @brief Select this stencil's DOF from the full matrix of DOF.
    /// @param X Full matrix of DOF (rowwise).
    /// @return This stencil's DOF.
    Eigen::VectorXd dof(Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the distance of the stencil.
    /// @param vertices Collision mesh vertices
    /// @return Squared distance of the stencil.
    virtual double
    compute_distance(Eigen::ConstRef<Eigen::MatrixXd> vertices) const = 0;

    /// @brief Compute the value of the GCP potential
    virtual double operator()(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const SmoothContactParameters& params) const = 0;

    /// @brief Compute the gradient of the GCP potential wrt. vertices involved
    virtual Vector<double, -1, ELEMENT_SIZE> gradient(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const SmoothContactParameters& params) const = 0;

    /// @brief Compute the Hessian of the GCP potential wrt. vertices involved
    virtual MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> hessian(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const SmoothContactParameters& params) const = 0;

    bool operator==(const SmoothCollision& other) const
    {
        return (
            primitive0 == other.primitive0 && primitive1 == other.primitive1);
    }

    bool operator!=(const SmoothCollision& other) const
    {
        return !(*this == other);
    }

    index_t operator[](int idx) const
    {
        if (idx == 0) {
            return primitive0;
        } else if (idx == 1) {
            return primitive1;
        } else {
            throw std::runtime_error("Invalid index in smooth_collision!");
        }
    }

    std::pair<index_t, index_t> get_hash() const
    {
        return std::make_pair(primitive0, primitive1);
    }

    double weight = 1;

protected:
    bool m_is_active = true;
    index_t primitive0, primitive1;
    double m_dhat;
    std::vector<index_t> m_vertex_ids;
};

/// @brief Templated class for various types of contact pairs
template <typename PrimitiveA, typename PrimitiveB>
class SmoothCollisionTemplate : public SmoothCollision {
public:
    using Super = SmoothCollision;
    /// @brief Distance type of the contact pair
    using DTYPE = typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type;
    /// @brief Number of points needed to compute the distance between two primitives
    static constexpr int N_CORE_POINTS =
        PrimitiveA::N_CORE_POINTS + PrimitiveB::N_CORE_POINTS;
    static constexpr int DIM = PrimitiveA::DIM;
    static constexpr int N_CORE_DOFS_A = PrimitiveA::N_CORE_POINTS * DIM;
    static constexpr int N_CORE_DOFS_B = PrimitiveB::N_CORE_POINTS * DIM;
    static constexpr int N_CORE_DOFS = N_CORE_POINTS * DIM;
    static constexpr int ELEMENT_SIZE = Super::ELEMENT_SIZE;

    SmoothCollisionTemplate(
        index_t primitive0,
        index_t primitive1,
        DTYPE dtype,
        const CollisionMesh& mesh,
        const SmoothContactParameters& params,
        const double dhat,
        const Eigen::MatrixXd& V);

    virtual ~SmoothCollisionTemplate() = default;

    std::string name() const override;

    int n_dofs() const override
    {
        return primitive_a->n_dofs() + primitive_b->n_dofs();
    }
    CollisionType type() const override;

    Vector<int, N_CORE_DOFS> get_core_indices() const;
    std::array<index_t, N_CORE_DOFS> core_vertex_ids() const;

    int num_vertices() const override
    {
        return primitive_a->n_vertices() + primitive_b->n_vertices();
    }

    template <typename T>
    Vector<T, N_CORE_DOFS> core_dof(const MatrixX<T>& X) const
    {
        return this->dof(X)(get_core_indices());
    }

    // ---- non distance type potential ----

    /// @brief Compute the GCP potential
    /// @param positions Vertex positions
    /// @param params GCP parameters
    /// @return GCP potential value
    double operator()(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const SmoothContactParameters& params) const override;

    /// @brief Compute the potential gradient wrt. positions
    /// @param positions Vertex positions
    /// @param params GCP parameters
    /// @return GCP potential gradient
    Vector<double, -1, ELEMENT_SIZE> gradient(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const SmoothContactParameters& params) const override;

    /// @brief Compute the potential Hessian wrt. positions
    /// @param positions Vertex positions
    /// @param params GCP parameters
    /// @return GCP potential Hessian
    MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> hessian(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const SmoothContactParameters& params) const override;

    // ---- distance ----

    /// @brief Compute the minimum squared distance between two primitives
    double
    compute_distance(Eigen::ConstRef<Eigen::MatrixXd> vertices) const override;

private:

    /// @brief The first primitive in the contact pair
    std::unique_ptr<PrimitiveA> primitive_a;
    /// @brief The second primitive in the contact pair
    std::unique_ptr<PrimitiveB> primitive_b;
};

} // namespace ipc