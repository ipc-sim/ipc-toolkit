#pragma once

#include <ipc/collisions/normal/normal_collision.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/smooth_contact/distance/primitive_distance.hpp>

namespace ipc {

enum class CollisionType {
    EdgeVertex,
    VertexVertex,
    FaceVertex,
    EdgeEdge,
};

class SmoothCollision {
public:
    constexpr static int ELEMENT_SIZE = 3 * MAX_VERT_3D;

    SmoothCollision(
        long primitive0_,
        long primitive1_,
        const double& dhat,
        const CollisionMesh& mesh)
        : primitive0(primitive0_)
        , primitive1(primitive1_)
        , dhat_(dhat)
    {
    }

    virtual ~SmoothCollision() = default;

    bool is_active() const { return is_active_; }

    double dhat() const { return dhat_; }

    virtual std::string name() const = 0;

    virtual int n_dofs() const = 0;

    virtual CollisionType type() const = 0;

    /// @brief Get the number of vertices in the collision stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the vertex IDs of the collision stencil.
    /// @return The vertex IDs of the collision stencil. Size is always 4, but elements i > num_vertices() are -1.
    std::vector<index_t> vertex_ids() const { return vertex_ids_; }

    /// @brief Get the vertex attributes of the collision stencil.
    /// @tparam T Type of the attributes
    /// @param vertices Vertex attributes
    /// @return The vertex positions of the collision stencil. Size is always 4, but elements i > num_vertices() are NaN.
    Eigen::MatrixXd vertices(Eigen::ConstRef<Eigen::MatrixXd> vertices) const
    {
        const int dim = vertices.cols();
        Eigen::MatrixXd stencil_vertices(vertex_ids_.size(), dim);
        for (int i = 0; i < vertex_ids_.size(); i++) {
            stencil_vertices.row(i) = vertices.row(vertex_ids_[i]);
        }

        return stencil_vertices;
    }

    /// @brief Select this stencil's DOF from the full matrix of DOF.
    /// @tparam T Type of the DOF
    /// @param X Full matrix of DOF (rowwise).
    /// @return This stencil's DOF.
    Eigen::VectorXd dof(Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Compute the distance of the stencil.
    /// @param vertices Collision mesh vertices
    /// @return Distance of the stencil.
    virtual double
    compute_distance(Eigen::ConstRef<Eigen::MatrixXd> vertices) const = 0;

    virtual double operator()(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const ParameterType& params) const = 0;

    virtual Vector<double, -1, ELEMENT_SIZE> gradient(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const ParameterType& params) const = 0;

    virtual MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> hessian(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const ParameterType& params) const = 0;

    bool operator==(const SmoothCollision& other) const
    {
        return (
            primitive0 == other.primitive0 && primitive1 == other.primitive1);
    }
    bool operator!=(const SmoothCollision& other) const
    {
        return !(*this == other);
    }
    const index_t& operator[](int idx) const
    {
        if (idx == 0) {
            return primitive0;
        }
        else if (idx == 1) {
            return primitive1;
        }
        else {
            throw std::runtime_error("Invalid index in smooth_collision!");
        }
    }

    std::pair<index_t, index_t> get_hash() const
    {
        return std::make_pair(primitive0, primitive1);
    }

    double weight = 1;

protected:
    bool is_active_ = true;
    index_t primitive0, primitive1;
    double dhat_;
    std::vector<index_t> vertex_ids_;
};

template <typename PrimitiveA, typename PrimitiveB>
class SmoothCollisionTemplate : public SmoothCollision {
public:
    using Super = SmoothCollision;
    using DTYPE = typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type;
    constexpr static int n_core_points =
        PrimitiveA::n_core_points + PrimitiveB::n_core_points;
    constexpr static int dim = PrimitiveA::dim;
    constexpr static int n_core_dofs_A = PrimitiveA::n_core_points * dim;
    constexpr static int n_core_dofs_B = PrimitiveB::n_core_points * dim;
    constexpr static int n_core_dofs = n_core_points * dim;
    constexpr static int ELEMENT_SIZE = Super::ELEMENT_SIZE;

    SmoothCollisionTemplate(
        index_t primitive0_,
        index_t primitive1_,
        DTYPE dtype,
        const CollisionMesh& mesh,
        const ParameterType& param,
        const double& dhat,
        const Eigen::MatrixXd& V);
    virtual ~SmoothCollisionTemplate() = default;

    std::string name() const override;

    inline int n_dofs() const override { return pA->n_dofs() + pB->n_dofs(); }
    CollisionType type() const override;

    Vector<int, n_core_dofs> get_core_indices() const;
    std::array<index_t, n_core_dofs> core_vertex_ids() const;

    inline int num_vertices() const override
    {
        return pA->n_vertices() + pB->n_vertices();
    }

    template <typename T>
    Vector<T, n_core_dofs> core_dof(const MatrixX<T>& X) const
    {
        return this->dof(X)(get_core_indices());
    }

    // ---- non distance type potential ----

    double operator()(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const ParameterType& params) const override;

    Vector<double, -1, ELEMENT_SIZE> gradient(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const ParameterType& params) const override;

    MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> hessian(
        Eigen::ConstRef<Vector<double, -1, ELEMENT_SIZE>> positions,
        const ParameterType& params) const override;

    // ---- distance ----

    double
    compute_distance(Eigen::ConstRef<Eigen::MatrixXd> vertices) const override;

private:
    std::unique_ptr<PrimitiveA> pA;
    std::unique_ptr<PrimitiveB> pB;
};

} // namespace ipc