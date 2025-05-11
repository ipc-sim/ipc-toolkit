#pragma once

#include <ipc/collisions/collision.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/smooth_contact/distance/primitive_distance.hpp>

namespace ipc {

enum class CollisionType
{
    EdgeVertex,
    VertexVertex,
    FaceVertex,
    EdgeEdge,
};

template <int max_vert> class SmoothCollision : public Collision<max_vert> {
protected:
    SmoothCollision(
        long primitive0_,
        long primitive1_,
        const double& dhat,
        const CollisionMesh& mesh)
        : primitive0(primitive0_)
        , primitive1(primitive1_)
        , dhat_(dhat)
    {
        vertices.fill(-1);
    }

public:
    constexpr static int max_size = Collision<max_vert>::max_size;
    virtual ~SmoothCollision() = default;

    bool is_active() const { return is_active_; }

    virtual int n_dofs() const = 0;
    virtual CollisionType type() const = 0;

    std::array<long, max_vert> vertex_ids(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const override
    {
        return vertices;
    }

    // ---- non distance type potential ----

    virtual double operator()(
        const Vector<double, -1, max_size>& positions,
        const ParameterType& params) const
    {
        return 0.;
    }

    virtual Vector<double, -1, max_size> gradient(
        const Vector<double, -1, max_size>& positions,
        const ParameterType& params) const
    {
        return Vector<double, -1, max_size>::Zero(positions.size());
    }

    virtual MatrixMax<double, max_size, max_size> hessian(
        const Vector<double, -1, max_size>& positions,
        const ParameterType& params) const
    {
        return MatrixMax<double, max_size, max_size>::Zero(
            positions.size(), positions.size());
    }

    // ---- distance ----

    Vector<double, -1, max_size> compute_distance_gradient(
        const Vector<double, -1, max_size>& positions) const override
    {
        return Vector<double, -1, max_size>::Zero(max_size);
    }

    MatrixMax<double, max_size, max_size> compute_distance_hessian(
        const Vector<double, -1, max_size>& positions) const override
    {
        return MatrixMax<double, max_size, max_size>::Zero(max_size, max_size);
    }

    bool operator==(const SmoothCollision& other) const
    {
        return (
            primitive0 == other.primitive0 && primitive1 == other.primitive1);
    }
    bool operator!=(const SmoothCollision& other) const
    {
        return !(*this == other);
    }
    const long& operator[](int idx) const
    {
        if (idx == 0)
            return primitive0;
        else if (idx == 1)
            return primitive1;
        else
            throw std::runtime_error("Invalid index in smooth_collision!");
    }

    std::pair<long, long> get_hash() const
    {
        return std::make_pair(primitive0, primitive1);
    }

    double get_dhat() const { return dhat_; }

    virtual std::string name() const { return ""; }

protected:
    bool is_active_ = true;
    long primitive0, primitive1;
    double dhat_;
    std::array<long, max_vert> vertices;
};

template <int max_vert, typename PrimitiveA, typename PrimitiveB>
class SmoothCollisionTemplate : public SmoothCollision<max_vert> {
public:
    using Super = SmoothCollision<max_vert>;
    using DTYPE = typename PrimitiveDistType<PrimitiveA, PrimitiveB>::type;
    constexpr static int n_core_points =
        PrimitiveA::n_core_points + PrimitiveB::n_core_points;
    constexpr static int dim = PrimitiveA::dim;
    constexpr static int n_core_dofs_A = PrimitiveA::n_core_points * dim;
    constexpr static int n_core_dofs_B = PrimitiveB::n_core_points * dim;
    constexpr static int n_core_dofs = n_core_points * dim;
    constexpr static int max_size = Collision<max_vert>::max_size;

    SmoothCollisionTemplate(
        long primitive0_,
        long primitive1_,
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
    std::array<long, n_core_dofs> core_vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const;

    inline int num_vertices() const override { return pA->n_vertices() + pB->n_vertices(); }

    template <typename T>
    Vector<T, n_core_dofs>
    core_dof(const MatrixX<T>& X,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const
    {
        return this->dof(X, edges, faces)(get_core_indices());
    }

    // ---- non distance type potential ----

    double operator()(
        const Vector<double, -1, max_size>& positions,
        const ParameterType& params) const override;

    Vector<double, -1, max_size> gradient(
        const Vector<double, -1, max_size>& positions,
        const ParameterType& params) const override;

    MatrixMax<double, max_size, max_size> hessian(
        const Vector<double, -1, max_size>& positions,
        const ParameterType& params) const override;

    // ---- distance ----

    double compute_distance(
        const Vector<double, -1, max_size>& positions) const override;

    Vector<double, -1, max_size> compute_distance_gradient(
        const Vector<double, -1, max_size>& positions) const override;

    MatrixMax<double, max_size, max_size> compute_distance_hessian(
        const Vector<double, -1, max_size>& positions) const override;

private:
    std::unique_ptr<PrimitiveA> pA;
    std::unique_ptr<PrimitiveB> pB;
};

} // namespace ipc