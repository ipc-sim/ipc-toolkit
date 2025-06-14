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

template <int max_vert>
class SmoothCollision : public CollisionStencil<max_vert> {
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
    constexpr static int max_size = max_vert * 3;
    virtual ~SmoothCollision() = default;

    bool is_active() const { return is_active_; }

    virtual int n_dofs() const = 0;
    virtual CollisionType type() const = 0;

    std::array<index_t, max_vert> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return vertices;
    }

    // ---- non distance type potential ----

    virtual double operator()(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions,
        const ParameterType& params) const = 0;

    virtual Vector<double, -1, max_size> gradient(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions,
        const ParameterType& params) const = 0;

    virtual Eigen::MatrixXd hessian(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions,
        const ParameterType& params) const = 0;

    // ---- distance ----

    Vector<double, -1, max_size> compute_distance_gradient(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions)
        const override
    {
        return Vector<double, -1, max_size>::Zero(max_size);
    }

    MatrixMax<double, max_size, max_size> compute_distance_hessian(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions)
        const override
    {
        return MatrixMax<double, max_size, max_size>::Zero(max_size, max_size);
    }

    VectorMax4d compute_coefficients(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions)
        const override
    {
        VectorMax4d coeffs(4);
        coeffs << 0, 0, 0, 0;
        return coeffs;
    }

    bool
    ccd(Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> vertices_t0,
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const override
    {
        return false;
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
    const index_t& operator[](int idx) const
    {
        if (idx == 0)
            return primitive0;
        else if (idx == 1)
            return primitive1;
        else
            throw std::runtime_error("Invalid index in smooth_collision!");
    }

    std::pair<index_t, index_t> get_hash() const
    {
        return std::make_pair(primitive0, primitive1);
    }

    double get_dhat() const { return dhat_; }

    virtual std::string name() const { return ""; }

    double weight = 1;

protected:
    bool is_active_ = true;
    index_t primitive0, primitive1;
    double dhat_;
    std::array<index_t, max_vert> vertices;
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
    constexpr static int max_size = SmoothCollision<max_vert>::max_size;

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
    std::array<index_t, n_core_dofs> core_vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const;

    inline int num_vertices() const override
    {
        return pA->n_vertices() + pB->n_vertices();
    }

    template <typename T>
    Vector<T, n_core_dofs> core_dof(
        const MatrixX<T>& X,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const
    {
        return this->dof(X, edges, faces)(get_core_indices());
    }

    // ---- non distance type potential ----

    double operator()(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions,
        const ParameterType& params) const override;

    Vector<double, -1, max_size> gradient(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions,
        const ParameterType& params) const override;

    Eigen::MatrixXd hessian(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions,
        const ParameterType& params) const override;

    // ---- distance ----

    double compute_distance(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions)
        const override;

    Vector<double, -1, max_size> compute_distance_gradient(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions)
        const override;

    MatrixMax<double, max_size, max_size> compute_distance_hessian(
        Eigen::ConstRef<Vector<double, -1, 3 * max_vert>> positions)
        const override;

private:
    std::unique_ptr<PrimitiveA> pA;
    std::unique_ptr<PrimitiveB> pB;
};

} // namespace ipc