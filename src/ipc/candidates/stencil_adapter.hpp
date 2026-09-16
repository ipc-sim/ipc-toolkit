#pragma once

#include <ipc/candidates/collision_stencil.hpp>

namespace ipc {

/// @brief Gives a candidate the runtime-polymorphic stencil interface.
///
/// The candidate types are deliberately not CollisionStencils: the broad phase
/// stores millions of them and a vptr would double their size. The collision
/// types are used through CollisionStencil&, so they derive from this instead,
/// which forwards every stencil operation to the candidate.
///
/// @tparam Candidate The candidate type to adapt.
template <class Candidate>
class StencilAdapter : public Candidate, public virtual CollisionStencil {
public:
    using Candidate::Candidate;

    /// @brief Adapt an existing candidate.
    StencilAdapter(const Candidate& candidate) : Candidate(candidate) { }

    // Both bases carry the mixin's operations (once through Candidate, once
    // through CollisionStencil), so every shared name has to be pinned to one
    // side or lookup is ambiguous. The polymorphic side is the one to keep.
    using CollisionStencil::compute_coefficients;
    using CollisionStencil::compute_distance;
    using CollisionStencil::compute_distance_gradient;
    using CollisionStencil::compute_distance_hessian;
    using CollisionStencil::compute_distance_vector;
    using CollisionStencil::compute_distance_vector_jacobian;
    using CollisionStencil::compute_normal;
    using CollisionStencil::compute_normal_jacobian;
    using CollisionStencil::contract_distance_vector_jacobian;
    using CollisionStencil::diag_distance_vector_outer;
    using CollisionStencil::diag_distance_vector_t_outer;
    using CollisionStencil::dim;
    using CollisionStencil::dof;
    using CollisionStencil::STENCIL_SIZE;
    using CollisionStencil::vertices;
    using CollisionStencil::write_ccd_query;

    int num_vertices() const override { return Candidate::num_vertices(); }

    std::array<index_t, CollisionStencil::STENCIL_SIZE> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return Candidate::vertex_ids(edges, faces);
    }

    double
    compute_distance(Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_distance(positions);
    }

    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_distance_gradient(positions);
    }

    MatrixMax12d compute_distance_hessian(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_distance_hessian(positions);
    }

    VectorMax4d
    compute_coefficients(Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_coefficients(positions);
    }

    bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const override
    {
        return Candidate::ccd(
            vertices_t0, vertices_t1, toi, min_distance, tmax,
            narrow_phase_ccd);
    }

protected:
    VectorMax3d compute_unnormalized_normal(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_unnormalized_normal(positions);
    }

    MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_unnormalized_normal_jacobian(positions);
    }
};

/// @brief A StencilAdapter for a candidate whose distance type can be known a priori.
///
/// The candidate's math takes the distance type as an argument; this supplies
/// it through a virtual so a collision can pin it to the value it cached when
/// the collision set was built, which is what keeps the classification fixed
/// across a solve.
///
/// @tparam Candidate The candidate type to adapt.
/// @tparam DType The candidate's distance type enum.
template <class Candidate, class DType>
class DTypeStencilAdapter : public Candidate, public virtual CollisionStencil {
public:
    using Candidate::Candidate;

    /// @brief Adapt an existing candidate.
    DTypeStencilAdapter(const Candidate& candidate) : Candidate(candidate) { }

    // Both bases carry the mixin's operations (once through Candidate, once
    // through CollisionStencil), so every shared name has to be pinned to one
    // side or lookup is ambiguous. The polymorphic side is the one to keep.
    using CollisionStencil::compute_coefficients;
    using CollisionStencil::compute_distance;
    using CollisionStencil::compute_distance_gradient;
    using CollisionStencil::compute_distance_hessian;
    using CollisionStencil::compute_distance_vector;
    using CollisionStencil::compute_distance_vector_jacobian;
    using CollisionStencil::compute_normal;
    using CollisionStencil::compute_normal_jacobian;
    using CollisionStencil::contract_distance_vector_jacobian;
    using CollisionStencil::diag_distance_vector_outer;
    using CollisionStencil::diag_distance_vector_t_outer;
    using CollisionStencil::dim;
    using CollisionStencil::dof;
    using CollisionStencil::STENCIL_SIZE;
    using CollisionStencil::vertices;
    using CollisionStencil::write_ccd_query;

    /// @brief The distance type known a priori; AUTO unless a collision pins it.
    virtual DType known_dtype() const { return DType::AUTO; }

    int num_vertices() const override { return Candidate::num_vertices(); }

    std::array<index_t, CollisionStencil::STENCIL_SIZE> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const override
    {
        return Candidate::vertex_ids(edges, faces);
    }

    double
    compute_distance(Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_distance(positions, known_dtype());
    }

    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_distance_gradient(positions, known_dtype());
    }

    MatrixMax12d compute_distance_hessian(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_distance_hessian(positions, known_dtype());
    }

    VectorMax4d
    compute_coefficients(Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_coefficients(positions, known_dtype());
    }

    bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const override
    {
        return Candidate::ccd(
            vertices_t0, vertices_t1, toi, min_distance, tmax,
            narrow_phase_ccd);
    }

protected:
    VectorMax3d compute_unnormalized_normal(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_unnormalized_normal(positions);
    }

    MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const override
    {
        return Candidate::compute_unnormalized_normal_jacobian(positions);
    }
};

} // namespace ipc
