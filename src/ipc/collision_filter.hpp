#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <functional>
#include <utility>

namespace ipc {

// ─────────────────────────────────────────────────────────────────────────────
// CollisionFilter
//
// A type-erased, composable predicate: operator()(size_t vi, size_t vj) -> bool
//
//   CollisionFilter f = make_vertex_patches_filter(patches);
//   CollisionFilter g = make_static_obstacle_filter(n_dynamic);
//   CollisionFilter h = f & g;  // both conditions must hold
//   if (h(i, j)) { ... }
//
// Filters are value types — cheap to copy (shared_ptr to immutable impl).
//   operator|  → union        (true if EITHER filter passes)
//   operator&  → intersection (true if BOTH filters pass)
//   operator!  → negation
// ─────────────────────────────────────────────────────────────────────────────

class CollisionFilter {
public:
    // ── Construction ─────────────────────────────────────────────────────────

    /// @brief Default filter: accept all pairs.
    CollisionFilter() : m_fn([](int, int) { return true; }) { }

    /// @brief Construct from any callable bool(size_t, size_t).
    /// @note Disabled when Fn is CollisionFilter itself to avoid shadowing
    ///       the copy constructor.
    template <
        typename Fn,
        typename = std::enable_if_t<
            std::is_invocable_r_v<bool, Fn, size_t, size_t>
            && !std::is_same_v<std::decay_t<Fn>, CollisionFilter>>>
    CollisionFilter(Fn&& fn) : m_fn(std::forward<Fn>(fn))
    {
    }

    // ── Call operator ────────────────────────────────────────────────────────

    /// @brief Test whether two vertices may collide.
    /// @param vi Index of the first vertex.
    /// @param vj Index of the second vertex.
    /// @return true if the pair should be considered for collision.
    bool operator()(size_t vi, size_t vj) const { return m_fn(vi, vj); }

    // ── Composition ──────────────────────────────────────────────────────────

    /// @brief Union: accept if EITHER filter passes.
    friend CollisionFilter operator|(CollisionFilter lhs, CollisionFilter rhs)
    {
        return CollisionFilter(
            [l = std::move(lhs.m_fn), r = std::move(rhs.m_fn)](int vi, int vj) {
                return l(vi, vj) || r(vi, vj);
            });
    }

    /// @brief Intersection: accept only if BOTH filters pass.
    friend CollisionFilter operator&(CollisionFilter lhs, CollisionFilter rhs)
    {
        return CollisionFilter(
            [l = std::move(lhs.m_fn), r = std::move(rhs.m_fn)](int vi, int vj) {
                return l(vi, vj) && r(vi, vj);
            });
    }

    /// @brief Negation: accept only if this filter rejects.
    CollisionFilter operator!() const
    {
        return CollisionFilter(
            [f = m_fn](int vi, int vj) { return !f(vi, vj); });
    }

    /// @brief Compound union assignment.
    CollisionFilter& operator|=(CollisionFilter rhs)
    {
        return *this = *this | std::move(rhs);
    }

    /// @brief Compound intersection assignment.
    CollisionFilter& operator&=(CollisionFilter rhs)
    {
        return *this = *this & std::move(rhs);
    }

private:
    std::function<bool(size_t, size_t)> m_fn;
};

// ─────────────────────────────────────────────────────────────────────────────
// Inline factory functions
// ─────────────────────────────────────────────────────────────────────────────

/// @brief Create a filter that only allows collisions between vertices in
///        different patches (e.g., different garment panels or bodies).
/// @param patch_ids Per-vertex patch label vector (one entry per vertex).
/// @return A CollisionFilter that blocks same-patch pairs.
inline CollisionFilter make_vertex_patches_filter(Eigen::VectorXi patch_ids)
{
    return CollisionFilter([ids = std::move(patch_ids)](size_t vi, size_t vj) {
        return ids[vi] != ids[vj];
    });
}

/// @brief Create a filter that prevents static obstacles from colliding with
///        each other. A vertex is considered "static" if its index is
///        >= n_dynamic. Pairs where both vertices are static are rejected.
/// @param n_dynamic Number of dynamic (simulated) vertices; static vertices
///        occupy indices [n_dynamic, n_verts).
/// @return A CollisionFilter that blocks static-static pairs.
inline CollisionFilter make_static_obstacle_filter(size_t n_dynamic)
{
    return CollisionFilter([n_dynamic](size_t vi, size_t vj) {
        return vi < n_dynamic || vj < n_dynamic;
    });
}

/// @brief Create a filter that prevents self-collisions within a connected
///        component of the face mesh. Two vertices in the same connected
///        component are blocked; cross-component pairs are allowed.
/// @param faces Face index matrix (#F × 3).
/// @return A CollisionFilter that blocks intra-component pairs.
/// @note Implemented in collision_filter.cpp (requires libigl internally).
CollisionFilter
make_connected_components_filter(Eigen::ConstRef<Eigen::MatrixXi> faces);

} // namespace ipc
