#pragma once

// Eigen ≥ 5 removed Eigen::internal::make_coherent (previously provided by
// unsupported/Eigen/AutoDiff), which MeshFEMCore's AutomaticDifferentiation.hh
// still references. This shim reimplements it with the Eigen 3.4 semantics:
// if exactly one of the two derivative vectors is empty, resize it to match
// the other and zero it.
//
// It must be included (or force-included via -include//FI) before any
// MeshFEMCore/MeshFEMSparse header. It intentionally has no includes of its
// own so it can be force-included into the MeshFEMSparse build (see
// cmake/recipes/meshfem_sparse.cmake).
//
// TODO: remove once MeshFEM supports Eigen 5 upstream.

// NOTE: EIGEN_WORLD_VERSION is 3 forever; Eigen 5 moved to semver with
// EIGEN_MAJOR_VERSION=5. If no Eigen header has been seen yet (the
// force-include case), assume Eigen 5 — this repository pins Eigen 5.
#if !defined(EIGEN_MAJOR_VERSION) || EIGEN_MAJOR_VERSION >= 5

namespace Eigen {
namespace internal {

    template <typename DerTypeA, typename DerTypeB>
    inline void make_coherent(const DerTypeA& a, const DerTypeB& b)
    {
        // Eigen 3.4's implementation const-casts too (the derivatives are
        // semantically mutable scratch space of the AutoDiffScalar pair).
        DerTypeA& a_ref = const_cast<DerTypeA&>(a); // NOLINT
        DerTypeB& b_ref = const_cast<DerTypeB&>(b); // NOLINT
        if (a_ref.size() == 0 && b_ref.size() != 0) {
            a_ref.resize(b_ref.size());
            a_ref.setZero();
        } else if (b_ref.size() == 0 && a_ref.size() != 0) {
            b_ref.resize(a_ref.size());
            b_ref.setZero();
        }
    }

} // namespace internal
} // namespace Eigen

#endif
