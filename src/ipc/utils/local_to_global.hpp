#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>
#include "MatrixCache.hpp"

namespace ipc {

template <typename DerivedLocalGrad, typename IDContainer, typename DerivedGrad>
void local_gradient_to_global_gradient(
    const Eigen::MatrixBase<DerivedLocalGrad>& local_grad,
    const IDContainer& ids,
    int dim,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        grad.segment(dim * ids[i], dim) += local_grad.segment(dim * i, dim);
    }
}

template <
    typename DerivedLocalGrad,
    typename IDContainer,
    typename Scalar = typename DerivedLocalGrad::Scalar>
void local_gradient_to_global_gradient(
    const Eigen::MatrixBase<DerivedLocalGrad>& local_grad,
    const IDContainer& ids,
    int dim,
    Eigen::SparseVector<Scalar>& grad)
{
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        for (int d = 0; d < dim; d++) {
            grad.coeffRef(dim * ids[i] + d) += local_grad(dim * i + d);
        }
    }
}

template <typename Derived, typename IDContainer>
void local_hessian_to_global_triplets(
    const Eigen::MatrixBase<Derived>& local_hessian,
    const IDContainer& ids,
    int dim,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    assert(local_hessian.rows() == local_hessian.cols());
    assert(local_hessian.rows() % dim == 0);
    const int n_verts = local_hessian.rows() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        for (int j = 0; j < n_verts; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    triplets.emplace_back(
                        dim * ids[i] + k, dim * ids[j] + l,
                        local_hessian(dim * i + k, dim * j + l));
                }
            }
        }
    }
}

class LocalThreadMatStorage
{
public:
    std::unique_ptr<MatrixCache> cache = nullptr;

    LocalThreadMatStorage() = delete;

    LocalThreadMatStorage(const int buffer_size, const int rows, const int cols)
    {
        init(buffer_size, rows, cols);
    }

    LocalThreadMatStorage(const int buffer_size, const MatrixCache &c)
    {
        init(buffer_size, c);
    }

    LocalThreadMatStorage(const LocalThreadMatStorage &other)
        : cache(other.cache->copy())
    {
    }

    LocalThreadMatStorage &operator=(const LocalThreadMatStorage &other)
    {
        assert(other.cache != nullptr);
        cache = other.cache->copy();
        return *this;
    }

    void init(const int buffer_size, const int rows, const int cols)
    {
        // assert(rows == cols);
        // cache = std::make_unique<DenseMatrixCache>();
        cache = std::make_unique<SparseMatrixCache>();
        cache->reserve(buffer_size);
        cache->init(rows, cols);
    }

    void init(const int buffer_size, const MatrixCache &c)
    {
        if (cache == nullptr)
            cache = c.copy();
        cache->reserve(buffer_size);
        cache->init(c);
    }
};

template <typename Derived, typename IDContainer>
void local_hessian_to_global_triplets(
    const Eigen::MatrixBase<Derived>& local_hessian,
    const IDContainer& ids,
    int dim,
    MatrixCache& triplets)
{
    assert(local_hessian.rows() == local_hessian.cols());
    assert(local_hessian.rows() % dim == 0);
    const int n_verts = local_hessian.rows() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        for (int j = 0; j < n_verts; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    const auto &val = local_hessian(dim * i + k, dim * j + l);
                    if (val != 0)
                        triplets.add_value(0, dim * ids[i] + k, dim * ids[j] + l, val);
                }
            }
        }
    }
}

} // namespace ipc
