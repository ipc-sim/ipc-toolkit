#pragma once

#include <ipc/utils/matrix_cache.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>

namespace ipc {

template <typename IDContainer>
void local_gradient_to_global_gradient(
    Eigen::ConstRef<Eigen::VectorXd> local_grad,
    const IDContainer& ids,
    const int dim,
    Eigen::Ref<Eigen::VectorXd> grad)
{
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        grad.segment(dim * ids[i], dim) += local_grad.segment(dim * i, dim);
    }
}

template <typename IDContainer>
void local_gradient_to_global_gradient(
    Eigen::ConstRef<Eigen::VectorXd> local_grad,
    const IDContainer& ids,
    const int dim,
    Eigen::SparseVector<double>& grad)
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

template <typename IDContainer>
void local_hessian_to_global_triplets(
    Eigen::ConstRef<Eigen::MatrixXd> local_hessian,
    const IDContainer& ids,
    const int dim,
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

template <typename Derived, typename IDContainer1, typename IDContainer2>
void local_jacobian_to_global_triplets(
    const Eigen::MatrixBase<Derived>& local_jacobian,
    const IDContainer1& row_ids,
    const IDContainer2& col_ids,
    int dim,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    // assert(local_jacobian.rows() == local_jacobian.cols());
    assert(local_jacobian.rows() % dim == 0);
    assert(local_jacobian.cols() % dim == 0);
    const int n_rows = local_jacobian.rows() / dim;
    const int n_cols = local_jacobian.cols() / dim;
    assert(row_ids.size() >= n_rows); // Can be extra ids
    assert(col_ids.size() >= n_cols); // Can be extra ids
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    triplets.emplace_back(
                        dim * row_ids[i] + k, dim * col_ids[j] + l,
                        local_jacobian(dim * i + k, dim * j + l));
                }
            }
        }
    }
}

class LocalThreadMatStorage {
public:
    std::unique_ptr<MatrixCache> cache = nullptr;

    LocalThreadMatStorage() = delete;

    LocalThreadMatStorage(const int buffer_size, const int rows, const int cols)
    {
        init(buffer_size, rows, cols);
    }

    LocalThreadMatStorage(const int buffer_size, const MatrixCache& c)
    {
        init(buffer_size, c);
    }

    LocalThreadMatStorage(const LocalThreadMatStorage& other)
        : cache(other.cache->copy())
    {
    }

    LocalThreadMatStorage& operator=(const LocalThreadMatStorage& other)
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

    void init(const int buffer_size, const MatrixCache& c)
    {
        if (cache == nullptr) {
            cache = c.copy();
        }
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
                    const auto& val = local_hessian(dim * i + k, dim * j + l);
                    if (val != 0) {
                        triplets.add_value(
                            0, dim * ids[i] + k, dim * ids[j] + l, val);
                    }
                }
            }
        }
    }
}

} // namespace ipc
