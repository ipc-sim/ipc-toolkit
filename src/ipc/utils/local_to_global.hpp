#pragma once

#include <ipc/config.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/matrix_cache.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>

namespace ipc {

/// @brief Add a local gradient to the global gradient, given the vertex ids and dimension.
/// @tparam IDContainer A container type that supports operator[] and size() to access vertex ids. Can be a std::vector, Eigen::VectorXi, etc.
/// @tparam GlobalOrder The storage order of the global gradient. Must be either Eigen::RowMajor or Eigen::ColMajor. If RowMajor, the global gradient is assumed to be ordered as [x0, y0, z0, x1, y1, z1, ...]. If ColMajor, the global gradient is assumed to be ordered as [x0, x1, ..., y0, y1, ..., z0, z1, ...].
/// @param local_grad The local gradient to be added to the global gradient. Must be of size n_verts * dim, where n_verts is the number of vertices in the local element and dim is the dimension of the problem.
/// @param ids The vertex ids corresponding to the local gradient. Must be of size at least n_verts, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @param dim The dimension of the problem (e.g., 2 for 2D, 3 for 3D).
/// @param grad The global gradient to which the local gradient will be added.
template <typename IDContainer, int GlobalOrder = VERTEX_DERIVATIVE_LAYOUT>
void local_gradient_to_global_gradient(
    Eigen::ConstRef<Eigen::VectorXd> local_grad,
    const IDContainer& ids,
    const int dim,
    Eigen::Ref<Eigen::VectorXd> grad)
{
    static_assert(
        GlobalOrder == Eigen::ColMajor || GlobalOrder == Eigen::RowMajor);
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    const int n_total_verts = grad.size() / dim;
    assert(grad.size() % dim == 0); // Ensure grad is properly sized
    for (int i = 0; i < n_verts; i++) {
        assert(ids[i] >= 0); // Ensure valid vertex id
        if constexpr (GlobalOrder == Eigen::RowMajor) {
            grad.segment(dim * ids[i], dim) += local_grad.segment(dim * i, dim);
        } else {
            for (int d = 0; d < dim; d++) {
                grad[n_total_verts * d + ids[i]] += local_grad[dim * i + d];
            }
        }
    }
}

/// @brief Add a local gradient to the global gradient, given the vertex ids and dimension.
/// @tparam IDContainer A container type that supports operator[] and size() to access vertex ids. Can be a std::vector, Eigen::VectorXi, etc.
/// @tparam GlobalOrder The storage order of the global gradient. Must be either Eigen::RowMajor or Eigen::ColMajor. If RowMajor, the global gradient is assumed to be ordered as [x0, y0, z0, x1, y1, z1, ...]. If ColMajor, the global gradient is assumed to be ordered as [x0, x1, ..., y0, y1, ..., z0, z1, ...].
/// @param local_grad The local gradient to be added to the global gradient. Must be of size n_verts * dim, where n_verts is the number of vertices in the local element and dim is the dimension of the problem.
/// @param ids The vertex ids corresponding to the local gradient. Must be of size at least n_verts, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @param dim The dimension of the problem (e.g., 2 for 2D, 3 for 3D).
/// @param grad The global gradient to which the local gradient will be added.
template <typename IDContainer, int GlobalOrder = VERTEX_DERIVATIVE_LAYOUT>
void local_gradient_to_global_gradient(
    Eigen::ConstRef<Eigen::VectorXd> local_grad,
    const IDContainer& ids,
    const int dim,
    Eigen::SparseVector<double>& grad)
{
    static_assert(
        GlobalOrder == Eigen::ColMajor || GlobalOrder == Eigen::RowMajor);
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    const int n_total_verts = grad.size() / dim;
    assert(grad.size() % dim == 0); // Ensure grad is properly sized
    for (int i = 0; i < n_verts; i++) {
        assert(ids[i] >= 0); // Ensure valid vertex id
        for (int d = 0; d < dim; d++) {
            if constexpr (GlobalOrder == Eigen::RowMajor) {
                grad.coeffRef(dim * ids[i] + d) += local_grad(dim * i + d);
            } else {
                grad.coeffRef(n_total_verts * d + ids[i]) +=
                    local_grad(dim * i + d);
            }
        }
    }
}

/// @brief Add a local Hessian to the global Hessian, given the vertex ids and dimension.
///
/// The local Hessian is added to the global Hessian in triplet form, so that it
/// can be efficiently assembled into a sparse matrix later.
///
/// @tparam IDContainer A container type that supports operator[] and size() to access vertex ids. Can be a std::vector, Eigen::VectorXi, etc.
/// @tparam GlobalOrder The storage order of the global Hessian. Must be either Eigen::RowMajor or Eigen::ColMajor. If RowMajor, the global Hessian is assumed to be ordered as [x0, y0, z0, x1, y1, z1, ...]. If ColMajor, the global Hessian is assumed to be ordered as [x0, x1, ..., y0, y1, ..., z0, z1, ...].
/// @param local_hessian The local Hessian to be added to the global Hessian. Must be of size (n_verts * dim) x (n_verts * dim), where n_verts is the number of vertices in the local element and dim is the dimension of the problem.
/// @param ids The vertex ids corresponding to the local gradient. Must be of size at least n_verts, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @param dim The dimension of the problem (e.g., 2 for 2D, 3 for 3D).
/// @param triplets The vector of triplets to which the local Hessian will be added.
/// @param n_total_verts The total number of vertices in the global mesh. Only required when GlobalOrder is ColMajor. Ignored when GlobalOrder is RowMajor.
template <typename IDContainer, int GlobalOrder = VERTEX_DERIVATIVE_LAYOUT>
void local_hessian_to_global_triplets(
    Eigen::ConstRef<Eigen::MatrixXd> local_hessian,
    const IDContainer& ids,
    const int dim,
    std::vector<Eigen::Triplet<double>>& triplets,
    const int n_total_verts = -1)
{
    static_assert(
        GlobalOrder == Eigen::ColMajor || GlobalOrder == Eigen::RowMajor);
    assert(local_hessian.rows() == local_hessian.cols());
    assert(local_hessian.rows() % dim == 0);
    const int n_verts = local_hessian.rows() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    if constexpr (GlobalOrder == Eigen::ColMajor) {
        if (n_total_verts <= 0) {
            log_and_throw_error(
                "Total number of vertices must be provided for ColMajor storage order.");
        }
    }
    for (int i = 0; i < n_verts; i++) {
        assert(ids[i] >= 0); // Ensure valid vertex id
        for (int j = 0; j < n_verts; j++) {
            assert(ids[j] >= 0); // Ensure valid vertex id
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    if constexpr (GlobalOrder == Eigen::RowMajor) {
                        triplets.emplace_back(
                            dim * ids[i] + k, dim * ids[j] + l,
                            local_hessian(dim * i + k, dim * j + l));
                    } else {
                        triplets.emplace_back(
                            n_total_verts * k + ids[i],
                            n_total_verts * l + ids[j],
                            local_hessian(dim * i + k, dim * j + l));
                    }
                }
            }
        }
    }
}

/// @brief Add a local Jacobian to the global Jacobian, given the vertex ids and dimension.
///
/// The local Jacobian is added to the global Jacobian in triplet form, so that
/// it can be efficiently assembled into a sparse matrix later.
///
/// @tparam Derived The derived type of the local Jacobian matrix. Must be a matrix expression that can be evaluated to an Eigen::MatrixXd. The local Jacobian must be of size (n_rows * dim) x (n_cols * dim), where n_rows and n_cols are the number of vertices in the local element for the row and column respectively, and dim is the dimension of the problem.
/// @tparam IDContainer1 A container type that supports operator[] and size() to access vertex ids for the row indices. Can be a std::vector, Eigen::VectorXi, etc. Must be of size at least n_rows, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @tparam IDContainer2 A container type that supports operator[] and size() to access vertex ids for the column indices. Can be a std::vector, Eigen::VectorXi, etc. Must be of size at least n_cols, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @tparam GlobalOrder The storage order of the global Jacobian. Must be either Eigen::RowMajor or Eigen::ColMajor. If RowMajor, the global Jacobian is assumed to be ordered as [x0, y0, z0, x1, y1, z1, ...]. If ColMajor, the global Jacobian is assumed to be ordered as [x0, x1, ..., y0, y1, ..., z0, z1, ...].
/// @param local_jacobian The local Jacobian to be added to the global Jacobian. Must be of size (n_rows * dim) x (n_cols * dim), where n_rows and n_cols are the number of vertices in the local element for the row and column respectively, and dim is the dimension of the problem.
/// @param row_ids The vertex ids corresponding to the row indices of the local Jacobian. Must be of size at least n_rows, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @param col_ids The vertex ids corresponding to the column indices of the local Jacobian. Must be of size at least n_cols, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @param dim The dimension of the problem (e.g., 2 for 2D, 3 for 3D).
/// @param triplets The vector of triplets to which the local Hessian will be added.
/// @param n_total_rows The total number of row vertices in the global mesh. Only required when GlobalOrder is ColMajor. Ignored when GlobalOrder is RowMajor.
/// @param n_total_cols The total number of column vertices in the global mesh. Only required when GlobalOrder is ColMajor. Ignored when GlobalOrder is RowMajor.
template <
    typename Derived,
    typename IDContainer1,
    typename IDContainer2,
    int GlobalOrder = VERTEX_DERIVATIVE_LAYOUT>
void local_jacobian_to_global_triplets(
    const Eigen::MatrixBase<Derived>& local_jacobian,
    const IDContainer1& row_ids,
    const IDContainer2& col_ids,
    int dim,
    std::vector<Eigen::Triplet<double>>& triplets,
    const int n_total_rows = -1,
    const int n_total_cols = -1)
{
    static_assert(
        GlobalOrder == Eigen::ColMajor || GlobalOrder == Eigen::RowMajor);
    // assert(local_jacobian.rows() == local_jacobian.cols());
    assert(local_jacobian.rows() % dim == 0);
    assert(local_jacobian.cols() % dim == 0);
    const int n_rows = local_jacobian.rows() / dim;
    const int n_cols = local_jacobian.cols() / dim;
    assert(row_ids.size() >= n_rows); // Can be extra ids
    assert(col_ids.size() >= n_cols); // Can be extra ids
    if constexpr (GlobalOrder == Eigen::ColMajor) {
        if (n_total_rows <= 0 || n_total_cols <= 0) {
            log_and_throw_error(
                "Total number of rows and columns must be provided for ColMajor storage order.");
        }
    }
    for (int i = 0; i < n_rows; i++) {
        assert(row_ids[i] >= 0); // Ensure valid vertex id
        for (int j = 0; j < n_cols; j++) {
            assert(col_ids[j] >= 0); // Ensure valid vertex id
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    if constexpr (GlobalOrder == Eigen::RowMajor) {
                        triplets.emplace_back(
                            dim * row_ids[i] + k, dim * col_ids[j] + l,
                            local_jacobian(dim * i + k, dim * j + l));
                    } else {
                        triplets.emplace_back(
                            n_total_rows * k + row_ids[i],
                            n_total_cols * l + col_ids[j],
                            local_jacobian(dim * i + k, dim * j + l));
                    }
                }
            }
        }
    }
}

/// @brief A helper class to store a local matrix cache for each thread.
/// This is useful for parallelizing the assembly of the global Hessian, where
/// each thread can have its own local cache to store the triplets before they
/// are merged into the global cache.
struct LocalThreadMatStorage {
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

/// @brief Add a local Hessian to the global Hessian, given the vertex ids and dimension.
///
/// The local Hessian is added to the global Hessian in triplet form, so that it
/// can be efficiently assembled into a sparse matrix later.
///
/// @tparam Derived The derived type of the local Hessian matrix. Must be a matrix expression that can be evaluated to an Eigen::MatrixXd. The local Hessian must be of size (n_verts * dim) x (n_verts * dim), where n_verts is the number of vertices in the local element and dim is the dimension of the problem.
/// @tparam IDContainer A container type that supports operator[] and size() to access vertex ids. Can be a std::vector, Eigen::VectorXi, etc. Must be of size at least n_verts, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @tparam GlobalOrder The storage order of the global gradient. Must be either Eigen::RowMajor or Eigen::ColMajor. If RowMajor, the global gradient is assumed to be ordered as [x0, y0, z0, x1, y1, z1, ...]. If ColMajor, the global gradient is assumed to be ordered as [x0, x1, ..., y0, y1, ..., z0, z1, ...].
/// @param local_hessian The local Hessian to be added to the global Hessian. Must be of size (n_verts * dim) x (n_verts * dim), where n_verts is the number of vertices in the local element and dim is the dimension of the problem.
/// @param ids The vertex ids corresponding to the local gradient. Must be of size at least n_verts, but can be larger (e.g., if the local element has fewer vertices than the number of ids).
/// @param dim The dimension of the problem (e.g., 2 for 2D, 3 for 3D).
/// @param triplets The MatrixCache to which the local Hessian will be added. The triplets will be added to the first matrix in the cache (i.e., matrix index 0). The cache should be initialized with the correct number of rows and columns for the global Hessian, and should have enough reserved space for the number of triplets that will be added.
/// @param n_total_verts The total number of vertices in the global mesh. Only required when GlobalOrder is ColMajor. Ignored when GlobalOrder is RowMajor.
template <
    typename Derived,
    typename IDContainer,
    int GlobalOrder = VERTEX_DERIVATIVE_LAYOUT>
void local_hessian_to_global_triplets(
    const Eigen::MatrixBase<Derived>& local_hessian,
    const IDContainer& ids,
    int dim,
    MatrixCache& triplets,
    const int n_total_verts = -1)
{
    static_assert(
        GlobalOrder == Eigen::ColMajor || GlobalOrder == Eigen::RowMajor);
    assert(local_hessian.rows() == local_hessian.cols());
    assert(local_hessian.rows() % dim == 0);
    const int n_verts = local_hessian.rows() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    if constexpr (GlobalOrder == Eigen::ColMajor) {
        if (n_total_verts <= 0) {
            log_and_throw_error(
                "Total number of vertices must be provided for ColMajor storage order.");
        }
    }
    for (int i = 0; i < n_verts; i++) {
        for (int j = 0; j < n_verts; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    const auto& val = local_hessian(dim * i + k, dim * j + l);
                    if (val != 0) {
                        if constexpr (GlobalOrder == Eigen::RowMajor) {
                            triplets.add_value(
                                0, dim * ids[i] + k, dim * ids[j] + l, val);
                        } else {
                            triplets.add_value(
                                0, n_total_verts * k + ids[i],
                                n_total_verts * l + ids[j], val);
                        }
                    }
                }
            }
        }
    }
}

} // namespace ipc
