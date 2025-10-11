#pragma once

#include "eigen_ext.hpp"

#include <memory>

namespace ipc {
/// abstract class used for caching
class MatrixCache {
public:
    MatrixCache() { }
    virtual ~MatrixCache() = default;

    virtual std::unique_ptr<MatrixCache> copy() const = 0;

    virtual void init(const size_t size) = 0;
    virtual void init(const size_t rows, const size_t cols) = 0;
    virtual void init(const MatrixCache& other) = 0;

    virtual void set_zero() = 0;

    virtual void reserve(const size_t size) = 0;
    virtual size_t entries_size() const = 0;
    virtual size_t capacity() const = 0;
    virtual size_t non_zeros() const = 0;
    virtual size_t triplet_count() const = 0;
    virtual bool is_sparse() const = 0;
    bool is_dense() const { return !is_sparse(); }

    virtual void
    add_value(const int e, const int i, const int j, const double value) = 0;
    virtual Eigen::SparseMatrix<double, Eigen::ColMajor>
    get_matrix(const bool compute_mapping = true) = 0;
    virtual void prune() = 0;

    virtual std::shared_ptr<MatrixCache>
    operator+(const MatrixCache& a) const = 0;
    virtual void operator+=(const MatrixCache& o) = 0;
};

class SparseMatrixCache : public MatrixCache {
public:
    // constructors (call init functions below)
    SparseMatrixCache() = default;
    SparseMatrixCache(const size_t size);
    SparseMatrixCache(const size_t rows, const size_t cols);
    SparseMatrixCache(const MatrixCache& other);
    SparseMatrixCache(
        const SparseMatrixCache& other, const bool copy_main_cache_ptr = false);

    std::unique_ptr<MatrixCache> copy() const override
    {
        // just copy main cache pointer
        return std::make_unique<SparseMatrixCache>(*this, true);
    }

    /// set matrix to be size x size
    void init(const size_t size) override;
    /// set matrix to be rows x cols
    void init(const size_t rows, const size_t cols) override;
    /// set matrix to be a matrix of all zeros with same size as other
    void init(const MatrixCache& other) override;
    /// set matrix to be a matrix of all zeros with same size as other
    /// (potentially with the same main cache)
    void init(
        const SparseMatrixCache& other, const bool copy_main_cache_ptr = false);

    /// set matrix values to zero
    /// modifies tmp_, m_mat, and values (setting all to zero)
    void set_zero() override;

    void reserve(const size_t size) override { m_entries.reserve(size); }

    size_t entries_size() const override { return m_entries.size(); }

    size_t capacity() const override { return m_entries.capacity(); }

    size_t non_zeros() const override
    {
        return m_mapping.empty() ? m_mat.nonZeros() : m_values.size();
    }

    size_t triplet_count() const override
    {
        return m_entries.size() + m_mat.nonZeros();
    }

    bool is_sparse() const override { return true; }

    size_t mapping_size() const { return m_mapping.size(); }

    /// e = element_index, i = global row_index, j = global column_index, value
    /// = value to add to matrix if the cache is yet to be constructed, save the
    /// row, column, and value to be added to the second cache
    ///     in this case, modifies_ entries_ and second_cache_entries_
    /// otherwise, save the value directly in the second cache
    ///     in this case, modfies values_
    void add_value(
        const int e, const int i, const int j, const double value) override;

    /// if the cache is yet to be constructed, save the
    /// cached (ordered) indices in inner_index_ and outer_index_
    /// then fill in map and second_cache_
    ///     in this case, modifies inner_index_, outer_index_, map, and
    ///     second_cache_ to reflect the matrix structure also empties
    ///     second_cache_entries and sets values_ to zero
    /// otherwise, update m_mat directly using the cached indices and values_
    ///     in this case, modifies m_mat and sets values_ to zero
    Eigen::SparseMatrix<double, Eigen::ColMajor>
    get_matrix(const bool compute_mapping = true) override;

    /// if caches have yet to be constructed, add the saved triplets to m_mat
    /// modifies tmp_ and m_mat, also sets entries_ to be empty after writing
    /// its values to m_mat
    void prune() override; ///< add saved entries to stored matrix

    std::shared_ptr<MatrixCache> operator+(const MatrixCache& a) const override;
    std::shared_ptr<MatrixCache> operator+(const SparseMatrixCache& a) const;
    void operator+=(const MatrixCache& o) override;
    void operator+=(const SparseMatrixCache& o);

    const Eigen::SparseMatrix<double, Eigen::ColMajor>& mat() const
    {
        return m_mat;
    }
    const std::vector<Eigen::Triplet<double>>& entries() const
    {
        return m_entries;
    }

private:
    const SparseMatrixCache* main_cache() const
    {
        return m_main_cache == nullptr ? this : m_main_cache;
    }

    const std::vector<std::vector<std::pair<int, size_t>>>& mapping() const
    {
        return main_cache()->m_mapping;
    }

    const std::vector<std::vector<int>>& second_cache() const
    {
        return main_cache()->m_second_cache;
    }

    size_t m_size = 0;

    Eigen::SparseMatrix<double, Eigen::ColMajor> m_tmp;

    Eigen::SparseMatrix<double, Eigen::ColMajor> m_mat;

    /// contains global matrix indices and corresponding value
    std::vector<Eigen::Triplet<double>> m_entries;

    /// maps row indices to column index/local index pairs
    std::vector<std::vector<std::pair<int, size_t>>> m_mapping;

    /// saves inner indices for sparse matrix
    std::vector<int> m_inner_index;

    /// saves inner indices for sparse matrix
    std::vector<int> m_outer_index;

    /// buffer for values (corresponds to inner/outer_index_  structure for
    /// sparse matrix)
    std::vector<double> m_values;

    const SparseMatrixCache* m_main_cache = nullptr;

    /// maps element index to local index
    std::vector<std::vector<int>> m_second_cache;

    /// maps element indices to global matrix indices
    std::vector<std::vector<std::pair<int, int>>> m_second_cache_entries;

    int m_current_e = -1;
    int m_current_e_index = -1;
};

class DenseMatrixCache : public MatrixCache {
public:
    DenseMatrixCache() { }
    DenseMatrixCache(const size_t size);
    DenseMatrixCache(const size_t rows, const size_t cols);
    DenseMatrixCache(const MatrixCache& other);
    DenseMatrixCache(const DenseMatrixCache& other);

    std::unique_ptr<MatrixCache> copy() const override
    {
        return std::make_unique<DenseMatrixCache>(*this);
    }

    void init(const size_t size) override;
    void init(const size_t rows, const size_t cols) override;
    void init(const MatrixCache& other) override;
    void init(const DenseMatrixCache& other);

    void set_zero() override;

    void reserve(const size_t size) override { }
    size_t entries_size() const override { return 0; }
    size_t capacity() const override { return m_mat.size(); }
    size_t non_zeros() const override { return m_mat.size(); }
    size_t triplet_count() const override { return non_zeros(); }
    bool is_sparse() const override { return false; }

    void add_value(
        const int e, const int i, const int j, const double value) override;
    Eigen::SparseMatrix<double, Eigen::ColMajor>
    get_matrix(const bool compute_mapping = true) override;
    void prune() override;

    std::shared_ptr<MatrixCache> operator+(const MatrixCache& a) const override;
    std::shared_ptr<MatrixCache> operator+(const DenseMatrixCache& a) const;
    void operator+=(const MatrixCache& o) override;
    void operator+=(const DenseMatrixCache& o);

    const Eigen::MatrixXd& mat() const { return m_mat; }

private:
    Eigen::MatrixXd m_mat;
};
} // namespace ipc
