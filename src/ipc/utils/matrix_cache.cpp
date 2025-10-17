#include "matrix_cache.hpp"

#include <ipc/utils/logger.hpp>
#include <ipc/utils/maybe_parallel_for.hpp>

namespace ipc {
SparseMatrixCache::SparseMatrixCache(const size_t size) { init(size); }

SparseMatrixCache::SparseMatrixCache(const size_t rows, const size_t cols)
{
    init(rows, cols);
}

SparseMatrixCache::SparseMatrixCache(const MatrixCache& other) { init(other); }

SparseMatrixCache::SparseMatrixCache(
    const SparseMatrixCache& other, const bool copy_main_cache_ptr)
{
    init(other, copy_main_cache_ptr);
}

void SparseMatrixCache::init(const size_t size)
{
    assert(mapping().empty() || m_size == size);

    m_size = size;
    m_tmp.resize(m_size, m_size);
    m_mat.resize(m_size, m_size);
    m_mat.setZero();
}

void SparseMatrixCache::init(const size_t rows, const size_t cols)
{
    assert(mapping().empty());

    m_size = rows == cols ? rows : 0;
    m_tmp.resize(rows, cols);
    m_mat.resize(rows, cols);
    m_mat.setZero();
}

void SparseMatrixCache::init(const MatrixCache& other)
{
    assert(this != &other);
    assert(&other == &dynamic_cast<const SparseMatrixCache&>(other));
    init(dynamic_cast<const SparseMatrixCache&>(other));
}

void SparseMatrixCache::init(
    const SparseMatrixCache& other, const bool copy_main_cache_ptr)
{
    assert(this != &other);
    if (copy_main_cache_ptr) {
        m_main_cache = other.m_main_cache;
    } else if (m_main_cache == nullptr) {
        m_main_cache = other.main_cache();
        // Only one level of cache
        assert(
            m_main_cache != this && m_main_cache != nullptr
            && m_main_cache->m_main_cache == nullptr);
    }
    m_size = other.m_size;

    m_values.resize(other.m_values.size());

    m_tmp.resize(other.m_mat.rows(), other.m_mat.cols());
    m_mat.resize(other.m_mat.rows(), other.m_mat.cols());
    m_mat.setZero();
    std::fill(m_values.begin(), m_values.end(), 0);
}

void SparseMatrixCache::set_zero()
{
    m_tmp.setZero();
    m_mat.setZero();

    std::fill(m_values.begin(), m_values.end(), 0);
}

void SparseMatrixCache::add_value(
    const int e, const int i, const int j, const double value)
{
    // caches have yet to be constructed (likely because the matrix has yet to
    // be fully assembled)
    if (mapping().empty()) {
        // save entry so it can be added to the matrix later
        m_entries.emplace_back(i, j, value);

        // save the index information so the cache can be built later
        if (m_second_cache_entries.size() <= e) {
            m_second_cache_entries.resize(e + 1);
        }
        m_second_cache_entries[e].emplace_back(i, j);
    } else {
        if (e != m_current_e) {
            m_current_e = e;
            m_current_e_index = 0;
        }

        // save entry directly to value buffer at the proper index
        m_values[second_cache()[e][m_current_e_index]] += value;
        m_current_e_index++;
    }
}

void SparseMatrixCache::prune()
{
    // caches have yet to be constructed (likely because the matrix has yet to
    // be fully assembled)
    if (mapping().empty()) {
        m_tmp.setFromTriplets(m_entries.begin(), m_entries.end());
        m_tmp.makeCompressed();
        m_mat += m_tmp;

        m_tmp.setZero();
        m_tmp.data().squeeze();
        m_mat.makeCompressed();

        m_entries.clear();

        m_mat.makeCompressed();
    }
}

Eigen::SparseMatrix<double, Eigen::ColMajor>
SparseMatrixCache::get_matrix(const bool compute_mapping)
{
    prune();

    // caches have yet to be constructed (likely because the matrix has yet to
    // be fully assembled)
    if (mapping().empty()) {
        if (compute_mapping && m_size > 0) {
            assert(m_main_cache == nullptr);

            m_values.resize(m_mat.nonZeros());
            m_inner_index.resize(m_mat.nonZeros());
            m_outer_index.resize(m_mat.rows() + 1);
            m_mapping.resize(m_mat.rows());

            // note: m_mat is column major
            const auto inn_ptr = m_mat.innerIndexPtr();
            const auto out_ptr = m_mat.outerIndexPtr();
            m_inner_index.assign(inn_ptr, inn_ptr + m_inner_index.size());
            m_outer_index.assign(out_ptr, out_ptr + m_outer_index.size());

            size_t index = 0;
            // loop over columns of the matrix
            for (size_t i = 0; i < m_mat.cols(); ++i) {
                const auto start = m_outer_index[i];
                const auto end = m_outer_index[i + 1];

                // loop over the nonzero elements of the given column
                for (size_t ii = start; ii < end; ++ii) {
                    // pick out current row
                    const auto j = m_inner_index[ii];
                    auto& map = m_mapping[j];
                    map.emplace_back(i, index);
                    ++index;
                }
            }

            logger().trace("Cache computed");

            m_second_cache.clear();
            m_second_cache.resize(m_second_cache_entries.size());
            // loop over each element
            for (int e = 0; e < m_second_cache_entries.size(); ++e) {
                // loop over each global index affected by the given element
                for (const auto& p : m_second_cache_entries[e]) {
                    const int i = p.first;
                    const int j = p.second;

                    // pick out column/sparse matrix index pairs for the given
                    // column
                    const auto& map = mapping()[i];
                    int local_index = -1;

                    // loop over column/sparse matrix index pairs
                    for (const auto& q : map) {
                        // match columns
                        if (q.first == j) {
                            assert(q.second < m_values.size());
                            local_index = q.second;
                            break;
                        }
                    }
                    assert(local_index >= 0);

                    // save the sparse matrix index used by this element
                    m_second_cache[e].emplace_back(local_index);
                }
            }

            m_second_cache_entries.resize(0);

            logger().trace("Second cache computed");
        }
    } else {
        assert(m_size > 0);
        const auto& outer_index = main_cache()->m_outer_index;
        const auto& inner_index = main_cache()->m_inner_index;
        // directly write the values to the matrix
        m_mat = Eigen::Map<const Eigen::SparseMatrix<double, Eigen::ColMajor>>(
            m_size, m_size, m_values.size(), outer_index.data(),
            inner_index.data(), m_values.data());

        m_current_e = -1;
        m_current_e_index = -1;
    }
    std::fill(m_values.begin(), m_values.end(), 0);
    return m_mat;
}

void SparseMatrixCache::operator+=(const MatrixCache& o)
{
    assert(&o == &dynamic_cast<const SparseMatrixCache&>(o));
    *this += dynamic_cast<const SparseMatrixCache&>(o);
}

void SparseMatrixCache::operator+=(const SparseMatrixCache& o)
{
    if (mapping().empty() || o.mapping().empty()) {
        m_mat += o.m_mat;

        const size_t this_e_size = m_second_cache_entries.size();
        const size_t o_e_size = o.m_second_cache_entries.size();

        m_second_cache_entries.resize(std::max(this_e_size, o_e_size));
        for (int e = 0; e < o_e_size; ++e) {
            assert(
                !m_second_cache_entries[e].empty()
                || !o.m_second_cache_entries[e].empty());
            m_second_cache_entries[e].insert(
                m_second_cache_entries[e].end(),
                o.m_second_cache_entries[e].begin(),
                o.m_second_cache_entries[e].end());
        }
    } else {
        const auto& outer_index = main_cache()->m_outer_index;
        const auto& inner_index = main_cache()->m_inner_index;
        const auto& oouter_index = o.main_cache()->m_outer_index;
        const auto& oinner_index = o.main_cache()->m_inner_index;
        assert(inner_index.size() == oinner_index.size());
        assert(outer_index.size() == oouter_index.size());
        assert(m_values.size() == o.m_values.size());

        maybe_parallel_for(
            o.m_values.size(), [&](int start, int end, int thread_id) {
                for (int i = start; i < end; ++i) {
                    m_values[i] += o.m_values[i];
                }
            });
    }
}

// ========================================================================

DenseMatrixCache::DenseMatrixCache(const size_t size)
{
    m_mat.setZero(size, size);
}

DenseMatrixCache::DenseMatrixCache(const size_t rows, const size_t cols)
{
    m_mat.setZero(rows, cols);
}

DenseMatrixCache::DenseMatrixCache(const MatrixCache& other) { init(other); }

DenseMatrixCache::DenseMatrixCache(const DenseMatrixCache& other)
{
    init(other);
}

void DenseMatrixCache::init(const size_t size) { m_mat.setZero(size, size); }

void DenseMatrixCache::init(const size_t rows, const size_t cols)
{
    m_mat.setZero(rows, cols);
}

void DenseMatrixCache::init(const MatrixCache& other)
{
    init(dynamic_cast<const DenseMatrixCache&>(other));
}

void DenseMatrixCache::init(const DenseMatrixCache& other)
{
    m_mat.setZero(other.m_mat.rows(), other.m_mat.cols());
}

void DenseMatrixCache::set_zero() { m_mat.setZero(); }

void DenseMatrixCache::add_value(
    const int e, const int i, const int j, const double value)
{
    m_mat(i, j) += value;
}

void DenseMatrixCache::prune() { }

Eigen::SparseMatrix<double, Eigen::ColMajor>
DenseMatrixCache::get_matrix(const bool compute_mapping)
{
    return m_mat.sparseView();
}

void DenseMatrixCache::operator+=(const MatrixCache& o)
{
    *this += dynamic_cast<const DenseMatrixCache&>(o);
}

void DenseMatrixCache::operator+=(const DenseMatrixCache& o)
{
    m_mat += o.m_mat;
}

} // namespace ipc