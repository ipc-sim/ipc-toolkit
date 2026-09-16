// Raw device memory for the ipc::cuda implementation files.
// This header is CUDA-only and must be included from .cu files exclusively.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/utils/cuda/device_utils.cuh>

#include <cuda_runtime.h>

#include <cstddef>
#include <utility>

namespace ipc::cuda {

/// @brief A device array that is never initialized, only ever grows, and is
/// released without throwing.
///
/// This replaces thrust::device_vector for per-frame scratch and output
/// buffers, for two reasons:
///
/// 1. device_vector::resize() value-initializes, launching a fill over every
///    appended element. For a buffer the next kernel writes from scratch that
///    is a full pass of wasted memory traffic, and for a candidate buffer sized
///    to its high-water mark it is the dominant per-call cost.
/// 2. device_vector's destructor throws if cudaFree fails. After a sticky
///    kernel fault every CUDA call fails, so unwinding past a live
///    device_vector calls std::terminate and the caller never sees the error
///    IPC_TOOLKIT_CUDA_CHECK raised. This buffer ignores cudaFree's result.
///
/// The capacity is a high-water mark: clear() and a smaller resize() keep the
/// allocation, so a per-frame rebuild of the same mesh allocates nothing.
/// Reallocating does NOT preserve the contents; every user here refills from
/// scratch.
template <typename T> class DeviceBuffer {
public:
    DeviceBuffer() = default;
    ~DeviceBuffer() { release(); }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& other) noexcept
        : m_data(std::exchange(other.m_data, nullptr))
        , m_size(std::exchange(other.m_size, size_t(0)))
        , m_capacity(std::exchange(other.m_capacity, size_t(0)))
    {
    }

    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept
    {
        if (this != &other) {
            release();
            m_data = std::exchange(other.m_data, nullptr);
            m_size = std::exchange(other.m_size, size_t(0));
            m_capacity = std::exchange(other.m_capacity, size_t(0));
        }
        return *this;
    }

    /// @brief Ensure room for @p n elements.
    /// @note Reallocating discards the contents; nothing is initialized.
    void reserve(const size_t n)
    {
        if (n <= m_capacity) {
            return;
        }
        release();
        IPC_TOOLKIT_CUDA_CHECK(
            cudaMalloc(reinterpret_cast<void**>(&m_data), n * sizeof(T)));
        m_capacity = n;
    }

    /// @brief Set the size to @p n, reserving as needed (see reserve()).
    void resize(const size_t n)
    {
        reserve(n);
        m_size = n;
    }

    /// @brief Set the size to zero, keeping the allocation.
    void clear() { m_size = 0; }

    /// @brief Free the allocation.
    /// @note Never throws: after a sticky device error cudaFree fails too, and
    /// a destructor cannot propagate that without terminating the program.
    void release() noexcept
    {
        if (m_data != nullptr) {
            static_cast<void>(cudaFree(m_data));
            m_data = nullptr;
        }
        m_size = 0;
        m_capacity = 0;
    }

    /// @brief Set every byte of the first size() elements to @p byte
    /// (asynchronous, on the default stream).
    void fill_bytes(const int byte)
    {
        if (m_size > 0) {
            IPC_TOOLKIT_CUDA_CHECK(
                cudaMemsetAsync(m_data, byte, m_size * sizeof(T)));
        }
    }

    /// @brief Zero the first size() elements (asynchronous).
    void zero() { fill_bytes(0); }

    /// @brief Resize to @p n and copy @p n elements from host memory.
    void upload(const T* host, const size_t n)
    {
        resize(n);
        if (n > 0) {
            IPC_TOOLKIT_CUDA_CHECK(cudaMemcpy(
                m_data, host, n * sizeof(T), cudaMemcpyHostToDevice));
        }
    }

    /// @brief Copy the first size() elements to host memory (synchronous).
    void download(T* host) const
    {
        if (m_size > 0) {
            IPC_TOOLKIT_CUDA_CHECK(cudaMemcpy(
                host, m_data, m_size * sizeof(T), cudaMemcpyDeviceToHost));
        }
    }

    T* data() { return m_data; }
    const T* data() const { return m_data; }
    size_t size() const { return m_size; }
    size_t capacity() const { return m_capacity; }
    bool empty() const { return m_size == 0; }

private:
    T* m_data = nullptr;
    size_t m_size = 0;
    size_t m_capacity = 0;
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
