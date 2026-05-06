#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_TRACY
#include <tracy/Tracy.hpp>
#endif

#include <string>

// Helper macro to stringify/paste after expansion
#define IPC_TOOLKIT_PROFILE_BLOCK_CONCAT_IMPL(a, b) a##b
#define IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(a, b)                                 \
    IPC_TOOLKIT_PROFILE_BLOCK_CONCAT_IMPL(a, b)

#ifdef IPC_TOOLKIT_WITH_PROFILER

// clang-format off
#include <spdlog/fmt/bundled/color.h>
#include <ipc/utils/logger.hpp>
// clang-format on

#include <nlohmann/json.hpp>

#include <chrono>
#include <iostream>
#include <thread>

#if defined(IPC_TOOLKIT_WITH_TRACY) && defined(TRACY_ENABLE)
// ZoneScoped + ZoneName supports both compile-time string literals and
// runtime std::string expressions, unlike ZoneScopedN which requires a
// constexpr const char*.
#define IPC_TOOLKIT_PROFILE_BLOCK(...)                                         \
    const std::string IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(                        \
        __ipc_zone_name_, __LINE__)(__VA_ARGS__);                              \
    ipc::ProfilePoint IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(                        \
        __ipc_profile_point_, __LINE__)(                                       \
        IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(__ipc_zone_name_, __LINE__));         \
    ZoneScoped;                                                                \
    ZoneName(                                                                  \
        IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(__ipc_zone_name_, __LINE__).c_str(),  \
        IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(__ipc_zone_name_, __LINE__).size())
#else
#define IPC_TOOLKIT_PROFILE_BLOCK(...)                                         \
    ipc::ProfilePoint IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(                        \
        __ipc_profile_point_, __COUNTER__)(__VA_ARGS__)
#endif

namespace ipc {

class Profiler {
public:
    Profiler();

    // ~Profiler() { print(); }

    /// @brief Clear all profiling data.
    void clear();

    /// @brief Start timing a new scope.
    /// @param name The name of the scope.
    void start(const std::string& name);

    /// @brief Stop timing the current scope and record the elapsed time.
    /// @param time_ms The elapsed time in milliseconds.
    void stop(const double time_ms);

    /// @brief Reset the profiler data and scopes.
    void reset();

    /// @brief Print the profiler data to the logger.
    void print() const;

    /// @brief Write the profiler data to an output stream in CSV format.
    void write_csv(std::ostream& os) const;

    /// @brief Write the profiler data to a CSV file.
    void write_csv(const std::string& filename) const;

    /// @brief Print the profiler data to standard output in CSV format.
    void print_csv() const { write_csv(std::cout); }

    /// @brief Access the profiling data as a JSON object.
    const nlohmann::json& data() const { return m_data; }

    /// @brief Access the profiling data as a JSON object.
    nlohmann::json& data() { return m_data; }

    bool is_recording_thread() const
    {
        return std::this_thread::get_id() == m_main_thread_id;
    }

protected:
    /// @brief The profiling data stored as a JSON object.
    nlohmann::json m_data;

    /// @brief The global scope pointer into the JSON data.
    nlohmann::json::json_pointer m_current_scope;

    /// @brief The thread that records data; calls from all other threads are
    ///        silently ignored, giving a single-thread estimate of block costs.
    std::thread::id m_main_thread_id;
};

Profiler& profiler();

class ChronoTimer {
public:
    void start() { m_start = std::chrono::high_resolution_clock::now(); }

    void stop() { m_end = std::chrono::high_resolution_clock::now(); }

    // NOLINTNEXTLINE(readability-identifier-naming)
    double getElapsedTimeInMilliSec() const
    {
        return std::chrono::duration<double, std::milli>(m_end - m_start)
            .count();
    }

private:
    std::chrono::high_resolution_clock::time_point m_start, m_end;
};

template <class Timer = ChronoTimer> class ProfilePoint {
public:
    ProfilePoint(const std::string& name) : ProfilePoint(profiler(), name) { }

    ProfilePoint(Profiler& p_profiler, const std::string& name)
        : m_profiler(p_profiler)
        , m_active(p_profiler.is_recording_thread())
    {
        if (m_active) {
            m_profiler.start(name);
            timer.start();
        }
    }

    ~ProfilePoint()
    {
        if (m_active) {
            timer.stop();
            m_profiler.stop(timer.getElapsedTimeInMilliSec());
        }
    }

protected:
    Profiler& m_profiler;
    Timer timer;
    bool m_active;
};

} // namespace ipc

#elif defined(IPC_TOOLKIT_WITH_TRACY) && defined(TRACY_ENABLE)

// Custom profiler disabled: Tracy zone only.
// ZoneScoped + ZoneName supports runtime strings, unlike ZoneScopedN.
#define IPC_TOOLKIT_PROFILE_BLOCK(...)                                         \
    const std::string IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(                        \
        __ipc_zone_name_, __LINE__)(__VA_ARGS__);                              \
    ZoneScoped;                                                                \
    ZoneName(                                                                  \
        IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(__ipc_zone_name_, __LINE__).c_str(),  \
        IPC_TOOLKIT_PROFILE_BLOCK_CONCAT(__ipc_zone_name_, __LINE__).size())

#else

// No profiling enabled: no-op.
#define IPC_TOOLKIT_PROFILE_BLOCK(...)

#endif