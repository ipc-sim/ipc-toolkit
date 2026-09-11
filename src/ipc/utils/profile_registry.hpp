#pragma once

#include <chrono>
#include <mutex>
#include <string>
#include <unordered_map>

namespace ipc {

/// @brief Thread-safe flat registry of named accumulating timers/counters.
///
/// Each registered name stores a running `total` and a `count`. For timings,
/// `total` is milliseconds; for counters, `total` is the summed value. The
/// registry is meant to replace noisy per-call log statements: callers push
/// samples via `add_time` / `add_value`, and the driver code (e.g. the
/// polyfem time-stepping loop) periodically flushes the state to a JSON file
/// via `dump_json`.
class ProfileRegistry {
public:
    struct Stat {
        double total = 0.0;
        std::size_t count = 0;
    };

    static ProfileRegistry& instance();

    /// @brief Add a timing sample (in milliseconds) to the named entry.
    void add_time(const std::string& name, double ms);

    /// @brief Add a numeric sample (e.g. a collision-set size) to the named
    /// entry. Total is the running sum, count is the number of samples.
    void add_value(const std::string& name, double value);

    /// @brief Clear all accumulated stats.
    void reset();

    /// @brief Overwrite `path` with the current registry contents as JSON.
    ///
    /// The output format is a flat object keyed by name:
    ///
    ///     {
    ///       "ho.broad_phase": { "total": 100342.5, "count": 4437, "mean": 22.6
    ///       }, "ho.collision_set.vertex_dicts": { "total": 10786200, "count":
    ///       4437, "mean": 2431.2 },
    ///       ...
    ///     }
    void dump_json(const std::string& path) const;

private:
    ProfileRegistry() = default;

    mutable std::mutex m_mutex;
    std::unordered_map<std::string, Stat> m_stats;
};

/// @brief RAII scope guard that accumulates elapsed wall time (ms) into the
/// registry when it goes out of scope.
class ScopedProfileTimer {
public:
    explicit ScopedProfileTimer(std::string name)
        : m_name(std::move(name))
        , m_start(std::chrono::high_resolution_clock::now())
    {
    }

    ~ScopedProfileTimer()
    {
        const auto end = std::chrono::high_resolution_clock::now();
        const double ms =
            std::chrono::duration<double, std::milli>(end - m_start).count();
        ProfileRegistry::instance().add_time(m_name, ms);
    }

    ScopedProfileTimer(const ScopedProfileTimer&) = delete;
    ScopedProfileTimer& operator=(const ScopedProfileTimer&) = delete;

private:
    std::string m_name;
    std::chrono::high_resolution_clock::time_point m_start;
};

#define IPC_PROFILE_CONCAT_IMPL(a, b) a##b
#define IPC_PROFILE_CONCAT(a, b) IPC_PROFILE_CONCAT_IMPL(a, b)

/// @brief Accumulate the wall time of the enclosing scope into the registry
/// under `name`.
#define IPC_PROFILE_SCOPE(name)                                                \
    ::ipc::ScopedProfileTimer IPC_PROFILE_CONCAT(                              \
        _ipc_profile_scope_, __LINE__)(name)

} // namespace ipc
