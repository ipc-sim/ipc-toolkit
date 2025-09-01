#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_PROFILER

// clang-format off
#include <spdlog/fmt/bundled/color.h>
#include <ipc/utils/logger.hpp>
// clang-format on

#include <igl/Timer.h>
#include <nlohmann/json.hpp>

#define IPC_TOOLKIT_PROFILE_BLOCK(...)                                         \
    ipc::ProfilePoint __ipc_profile_point(__VA_ARGS__)

namespace ipc {

class Profiler {
public:
    Profiler() = default;

    // ~Profiler() { print(); }

    void clear() { m_data.clear(); }

    void start(const std::string& name)
    {
        current_scope.push_back(name);
        if (!m_data.contains(current_scope)) {
            m_data[current_scope] = {
                { "time_ms", 0 },
                { "count", 0 },
            };
        }
    }

    void stop(const double time_ms)
    {
        const static std::string log_fmt_text = fmt::format(
            "[{}] {{}} {{:.6f}} ms",
            fmt::format(fmt::fg(fmt::terminal_color::magenta), "timing"));

        logger().trace(log_fmt_text, current_scope.to_string(), time_ms);

        assert(m_data.contains(current_scope));
        assert(m_data.at(current_scope).contains("time_ms"));
        assert(m_data.at(current_scope).contains("count"));
        m_data[current_scope]["time_ms"] =
            m_data[current_scope]["time_ms"].get<double>() + time_ms;
        m_data[current_scope]["count"] =
            m_data[current_scope]["count"].get<size_t>() + 1;
        current_scope.pop_back();
    }

    void print() const
    {
        logger().trace(
            "[{}] profiler: {}",
            fmt::format(fmt::fg(fmt::terminal_color::magenta), "timing"),
            m_data.dump(2));
    }

    const nlohmann::json& data() const { return m_data; }
    nlohmann::json& data() { return m_data; }

protected:
    nlohmann::json m_data;
    nlohmann::json::json_pointer current_scope;
};

Profiler& profiler();

template <class Timer = igl::Timer> class ProfilePoint {
public:
    ProfilePoint(const std::string& name) : ProfilePoint(profiler(), name) { }

    ProfilePoint(Profiler& p_profiler, const std::string& name)
        : m_profiler(p_profiler)
    {
        m_profiler.start(name);
        timer.start();
    }

    ~ProfilePoint()
    {
        timer.stop();
        m_profiler.stop(timer.getElapsedTimeInMilliSec());
    }

protected:
    Profiler& m_profiler;
    Timer timer;
};

} // namespace ipc

#else

#define IPC_TOOLKIT_PROFILE_BLOCK(...)

#endif