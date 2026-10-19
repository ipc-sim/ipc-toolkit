#pragma once

// clang-format off
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
// clang-format on

namespace ipc {

/// Retrieves the current logger.
/// @return A const reference to the logger object.
spdlog::logger& logger();

/// Setup a logger object. Calling this function with other function is not
/// thread-safe.
/// @param logger New logger object to be used.
void set_logger(std::shared_ptr<spdlog::logger> logger);

template <typename Exception = std::runtime_error>
[[noreturn]] void log_and_throw_error(const std::string& msg)
{
    logger().error(msg);
    throw Exception(msg);
}

template <typename... Args>
[[noreturn]] void
log_and_throw_error(const std::string& msg, const Args&... args)
{
    log_and_throw_error<std::runtime_error>(
        fmt::format(fmt::runtime(msg), args...));
}

template <typename Exception, typename... Args>
[[noreturn]] void
log_and_throw_error(const std::string& msg, const Args&... args)
{
    log_and_throw_error<Exception>(fmt::format(fmt::runtime(msg), args...));
}

} // namespace ipc
