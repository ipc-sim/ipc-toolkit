#pragma once

#ifdef IPC_TOOLKIT_WITH_LOGGER

// clang-format off
#include <spdlog/fmt/bundled/ranges.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
// clang-format on

/// Helper macro that makes it easier to disable logging completly at compile
/// time (e.g., IPC_LOG(info("A useful info log"))).
#define IPC_LOG(message) ipc::logger().message

namespace ipc {

/// Retrieves the current logger.
/// @return A const reference to the logger object.
spdlog::logger& logger();

/// Setup a logger object. Calling this function with other function is not
/// thread-safe.
/// @param[in] logger New logger object to be used.
void set_logger(std::shared_ptr<spdlog::logger> logger);

} // namespace ipc

#else

#include <fmt/format.h>

// Dummy macros to disable logging completly.
#define IPC_LOG(message)

#endif
