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

} // namespace ipc
