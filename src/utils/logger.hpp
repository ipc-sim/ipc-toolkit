#pragma once

#ifdef IPC_TOOLKIT_WITH_LOGGER

#include <spdlog/spdlog.h>
#include <spdlog/async.h>
#include <spdlog/fmt/bundled/ranges.h>

#include <ostream>

#define IPC_LOG(message) logger().message

namespace ipc {

struct Logger {
    static std::shared_ptr<spdlog::async_logger> _logger;

    // By default, write to stdout, but don't write to any file
    static void init(
        bool use_cout = true,
        const std::string& filename = "",
        bool truncate = true);
    static void init(std::ostream& os);
    static void init(std::vector<spdlog::sink_ptr>& sinks);
};

// Retrieve current logger, or create one if not available
inline spdlog::async_logger& logger()
{
    if (!Logger::_logger) {
        Logger::init();
    }
    return *Logger::_logger;
}

} // namespace ipc

#else

#include <fmt/format.h>

#define IPC_LOG(message)

#endif
