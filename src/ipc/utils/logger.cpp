#include "logger.hpp"

#include <spdlog/sinks/stdout_color_sinks.h>

#include <sstream>

namespace ipc {

namespace {

    // Custom logger instance defined by the user, if any
    std::shared_ptr<spdlog::logger>& get_shared_logger()
    {
        static std::shared_ptr<spdlog::logger> logger;
        return logger;
    }

} // namespace

// Retrieve current logger
spdlog::logger& logger()
{
    if (get_shared_logger()) {
        return *get_shared_logger();
    } else {
        // Create the logger manually instead of using spdlog's factory
        // methods (_st and _mt functions): the factories register the logger
        // in spdlog's global registry, whose name-uniqueness check fails when
        // several shared libraries each embed a copy of the toolkit (e.g.,
        // two Python modules). See
        // https://github.com/gabime/spdlog/wiki/2.-Creating-loggers
        static auto default_logger = std::make_shared<spdlog::logger>(
            "ipctk",
            std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
        return *default_logger;
    }
}

// Use a custom logger
void set_logger(std::shared_ptr<spdlog::logger> logger)
{
    get_shared_logger() = std::move(logger);
}

void log_and_throw_error(const std::string& msg)
{
    logger().error(msg);
    throw std::runtime_error(msg);
}

} // namespace ipc
