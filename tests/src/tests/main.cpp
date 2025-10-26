// Catch2 Documentation: https://github.com/catchorg/Catch2/tree/master/docs

#include <catch2/catch_session.hpp>

#include <ipc/utils/logger.hpp>

#include <tbb/global_control.h>
#include <tbb/info.h>

Catch::Clara::ParserResult parse_log_level(int const d, int& log_level)
{
    if (d < 0 || d > spdlog::level::off) {
        return Catch::Clara::ParserResult::runtimeError(
            "Log level must be between 0 and 6");
    } else {
        log_level = d;
        return Catch::Clara::ParserResult::ok(
            Catch::Clara::ParseResultType::Matched);
    }
}

Catch::Clara::ParserResult parse_num_threads(int const d, int& num_threads)
{
    if (num_threads <= 0) {
        num_threads = tbb::info::default_concurrency();
    } else if (num_threads > tbb::info::default_concurrency()) {
        ipc::logger().warn(
            "Attempting to use more threads than available ({:d} > {:d})!",
            num_threads, tbb::info::default_concurrency());
        num_threads = tbb::info::default_concurrency();
    } else {
        num_threads = d;
    }
    return Catch::Clara::ParserResult::ok(
        Catch::Clara::ParseResultType::Matched);
}

int main(int argc, char** argv)
{
    Catch::Session session; // There must be exactly one instance

    // Build a new parser on top of Catch's
    using namespace Catch::Clara;
    auto cli = session.cli();

    int log_level = spdlog::level::warn;
    cli |=
        Opt([&log_level](int const d) { return parse_log_level(d, log_level); },
            "log_level")["--log"]["--logger-level"](
            "logger verbosity level int (0-6)");

    int num_threads = tbb::info::default_concurrency();
    cli |=
        Opt([&num_threads](
                int const d) { return parse_num_threads(d, num_threads); },
            "num_threads")["--nthreads"]["--num-threads"](
            "maximum number of threads to use");

    session.cli(cli);

    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) // Indicates a command line error
        return returnCode;

    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
    ipc::logger().set_level(static_cast<spdlog::level::level_enum>(log_level));

    tbb::global_control thread_limiter(
        tbb::global_control::max_allowed_parallelism, num_threads);

    return session.run();
}
