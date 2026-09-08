#include "distance_type.hpp"

#include <ipc/utils/logger.hpp>

#include <spdlog/spdlog.h>

#include <stdexcept>

// Host-only definitions for the error reporting helpers declared in
// distance_type.hpp.
//
// These live apart from distance_type.cpp because that file is compiled as a
// CUDA translation unit when CUDA is enabled, and nvcc cannot parse spdlog.
// Keeping the logger and the exceptions here means the shared device sources
// never have to include either one. The inline wrappers in the header pick
// between these and a device-side trap.
namespace ipc::detail {

void warn_degenerate_point_edge_host() noexcept
{
    logger().warn("Degenerate edge in point_edge_distance_type!");
}

void throw_invalid_distance_type_host(const char* function)
{
    throw std::invalid_argument(
        fmt::format("{}: invalid distance type", function));
}

void throw_auto_requires_explicit_dtype_host(const char* function)
{
    throw std::invalid_argument(
        fmt::format(
            "{}: an explicit distance type is required for non-floating-point "
            "scalars; resolving AUTO means comparing single ordered values, "
            "which an autodiff, SIMD batch, or interval scalar does not "
            "provide",
            function));
}

} // namespace ipc::detail
