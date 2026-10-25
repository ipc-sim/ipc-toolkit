#pragma once

// We include <cmath> first so M_PI is visible whenever the platform provides
// it, regardless of what the including file pulled in before us.
#include <cmath> // IWYU pragma: keep (provides M_PI where available)

namespace ipc {

// M_PI is not standard C++ (MSVC only defines it when _USE_MATH_DEFINES is
// set before <cmath>), so we use it when available and fall back to our own
// literal otherwise, following libigl's igl/PI.h. We avoid EIGEN_PI because it
// is a long double on GCC, and narrowing it trips -Werror=float-conversion.
#ifdef M_PI
/// @brief π
constexpr double PI = M_PI;
#else
/// @brief π
constexpr double PI = 3.1415926535897932384626433832795;
#endif

} // namespace ipc
