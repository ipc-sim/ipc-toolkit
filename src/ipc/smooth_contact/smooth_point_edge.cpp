#include <ipc/utils/quadrature.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include "smooth_point_edge.hpp"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <stdexcept>
#include <iostream>

DECLARE_DIFFSCALAR_BASE();

namespace ipc {

}