#include "interval.hpp"

namespace ipc {

filib::Interval norm(const Eigen::Ref<const VectorMax3I>& x)
{
    filib::Interval sqr_norm(0);
    for (int i = 0; i < x.size(); i++) {
        sqr_norm += sqr(x[i]);
    }
    return sqrt(sqr_norm);
}

} // namespace ipc