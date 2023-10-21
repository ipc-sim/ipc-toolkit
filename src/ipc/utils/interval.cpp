#include "interval.hpp"

namespace ipc {

filib::Interval squared_norm(const Eigen::Ref<const VectorMax3I>& x)
{
    filib::Interval sqr_norm(0);
    for (int i = 0; i < x.size(); i++) {
        sqr_norm += sqr(x[i]);
    }
    return sqr_norm;
}

filib::Interval norm(const Eigen::Ref<const VectorMax3I>& x)
{
    return sqrt(squared_norm(x));
}

} // namespace ipc