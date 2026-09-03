#include "area.hpp"

#include <cmath>

namespace ipc::autogen {

// dA is (9×1) flattened in column-major order
template <typename T>
void triangle_area_gradient(
    T t0_x,
    T t0_y,
    T t0_z,
    T t1_x,
    T t1_y,
    T t1_z,
    T t2_x,
    T t2_y,
    T t2_z,
    T dA[9])
{
    const T t0 = -t2_y;
    const T t1 = t0 + t1_y;
    const T t2 = t0_x - t1_x;
    const T t3 = t0 + t0_y;
    const T t4 = -t2_x;
    const T t5 = t0_x + t4;
    const T t6 = t0_y - t1_y;
    const T t7 = t2 * t3 - t5 * t6;
    const T t8 = -t2_z;
    const T t9 = t1_z + t8;
    const T t10 = t0_z + t8;
    const T t11 = t0_z - t1_z;
    const T t12 = t10 * t2 - t11 * t5;
    const T t13 = t10 * t6 - t11 * t3;
    const T t14 = T(0.5) / std::sqrt(t12 * t12 + t13 * t13 + t7 * t7);
    const T t15 = t1_x + t4;
    dA[0] = t14 * (t1 * t7 + t12 * t9);
    dA[1] = -t14 * (-t13 * t9 + t15 * t7);
    dA[2] = -t14 * (t1 * t13 + t12 * t15);
    dA[3] = -t14 * (t10 * t12 + t3 * t7);
    dA[4] = t14 * (-t10 * t13 + t5 * t7);
    dA[5] = t14 * (t12 * t5 + t13 * t3);
    dA[6] = t14 * (t11 * t12 + t6 * t7);
    dA[7] = -t14 * (-t11 * t13 + t2 * t7);
    dA[8] = -t14 * (t12 * t2 + t13 * t6);
}

#define IPC_INSTANTIATE_AREA_AUTOGEN(T)                                        \
    template void triangle_area_gradient<T>(T, T, T, T, T, T, T, T, T, T[9])

IPC_INSTANTIATE_AREA_AUTOGEN(float);
IPC_INSTANTIATE_AREA_AUTOGEN(double);
#undef IPC_INSTANTIATE_AREA_AUTOGEN

} // namespace ipc::autogen
