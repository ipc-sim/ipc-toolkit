#pragma once

#include <cstdint>
#include <vector>

namespace ipc {

enum class HeavisideType : uint8_t { ZERO = 0, ONE = 1, VARIANT = 2 };

struct OrientationTypes {

    static HeavisideType
    compute_type(const double val, const double alpha, const double beta);

    int size() const { return m_size; }
    void set_size(const int size);
    const HeavisideType& tangent_type(const int i) const
    {
        return tangent_types[i];
    }
    const HeavisideType& normal_type(const int i) const
    {
        return normal_types[i];
    }
    HeavisideType& tangent_type(const int i) { return tangent_types[i]; }
    HeavisideType& normal_type(const int i) { return normal_types[i]; }

    int m_size = 0;
    std::vector<HeavisideType> tangent_types, normal_types;
};

} // namespace ipc
