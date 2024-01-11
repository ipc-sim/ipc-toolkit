#pragma once

#include "unordered_map_and_set.hpp"
#include <array>
#include <ipc/config.hpp>

namespace ipc {

class unordered_tuple {
public:
    unordered_tuple(const long &_a, const long &_b)
    {
        x = {{std::min(_a, _b), std::max(_a, _b)}};
    }

    const std::array<long, 2> &get() const { return x; }

    long operator[](const int &idx) const
    {
        return x[idx];
    }

    template <typename H>
    friend H AbslHashValue(H h, const unordered_tuple& other)
    {
        return H::combine(std::move(h), other.x[0], other.x[1]);
    }

    bool operator==(const unordered_tuple& other) const
    {
        return x[0] == other.x[0] && x[1] == other.x[1];
    }
    bool operator!=(const unordered_tuple& other) const
    {
        return !(*this == other);
    }
private:
    std::array<long, 2> x;
};

}