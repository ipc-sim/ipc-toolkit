#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION

#include <gmp.h>
#include <iostream>

namespace ipc {

// https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/rational-gmp_8cpp-source.html
class Rational {
public:
    mpq_t value;
    void canonicalize() { mpq_canonicalize(value); }
    int get_sign() const { return mpq_sgn(value); }

    Rational()
    {
        mpq_init(value);
        mpq_set_d(value, 0);
    }

    Rational(double d)
    {
        mpq_init(value);
        mpq_set_d(value, d);
        // canonicalize();
    }

    Rational(const mpq_t& v_)
    {
        mpq_init(value);
        mpq_set(value, v_);
        // canonicalize();
    }

    Rational(const Rational& other)
    {
        mpq_init(value);
        mpq_set(value, other.value);
    }

    ~Rational() { mpq_clear(value); }

    friend Rational operator-(const Rational& v)
    {
        Rational r_out;
        mpq_neg(r_out.value, v.value);
        return r_out;
    }

    friend Rational operator+(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_add(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Rational operator-(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_sub(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Rational operator*(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_mul(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Rational operator/(const Rational& x, const Rational& y)
    {
        Rational r_out;
        mpq_div(r_out.value, x.value, y.value);
        return r_out;
    }

    Rational& operator=(const Rational& x)
    {
        if (this == &x)
            return *this;
        mpq_set(value, x.value);
        return *this;
    }

    Rational& operator=(const double x)
    {
        mpq_set_d(value, x);
        // canonicalize();
        return *this;
    }

    friend bool operator<(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) < 0;
    }

    friend bool operator>(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) > 0;
    }

    friend bool operator<=(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) <= 0;
    }

    friend bool operator>=(const Rational& r, const Rational& r1)
    {
        return mpq_cmp(r.value, r1.value) >= 0;
    }

    friend bool operator==(const Rational& r, const Rational& r1)
    {
        return mpq_equal(r.value, r1.value);
    }

    friend bool operator!=(const Rational& r, const Rational& r1)
    {
        return !mpq_equal(r.value, r1.value);
    }

    double to_double() { return mpq_get_d(value); }

    friend std::ostream& operator<<(std::ostream& os, const Rational& r)
    {
        os << mpq_get_d(r.value);
        return os;
    }
};

} // namespace ipc

#endif
