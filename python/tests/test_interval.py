import math
import numpy as np

import find_ipctk
from ipctk import filib
from ipctk.filib import Interval


def np_sqr(x): return x ** 2
def np_exp10(x): return 10 ** x
def np_cot(x): return 1 / np.tan(x)
def np_arccot(x): return np.pi / 2 - np.arctan(x)
def np_coth(x): return 1 / np.tanh(x)
def np_arccoth(x): return np.arctanh(1/x)


def test_interval():
    x, y = 2, 0.5
    xi, yi = Interval(x), Interval(y)

    assert x + y in xi + yi
    assert x + 2 in xi + 2
    assert 2 + x in 2 + xi

    assert x - y in xi - yi
    assert x - 2 in xi - 2
    assert 2 - x in 2 - xi

    assert x * y in xi * yi
    assert x * 2 in xi * 2
    assert 2 * x in 2 * xi

    assert x / y in xi / yi
    assert x / 2 in xi / 2
    assert 2 / x in 2 / xi

    tmp = Interval(xi.INF, xi.SUP)
    tmp += Interval(1)
    assert x + 1 in tmp
    tmp += 1
    assert x + 2 in tmp
    tmp -= Interval(1)
    assert x + 1 in tmp
    tmp -= 1
    assert x in tmp
    tmp *= Interval(2)
    assert 2 * x in tmp
    tmp *= 2
    assert 4 * x in tmp
    tmp /= Interval(2)
    assert 2 * x in tmp
    tmp /= 2
    assert x in tmp

    assert +x == +xi
    assert -x == -xi

    assert xi == xi
    assert xi == x
    assert x == xi

    assert xi != yi
    assert xi != y
    assert y != xi

    assert yi < xi
    assert 0 < xi
    assert yi < x

    assert yi <= xi
    assert 0 <= xi
    assert yi <= x

    assert xi > yi
    assert xi > y
    assert x > yi

    assert xi >= yi
    assert xi >= y
    assert x >= yi

    assert Interval(3, 2).empty()

    assert not (xi | yi).empty()
    assert (xi & yi).empty()
    assert filib.disjoint(xi, yi)

    assert x in xi
    assert x not in yi
    assert xi in xi.blow(1e-16)
    assert xi not in yi

    assert filib.max(xi, yi) == xi
    assert filib.max(xi, y) == xi
    assert filib.max(y, xi) == xi

    assert filib.min(xi, yi) == yi
    assert filib.min(yi, y) == yi
    assert filib.min(y, yi) == yi

    assert np_sqr(x) in filib.sqr(xi)
    assert np_sqr(x) in np_sqr(xi)
    assert np.sqrt(x) in np.sqrt(xi)

    assert np.exp(x) in np.exp(xi)
    assert np.e**x in np.e**xi
    assert np.exp2(x) in np.exp2(xi)
    assert 2**x in 2**xi
    assert np_exp10(x) in filib.exp10(xi)
    assert np_exp10(x) in np_exp10(xi)
    assert np.expm1(x) in np.expm1(xi)

    assert np.log(x) in np.log(xi)
    assert np.log2(x) in np.log2(xi)
    assert np.log10(x) in np.log10(xi)
    assert np.log1p(x) in np.log1p(xi)

    assert np.sin(x) in np.sin(xi)
    assert np.cos(x) in np.cos(xi)
    assert np.tan(x) in np.tan(xi)
    assert np_cot(x) in filib.cot(xi)
    assert np_cot(x) in np_cot(xi)

    assert np.arcsin(y) in np.arcsin(yi)
    assert np.arccos(y) in np.arccos(yi)
    assert np.arctan(y) in np.arctan(yi)
    assert np_arccot(y) in filib.acot(yi)
    assert np_arccot(y) in np_arccot(yi)

    assert np.sinh(x) in np.sinh(xi)
    assert np.cosh(x) in np.cosh(xi)
    assert np.tanh(x) in np.tanh(xi)
    assert np_coth(x) in filib.coth(xi)

    assert np.arcsinh(y) in np.arcsinh(yi)
    assert np.arccosh(x) in np.arccosh(xi)
    assert np.arctanh(y) in np.arctanh(yi)
    assert np_arccoth(x) in filib.acoth(xi)
    assert np_arccoth(x) in np_arccoth(xi)

    assert math.erf(x) in filib.erf(xi)
    assert math.erfc(x) in filib.erfc(xi)

    a = np.array([Interval(-1, 1), Interval(-1, 1), Interval(-1, 1)])
    assert np.linalg.norm(a).INF == 0
