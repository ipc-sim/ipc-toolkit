#include <common.hpp>
#include <pybind11/numpy.h>

#include <ipc/utils/interval.hpp>

namespace py = pybind11;
using namespace filib;

// ============================================================================

#define RETURN_INTERVAL(func)                                                  \
    [](const Interval& i) -> Interval { return func(i); }

#define OP_II(op)                                                              \
    [](const Interval& a, const Interval& b) -> Interval { return op; }
#define OP_ID(op) [](const Interval& a, double b) -> Interval { return op; }
#define OP_DI(op) [](const Interval& b, double a) -> Interval { return op; }

#define DEF_OP(name, op) def("__" name "__", op, py::is_operator())
#define DEF_BIN_OP(name, op)                                                   \
    DEF_OP(name, OP_II(op)).DEF_OP(name, OP_ID(op)).DEF_OP("r" name, OP_DI(op))
#define DEF_UN_OP(name, op)                                                    \
    DEF_OP(name, [](const Interval& i) -> Interval { return op; })

#define DEF_IOP(name, RHSType, op)                                             \
    def(                                                                       \
        "__i" name "__",                                                       \
        [](Interval& self, RHSType rhs) -> Interval& {                         \
            op;                                                                \
            return self;                                                       \
        },                                                                     \
        py::is_operator())

// ============================================================================

void define_interval(py::module_& m)
{
#ifdef IPC_TOOLKIT_WITH_FILIB
    auto m_filib = m.def_submodule("filib", "Fast Interval Library");

    py::class_<Interval>(m_filib, "Interval")
        .def(py::init())
        .def(py::init<double>(), py::arg("x"))
        .def(py::init<double, double>(), py::arg("x"), py::arg("y"))
        .def(
            "__str__",
            [](const Interval& self) {
                std::stringstream ss;
                ss << self;
                return ss.str();
            })
        .def(
            "__repr__",
            [](const Interval& self) {
                std::stringstream ss;
                ss << "Interval(" << self.INF << ", " << self.SUP << ")";
                return ss.str();
            })
        .DEF_BIN_OP("add", a + b)
        .DEF_UN_OP("pos", +i)
        .DEF_IOP("add", const Interval&, self += rhs)
        .DEF_IOP("add", double, self += rhs)
        .DEF_BIN_OP("sub", a - b)
        .DEF_UN_OP("neg", -i)
        .DEF_IOP("sub", const Interval&, self -= rhs)
        .DEF_IOP("sub", double, self -= rhs)
        // .DEF_BIN_OP("mul", a * b)
        .DEF_OP(
            "mul",
            [](const Interval& a, const Interval& b) -> Interval {
                if (a == b) {
                    return sqr(a); // More accurate
                }
                return a * b;
            })
        .DEF_OP("mul", OP_ID(a * b))
        .DEF_OP("rmul", OP_DI(a * b))
        .DEF_IOP("mul", const Interval&, self *= rhs)
        .DEF_IOP("mul", double, self *= rhs)
        .DEF_BIN_OP("truediv", a / b)
        .DEF_IOP("truediv", const Interval&, self /= rhs)
        .DEF_IOP("truediv", double, self /= rhs)
        .def(py::self == py::self)
        .def(py::self == double())
        .def(double() == py::self)
        .def(py::self != py::self)
        .def(py::self != double())
        .def(double() != py::self)
        .def(py::self < py::self)
        .def(py::self < double())
        .def(double() < py::self)
        .def(py::self <= py::self)
        .def(py::self <= double())
        .def(double() <= py::self)
        .def(py::self > py::self)
        .def(py::self > double())
        .def(double() > py::self)
        .def(py::self >= py::self)
        .def(py::self >= double())
        .def(double() >= py::self)
        .def(
            "__contains__",
            [](const Interval& self, double x) { return in(x, self); },
            py::is_operator())
        .def(
            "__contains__",
            [](const Interval& self, const Interval& i) { return in(i, self); },
            py::is_operator())
        .def("empty", [](const Interval& self) { return empty(self); })
        .DEF_OP("or", OP_II(a | b))
        .DEF_OP("and", OP_II(a & b))
        .def("mid", [](const Interval& self) { return mid(self); })
        .def("diam", [](const Interval& self) { return diam(self); })
        .def("drel", [](const Interval& self) { return drel(self); })
        .def(
            "blow",
            [](const Interval& self, double eps) -> Interval {
                return blow(self, eps);
            })
        .def("sqr", RETURN_INTERVAL(sqr))
        .def(
            "__pow__",
            [](const Interval& self, double x) -> Interval {
                if (x != 2.0) {
                    throw std::invalid_argument(
                        "Only x == 2.0 is supported for Interval**double");
                }
                return sqr(self);
            },
            py::is_operator())
        .def("sqrt", RETURN_INTERVAL(sqrt))
        .def("exp", RETURN_INTERVAL(exp))
        .def("exp2", RETURN_INTERVAL(exp2))
        .def("exp10", RETURN_INTERVAL(exp10))
        .def(
            "__rpow__",
            [](const Interval& self, double b) -> Interval {
                if (b == 2.0) {
                    return sqr(self);
                } else if (b == exp(1.0)) {
                    return exp(self);
                } else if (b == 10.0) {
                    return exp10(self);
                } else {
                    throw std::invalid_argument(
                        "Only b in [2.0, e, 10.0] are supported for double**Interval");
                }
            },
            py::is_operator())
        .def("expm1", RETURN_INTERVAL(expm1))
        .def("log", RETURN_INTERVAL(log))
        .def("log2", RETURN_INTERVAL(log2))
        .def("log10", RETURN_INTERVAL(log10))
        .def("log1p", RETURN_INTERVAL(log1p))
        .def("sin", RETURN_INTERVAL(sin))
        .def("cos", RETURN_INTERVAL(cos))
        .def("cot", RETURN_INTERVAL(cot))
        .def("tan", RETURN_INTERVAL(tan))
        .def("asin", RETURN_INTERVAL(asin))
        .def("acos", RETURN_INTERVAL(acos))
        .def("atan", RETURN_INTERVAL(atan))
        .def("acot", RETURN_INTERVAL(acot))
        .def("arcsin", RETURN_INTERVAL(asin))
        .def("arccos", RETURN_INTERVAL(acos))
        .def("arctan", RETURN_INTERVAL(atan))
        .def("arccot", RETURN_INTERVAL(acot))
        .def("sinh", RETURN_INTERVAL(sinh))
        .def("cosh", RETURN_INTERVAL(cosh))
        .def("tanh", RETURN_INTERVAL(tanh))
        .def("coth", RETURN_INTERVAL(coth))
        .def("asinh", RETURN_INTERVAL(asinh))
        .def("acosh", RETURN_INTERVAL(acosh))
        .def("atanh", RETURN_INTERVAL(atanh))
        .def("acoth", RETURN_INTERVAL(acoth))
        .def("arcsinh", RETURN_INTERVAL(asinh))
        .def("arccosh", RETURN_INTERVAL(acosh))
        .def("arctanh", RETURN_INTERVAL(atanh))
        .def("arccoth", RETURN_INTERVAL(acoth))
        .def("erf", RETURN_INTERVAL(erf))
        .def("erfc", RETURN_INTERVAL(erfc))
        .def_readwrite("INF", &interval::INF)
        .def_readwrite("SUP", &interval::SUP);

    m_filib.def("disjoint", [](const Interval& i, const Interval& j) {
        return disjoint(i, j);
    });
    m_filib.def("max", [](const Interval& i, const Interval& j) -> Interval {
        return max(i, j);
    });
    m_filib.def("max", [](const Interval& i, double j) -> Interval {
        return max(i, j);
    });
    m_filib.def("max", [](double i, const Interval& j) -> Interval {
        return max(i, j);
    });
    m_filib.def("min", [](const Interval& i, const Interval& j) -> Interval {
        return min(i, j);
    });
    m_filib.def("min", [](const Interval& i, double j) -> Interval {
        return min(i, j);
    });
    m_filib.def("min", [](double i, const Interval& j) -> Interval {
        return min(i, j);
    });
    m_filib.def("sqr", RETURN_INTERVAL(sqr));
    m_filib.def("sqrt", RETURN_INTERVAL(sqrt));
    m_filib.def("exp", RETURN_INTERVAL(exp));
    m_filib.def("exp2", RETURN_INTERVAL(exp2));
    m_filib.def("exp10", RETURN_INTERVAL(exp10));
    m_filib.def("expm1", RETURN_INTERVAL(expm1));
    m_filib.def("log", RETURN_INTERVAL(log));
    m_filib.def("log2", RETURN_INTERVAL(log2));
    m_filib.def("log10", RETURN_INTERVAL(log10));
    m_filib.def("log1p", RETURN_INTERVAL(log1p));
    m_filib.def("sin", RETURN_INTERVAL(sin));
    m_filib.def("cos", RETURN_INTERVAL(cos));
    m_filib.def("cot", RETURN_INTERVAL(cot));
    m_filib.def("tan", RETURN_INTERVAL(tan));
    m_filib.def("asin", RETURN_INTERVAL(asin));
    m_filib.def("acos", RETURN_INTERVAL(acos));
    m_filib.def("atan", RETURN_INTERVAL(atan));
    m_filib.def("acot", RETURN_INTERVAL(acot));
    m_filib.def("sinh", RETURN_INTERVAL(sinh));
    m_filib.def("cosh", RETURN_INTERVAL(cosh));
    m_filib.def("coth", RETURN_INTERVAL(coth));
    m_filib.def("tanh", RETURN_INTERVAL(tanh));
    m_filib.def("asinh", RETURN_INTERVAL(asinh));
    m_filib.def("acosh", RETURN_INTERVAL(acosh));
    m_filib.def("acoth", RETURN_INTERVAL(acoth));
    m_filib.def("atanh", RETURN_INTERVAL(atanh));
    m_filib.def("erf", RETURN_INTERVAL(erf));
    m_filib.def("erfc", RETURN_INTERVAL(erfc));

    PYBIND11_NUMPY_DTYPE(Interval, INF, SUP);
#endif
}