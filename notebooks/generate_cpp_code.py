import numpy as np
import sympy
from sympy import Matrix, MatrixSymbol
from sympy.printing import ccode
import subprocess

from utils import jacobian


def generate_code(expr, out_var_name=None):
    CSE_results = sympy.cse(
        expr, sympy.numbered_symbols("t"), optimizations='basic')
    lines = []
    for helper in CSE_results[0]:
        if isinstance(helper[1], MatrixSymbol):
            lines.append(f'const auto {helper[0]}[{helper[1].size}];')
            lines.append(sympy.cxxcode(helper[1], helper[0]))
        else:
            lines.append(
                f'const auto {sympy.cxxcode(helper[0])} = {sympy.cxxcode(helper[1])};')

    if out_var_name != None:
        for i, result in enumerate(CSE_results[1]):
            lines.append(sympy.cxxcode(result, out_var_name))
    else:
        for i, result in enumerate(CSE_results[1]):
            lines.append(f"return {sympy.cxxcode(result)};")

    return '\n'.join(lines)


class CXXFunctionGenerator:
    def __init__(self, expr, params, name):
        self.expr = expr
        self.params = params
        self.name = name
        self.comment = ""
        self.out_param = None

    def signature(self):
        params = ", ".join(f"double {ccode(var)}" for var in self.params)
        if self.out_param != None:
            params += f", double {self.out_param}"
        ret_type = "double" if self.out_param == None else "void"
        return f"{self.comment}\n{ret_type} {self.name}({params});"

    def __call__(self):
        signature = self.signature()[:-1]  # remove semicolon
        return f"""
{signature}{{
{generate_code(self.expr, None if self.out_param ==
               None else self.out_param.split("[")[0])}
}}
"""


class CXXGradientGenerator(CXXFunctionGenerator):
    def __init__(self, expr, params, name, out_param_name="grad"):
        super().__init__(expr, params, name)
        self.expr = self.expr.diff(Matrix(params))
        self.out_param = f"{out_param_name}[{np.prod(self.expr.shape)}]"


class CXXJacobianGenerator(CXXFunctionGenerator):
    def __init__(self, expr, params, name, out_param_name="J"):
        super().__init__(expr, params, name)
        self.expr = jacobian(expr, params)
        J_shape = self.expr.shape
        self.expr = Matrix(np.array(self.expr).flatten(order="F"))
        self.out_param = f"{out_param_name}[{np.prod(self.expr.shape)}]"
        self.comment = f"// {out_param_name} is ({J_shape[0]}×{J_shape[1]}) flattened in column-major order"


class CXXHessianGenerator(CXXJacobianGenerator):
    def __init__(self, expr, params, name, out_param_name="hess"):
        super().__init__(expr.diff(Matrix(params)), params, name, out_param_name)


def generate_hpp_file(code_generators, file_name, transformer=None):
    newline = "\n"
    with open(file_name, 'w') as f:
        f.write(f"""\
#include <{file_name[:-4]}.hpp>

namespace ipc{{
namespace autogen{{
    {newline.join(code_generator.signature()
                  for code_generator in code_generators)}
}}
}}
""")
    subprocess.run(["clang-format", str(file_name), "-i"])


def generate_cpp_file(code_generators, file_name):
    newline = "\n"
    with open(file_name, 'w') as f:
        f.write(f"""\
#include <{file_name[:-4]}.hpp>

namespace ipc{{
namespace autogen{{
    {newline.join(code_generator() for code_generator in code_generators)}
}}
}}
""")
    subprocess.run(["clang-format", str(file_name), "-i"])
