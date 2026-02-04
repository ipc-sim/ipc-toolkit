import numpy as np
import sympy


def sq_norm(x):
    """Compute the squared norm of a vector x."""
    return x.dot(x)


def norm(x):
    """Compute the norm of a vector x."""
    return sympy.sqrt(sq_norm(x))


def normalize(x):
    """Normalize a vector x."""
    return x / norm(x)


def jacobian(F, x):
    """Compute the Jacobian of a matrix-valued function F with respect to vector x."""
    J = np.empty((F.size, x.size), dtype=object)
    for xi in range(x.size):
        # Flatten column-major order
        J[:, xi] = np.array(sympy.Matrix(F).diff(x[xi])).reshape(F.size, order="F")
    return J
