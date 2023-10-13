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
    J = np.empty((x.size * F.shape[0], F.shape[1]), dtype=object)
    for xi in range(x.size):
        for Fi in range(F.shape[0]):
            for Fj in range(F.shape[1]):
                J[xi * F.shape[0] + Fi, Fj] = F[Fi, Fj].diff(x[xi])
    return J
