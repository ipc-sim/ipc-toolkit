"""PyTorch autograd wrappers for IPC Toolkit.

This is a thin, pure-Python layer over the ``ipctk`` bindings. We don't
re-derive anything here: the forward pass calls the library's value, the
backward pass calls its analytic gradient, and the double backward calls its
analytic Hessian. That means the derivatives you get through autograd are the
same ones the C++ library uses, and second-order tools like
``torch.autograd.functional.hessian`` work out of the box.

Example:
    >>> import ipctk, ipctk.torch, torch
    >>> mesh = ipctk.CollisionMesh(V, E, F)
    >>> X = torch.tensor(V, requires_grad=True)
    >>> collisions = ipctk.NormalCollisions()
    >>> collisions.build(mesh, X.detach().numpy(), dhat)
    >>> B = ipctk.torch.Potential(ipctk.BarrierPotential(dhat, kappa))
    >>> B(collisions, mesh, X).backward()  # X.grad is the barrier gradient

Two things are treated as constants, so they never receive a gradient:

* Topology. The ``CollisionMesh`` (edges, faces, DOF maps) is discrete.
* The active set. We build ``collisions`` outside of autograd and hold it fixed
  across forward and backward. This is the same assumption a Newton solver
  makes within one iteration. The barrier goes to zero smoothly at ``dhat``, so
  a pair entering or leaving the set does not cause a jump in the energy, but
  you should rebuild the collisions whenever the positions change.

All computation happens in float64 on the CPU. We copy CUDA tensors to the host
and cast float32 tensors to float64 on the way in, and cast the results back to
the input's dtype and device on the way out. For float64 CPU tensors the
conversion to NumPy is zero-copy.

The library assembles Hessians in parallel, so second derivatives of a
``Potential`` can differ at round-off level between two identical calls. When
running ``torch.autograd.gradgradcheck``, pass a small ``nondet_tol``.
"""

import numpy as np

try:
    import torch
except ImportError as e:
    raise ImportError(
        "ipctk.torch requires PyTorch. Install it with `pip install torch` "
        "or `pip install ipctk[torch]`."
    ) from e
from torch.autograd.function import once_differentiable

import ipctk

__all__ = [
    "Potential",
    "point_point_distance",
    "point_edge_distance",
    "edge_edge_distance",
    "point_triangle_distance",
]


def _to_numpy(x):
    return x.detach().cpu().numpy().astype(np.float64, copy=False)


def _like(array, x):
    """Convert a NumPy array to a tensor with the dtype and device of x."""
    # We go through np.asarray first because the bindings return plain Python
    # floats for scalar values, and torch.as_tensor would round those to the
    # default float32 before we ever get a chance to cast them.
    return torch.as_tensor(np.asarray(array), dtype=x.dtype, device=x.device)


# -----------------------------------------------------------------------------
# Potentials


class _PotentialValue(torch.autograd.Function):
    @staticmethod
    def forward(ctx, X, potential, collisions, mesh):
        ctx.save_for_backward(X)
        ctx.args = (potential, collisions, mesh)
        return _like(potential(collisions, mesh, _to_numpy(X)), X)

    @staticmethod
    def backward(ctx, grad_out):
        (X,) = ctx.saved_tensors
        # We route the gradient through a second Function instead of computing
        # it inline. When autograd is asked for higher-order derivatives
        # (create_graph=True), this records a node whose own backward is the
        # Hessian-vector product.
        grad = _PotentialGradient.apply(X, *ctx.args)
        return grad_out * grad, None, None, None


class _PotentialGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, X, potential, collisions, mesh):
        ctx.save_for_backward(X)
        ctx.args = (potential, collisions, mesh)
        # The library flattens X row by row (x0, y0, z0, x1, ...), which
        # matches NumPy's default C-order reshape.
        grad = potential.gradient(collisions, mesh, _to_numpy(X))
        return _like(grad.reshape(X.shape), X)

    @staticmethod
    @once_differentiable
    def backward(ctx, grad_out):
        (X,) = ctx.saved_tensors
        potential, collisions, mesh = ctx.args
        H = potential.hessian(collisions, mesh, _to_numpy(X))
        # H is symmetric, so the vector-Jacobian product H^T v is just H v.
        Hv = H @ _to_numpy(grad_out).reshape(-1)
        return _like(Hv.reshape(X.shape), X), None, None, None


class Potential:
    """Make an ipctk potential differentiable with PyTorch autograd.

    This works with any potential that shares the ``ipctk`` potential
    interface, i.e., ``BarrierPotential``, ``NormalAdhesionPotential``,
    ``FrictionPotential``, and ``TangentialAdhesionPotential``.

    Parameters:
        potential: The ipctk potential to wrap.
    """

    def __init__(self, potential):
        self.potential = potential

    def __call__(self, collisions, mesh, X):
        """Compute the potential as a differentiable scalar tensor.

        Parameters:
            collisions: The (fixed) set of collisions.
            mesh: The collision mesh.
            X: Tensor of shape (num_vertices, dim) holding the degrees of
                freedom of the collision mesh (e.g., vertices or velocities).

        Returns:
            A 0-d tensor with the dtype and device of X.
        """
        return _PotentialValue.apply(X, self.potential, collisions, mesh)

    def hessian(
        self,
        collisions,
        mesh,
        X,
        project_hessian_to_psd=ipctk.PSDProjectionMethod.NONE,
    ):
        """Assemble the Hessian as a sparse COO tensor.

        Autograd can recover the Hessian too, but it does so one dense
        Hessian-vector product at a time. A Newton solver wants the sparse
        matrix directly, optionally projected to be positive semi-definite.

        Returns:
            A sparse tensor of shape (X.numel(), X.numel()). This is not part
            of the autograd graph.
        """
        H = self.potential.hessian(
            collisions, mesh, _to_numpy(X), project_hessian_to_psd
        ).tocoo()
        indices = np.vstack([H.row, H.col]).astype(np.int64)
        return torch.sparse_coo_tensor(
            torch.from_numpy(indices),
            torch.from_numpy(H.data),
            H.shape,
            dtype=X.dtype,
            device=X.device,
        ).coalesce()


# -----------------------------------------------------------------------------
# Distances


class _DistanceFunctions:
    """The value, gradient, and Hessian bindings for one distance type."""

    def __init__(self, name):
        self.value = getattr(ipctk, f"{name}_distance")
        self.gradient = getattr(ipctk, f"{name}_distance_gradient")
        self.hessian = getattr(ipctk, f"{name}_distance_hessian")


def _stack_points(points):
    """Stack k tensors of shape (..., dim) into an array of (B, k, dim)."""
    P = np.stack([_to_numpy(p) for p in points], axis=-2)
    return P.reshape(-1, *P.shape[-2:])


class _DistanceValue(torch.autograd.Function):
    @staticmethod
    def forward(ctx, fns, *points):
        ctx.save_for_backward(*points)
        ctx.fns = fns
        # The bindings take a single stencil, so we loop over the batch. This
        # keeps the prototype simple; a batched binding would remove the
        # per-element Python overhead.
        d = np.array([fns.value(*stencil) for stencil in _stack_points(points)])
        return _like(d.reshape(points[0].shape[:-1]), points[0])

    @staticmethod
    def backward(ctx, grad_out):
        points = ctx.saved_tensors
        grad = _DistanceGradient.apply(ctx.fns, *points)  # (..., k, dim)
        grad = grad_out[..., None, None] * grad
        return None, *grad.unbind(-2)


class _DistanceGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, fns, *points):
        ctx.save_for_backward(*points)
        ctx.fns = fns
        P = _stack_points(points)
        g = np.stack([fns.gradient(*stencil) for stencil in P])
        shape = (*points[0].shape[:-1], *P.shape[-2:])
        return _like(g.reshape(shape), points[0])

    @staticmethod
    @once_differentiable
    def backward(ctx, grad_out):
        points = ctx.saved_tensors
        P = _stack_points(points)
        v = _to_numpy(grad_out).reshape(P.shape[0], -1)
        Hv = np.stack(
            [ctx.fns.hessian(*stencil) @ vi for stencil, vi in zip(P, v)]
        )
        Hv = _like(Hv.reshape(grad_out.shape), points[0])
        return None, *Hv.unbind(-2)


def _distance(fns, *points):
    # We broadcast here, before entering the Function, so that autograd itself
    # sums the gradient over any broadcast dimensions. For example, the
    # gradient of one triangle against a batch of points accumulates
    # correctly without any reduction logic in our backward.
    return _DistanceValue.apply(fns, *torch.broadcast_tensors(*points))


_POINT_POINT = _DistanceFunctions("point_point")
_POINT_EDGE = _DistanceFunctions("point_edge")
_EDGE_EDGE = _DistanceFunctions("edge_edge")
_POINT_TRIANGLE = _DistanceFunctions("point_triangle")


def point_point_distance(p0, p1):
    """Squared distance between two points.

    Inputs have shape (..., dim) and are broadcast against each other.
    """
    return _distance(_POINT_POINT, p0, p1)


def point_edge_distance(p, e0, e1):
    """Squared distance between a point and an edge.

    Inputs have shape (..., dim) and are broadcast against each other.
    """
    return _distance(_POINT_EDGE, p, e0, e1)


def edge_edge_distance(ea0, ea1, eb0, eb1):
    """Squared distance between two 3D edges.

    Inputs have shape (..., 3) and are broadcast against each other.

    Note:
        This is not smooth when the edges are parallel. IPC handles that case
        with a mollifier inside the barrier potential, so prefer
        ``Potential`` over this function when building contact energies.
    """
    return _distance(_EDGE_EDGE, ea0, ea1, eb0, eb1)


def point_triangle_distance(p, t0, t1, t2):
    """Squared distance between a point and a 3D triangle.

    Inputs have shape (..., 3) and are broadcast against each other.
    """
    return _distance(_POINT_TRIANGLE, p, t0, t1, t2)
