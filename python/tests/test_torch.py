"""Tests for the PyTorch autograd wrappers in ipctk.torch.

PyTorch is an optional dependency, so every test skips when it is missing.
"""

import unittest

import numpy as np
from find_ipctk import ipctk
from utils import load_mesh

try:
    import torch

    import ipctk.torch
except ImportError:
    torch = None

# two-cubes-close.ply has a ~0.069 gap between the cubes, so dhat=0.1 activates
# a few hundred collisions while staying small enough for finite differences.
DHAT = 0.1


def require_torch():
    if torch is None:
        raise unittest.SkipTest("PyTorch is not installed")


def two_cubes():
    V, E, F = load_mesh("two-cubes-close.ply")
    mesh = ipctk.CollisionMesh(V, E, F)
    collisions = ipctk.NormalCollisions()
    collisions.build(mesh, mesh.rest_positions, DHAT)
    X = torch.tensor(mesh.rest_positions, requires_grad=True)
    return mesh, collisions, X


def potentials():
    yield ipctk.torch.Potential(ipctk.BarrierPotential(DHAT, 1.0))
    yield ipctk.torch.Potential(
        ipctk.NormalAdhesionPotential(DHAT / 2, DHAT, 1.0, 0.01))


def distances():
    yield ipctk.torch.point_point_distance, 2
    yield ipctk.torch.point_edge_distance, 3
    yield ipctk.torch.edge_edge_distance, 4
    yield ipctk.torch.point_triangle_distance, 4


def test_potential_derivatives_match_finite_differences():
    require_torch()
    mesh, collisions, X = two_cubes()
    for P in potentials():
        def f(X):
            return P(collisions, mesh, X)

        # The library assembles the Hessian in parallel, so the summation
        # order (and hence the last few bits) can differ between two calls.
        # We allow round-off level nondeterminism instead of requiring
        # bitwise-identical double backward passes.
        assert torch.autograd.gradcheck(f, (X,), nondet_tol=1e-12)
        assert torch.autograd.gradgradcheck(f, (X,), nondet_tol=1e-12)


def test_sparse_hessian_matches_autograd_hessian():
    require_torch()
    mesh, collisions, X = two_cubes()
    P = ipctk.torch.Potential(ipctk.BarrierPotential(DHAT, 1.0))

    H = P.hessian(collisions, mesh, X)
    H_autograd = torch.autograd.functional.hessian(
        lambda x: P(collisions, mesh, x.reshape(X.shape)),
        X.detach().reshape(-1))

    assert H.shape == (X.numel(), X.numel())
    assert torch.allclose(H.to_dense(), H_autograd)


def test_barrier_gradient_descent_separates_surfaces():
    require_torch()
    mesh, collisions, X = two_cubes()
    P = ipctk.torch.Potential(ipctk.BarrierPotential(DHAT, 1.0))

    P(collisions, mesh, X).backward()
    X_new = (X - 1e-4 * X.grad).detach().numpy()

    d_before = collisions.compute_minimum_distance(mesh, mesh.rest_positions)
    d_after = collisions.compute_minimum_distance(mesh, X_new)
    assert d_after > d_before


def test_distance_derivatives_match_finite_differences():
    require_torch()
    rng = np.random.default_rng(0)
    for distance, num_points in distances():
        points = tuple(
            torch.tensor(rng.random((4, 3)), requires_grad=True)
            for _ in range(num_points))

        assert distance(*points).shape == (4,)
        assert torch.autograd.gradcheck(distance, points)
        assert torch.autograd.gradgradcheck(distance, points)


def test_distance_broadcasts_shared_triangle_over_points():
    require_torch()
    rng = np.random.default_rng(1)
    p = torch.tensor(rng.random((5, 3)))
    t = [torch.tensor(rng.random(3), requires_grad=True) for _ in range(3)]

    d = ipctk.torch.point_triangle_distance(p, *t)
    d.sum().backward()

    expected = [ipctk.point_triangle_distance(pi, *(ti.detach().numpy() for ti in t))
                for pi in p.numpy()]
    assert np.allclose(d.detach().numpy(), expected)
    # Each triangle vertex's gradient accumulates over every point.
    expected_grad = sum(
        ipctk.point_triangle_distance_gradient(
            pi, *(ti.detach().numpy() for ti in t))
        for pi in p.numpy()).reshape(4, 3)
    for i, ti in enumerate(t):
        assert np.allclose(ti.grad.numpy(), expected_grad[i + 1])


def test_float32_inputs_stay_float32():
    require_torch()
    p0 = torch.tensor([0.0, 0.0, 0.0], requires_grad=True)
    p1 = torch.tensor([1.0, 2.0, 2.0], requires_grad=True)
    assert p0.dtype == torch.float32

    d = ipctk.torch.point_point_distance(p0, p1)
    d.backward()

    assert d.dtype == torch.float32
    assert d.item() == 9.0  # squared distance
    assert p0.grad.dtype == torch.float32
    assert torch.allclose(p0.grad, 2 * (p0 - p1).detach())
