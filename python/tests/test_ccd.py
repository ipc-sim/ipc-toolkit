import numpy as np
from utils import *
import ipctk
from ipctk.filib import Interval


def test_ccd():
    V0, E, F = load_mesh("two-cubes-close.obj")
    V1, E, F = load_mesh("two-cubes-intersecting.obj")

    mesh = ipctk.CollisionMesh(V0, E, F)

    assert not ipctk.has_intersections(mesh, V0)
    assert ipctk.has_intersections(mesh, V1)

    assert not ipctk.is_step_collision_free(mesh, V0, V1)

    toi = ipctk.compute_collision_free_stepsize(mesh, V0, V1)
    assert 0 < toi < 1

    assert ipctk.is_step_collision_free(mesh, V0, (V1 - V0) * toi + V0)


def test_nonlinear_ccd():
    class LinearTrajectory(ipctk.NonlinearTrajectory):
        def __init__(self, p0, p1):
            ipctk.NonlinearTrajectory.__init__(self)
            self.p0 = p0
            self.p1 = p1

        def __call__(self, t: float):
            return (self.p1 - self.p0) * t + self.p0

        def max_distance_from_linear(self, t0: float, t1: float):
            return 0

    p0 = LinearTrajectory(np.array([-1, 0, 0]), np.array([1, 0, 0]))
    p1 = LinearTrajectory(np.array([0, -1, 0]), np.array([0, 1, 0]))

    assert np.all(p0(0) == np.array([-1, 0, 0]))
    assert np.all(p0(1) == np.array([1, 0, 0]))
    assert np.all(p1(0) == np.array([0, -1, 0]))
    assert np.all(p1(1) == np.array([0, 1, 0]))

    assert np.all(p0(0.5) == np.array([0, 0, 0]))
    assert np.all(p1(0.5) == np.array([0, 0, 0]))

    assert p0.max_distance_from_linear(0, 1) == 0
    assert p1.max_distance_from_linear(0, 1) == 0

    is_colliding, toi = ipctk.point_point_nonlinear_ccd(p0, p1)

    assert is_colliding
    assert 0 < toi < 1
    assert np.abs(toi - 0.5) < 1e-6

    class IntervalLinearTrajectory(ipctk.IntervalNonlinearTrajectory):
        def __init__(self, p0, p1):
            ipctk.IntervalNonlinearTrajectory.__init__(self)
            self.p0 = p0
            self.p1 = p1

        def __call__(self, t):
            r = ((self.p1 - self.p0) * t + self.p0)
            if isinstance(t, Interval):
                r = np.array([(ri.INF, ri.SUP) for ri in r], dtype="f8,f8")
            return r

    p0 = IntervalLinearTrajectory(np.array([-1, 0, 0]), np.array([1, 0, 0]))
    p1 = IntervalLinearTrajectory(np.array([0, -1, 0]), np.array([0, 1, 0]))

    is_colliding, toi = ipctk.point_point_nonlinear_ccd(p0, p1)

    assert is_colliding
    assert 0 < toi < 1
    assert np.abs(toi - 0.5) < 1e-2
