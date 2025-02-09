import numpy as np

from find_ipctk import ipctk
from ipctk.filib import Interval

from utils import load_mesh


def test_ccd():
    V0, E, F = load_mesh("two-cubes-close.ply")
    V1, E, F = load_mesh("two-cubes-intersecting.ply")

    mesh = ipctk.CollisionMesh(V0, E, F)

    assert not ipctk.has_intersections(mesh, V0)
    assert ipctk.has_intersections(mesh, V1)

    assert not ipctk.is_step_collision_free(mesh, V0, V1)

    toi = ipctk.compute_collision_free_stepsize(mesh, V0, V1)
    assert 0 <= toi <= 1

    assert ipctk.is_step_collision_free(mesh, V0, (V1 - V0) * toi + V0)

def test_custom_ccd():
    class DumbNarrowPhaseCCD(ipctk.NarrowPhaseCCD):
        def __init__(self):
            ipctk.NarrowPhaseCCD.__init__(self)
        def point_point_ccd(self, p0_t0, p1_t0, p0_t1, p1_t1, min_distance=0.0, tmax=1.0):
            return True, 0.0
        def point_edge_ccd(self, p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, min_distance=0.0, tmax=1.0):
            return True, 0.0
        def point_triangle_ccd(self, p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, min_distance=0.0, tmax=1.0):
            return True, 0.0
        def edge_edge_ccd(self, ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, min_distance=0.0, tmax=1.0):
            return True, 0.0

    V0, E, F = load_mesh("two-cubes-close.ply")
    V1, E, F = load_mesh("two-cubes-intersecting.ply")

    mesh = ipctk.CollisionMesh(V0, E, F)

    ipctk.set_num_threads(1)

    assert not ipctk.is_step_collision_free(
        mesh, V0, V1, narrow_phase_ccd=DumbNarrowPhaseCCD())

    toi = ipctk.compute_collision_free_stepsize(
        mesh, V0, V1, narrow_phase_ccd=DumbNarrowPhaseCCD())
    assert 0 <= toi <= 1

def test_custom_broad_phase():
    class DumbBroadPhase(ipctk.BroadPhase):
        def __init__(self):
            ipctk.BroadPhase.__init__(self)
        def name(self):
            return "DumbBroadPhase"
        def detect_vertex_vertex_candidates(self):
            return []
        def detect_edge_vertex_candidates(self):
            return []
        def detect_edge_edge_candidates(self):
            return []
        def detect_face_vertex_candidates(self):
            return []
        def detect_edge_face_candidates(self):
            return []
        def detect_face_face_candidates(self):
            return []

    V0, E, F = load_mesh("two-cubes-close.ply")
    V1, E, F = load_mesh("two-cubes-intersecting.ply")

    mesh = ipctk.CollisionMesh(V0, E, F)

    ipctk.set_num_threads(1)

    assert ipctk.is_step_collision_free(mesh, V0, V1, broad_phase=DumbBroadPhase())

    toi = ipctk.compute_collision_free_stepsize(mesh, V0, V1, broad_phase=DumbBroadPhase())
    assert toi == 1

    assert not ipctk.has_intersections(mesh, V0, broad_phase=DumbBroadPhase())

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

    # class IntervalLinearTrajectory(ipctk.IntervalNonlinearTrajectory):
    #     def __init__(self, p0, p1):
    #         ipctk.IntervalNonlinearTrajectory.__init__(self)
    #         self.p0 = p0
    #         self.p1 = p1

    #     def __call__(self, t):
    #         r = ((self.p1 - self.p0) * t + self.p0)
    #         if isinstance(t, Interval):
    #             r = np.array([(ri.INF, ri.SUP) for ri in r], dtype="f8,f8")
    #         return r

    # p0 = IntervalLinearTrajectory(np.array([-1, 0, 0]), np.array([1, 0, 0]))
    # p1 = IntervalLinearTrajectory(np.array([0, -1, 0]), np.array([0, 1, 0]))

    # is_colliding, toi = ipctk.point_point_nonlinear_ccd(p0, p1)

    # assert is_colliding
    # assert 0 < toi < 1
    # assert np.abs(toi - 0.5) < 1e-2

    # BEGIN_RIGID_2D_TRAJECTORY
    class Rigid2DTrajectory(ipctk.NonlinearTrajectory):
        def __init__(self, position, translation, delta_translation, rotation, delta_rotation):
            ipctk.NonlinearTrajectory.__init__(self)
            self.position = position
            self.translation = translation
            self.delta_translation = delta_translation
            self.rotation = rotation
            self.delta_rotation = delta_rotation

        def __call__(self, t):
            theta = self.rotation + t * self.delta_rotation
            R = np.array([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta), np.cos(theta)]])
            return R @ self.position + self.translation + t * self.delta_translation

        def max_distance_from_linear(self, t0, t1):
            if self.delta_rotation * (t1 - t0) >= 2 * np.pi:
                # This is the most conservative estimate
                return 2 * np.linalg.norm(self.position)  # 2 * radius
            p_t0 = self(t0)
            p_t1 = self(t1)
            return np.linalg.norm(self((t0 + t1) / 2) - ((p_t1 - p_t0) * 0.5 + p_t0))
    # END_RIGID_2D_TRAJECTORY

    # BEGIN_TEST_RIGID_2D_TRAJECTORY
    p = Rigid2DTrajectory(
        np.array([0, 0.5]), np.zeros(2), np.zeros(2), 0, 0)
    e0 = Rigid2DTrajectory(
        np.array([-1, 0]), np.zeros(2), np.zeros(2), 0, np.pi)
    e1 = Rigid2DTrajectory(
        np.array([1, 0]), np.zeros(2), np.zeros(2), 0, np.pi)

    # increase the conservative_rescaling from 0.8 to 0.9 to get a more accurate estimate
    collision, toi = ipctk.point_edge_nonlinear_ccd(
        p, e0, e1, conservative_rescaling=0.9)

    assert collision
    assert 0.49 <= toi <= 0.5  # conservative estimate
    # END_TEST_RIGID_2D_TRAJECTORY
