"""Tests for the potential bindings.

Covers the Python-side input validation (which the C++ API only enforces with
assert(), compiled out under NDEBUG) and the smooth-contact/GCP bindings.

Uses unittest.TestCase so that assertRaises is available under both nose2 (the
runner used in CI) and pytest, without adding a pytest dependency.
"""

import unittest

import numpy as np
from find_ipctk import ipctk
from utils import load_mesh

# two-cubes-close.ply has a ~0.069 gap between the cubes, so dhat=0.1 activates
# a few hundred collisions while staying fast (~1 ms to build).
DHAT = 0.1
STIFFNESS = 1.0


def two_cubes():
    V, E, F = load_mesh("two-cubes-close.ply")
    mesh = ipctk.CollisionMesh(V, E, F)
    return mesh, mesh.rest_positions


def normal_collisions(mesh, vertices, dhat=DHAT):
    collisions = ipctk.NormalCollisions()
    collisions.build(mesh, vertices, dhat)
    return collisions


class TestBarrierPotentialValidation(unittest.TestCase):
    """The C++ ctor/setters assert dhat > 0, stiffness > 0, barrier != nullptr.

    assert() is compiled out under NDEBUG, so the bindings must validate and
    raise ValueError instead of admitting undefined behavior in a release build.
    """

    INVALID = [0.0, -1.0, -1e-3, float("nan")]

    def test_ctor_rejects_invalid_dhat(self):
        for dhat in self.INVALID:
            with self.subTest(dhat=dhat):
                with self.assertRaises(ValueError):
                    ipctk.BarrierPotential(dhat, STIFFNESS)

    def test_ctor_rejects_invalid_stiffness(self):
        for stiffness in self.INVALID:
            with self.subTest(stiffness=stiffness):
                with self.assertRaises(ValueError):
                    ipctk.BarrierPotential(DHAT, stiffness)

    def test_ctor_rejects_none_barrier(self):
        with self.assertRaises(ValueError):
            ipctk.BarrierPotential(None, DHAT, STIFFNESS)

    def test_barrier_ctor_rejects_invalid_dhat(self):
        for dhat in self.INVALID:
            with self.subTest(dhat=dhat):
                with self.assertRaises(ValueError):
                    ipctk.BarrierPotential(
                        ipctk.ClampedLogBarrier(), dhat, STIFFNESS)

    def test_dhat_setter_rejects_invalid(self):
        for dhat in self.INVALID:
            with self.subTest(dhat=dhat):
                B = ipctk.BarrierPotential(DHAT, STIFFNESS)
                with self.assertRaises(ValueError):
                    B.dhat = dhat

    def test_stiffness_setter_rejects_invalid(self):
        for stiffness in self.INVALID:
            with self.subTest(stiffness=stiffness):
                B = ipctk.BarrierPotential(DHAT, STIFFNESS)
                with self.assertRaises(ValueError):
                    B.stiffness = stiffness

    def test_barrier_setter_rejects_none(self):
        B = ipctk.BarrierPotential(DHAT, STIFFNESS)
        with self.assertRaises(ValueError):
            B.barrier = None

    def test_state_unchanged_after_rejected_set(self):
        """A rejected assignment must not partially apply."""
        B = ipctk.BarrierPotential(DHAT, STIFFNESS)
        for attr, bad in (("dhat", 0.0), ("stiffness", -1.0)):
            with self.subTest(attr=attr):
                with self.assertRaises(ValueError):
                    setattr(B, attr, bad)
        self.assertEqual(B.dhat, DHAT)
        self.assertEqual(B.stiffness, STIFFNESS)

    def test_valid_values_accepted(self):
        B = ipctk.BarrierPotential(DHAT, STIFFNESS)
        B.dhat = 2e-3
        self.assertEqual(B.dhat, 2e-3)
        B.stiffness = 5.0
        self.assertEqual(B.stiffness, 5.0)
        B.barrier = ipctk.ClampedLogBarrier()
        # use_physical_barrier has no precondition
        B.use_physical_barrier = True
        self.assertTrue(B.use_physical_barrier)
        B.use_physical_barrier = False
        self.assertFalse(B.use_physical_barrier)
        # Tiny but positive must be allowed; only <= 0 and NaN are rejected.
        ipctk.BarrierPotential(1e-300, 1e-300)
        ipctk.BarrierPotential(ipctk.ClampedLogBarrier(), DHAT, STIFFNESS)


class TestBarrierPotentialProperties(unittest.TestCase):
    """The stiffness/use_physical_barrier properties must reach the evaluation
    path, not merely round-trip through a stored field."""

    @classmethod
    def setUpClass(cls):
        cls.mesh, cls.vertices = two_cubes()
        cls.collisions = normal_collisions(cls.mesh, cls.vertices)
        assert len(cls.collisions) > 0, "fixture produced no collisions"

    def test_ctor_roundtrip(self):
        B = ipctk.BarrierPotential(DHAT, 2.5, use_physical_barrier=True)
        self.assertEqual(B.dhat, DHAT)
        self.assertEqual(B.stiffness, 2.5)
        self.assertTrue(B.use_physical_barrier)

    def test_stiffness_scales_potential_linearly(self):
        """kappa multiplies the barrier potential, so tripling it must triple
        the value. Guards against a setter that stores but is never read."""
        B1 = ipctk.BarrierPotential(DHAT, 1.0)
        B3 = ipctk.BarrierPotential(DHAT, 1.0)
        B3.stiffness = 3.0

        p1 = B1(self.collisions, self.mesh, self.vertices)
        p3 = B3(self.collisions, self.mesh, self.vertices)
        self.assertGreater(p1, 0.0)
        # rtol, not exact: the sum is a parallel reduction, so the operand
        # order (and thus rounding) is not guaranteed to be reproducible.
        np.testing.assert_allclose(p3, 3.0 * p1, rtol=1e-9)

        g1 = B1.gradient(self.collisions, self.mesh, self.vertices)
        g3 = B3.gradient(self.collisions, self.mesh, self.vertices)
        np.testing.assert_allclose(g3, 3.0 * g1, rtol=1e-9)

    def test_use_physical_barrier_setter_matches_ctor(self):
        """Setting the property must be equivalent to passing the ctor kwarg."""
        for flag in (False, True):
            with self.subTest(use_physical_barrier=flag):
                via_ctor = ipctk.BarrierPotential(
                    DHAT, 2.5, use_physical_barrier=flag)
                via_setter = ipctk.BarrierPotential(DHAT, 1.0)
                via_setter.stiffness = 2.5
                via_setter.use_physical_barrier = flag

                args = (self.collisions, self.mesh, self.vertices)
                np.testing.assert_allclose(
                    via_setter(*args), via_ctor(*args), rtol=1e-9)
                np.testing.assert_allclose(
                    via_setter.gradient(*args), via_ctor.gradient(*args),
                    rtol=1e-9)
                np.testing.assert_allclose(
                    via_setter.hessian(*args).todense(),
                    via_ctor.hessian(*args).todense(), rtol=1e-9)

    def test_use_physical_barrier_changes_result(self):
        """The flag must actually alter the potential, otherwise the two
        branches of the test above would agree trivially."""
        off = ipctk.BarrierPotential(DHAT, STIFFNESS, False)
        on = ipctk.BarrierPotential(DHAT, STIFFNESS, True)
        args = (self.collisions, self.mesh, self.vertices)
        self.assertNotEqual(off(*args), on(*args))


class TestSmoothContactParameters(unittest.TestCase):
    def test_adaptive_dhat_ratio_default(self):
        params = ipctk.SmoothContactParameters(DHAT, 0.5, 0.0, 0.1, 0.0, 2)
        self.assertEqual(params.adaptive_dhat_ratio, 0.5)

    def test_adaptive_dhat_ratio_roundtrip(self):
        params = ipctk.SmoothContactParameters(DHAT, 0.5, 0.0, 0.1, 0.0, 2)
        for ratio in (0.1, 0.25, 0.9):
            params.adaptive_dhat_ratio = ratio
            self.assertEqual(params.adaptive_dhat_ratio, ratio)

    def test_adaptive_dhat_ratio_affects_adaptive_dhat(self):
        """A larger ratio yields larger per-element dhat, so a deformed
        configuration activates more collisions. Guards against the property
        being stored but never reaching compute_adaptive_dhat()."""
        mesh, rest = two_cubes()

        # Move the right cube toward the left one so the deformed state is
        # close enough to activate, but only for large enough adaptive dhat.
        deformed = rest.copy()
        deformed[rest[:, 0] > 0.96, 0] -= 0.04

        counts = []
        for ratio in (0.1, 0.5, 0.9):
            params = ipctk.SmoothContactParameters(DHAT, 0.5, 0.0, 0.1, 0.0, 2)
            params.adaptive_dhat_ratio = ratio
            collisions = ipctk.SmoothCollisions()
            collisions.compute_adaptive_dhat(mesh, rest, params)
            collisions.build(mesh, deformed, params, True)
            counts.append(len(collisions))

        self.assertEqual(counts[0], 0, f"expected no activation, got {counts}")
        self.assertLess(counts[0], counts[1], f"not monotonic: {counts}")
        self.assertLess(counts[1], counts[2], f"not monotonic: {counts}")


class TestSmoothCollisionsAdaptiveDhat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh, cls.rest = two_cubes()
        cls.params = ipctk.SmoothContactParameters(
            DHAT, 0.5, 0.0, 0.1, 0.0, 2)
        cls.potential = ipctk.SmoothContactPotential(cls.params)

    def _build(self, use_adaptive_dhat):
        collisions = ipctk.SmoothCollisions()
        if use_adaptive_dhat:
            collisions.compute_adaptive_dhat(self.mesh, self.rest, self.params)
        collisions.build(self.mesh, self.rest, self.params, use_adaptive_dhat)
        return collisions

    def test_non_adaptive_has_spurious_rest_forces(self):
        """Baseline: without adaptive dhat, this mesh/dhat pair produces a
        nonzero potential at rest. This is what adaptive dhat exists to fix, so
        if it ever becomes zero the test below stops proving anything."""
        collisions = self._build(use_adaptive_dhat=False)
        potential = self.potential(collisions, self.mesh, self.rest)
        gradient = self.potential.gradient(collisions, self.mesh, self.rest)
        self.assertGreater(potential, 0.0)
        self.assertGreater(np.abs(gradient).max(), 0.0)

    def test_adaptive_dhat_eliminates_spurious_rest_forces(self):
        """GCP's 'no spurious forces' guarantee: with adaptive dhat computed
        from the rest configuration, the potential and its gradient are exactly
        zero in that configuration."""
        collisions = self._build(use_adaptive_dhat=True)
        potential = self.potential(collisions, self.mesh, self.rest)
        gradient = self.potential.gradient(collisions, self.mesh, self.rest)
        self.assertEqual(potential, 0.0)
        np.testing.assert_array_equal(gradient, np.zeros_like(gradient))

    def test_broad_phase_argument_accepted(self):
        collisions = ipctk.SmoothCollisions()
        collisions.compute_adaptive_dhat(
            self.mesh, self.rest, self.params, ipctk.LBVH())


class TestSmoothContactPotentialNaming(unittest.TestCase):
    """The Python class was previously exposed as "SmoothPotential", which did
    not match the C++ name. Guard the rename in both directions."""

    def test_matches_cpp_name(self):
        self.assertTrue(hasattr(ipctk, "SmoothContactPotential"))

    def test_old_name_removed(self):
        self.assertFalse(hasattr(ipctk, "SmoothPotential"))


if __name__ == "__main__":
    unittest.main()
