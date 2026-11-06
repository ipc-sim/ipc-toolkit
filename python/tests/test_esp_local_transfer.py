"""Python access to the experimental local-transfer formulation."""
import unittest
import numpy as np
from find_ipctk import ipctk


class TestESPLocalTransfer(unittest.TestCase):
    def test_opt_in_and_derivative(self):
        V = np.array([[-.01, -.0003, .0006], [.01, .0002, .0008],
                      [-.005, -.021, 0], [-.0047, .019, 0]])
        mesh = ipctk.CollisionMesh(V, np.array([[0, 1], [2, 3]], dtype=int))
        mesh.can_collide = ipctk.make_vertex_patches_filter(np.array([0, 0, 1, 1]))
        mesh.init_adjacencies()
        params = ipctk.ESPParameters(.03, .05)
        self.assertFalse(params.local_quadrature_transfer)
        params.local_quadrature_transfer = True
        potential = ipctk.ESPPotential(params, use_near_far=True)
        contacts = ipctk.ESPCollisions()
        contacts.build(mesh, V, params)
        gradient = potential.gradient(contacts, mesh, V)
        self.assertTrue(np.isfinite(gradient).all())
        h = 1e-9
        values = []
        for sign in [-1, 1]:
            X = V.copy()
            X[0, 0] += sign*h
            contacts.build(mesh, X, params)
            values.append(potential(contacts, mesh, X))
        self.assertAlmostEqual((values[1]-values[0])/(2*h), gradient[0], delta=1e-5)
        with self.assertRaises(ValueError):
            ipctk.ESPPotential(params, use_near_far=False)

    def test_edge_edge_options(self):
        params = ipctk.ESPParameters(.03, .05)
        self.assertEqual(params.ee_mollifier.type,
                         ipctk.EEMollifierType.SMOOTHSTEP)
        self.assertEqual(params.ee_mollifier.rel_eps, 1.0)
        self.assertTrue(params.ee_area_scaling)
        self.assertEqual(params.ee_area_ref, 0)
        self.assertEqual(params.ee_weight, 1.5)
        self.assertEqual(params.ee_mollifier_power, 0)
        self.assertEqual(params.ee_mollifier_exponent, 1)  # log barrier
        params.barrier = ipctk.InversePowerBarrier(2)
        self.assertEqual(params.ee_mollifier_exponent, 1)  # smoothstep
        mollifier = ipctk.EEMollifierConfig()
        mollifier.type = ipctk.EEMollifierType.GCP
        params.ee_mollifier = mollifier
        self.assertEqual(params.ee_mollifier.type, ipctk.EEMollifierType.GCP)
        self.assertEqual(params.ee_mollifier_exponent, 3)  # GCP, p = 2
        params.ee_mollifier_power = 1
        self.assertEqual(params.ee_mollifier_exponent, 1)

    def test_face_rules(self):
        # Two triangles facing each other, 4 mm apart.
        V = np.array([[0, 0, 0], [.1, 0, 0], [0, .1, 0],
                      [.02, .01, .004], [.09, .02, .006], [.01, .08, .005]])
        F = np.array([[0, 1, 2], [3, 5, 4]])
        E = np.array([[0, 1], [1, 2], [2, 0], [3, 5], [5, 4], [4, 3]])
        mesh = ipctk.CollisionMesh(V, E, F)
        mesh.init_adjacencies()

        def energy(params):
            contacts = ipctk.ESPCollisions()
            contacts.build(mesh, V, params)
            return ipctk.ESPPotential(params, use_near_far=True)(contacts, mesh, V)

        default = ipctk.ESPParameters(.01, 1.0)
        explicit = ipctk.ESPParameters(.01, 1.0)
        # The vertex rule with unnormalized weights: same energy as the default.
        explicit.set_quad_rule(np.eye(3), np.array([2., 2., 2.]))
        self.assertGreater(energy(default), 0)
        self.assertAlmostEqual(energy(explicit) / energy(default), 1, places=12)
        lattice = ipctk.ESPParameters(.01, 1.0)
        lattice.set_lattice_rule(1)
        self.assertAlmostEqual(energy(lattice) / energy(default), 1, places=12)
        with self.assertRaises(ValueError):
            explicit.set_quad_rule(np.zeros((0, 3)), np.zeros(0))
        with self.assertRaises(ValueError):
            explicit.set_quad_rule(np.eye(3), np.array([1., -1., 1.]))
