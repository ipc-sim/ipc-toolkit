"""Python access to the experimental, frictionless local-transfer formulation."""
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
