import numpy as np

import find_ipctk
import ipctk

import utils


def check_ipc_derivatives(broad_phase_method, use_convergent_formulation, mesh_name, dhat, all_vertices_on_surface):
    vertices, edges, faces = utils.load_mesh(mesh_name)

    if all_vertices_on_surface:
        mesh = ipctk.CollisionMesh(vertices, edges, faces)
    else:
        mesh = ipctk.CollisionMesh.build_from_full_mesh(vertices, edges, faces)
        vertices = mesh.vertices(vertices)

    collisions = ipctk.NormalCollisions()
    collisions.use_area_weighting = use_convergent_formulation
    collisions.use_improved_max_approximator = use_convergent_formulation
    collisions.build(mesh, vertices, dhat,
                     broad_phase_method=broad_phase_method)
    assert len(collisions) > 0

    B = ipctk.BarrierPotential(
        dhat, use_physical_barrier=use_convergent_formulation)

    grad_b = B.gradient(collisions, mesh, vertices)
    fgrad_b = utils.finite_gradient(
        vertices.flatten(), lambda x: B(collisions, mesh, x.reshape(vertices.shape)))

    assert np.linalg.norm(grad_b) > 1e-8
    assert np.allclose(grad_b, fgrad_b)

    hess_b = B.hessian(collisions, mesh, vertices).toarray()
    fhess_b = utils.finite_jacobian(
        vertices.flatten(), lambda x: B.gradient(collisions, mesh, x.reshape(vertices.shape)))

    assert np.linalg.norm(hess_b) > 1e-8
    assert np.allclose(hess_b, fhess_b, atol=1e-5)


def test_ipc():
    for method in utils.broad_phase_methods():
        for use_convergent_formulation in (True, False):
            yield check_ipc_derivatives, method, use_convergent_formulation, "cube.ply", np.sqrt(2.0), True
            yield check_ipc_derivatives, method, use_convergent_formulation, "two-cubes-far.ply", 1e-1, False
            yield check_ipc_derivatives, method, use_convergent_formulation, "two-cubes-close.ply", 1e-1, False
            yield check_ipc_derivatives, method, use_convergent_formulation, "bunny.ply", 5e-3, True
