import numpy as np
from utils import *
import ipctk


def check_ipc_derivatives(broad_phase_method, use_convergent_formulation, mesh_name, dhat, all_vertices_on_surface):
    V, E, F = load_mesh(mesh_name)

    if all_vertices_on_surface:
        mesh = ipctk.CollisionMesh(V, E, F)
    else:
        mesh = ipctk.CollisionMesh.build_from_full_mesh(V, E, F)
        V = mesh.vertices(V)

    collision_constraints = ipctk.CollisionConstraints()
    collision_constraints.use_convergent_formulation = use_convergent_formulation
    collision_constraints.build(
        mesh, V, dhat, broad_phase_method=broad_phase_method)
    assert len(collision_constraints) > 0

    grad_b = collision_constraints.compute_potential_gradient(mesh, V, dhat)
    fgrad_b = finite_gradient(
        V.flatten(), lambda x: collision_constraints.compute_potential(mesh, x.reshape(V.shape), dhat))

    assert np.linalg.norm(grad_b) > 1e-8
    assert np.allclose(grad_b, fgrad_b)

    hess_b = collision_constraints.compute_potential_hessian(mesh, V, dhat).A
    fhess_b = finite_jacobian(
        V.flatten(), lambda x: collision_constraints.compute_potential_gradient(mesh, x.reshape(V.shape), dhat))

    assert np.linalg.norm(hess_b) > 1e-8
    assert np.allclose(hess_b, fhess_b, atol=1e-5)


def test_ipc():
    for method in broad_phase_methods():
        for use_convergent_formulation in (True, False):
            yield check_ipc_derivatives, method, use_convergent_formulation, "cube.obj", np.sqrt(2.0), True
            yield check_ipc_derivatives, method, use_convergent_formulation, "two-cubes-far.obj", 1e-1, False
            yield check_ipc_derivatives, method, use_convergent_formulation, "two-cubes-close.obj", 1e-1, False
            yield check_ipc_derivatives, method, use_convergent_formulation, "bunny.obj", 5e-3, True
