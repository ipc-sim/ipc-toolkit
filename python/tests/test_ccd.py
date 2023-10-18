import numpy as np
from utils import *
import ipctk


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
