import numpy as np

from find_ipctk import ipctk
from ipctk.filib import Interval

from utils import load_mesh


def test_candidates():
    V0, E, F = load_mesh("two-cubes-close.ply")
    V1, E, F = load_mesh("two-cubes-intersecting.ply")

    mesh = ipctk.CollisionMesh(V0, E, F)

    candidates = ipctk.Candidates()
    candidates.build(mesh, V0, V1)

    assert len(candidates) > 0, "No candidates generated."