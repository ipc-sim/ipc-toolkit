import numpy as np

from find_ipctk import ipctk


def two_close_triangles():
    """Two parallel triangles 0.05 apart, so an inflation radius of 0.1 makes
    every vertex of one overlap the face of the other."""
    V = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0],
        [0, 0, 0.05], [1, 0, 0.05], [0, 1, 0.05],
    ], dtype=float)
    E = np.array([[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3]])
    F = np.array([[0, 1, 2], [3, 4, 5]])
    return V, E, F


def test_cuda_is_available_predicts_gpu_broad_phase():
    available = ipctk.cuda.is_available()
    assert isinstance(available, bool)

    if not hasattr(ipctk.cuda, "LBVH"):
        # A CPU-only build has no GPU classes, so it can never have CUDA.
        assert not available
        return

    V, E, F = two_close_triangles()
    lbvh = ipctk.cuda.LBVH()
    if available:
        lbvh.build(V, E, F, inflation_radius=0.1)
        assert len(lbvh.detect_collision_candidates()) > 0
    else:
        # Without a device the first CUDA call fails, and we surface it as a
        # Python exception rather than crashing the interpreter.
        try:
            lbvh.build(V, E, F, inflation_radius=0.1)
        except RuntimeError:
            pass
        else:
            assert False, "LBVH.build should raise without a CUDA device"
