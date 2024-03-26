import numpy as np

import find_ipctk
import ipctk


def test_segment_segment_intersect():
    assert ipctk.segment_segment_intersect(
        np.array([-1, 0]), np.array([1, 0]), np.array([0, -1]), np.array([0, 1]))
