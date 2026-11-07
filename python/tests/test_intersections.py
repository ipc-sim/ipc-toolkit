import numpy as np

from find_ipctk import ipctk


def test_segment_segment_intersect():
    assert ipctk.segment_segment_intersect(
        np.array([-1, 0]), np.array([1, 0]), np.array([0, -1]), np.array([0, 1]))

def test_edge_triangle_intersection():
    t0 = np.array([-1.0, -1.0, 0.0])
    t1 = np.array([1.0, -1.0, 0.0])
    t2 = np.array([0.0, 1.0, 0.0])

    # Vertical edge through the interior.
    e0 = np.array([0.1, -0.3, -2.0])
    e1 = np.array([0.1, -0.3, 1.0])

    intersects, u, v, t = ipctk.edge_triangle_intersection(e0, e1, t0, t1, t2)
    assert intersects
    assert u >= 0 and v >= 0 and u + v <= 1 + 1e-12
    assert 0 <= t <= 1

    # The two parameterizations of the hit point must agree.
    p_tri = t0 + u * (t1 - t0) + v * (t2 - t0)
    p_edge = e0 + t * (e1 - e0)
    assert np.linalg.norm(p_tri - p_edge) < 1e-12

    # Consistent with the boolean-only predicate.
    assert intersects == ipctk.is_edge_intersecting_triangle(
        e0, e1, t0, t1, t2)


def test_edge_triangle_intersection_miss():
    t0 = np.array([-1.0, -1.0, 0.0])
    t1 = np.array([1.0, -1.0, 0.0])
    t2 = np.array([0.0, 1.0, 0.0])

    # Edge entirely on one side of the triangle's plane.
    intersects, u, v, t = ipctk.edge_triangle_intersection(
        np.array([0.0, 0.0, 0.5]), np.array([0.0, 0.0, 1.5]), t0, t1, t2)
    assert not intersects
    assert np.isnan([u, v, t]).all()

    # Edge crosses the plane, but outside the triangle.
    intersects, u, v, t = ipctk.edge_triangle_intersection(
        np.array([5.0, 5.0, -1.0]), np.array([5.0, 5.0, 1.0]), t0, t1, t2)
    assert not intersects
    assert np.isnan([u, v, t]).all()
