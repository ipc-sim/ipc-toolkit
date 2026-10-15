import time

import find_ipctk
import numpy as np
import scipy
from ipctk import (
    CollisionMesh,
    NormalCollisions,
    make_sparse_filter,
    make_vertex_patches_filter,
)


def test_collision_mesh():
    V = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    E = np.array([[0, 1], [1, 3], [3, 2], [2, 0]], dtype=int)

    W = scipy.sparse.lil_matrix((4, 3), dtype=float)
    W[0, 0] = W[1, 1] = W[2, 2] = W[3, 1] = W[3, 2] = 1

    mesh = CollisionMesh(V, E, displacement_map=W)

    U = np.array([[0, 0], [1, 1], [0, 0]], dtype=float)

    Uc = mesh.map_displacements(U)
    expected_Uc = np.array([[0, 0], [1, 1], [0, 0], [1, 1]], dtype=float)
    assert (Uc == expected_Uc).all()

    Vc = mesh.displace_vertices(U)
    expected_Vc = np.array([[0, 0], [2, 1], [0, 1], [2, 2]], dtype=float)
    assert (Vc == expected_Vc).all()

    g = np.array([1, 1, -1, 1, 1, -1, -1, -1], dtype=float)
    gf = mesh.to_full_dof(g)
    expected_gf = np.array([1, 1, -2, 0, 0, -2], dtype=float)
    assert (gf == expected_gf).all()

    H = np.eye(8)
    Hf = mesh.to_full_dof(scipy.sparse.csc_matrix(H)).toarray()
    expected_Hf = np.array([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 2, 0, 1, 0],
        [0, 0, 0, 2, 0, 1],
        [0, 0, 1, 0, 2, 0],
        [0, 0, 0, 1, 0, 2],
    ], dtype=float)
    assert (Hf == expected_Hf).all()


def check_faces_to_edges(E, expected_F2E):
    F = np.array([[0, 1, 2]])
    assert (CollisionMesh.construct_faces_to_edges(F, E) == expected_F2E).all()


def test_faces_to_edges():
    yield check_faces_to_edges, np.array([[0, 1], [1, 2], [2, 0]]), np.array([0, 1, 2])
    yield check_faces_to_edges, np.array([[2, 0], [2, 1], [1, 0]]), np.array([2, 1, 0])
    yield check_faces_to_edges, np.array([[0, 1], [2, 0], [2, 1]]), np.array([0, 2, 1])
    # Shouldnt work
    try:
        check_faces_to_edges(np.array([[0, 1], [1, 2], [0, 3]]), None)
        assert False
    except RuntimeError as e:
        assert str(e) == "Unable to find edge!"
    except:
        assert False


def test_codim_points_collision_mesh():
    V = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    mesh = CollisionMesh(V)
    expected_codim_vertices = np.array([0, 1, 2, 3], dtype=int)
    assert (mesh.codim_vertices == expected_codim_vertices).all()


def test_collision_mesh_no_faces():
    # Based on https://github.com/ipc-sim/ipc-toolkit/issues/103
    points = np.array([[-0.5, 0.5], [0.5, 0.5], [0.5, -0.5]])
    edges = np.array([[0, 1], [1, 2], [2, 0]])
    mesh = CollisionMesh(points, edges)
    assert mesh.faces.size == 0


def test_collision_mesh_does_not_segfault():
    # Based on https://github.com/ipc-sim/ipc-toolkit/issues/102
    points = np.array([[-0.5, 0.5], [0.5, 0.5], [0.5, -0.5]])
    edges = np.array([[0, 1], [1, 2], [2, 0]])
    faces = np.array([[0, 1, 2]])
    mesh = CollisionMesh(points, edges, faces)
    assert (mesh.rest_positions == points).all()
    assert (mesh.edges == edges).all()
    assert (mesh.faces == faces).all()


def test_can_collide():
    V = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    E = np.array([[0, 1], [1, 3], [3, 2], [2, 0]], dtype=int)

    V = np.vstack([V, V + 1.001])
    E = np.vstack([E, E + V.shape[0] // 2])

    mesh = CollisionMesh(V, E)

    def default_can_collide(i, j):
        return True

    patches = np.concatenate([np.zeros(4, dtype=int), np.ones(4, dtype=int)])
    print(patches.size)

    def patches_can_collide(i, j):
        return patches[i] != patches[j]

    dict_can_collide = {}
    for i in range(V.shape[0]):
        for j in range(V.shape[0]):
            if i < j and not patches_can_collide(i, j):
                dict_can_collide[(i, j)] = False

    can_collides = [
        default_can_collide,
        patches_can_collide,
        make_sparse_filter(dict_can_collide, True),
        make_vertex_patches_filter(patches),
    ]

    for can_collide in can_collides:
        if can_collide != default_can_collide:
            mesh.can_collide = can_collide

        for i in range(V.shape[0]):
            for j in range(V.shape[0]):
                assert mesh.can_collide(i, j) == can_collide(i, j)


def test_face_normals():
    # Two triangles of a unit square in the z = 0 plane.
    V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    F = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    E = np.array([[0, 1], [1, 2], [2, 0], [2, 3], [3, 0]], dtype=int)

    mesh = CollisionMesh(V, E, F)

    N = mesh.face_normals(V)
    assert N.shape == (mesh.num_faces, 3)
    assert np.allclose(N, np.array([[0, 0, 1], [0, 0, 1]], dtype=float))
    assert np.allclose(np.linalg.norm(N, axis=1), 1.0)

    # Rotate 90° about the x-axis: +z normal ⇒ −y normal.
    V_rot = np.column_stack((V[:, 0], -V[:, 2], V[:, 1]))
    N_rot = mesh.face_normals(V_rot)
    assert np.allclose(N_rot, np.array([[0, -1, 0], [0, -1, 0]], dtype=float))


def test_face_normals_2d_raises():
    V = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    E = np.array([[0, 1], [1, 3], [3, 2], [2, 0]], dtype=int)

    mesh = CollisionMesh(V, E)

    try:
        mesh.face_normals(V)
        assert False
    except ValueError:
        pass


def _two_segments():
    V = np.array([[0, 0], [1, 0], [0, 0.1], [1, 0.1]], dtype=float)
    E = np.array([[0, 1], [2, 3]], dtype=int)
    return V, E


def test_custom_areas_become_collision_weights():
    V, E = _two_segments()
    scale = 2 * np.pi * 3.0

    default_mesh = CollisionMesh(V, E)
    default_collisions = NormalCollisions()
    default_collisions.use_area_weighting = True
    default_collisions.build(default_mesh, V, dhat=0.2)

    mesh = CollisionMesh(V, E)
    mesh.vertex_areas = scale * mesh.vertex_areas
    mesh.edge_areas = scale * mesh.edge_areas
    assert np.allclose(mesh.vertex_areas, scale * default_mesh.vertex_areas)
    assert np.allclose(mesh.edge_areas, scale * default_mesh.edge_areas)

    collisions = NormalCollisions()
    collisions.use_area_weighting = True
    collisions.build(mesh, V, dhat=0.2)

    assert len(collisions) > 0
    assert len(collisions) == len(default_collisions)
    for i in range(len(collisions)):
        assert np.isclose(
            collisions[i].weight, scale * default_collisions[i].weight)


def test_invalid_custom_areas_are_rejected():
    V, E = _two_segments()
    mesh = CollisionMesh(V, E)
    original = mesh.vertex_areas.copy()

    for bad in [np.ones(3), np.full(4, -1.0), np.full(4, np.nan), np.full(4, np.inf)]:
        try:
            mesh.vertex_areas = bad
            assert False
        except (ValueError, RuntimeError):
            pass
    try:
        mesh.edge_areas = np.ones(3)
        assert False
    except (ValueError, RuntimeError):
        pass

    assert np.array_equal(mesh.vertex_areas, original)


def test_custom_areas_and_area_jacobians_are_exclusive():
    V, E = _two_segments()

    # Setting after initializing the Jacobians
    mesh = CollisionMesh(V, E)
    mesh.init_area_jacobians()
    try:
        mesh.vertex_areas = 2 * mesh.vertex_areas
        assert False
    except (ValueError, RuntimeError):
        pass

    # Initializing after setting
    mesh = CollisionMesh(V, E)
    mesh.vertex_areas = 2 * mesh.vertex_areas
    try:
        mesh.init_area_jacobians()
        assert False
    except (ValueError, RuntimeError):
        pass


def test_area_getters_are_read_only_views():
    # In-place edits must fail loudly instead of bypassing the validating setter.
    V, E = _two_segments()
    mesh = CollisionMesh(V, E)
    try:
        mesh.vertex_areas[0] = 5.0
        assert False
    except ValueError:
        pass


def test_zero_custom_areas_drop_their_collisions():
    # E.g., axisymmetric primitives on the symmetry axis have no ring area.
    V, E = _two_segments()
    mesh = CollisionMesh(V, E)
    mesh.vertex_areas = np.zeros(4)
    mesh.edge_areas = np.zeros(2)

    collisions = NormalCollisions()
    collisions.use_area_weighting = True
    collisions.build(mesh, V, dhat=0.2)

    # Zero-weight collisions are dropped, so no contact force is produced.
    assert len(collisions) == 0
