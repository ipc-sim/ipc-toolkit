import numpy as np
import scipy

import find_ipctk import ipctk.CollisionMesh


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
