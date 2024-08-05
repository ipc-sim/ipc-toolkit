import pathlib
import numpy as np
import meshio

import find_ipctk
import ipctk


def download_test_data_if_needed(directory):
    if directory.exists():
        return

    # Clone the test data repository
    import subprocess
    subprocess.run([
        'git', 'clone', 'https://github.com/ipc-sim/ipc-toolkit-tests-data',
        str(directory)
    ])


def test_data_dir():
    _test_data_dir = pathlib.Path(__file__).parents[2] / 'tests' / 'data'
    download_test_data_if_needed(_test_data_dir)
    return _test_data_dir


def load_mesh(mesh_name):
    mesh = meshio.read(test_data_dir() / mesh_name)
    return mesh.points, ipctk.edges(mesh.cells_dict['triangle']), mesh.cells_dict['triangle']


def broad_phase_methods():
    yield ipctk.BroadPhaseMethod.BRUTE_FORCE
    yield ipctk.BroadPhaseMethod.HASH_GRID
    yield ipctk.BroadPhaseMethod.SPATIAL_HASH
    yield ipctk.BroadPhaseMethod.BOUNDING_VOLUME_HIERARCHY
    yield ipctk.BroadPhaseMethod.SWEEP_AND_PRUNE


def finite_jacobian(x, f, h=1e-8):
    external_coeffs = [1, -1]
    internal_coeffs = [1, -1]
    assert len(external_coeffs) == len(internal_coeffs)
    inner_steps = len(internal_coeffs)
    denom = 2 * h

    jac = np.zeros((f(x).shape[0], x.shape[0]))

    x_mutable = x.copy()
    for i in range(x.shape[0]):
        for ci in range(inner_steps):
            x_mutable[i] += internal_coeffs[ci] * h
            jac[:, i] += external_coeffs[ci] * f(x_mutable)
            x_mutable[i] = x[i]
        jac[:, i] /= denom

    return jac


def finite_gradient(x, f, h=1e-8):
    return finite_jacobian(x, lambda x: np.array([f(x)]), h)
