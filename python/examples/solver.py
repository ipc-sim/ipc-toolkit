import numpy as np
import scipy
import meshio
import ipctk

mesh = meshio.read("../tests/data/two-cubes-close.ply")
vertices = mesh.points
faces = mesh.cells_dict["triangle"]


class Inertia:
    """Simple inertia with prescribed target positions."""

    def __init__(self):
        import igl
        mid_x = (vertices[:, 0].min() + vertices[:, 0].max()) / 2
        self.xhat = vertices.copy()
        self.xhat[vertices[:, 0] < mid_x, 0] += 0.2
        self.xhat[vertices[:, 0] > mid_x, 0] -= 0.2
        self.mass = igl.massmatrix(
            vertices, faces, igl.MASSMATRIX_TYPE_VORONOI)
        self.mass = np.repeat(self.mass.diagonal(), 3)
        self.mass = scipy.sparse.diags(self.mass)

    def __call__(self, x):
        dx = (x - self.xhat).flatten()
        return dx.transpose() @ self.mass @ dx

    def gradient(self, x):
        return self.mass @ (x - self.xhat).flatten()

    def hessian(self, x):
        return self.mass


inertia = Inertia()

# -----------------------------------------------------------------------------

# Contact parameters
dhat = 0.001
kappa = 100 * inertia.mass.toarray().mean()

# Build the collision mesh
collision_mesh = ipctk.CollisionMesh(vertices, ipctk.edges(faces), faces)

# Build the collision set
C = ipctk.NormalCollisions()
C.use_area_weighting = True
C.use_improved_max_approximator = True
C.build(collision_mesh, vertices, dhat)

# Create a barrier potential
B = ipctk.BarrierPotential(dhat, use_physical_barrier=True)

# Initial energy and gradient
prev_energy = (
    inertia(vertices)
    + kappa * B(C, collision_mesh, vertices)
    # + ...
)
grad = (
    inertia.gradient(vertices)
    + kappa * B.gradient(C, collision_mesh, vertices)
    # + ...
)

iter = 0
while np.linalg.norm(grad) > 1e-3 and iter < 1000:
    meshio.Mesh(vertices, {"triangle": faces}).write(
        f"solver-step-{iter:03d}.ply")

    # Compute the Hessian
    hess = (
        inertia.hessian(vertices)
        + kappa * B.hessian(
            C, collision_mesh, vertices,
            project_hessian_to_psd=ipctk.PSDProjectionMethod.CLAMP)
        # + ...
    )

    # Solve for the step direction dx = -hess⁻¹ grad
    dx = scipy.sparse.linalg.spsolve(hess, -grad)
    assert grad.dot(dx) < 0

    next_vertices = vertices + dx.reshape(-1, 3, order="C")

    # ------------------------------------------------------------------------------
    # Perform a line search to find a step size that decreases the barrier potential

    # Candidates for continuous collision detection
    candidates = ipctk.Candidates()
    candidates.build(collision_mesh, vertices, next_vertices, dhat)

    # Initial step size for line search (ensure no collisions)
    alpha = candidates.compute_collision_free_stepsize(
        collision_mesh, vertices, next_vertices)

    # Backtrack until the barrier potential decreases
    condition = True
    ls_iter = 0
    while condition:
        next_vertices = vertices + alpha * dx.reshape(-1, 3, order="C")
        C.build(candidates, collision_mesh, next_vertices, dhat)
        curr_energy = (
            inertia(next_vertices)
            + kappa * B(C, collision_mesh, next_vertices)
            # + ...
        )
        condition = curr_energy > prev_energy and ls_iter < 100
        if condition:
            alpha *= 0.5
            ls_iter += 1

    vertices = next_vertices

    prev_energy = curr_energy
    grad = (
        inertia.gradient(vertices)
        + kappa * B.gradient(C, collision_mesh, vertices)
        # + ...
    )

    print("||grad||:", np.linalg.norm(grad))

    iter += 1

meshio.Mesh(vertices, {"triangle": faces}).write(f"solver-step-{iter:03d}.ply")
