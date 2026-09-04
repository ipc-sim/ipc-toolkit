"""Double pendulum via Affine Body Dynamics joint constraints [Chen et al. 2022].

Two rods released from horizontal: the upper rod is pinned to a fixed point at
its left end, and the lower rod's left end is glued to the upper rod's free end
with a point connection. Under gravity the linkage swings chaotically.

The dynamics are solved at a small timestep (dt = 1/600) for accuracy, but the
viewer only redraws every SUBSTEPS = 10 steps, so it still runs at ~60 fps.

Joint constraints are linear equalities in the affine DOFs, so they require the
affine dynamics model (BodyDynamics.AFFINE). The two rods overlap at the shared
joint, so collision between them is disabled with a vertex-patch filter (there
is no collision filtering between jointed bodies otherwise).
"""

import numpy as np
import polyscope as ps
from find_ipctk import ipctk
from polyscope import imgui

ipctk.set_logger_level(ipctk.warn)


def box(half_extents):
    """Axis-aligned box centered at the origin with outward-facing normals."""
    V = (
        np.array(
            [
                [-1, -1, -1],
                [1, -1, -1],
                [1, 1, -1],
                [-1, 1, -1],
                [-1, -1, 1],
                [1, -1, 1],
                [1, 1, 1],
                [-1, 1, 1],
            ],
            dtype=float,
        )
        * half_extents
    )
    F = np.array(
        [
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [2, 3, 7],
            [2, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ],
        dtype=np.int32,
    )
    return V, ipctk.edges(F), F


# --- Scene: two rods forming a double pendulum --------------------------------

half_length = 0.6  # each rod spans x in [-half_length, half_length]
y0 = 2.5  # height of the pivot

anchor = np.array([0.0, y0, 0.0])  # fixed pivot (upper rod's left end)
joint = np.array([2 * half_length, y0, 0.0])  # upper's right end = lower's left

V_rod, E_rod, F_rod = box([half_length, 0.06, 0.06])
initial_poses = ipctk.Poses(
    [
        # Upper rod: horizontal, left end at the pivot.
        ipctk.Pose(np.array([half_length, y0, 0.0]), np.zeros(3)),
        # Lower rod: horizontal, left end at the joint.
        ipctk.Pose(np.array([3 * half_length, y0, 0.0]), np.zeros(3)),
    ]
)

bodies = ipctk.RigidBodies(
    rest_positions=[V_rod, V_rod],
    edges=[E_rod, E_rod],
    faces=[F_rod, F_rod],
    densities=np.full(2, 1000.0),
    initial_poses=initial_poses,
)

# The rods meet at the joint, so disable collision between them: putting all
# vertices in one patch blocks every same-patch (i.e. every) collision pair.
bodies.can_collide = ipctk.make_vertex_patches_filter(
    np.zeros(2 * V_rod.shape[0], dtype=np.int32)
)

joints = ipctk.JointConstraints(bodies, initial_poses)
joints.add_fixed_point(0, anchor)  # pin the upper rod
joints.add_point_connection(0, 1, joint)  # link lower rod to upper rod

settings = ipctk.demo.Simulator.Settings()
settings.body_dynamics = ipctk.demo.Simulator.BodyDynamics.AFFINE  # joints require affine DOFs

# Use second-order BDF for more accurate dynamics (order 1 is implicit Euler).
# The effective order ramps up from 1 over the first two steps.
settings.bdf_order = 2

# Tighten the Newton convergence tolerance for accurate (near energy-conserving)
# dynamics. The default (1e-2) terminates each Newton solve while the velocity
# is still under-resolved; that small per-step error is dissipative, and since
# the tolerance is a fixed velocity independent of dt, it accumulates over more
# steps at a smaller dt -- so the default actually damps *more* at dt=1/600 than
# at dt=1/60. At 1e-4 the swing amplitude is conserved and, as expected, the
# smaller timestep is the less-damped one.

# Solve at a small timestep for accuracy, but only redraw every SUBSTEPS steps
# so the viewer still runs at ~60 fps (600 steps/s ÷ 10 = 60 frames/s).
DT = 1 / 60.0
SUBSTEPS = 10
DT /= SUBSTEPS

sim = ipctk.demo.Simulator(
    bodies=bodies,
    initial_poses=initial_poses,
    joints=joints,
    dt=DT,
    settings=settings,
)


# build_from_meshes folds each body's principal-axis rotation into its rest
# frame, so track material points by their initial world position instead of
# assuming a body-frame axis.
lower_tip_rest = joint + np.array([2 * half_length, 0.0, 0.0])


def track(body, world_point_at_rest):
    """World position at the latest step of the material point that was at
    ``world_point_at_rest`` in the initial configuration."""
    p0 = initial_poses[body]
    R0 = np.array(p0.rotation_matrix())
    material = R0.T @ (np.asarray(world_point_at_rest) - np.array(p0.position))
    # sim.rigid_poses is the current frame only (O(num_bodies)); rigid_pose_history
    # would copy the entire growing history on every call.
    pose = sim.rigid_poses[body]
    return pose.transform_vertices(material.reshape(1, 3))[0]


def joint_pos():
    """Current world position of the moving joint (upper rod's free end)."""
    return track(0, joint)


def lower_tip():
    """World position of the lower rod's free (swinging) end."""
    return track(1, lower_tip_rest)


# --- Viewer -------------------------------------------------------------------

ps.init()
ps.set_give_focus_on_show(True)
ps.look_at((0.0, y0, 6.0), (2 * half_length, y0 - 1.0, 0.0))

ps_mesh = ps.register_surface_mesh(
    "double pendulum", bodies.vertices(initial_poses), bodies.faces
)
ps_pivots = ps.register_point_cloud("pivots", np.vstack([anchor, joint]))
ps_pivots.set_color((0.9, 0.3, 0.1))
ps_pivots.set_radius(0.012)

# A fading trail of the swinging tip makes the chaotic motion visible. One
# point is recorded per solver step (SUBSTEPS per frame), so ~2 s of history.
TRAIL_LEN = 120 * SUBSTEPS
trail = np.tile(lower_tip(), (TRAIL_LEN, 1))
ps_trail = ps.register_point_cloud("tip trail", trail)
ps_trail.set_color((0.1, 0.5, 0.9))
ps_trail.set_radius(0.004)
ps_trail.add_scalar_quantity(
    "age", np.linspace(0.0, 1.0, TRAIL_LEN), enabled=True, cmap="blues"
)

playing = False


def push_trail(point):
    """Append a tip position to the rolling trail buffer (recorded per step so
    the fine sub-frame motion is captured)."""
    global trail
    trail = np.roll(trail, -1, axis=0)
    trail[-1] = point


def render():
    # sim.rigid_poses is the current frame converted to the Poses type expected
    # by bodies.vertices (O(num_bodies), unlike rigid_pose_history).
    ps_mesh.update_vertex_positions(bodies.vertices(sim.rigid_poses))
    ps_pivots.update_point_positions(np.vstack([anchor, joint_pos()]))
    ps_trail.update_point_positions(trail)


def advance():
    """Take SUBSTEPS simulation steps, recording the tip each step. Returns
    False if a step failed to converge."""
    for _ in range(SUBSTEPS):
        if not sim.step():
            return False
        push_trail(lower_tip())
    return True


def callback():
    global playing, trail
    if imgui.Button("Play" if not playing else "Pause") or imgui.IsKeyPressed(
        imgui.ImGuiKey_Space
    ):
        playing = not playing
    imgui.SameLine()
    if imgui.Button("Step") or playing:
        if not advance():
            playing = False
        render()
    imgui.SameLine()
    if imgui.Button("Reset"):
        sim.reset()
        trail = np.tile(lower_tip(), (TRAIL_LEN, 1))
        render()


ps.set_user_callback(callback)
ps.show()
