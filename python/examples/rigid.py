"""Interactive rigid body simulation demo.

Usage:
    python rigid.py                            # built-in bunnies-in-a-bowl scene
    python rigid.py scene.glb [settings.py]    # glTF geometry + Python settings

The glTF file only provides geometry: each mesh node becomes a body named after
its node, with a default density of 1000. Everything physical comes from the
optional settings file, a plain Python module that may define:

    dt              Time step in seconds (default 1/60).
    convert_planes  Turn two-triangle squares into infinite planes (default True).
    is_dof_fixed    {body name: [x, y, z, rx, ry, rz] booleans}. These are read
                    while loading because fixing rotational DOFs changes how a
                    body's inertia is stored, so they cannot change afterwards.
    configure(scene)
                    Called after loading and before building the simulator. It
                    edits the real ipctk objects through the Scene below, so
                    anything the bindings expose can be set (see scenes/).
"""

import pathlib
import runpy
import sys

import meshio
import numpy as np
import polyscope as ps
from find_ipctk import ipctk
from polyscope import imgui
from scipy.spatial.transform import Rotation

ipctk.set_logger_level(ipctk.debug)


class Scene:
    """The state a settings file's configure(scene) may change.

    Attributes:
        bodies: The ipctk.RigidBodies (set density, type, external_force, ...).
        poses: The initial poses (aligned with bodies).
        settings: The ipctk.demo.Simulator.Settings used to build the simulator.
            Dict-valued fields (solver_params, linear_solver_params) return a
            copy, so assign a whole dict back rather than editing in place.
        joints: Optional ipctk.JointConstraints (affine dynamics only).
    """

    def __init__(self, bodies, poses):
        self.bodies = bodies
        self.poses = poses
        self.settings = ipctk.demo.Simulator.Settings()
        self.settings.body_dynamics = ipctk.demo.Simulator.BodyDynamics.RIGID
        self.settings.bdf_order = 1
        self.joints = None
        self._velocities = None
        self._drivers = {}

    @property
    def names(self):
        """The body names, in body order."""
        return [self.bodies.body_name(i) for i in range(len(self.bodies))]

    def index(self, name):
        """The index of the body called ``name``."""
        i = self.bodies.body_index(name)
        if i is None:
            raise KeyError(f"no body named {name!r}; the bodies are {self.names}")
        return i

    def body(self, name):
        """The ipctk.RigidBody called ``name`` (edits apply to the scene)."""
        return self.bodies[self.index(name)]

    def velocity(self, name):
        """The initial velocity of body ``name`` as a Pose (position is the
        linear velocity, rotation the angular velocity, both in world space).
        Edits apply to the scene; bodies start at rest unless set."""
        if self._velocities is None:
            self._velocities = ipctk.Poses(
                [
                    ipctk.Pose(np.zeros(3), np.zeros(3))
                    for _ in range(len(self.bodies))
                ]
            )
        return self._velocities[self.index(name)]

    def drive(self, name, driver):
        """Attach an ipctk.demo.KinematicDriver to the KINEMATIC body ``name``."""
        self._drivers[self.index(name)] = driver

    def build_simulator(self, dt):
        if self.joints is not None and self._velocities is not None:
            raise ValueError("the Simulator takes either joints or initial velocities")
        kwargs = dict(
            bodies=self.bodies,
            initial_poses=self.poses,
            dt=dt,
            settings=self.settings,
        )
        if self.joints is not None:
            kwargs["joints"] = self.joints
        elif self._velocities is not None:
            kwargs["initial_velocities"] = self._velocities
        sim = ipctk.demo.Simulator(**kwargs)
        # Drivers attach after construction, once the body types are final.
        for i, driver in self._drivers.items():
            sim.set_kinematic_driver(i, driver)
        return sim


config = {}
if len(sys.argv) > 2:
    # Run the settings file as a module and keep its globals.
    config = runpy.run_path(sys.argv[2], run_name="__ipctk_settings__")
dt = config.get("dt", 1 / 60.0)

if len(sys.argv) > 1:
    if not hasattr(ipctk, "read_gltf"):
        sys.exit(
            "This ipctk build has no glTF support. Reconfigure it with "
            "-DIPC_TOOLKIT_WITH_GLTF=ON and rebuild."
        )
    bodies, initial_poses = ipctk.read_gltf(
        sys.argv[1],
        convert_planes=config.get("convert_planes", True),
        is_dof_fixed=config.get("is_dof_fixed", {}),
    )
else:
    mesh_names = [
        "bunny (lowpoly).ply",
        "bunny (lowpoly).ply",
        "bowl.ply",
    ]

    rest_positions = []
    edges = []
    faces = []
    for mesh_name in mesh_names:
        mesh = meshio.read(
            pathlib.Path(__file__).parents[2] / "tests" / "data" / mesh_name
        )
        rest_positions.append(mesh.points.astype("float64"))
        edges.append(ipctk.edges(mesh.cells_dict["triangle"]))
        faces.append(mesh.cells_dict["triangle"])

    initial_poses = ipctk.Poses(
        [
            ipctk.Pose(np.array([1.0, 1.5, 0]), np.random.random(3)),
            ipctk.Pose(np.array([-1.0, 2.0, 0.0]), np.array([0.0, np.pi / 4, 0.0])),
            ipctk.Pose(np.array([0.0, 1.1, 0.0]), np.zeros(3)),
        ]
    )

    bodies = ipctk.RigidBodies(
        rest_positions=rest_positions,
        edges=edges,
        faces=faces,
        densities=np.full(len(mesh_names), 1000.0),
        initial_poses=initial_poses,
    )

    # Add a ground plane at y = 0 for the bodies to rest on. `planes` is a
    # read/write property that returns a copy, so assign the whole list.
    bodies.planes = [ipctk.Hyperplane(np.array([0.0, 1.0, 0.0]), np.zeros(3))]

scene = Scene(bodies, initial_poses)
if "configure" in config:
    config["configure"](scene)

ps.init()

if len(bodies.planes) > 0:
    ps.set_ground_plane_height_mode("manual")
    ps.set_ground_plane_height(bodies.planes[0].origin()[1])
else:
    ps.set_ground_plane_mode("none")

ps.set_give_focus_on_show(True)

ps_mesh = ps.register_surface_mesh(
    "rigid body",
    bodies.vertices(initial_poses),
    bodies.faces,
)

# The bodies are merged into a single surface mesh, so map each face back to the
# body that owns it (via any of its vertices) to resolve picks to a body.
face_body = np.array([bodies.vertex_to_body(int(v)) for v in bodies.faces[:, 0]])

DYNAMIC_COLOR = np.array([0.85, 0.85, 0.9])
STATIC_COLOR = np.array([0.85, 0.3, 0.3])


def refresh_body_colors():
    """Color each body's faces by its type (static bodies stand out)."""
    body_colors = np.array(
        [
            STATIC_COLOR
            if bodies[i].type == ipctk.RigidBody.Type.STATIC
            else DYNAMIC_COLOR
            for i in range(len(bodies))
        ]
    )
    ps_mesh.add_color_quantity(
        "body type", body_colors[face_body], defined_on="faces", enabled=True
    )


refresh_body_colors()

ps_com = ps.register_point_cloud(
    "rigid body coms",
    np.vstack([initial_poses[i].position.reshape(-1, 3) for i in range(len(bodies))]),
)

for d in range(3):
    dim = np.zeros((len(bodies), 3))
    for i in range(len(bodies)):
        dim[i, d] = 1  # 100 * bodies[i].moment_of_inertia[d]
        R = Rotation.from_rotvec(initial_poses[i].rotation.copy()).as_matrix()
        dim[i, :] = R @ dim[i, :]
    ps_com.add_vector_quantity(
        "xyz"[d],
        dim,
        enabled=True,
        vectortype="standard",
    )

sim = scene.build_simulator(dt)

playing = False

# Screen position where the right mouse button was pressed (to tell a click,
# which converts a body to static, from a drag, which pans the camera).
right_click_press = None


def _mouse_screen_coords():
    p = imgui.GetIO().MousePos
    try:
        return (float(p.x), float(p.y))
    except AttributeError:
        return (float(p[0]), float(p[1]))


def make_body_static_under_cursor(screen_coords):
    """Convert the body under the cursor to a STATIC (immovable) body."""
    pick = ps.pick(screen_coords=screen_coords)
    if not pick.is_hit or pick.structure_name != "rigid body":
        return

    # Resolve the picked element to one of the merged mesh's vertices, then to
    # the body that owns it.
    element_type = pick.structure_data.get("element_type")
    index = pick.structure_data.get("index")
    if index is None:
        return
    if element_type == "face":
        vertex = int(bodies.faces[index][0])
    elif element_type == "vertex":
        vertex = int(index)
    elif element_type == "edge":
        vertex = int(bodies.edges[index][0])
    else:
        return

    body = bodies.vertex_to_body(vertex)
    if bodies[body].type != ipctk.RigidBody.Type.STATIC:
        # The simulator rebuilds its DOF mask from the body types every step,
        # so the body freezes at its current pose on the next step.
        bodies[body].convert_to_static()
        print(f"Body {body} converted to STATIC")
        refresh_body_colors()


def update_mesh():
    # sim.rigid_poses / sim.poses return only the current frame (O(num_bodies));
    # the *_history properties copy the entire growing history on each access.
    poses = sim.rigid_poses
    affine_poses = sim.poses
    ps_mesh.update_vertex_positions(bodies.vertices(poses))
    ps_com.update_point_positions(
        np.vstack([poses[i].position.reshape(-1, 3) for i in range(len(bodies))])
    )
    for d in range(3):
        dim = np.zeros((len(bodies), 3))
        for i in range(len(bodies)):
            dim[i, d] = 1  # 100 * bodies[i].moment_of_inertia[d]
            R = np.array(affine_poses[i].rotation)
            dim[i, :] = R @ dim[i, :]
        ps_com.add_vector_quantity(
            "xyz"[d],
            dim,
            enabled=True,
            vectortype="standard",
        )


def callback():
    global playing, right_click_press
    if imgui.Button("Play" if not playing else "Pause") or imgui.IsKeyPressed(
        imgui.ImGuiKey_Space
    ):
        playing = not playing
    imgui.SameLine()
    if imgui.Button("Step") or playing:
        if not sim.step():
            playing = False
        update_mesh()
    imgui.SameLine()
    if imgui.Button("Reset"):
        sim.reset()
        update_mesh()
        # reset() restores the original body types, so update the colors too.
        refresh_body_colors()
    imgui.SameLine()
    if imgui.Button("Save"):
        ipctk.write_gltf(
            "output.glb",
            bodies,
            sim.rigid_pose_history,
            dt,
        )

    imgui.TextUnformatted("Right-click a body to make it static.")

    # Right-click a body to make it static. A right *drag* still pans the
    # camera, so only convert on a click that barely moved.
    io = imgui.GetIO()
    if imgui.IsMouseClicked(imgui.ImGuiMouseButton_Right):
        right_click_press = _mouse_screen_coords()
    if (
        imgui.IsMouseReleased(imgui.ImGuiMouseButton_Right)
        and right_click_press is not None
    ):
        release = _mouse_screen_coords()
        moved_sq = (release[0] - right_click_press[0]) ** 2 + (
            release[1] - right_click_press[1]
        ) ** 2
        if moved_sq < 25.0 and not io.WantCaptureMouse:
            make_body_static_under_cursor(right_click_press)
        right_click_press = None


ps.set_user_callback(callback)

ps.show()
