import pathlib
import sys

import meshio
import numpy as np
import polyscope as ps
from find_ipctk import ipctk
from polyscope import imgui
from scipy.spatial.transform import Rotation

ipctk.set_logger_level(ipctk.debug)

if len(sys.argv) > 1:
    bodies, initial_poses = ipctk.read_gltf(sys.argv[1], True)
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

# for pose in initial_poses:
#     print(pose)

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

settings = ipctk.demo.Simulator.Settings()
settings.body_dynamics = ipctk.demo.Simulator.BodyDynamics.RIGID
settings.bdf_order = 2

sim = ipctk.demo.Simulator(
    bodies=bodies,
    initial_poses=initial_poses,
    dt=1 / 60.0,
    settings=settings,
)

playing = False


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
    global playing
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
    imgui.SameLine()
    if imgui.Button("Save"):
        ipctk.write_gltf(
            "output.glb",
            bodies,
            sim.rigid_pose_history,
            1 / 60.0,
        )


ps.set_user_callback(callback)

ps.show()
