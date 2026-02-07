import pathlib

import meshio
import numpy as np
import polyscope as ps
from find_ipctk import ipctk
from polyscope import imgui
from scipy.spatial.transform import Rotation

mesh_names = [
    "bunny (lowpoly).ply",
    "bunny (lowpoly).ply",
    "bowl.ply",
]

rest_positions = []
edges = []
faces = []
for mesh_name in mesh_names:
    mesh = meshio.read(pathlib.Path(__file__).parents[2] / "tests" / "data" / mesh_name)
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
    densities=np.ones(len(mesh_names)),
    initial_poses=initial_poses,
)

# for pose in initial_poses:
#     print(pose)

ps.init()

ps.set_ground_plane_height_mode("manual")
ps.set_ground_plane_height(0.0)

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
        dim[i, d] = 100 * bodies[i].moment_of_inertia[d]
        R = Rotation.from_rotvec(initial_poses[i].rotation.copy()).as_matrix()
        dim[i, :] = R @ dim[i, :]
    ps_com.add_vector_quantity(
        "xyz"[d],
        dim,
        enabled=True,
        vectortype="ambient",
    )

sim = ipctk.Simulator(
    bodies=bodies,
    initial_poses=initial_poses,
    dt=1 / 60.0,
)

playing = False


def update_mesh():
    ps_mesh.update_vertex_positions(bodies.vertices(sim.pose_history[-1]))
    ps_com.update_point_positions(
        np.vstack(
            [
                sim.pose_history[-1][i].position.reshape(-1, 3)
                for i in range(len(bodies))
            ]
        )
    )
    for d in range(3):
        dim = np.zeros((len(bodies), 3))
        for i in range(len(bodies)):
            dim[i, d] = 100 * bodies[i].moment_of_inertia[d]
            R = Rotation.from_rotvec(
                sim.pose_history[-1][i].rotation.copy()
            ).as_matrix()
            dim[i, :] = R @ dim[i, :]
        ps_com.add_vector_quantity(
            "xyz"[d],
            dim,
            enabled=True,
            vectortype="ambient",
        )


def callback():
    global playing
    if imgui.Button("Play" if not playing else "Pause") or imgui.IsKeyPressed(
        imgui.ImGuiKey_Space
    ):
        playing = not playing
    imgui.SameLine()
    if imgui.Button("Step") or playing:
        sim.step()
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
            sim.pose_history,
            1 / 60.0,
        )


ps.set_user_callback(callback)

ps.show()
