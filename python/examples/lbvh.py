from find_ipctk import ipctk
import meshio
import polyscope as ps
from polyscope import imgui
import numpy as np

import pathlib

mesh = meshio.read(pathlib.Path(
    __file__).parents[2] / "tests/data/puffer-ball/20.ply")

lbvh = ipctk.LBVH()
lbvh.build(mesh.points, np.array([], dtype=int), mesh.cells_dict["triangle"])

ps.init()

ps.set_give_focus_on_show(True)

ps_mesh = ps.register_surface_mesh(
    "bunny",
    mesh.points,
    mesh.cells_dict["triangle"]
)

nodes = lbvh.vertex_nodes


def traverse_lbvh(node, max_depth):
    if node.is_inner and max_depth > 0:
        V_left, E_left = traverse_lbvh(nodes[node.left], max_depth - 1)
        V_right, E_right = traverse_lbvh(nodes[node.right], max_depth - 1)
        return np.vstack([V_left, V_right]), np.vstack([E_left, E_right + V_left.shape[0]])

    E = np.array([
        [0, 1],
        [0, 2],
        [0, 3],
        [1, 5],
        [1, 4],
        [2, 4],
        [2, 6],
        [3, 5],
        [3, 6],
        [7, 4],
        [7, 5],
        [7, 6],
    ])
    V = np.array([
        node.aabb_min,
        [node.aabb_min[0], node.aabb_min[1], node.aabb_max[2]],
        [node.aabb_min[0], node.aabb_max[1], node.aabb_min[2]],
        [node.aabb_max[0], node.aabb_min[1], node.aabb_min[2]],
        [node.aabb_min[0], node.aabb_max[1], node.aabb_max[2]],
        [node.aabb_max[0], node.aabb_min[1], node.aabb_max[2]],
        [node.aabb_max[0], node.aabb_max[1], node.aabb_min[2]],
        node.aabb_max
    ])
    return V, E


max_depth = 0
bvh_nodes, bvh_edges = traverse_lbvh(nodes[0], max_depth=max_depth)

ps.register_curve_network("bvh", bvh_nodes, bvh_edges)


def foo():
    global max_depth
    changed, max_depth = imgui.SliderInt("max depth", max_depth, 0, 20)
    if changed:
        bvh_nodes, bvh_edges = traverse_lbvh(nodes[0], max_depth=max_depth)
        ps.register_curve_network("bvh", bvh_nodes, bvh_edges)


ps.set_user_callback(foo)

ps.show()
