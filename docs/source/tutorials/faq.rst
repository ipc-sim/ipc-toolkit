.. _faq:

Frequently Asked Questions
==========================

.. role:: cpp(code)
   :language: c++
.. role:: cmake(code)
   :language: cmake
.. role:: python(code)
   :language: python

How do I include IPC Toolkit in my project?
-------------------------------------------

If you are using CMake, the public include directory is added to the `ipc::toolkit` cmake target which means that any lib/bin that includes `ipc::toolkit` as a dependency also adds those include directories too.

If you are not using CMake, the include path is ``src``.

Files are included with the prefix :cpp:`#include <ipc/...>` in C++ and :python:`import ipctk` in Python.

How do I determine which edges intersect?
-----------------------------------------

We do not provide an edge-edge intersection function in 3D, but you can approximate it by computing the distance between the two edges and checking if it is less than a threshold.

Inside any :cpp:`BroadPhase` classes, the function :cpp:`detect_edge_edge_candidates` determines which edges intersect based on their bounding boxes. You can then use :cpp:`edge_edge_distance` to check if they approximatly intersect by computing the distance.

How do I build the edge matrix from the face matrix?
----------------------------------------------------

To build the edge matrix you can use :cpp:`igl::edges(faces, edges);` in C++ or :python:`ipctk.edges(faces)` in Python.

Is there a way to ignore select collisions?
-------------------------------------------

Yes. Both :cpp:`CollisionMesh::can_collide` and :cpp:`BroadPhase::can_vertices_collide` are :cpp:`CollisionFilter` objects — composable predicates of the form :cpp:`bool(size_t vi, size_t vj)` that control which vertex pairs enter the collision pipeline.

When building candidates through :cpp:`Candidates::build`, the broad phase automatically inherits :cpp:`CollisionMesh::can_collide`, so you only need to set it once on the mesh.

**How primitive-level checks work**

:cpp:`BroadPhase` expands vertex-level decisions to primitive pairs. For example, checking whether vertex ``vi`` can collide with triangle ``f = (vj, vk, vl)`` evaluates:

.. code-block::

    can_face_vertex_collide(f, vi) :=
        can_vertices_collide(vj, vi)
        && can_vertices_collide(vk, vi)
        && can_vertices_collide(vl, vi)

This means the filter acts on the one-ring of a vertex rather than a single primitive pair. For finer control you can subclass :cpp:`BroadPhase` and override the virtual methods:

.. code-block:: c++

    virtual bool can_edge_vertex_collide(size_t ei, size_t vi) const;
    virtual bool can_edges_collide(size_t eai, size_t ebi) const;
    virtual bool can_face_vertex_collide(size_t fi, size_t vi) const;
    virtual bool can_edge_face_collide(size_t ei, size_t fi) const;
    virtual bool can_faces_collide(size_t fai, size_t fbi) const;

:cpp:`CollisionFilter` wraps any :cpp:`bool(size_t, size_t)` callable and supports logical composition via ``|`` (union), ``&`` (intersection), and ``!``/``~`` (negation).

The available factory functions are:

- :cpp:`make_connected_components_filter(faces)` — blocks pairs within the same connected component (prevents self-collision).
- :cpp:`make_static_obstacle_filter(n_dynamic)` — blocks static-vs-static pairs; vertices with index ``>= n_dynamic`` are considered static.
- :cpp:`make_vertex_patches_filter(patch_ids)` — blocks pairs that share the same integer patch label.
- :cpp:`make_sparse_filter(explicit_values, default_value)` — sparse explicit overrides with a fallback default (Python only).

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/collision_filter.hpp>

            // Built-in factory functions
            CollisionFilter no_self   = make_connected_components_filter(mesh.faces());
            CollisionFilter no_static = make_static_obstacle_filter(n_dynamic_verts);
            CollisionFilter by_patch  = make_vertex_patches_filter(patch_ids);

            // Combine with | (union), & (intersection), ! (negation)
            mesh.can_collide = no_self & no_static;

            // Or construct directly from a lambda
            mesh.can_collide = CollisionFilter([&](size_t vi, size_t vj) {
                return group[vi] != group[vj];
            });

    .. md-tab-item:: Python

        .. code-block:: python

            import numpy as np
            from ipctk import (
                CollisionFilter,
                make_connected_components_filter,
                make_sparse_filter,
                make_static_obstacle_filter,
                make_vertex_patches_filter,
            )

            # Built-in factory functions
            no_self   = make_connected_components_filter(mesh.faces)
            no_static = make_static_obstacle_filter(n_dynamic)
            by_patch  = make_vertex_patches_filter(np.array([0, 0, 1, 1], dtype=np.int32))

            # Composition: &, |, ~
            mesh.can_collide = no_self & no_static

            # Sparse explicit overrides (e.g. always allow pair (2, 5))
            mesh.can_collide = make_sparse_filter({(2, 5): True}, default_value=False)

            # Arbitrary callable (slower — prefer factory functions for large meshes)
            mesh.can_collide = CollisionFilter(lambda i, j: group[i] != group[j])

My question is not answered here. What should I do?
---------------------------------------------------

Please open an Q&A discussion post on `GitHub <https://github.com/ipc-sim/ipc-toolkit/discussions>`_ and we will do our best to help you.