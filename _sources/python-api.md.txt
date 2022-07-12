# API

% .. automodule:: ipctk
% :members:
% :undoc-members:
% :show-inheritance:

## Main Functions

```{eval-rst}
.. autofunction:: ipctk.construct_constraint_set
```

```{eval-rst}
.. autofunction:: ipctk.compute_barrier_potential
```

```{eval-rst}
.. autofunction:: ipctk.compute_barrier_potential_gradient
```

```{eval-rst}
.. autofunction:: ipctk.compute_barrier_potential_hessian
```

## Collision Constraints

```{eval-rst}
.. autoclass:: ipctk.CollisionConstraint
    :members:
    :show-inheritance:
```

```{eval-rst}
.. autoclass:: ipctk.VertexVertexConstraint
    :members:
    :show-inheritance:
```

```{eval-rst}
.. autoclass:: ipctk.EdgeVertexConstraint
    :members:
    :show-inheritance:
```

```{eval-rst}
.. autoclass:: ipctk.EdgeEdgeConstraint
    :members:
    :show-inheritance:
```

```{eval-rst}
.. autoclass:: ipctk.FaceVertexConstraint
    :members:
    :show-inheritance:
```

```{eval-rst}
.. autoclass:: ipctk.Constraints
    :members:
```

## Barrier

```{eval-rst}
.. autofunction:: ipctk.barrier
```

```{eval-rst}
.. autofunction:: ipctk.barrier_gradient
```

```{eval-rst}
.. autofunction:: ipctk.barrier_hessian
```

## Distance

### Distance Type

```{eval-rst}
.. autoclass:: ipctk.PointEdgeDistanceType
```

```{eval-rst}
.. autoclass:: ipctk.EdgeEdgeDistanceType
```

```{eval-rst}
.. autoclass:: ipctk.PointTriangleDistanceType
```

```{eval-rst}
.. autofunction:: ipctk.point_edge_distance_type
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_distance_type
```

```{eval-rst}
.. autofunction:: ipctk.point_triangle_distance_type
```

### Edge-Edge Mollifier

```{eval-rst}
.. autofunction:: ipctk.edge_edge_mollifier_threshold
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_cross_squarednorm
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_cross_squarednorm_gradient
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_cross_squarednorm_hessian
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_mollifier
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_mollifier_gradient
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_mollifier_hessian
```

### Edge-Edge

```{eval-rst}
.. autofunction:: ipctk.edge_edge_distance
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_distance_gradient
```

```{eval-rst}
.. autofunction:: ipctk.edge_edge_distance_hessian
```

### Line-Line

```{eval-rst}
.. autofunction:: ipctk.line_line_distance
```

```{eval-rst}
.. autofunction:: ipctk.line_line_distance_gradient
```

```{eval-rst}
.. autofunction:: ipctk.line_line_distance_hessian
```

### Point-Edge

```{eval-rst}
.. autofunction:: ipctk.point_edge_distance
```

```{eval-rst}
.. autofunction:: ipctk.point_edge_distance_gradient
```

```{eval-rst}
.. autofunction:: ipctk.point_edge_distance_hessian
```

### Point-Line

```{eval-rst}
.. autofunction:: ipctk.point_line_distance
```

```{eval-rst}
.. autofunction:: ipctk.point_line_distance_gradient
```

```{eval-rst}
.. autofunction:: ipctk.point_line_distance_hessian
```

### Point-Plane

```{eval-rst}
.. autofunction:: ipctk.point_plane_distance
```

```{eval-rst}
.. autofunction:: ipctk.point_plane_distance_gradient
```

```{eval-rst}
.. autofunction:: ipctk.point_plane_distance_hessian
```

### Point-Point

```{eval-rst}
.. autofunction:: ipctk.point_point_distance
```

```{eval-rst}
.. autofunction:: ipctk.point_point_distance_gradient
```

```{eval-rst}
.. autofunction:: ipctk.point_point_distance_hessian
```

### Point-Triangle

```{eval-rst}
.. autofunction:: ipctk.point_triangle_distance
```

```{eval-rst}
.. autofunction:: ipctk.point_triangle_distance_gradient
```

```{eval-rst}
.. autofunction:: ipctk.point_triangle_distance_hessian
```

## CCD

### Broad-Phase

```{eval-rst}
.. autoclass:: ipctk.BroadPhaseMethod
```

#### Candidates

```{eval-rst}
.. autoclass:: ipctk.VertexVertexCandidate
    :members:
```

```{eval-rst}
.. autoclass:: ipctk.EdgeVertexCandidate
    :members:
```

```{eval-rst}
.. autoclass:: ipctk.EdgeEdgeCandidate
    :members:
```

```{eval-rst}
.. autoclass:: ipctk.EdgeFaceCandidate
    :members:
```

```{eval-rst}
.. autoclass:: ipctk.FaceVertexCandidate
    :members:
```

```{eval-rst}
.. autoclass:: ipctk.Candidates
```

### Narrow-Phase

## Utils

```{eval-rst}
.. autofunction:: ipctk.has_intersections
```

### Logger

```{eval-rst}
.. autoclass:: ipctk.LoggerLevel
```

```{eval-rst}
.. autofunction:: ipctk.set_logger_level
```

### Multi-Threading

```{eval-rst}
.. autofunction:: ipctk.get_num_threads
```

```{eval-rst}
.. autofunction:: ipctk.set_num_threads
```
