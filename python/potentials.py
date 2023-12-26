from mymath import *
import torch
import torch.nn.functional as F

# given some uv sampled on the edge (e0, e1), compute the integrand inside the Psi(x) for point p outside of edge
# alpha > 0, the smaller the less smooth
# eps > 0, the smaller the less smooth
# r >= 1, the larger the less smooth
def point_edge_potential_exact_pointwise(points, e0, e1, eps, r, alpha, uv):
    n_pts = points.shape[0]
    n_samples = uv.shape[0]
    dim = points.shape[1]

    tangent = F.normalize(e1 - e0, p=2, dim=None)
    diff = (e0.view(1, -1) - points).view(-1, 2, 1) + torch.outer(e1 - e0, uv).view(1, 2, -1)
    dists = torch.norm(diff, dim=1)
    barrier_val = inv_barrier(dists, eps, r).view(n_pts, n_samples)
    Phi_val = spline3(torch.tensordot(tangent, diff, dims=[[0],[1]]).view(n_pts, n_samples) / dists * (2 / alpha))
    return barrier_val * Phi_val

# compute integral of point_edge_potential_exact_pointwise() using composite trapezoidal quadrature rule
def point_edge_potential_exact(points, e0, e1, eps, r, alpha, n_samples):
    uv = torch.linspace(0, 1, n_samples)
    integrand = point_edge_potential_exact_pointwise(points, e0, e1, eps, r, alpha, uv)
    weights = torch.ones((uv.shape[0])) * torch.norm(e1 - e0) / (n_samples - 1)
    weights[0] /= 2
    weights[-1] /= 2
    return torch.tensordot(integrand, weights, dims=[[1],[0]]).view(-1)

# discrete version of point edge potential Psi(x) using the smoothed closest point
# points: N x d
# smooth_factor >= 1, no smoothing at all if = 1
# alpha > 0, the smaller the less smooth
# eps > 0, the smaller the less smooth
# r >= 1, the larger the less smooth
def point_edge_potential_discrete(points, e0, e1, eps, r, alpha, smooth_factor):
    tangent = e1 - e0
    length = torch.norm(tangent)
    normal = torch.tensor([-tangent[1], tangent[0]]) / length

    # uv of point projected onto the line
    s = ((points - e0.view(1, -1)) @ (tangent / length**2)).view(-1)
    L = L_ns(s)**smooth_factor
    # smooth closest point
    y = e0 + torch.outer((s - L).view(-1), tangent)

    Phi = (F.normalize(y - points, p=2, dim=1) @ tangent.view(-1,1)).view(-1) ** 2 / length**2
    dist_sqr = ((points - e0) @ normal.view(-1,1)).view(-1) **2 + length**2 * L**2
    
    return torch.norm(e1 - e0) * spline3(Phi * (2 / alpha)) * inv_barrier(dist_sqr, eps, r)

# compute the convergent IPC potential wrt. a closed curve prescribed by vertices on some spatial points
# points: N x d
# vertices: M x d
def convergent_ipc_potential(points, vertices, dhat):
    n_edges = vertices.shape[0] - 1
    values = torch.zeros((points.shape[0]))
    for e in range(n_edges):
        tangent = vertices[e+1] - vertices[e]
        length = torch.norm(tangent)

        s = ((points - vertices[e].view(1, -1)) @ (tangent / length**2)).view(-1)
        y = vertices[e] + torch.outer((s - L_ns(s)).view(-1), tangent)

        dist = torch.norm(points - y.view(-1, points.shape[1]), dim=1).view(-1)
        values += log_barrier(dist, dhat)
    
    for v in range(n_edges+1):
        dist = torch.norm(points - vertices[v, :], dim=1).view(-1)
        values -= log_barrier(dist, dhat) * (0.5 if v in [0, n_edges] else 1)

    return values