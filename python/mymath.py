import torch
import numpy as np
import matplotlib.pyplot as plt

def soft_min(x):
    return (torch.dot(torch.exp(-x), x) / torch.sum(torch.exp(-x))).item()

def soft_max(x):
    return (torch.dot(torch.exp(x), x) / torch.sum(torch.exp(x))).item()

# x: torch tensor
def spline3(x):
    branches = (lambda x: 0 * x, lambda x: (x+2)**3 / 6, lambda x: 2 / 3 - x**2 * (1 - torch.abs(x) / 2), lambda x: -(x-2)**3 / 6, lambda x: 0 * x)
    bounds = [-2,-1,1,2]  # boundaries between branches
    y = torch.zeros_like(x)
    for i in range(1, len(bounds)):
        mask = torch.all(torch.stack((x < bounds[i], x >= bounds[i-1])), dim=0)
        y[mask] = branches[i](x[mask])
    return y

# x: torch tensor
# eps: float, size of support
# r: float, > 1
def inv_barrier(x, eps, r):
    return spline3(x * (2 / eps)) / torch.abs(x)**r

def L_ns(x):
    branches = (lambda x: x, lambda x: x-1)
    y = torch.zeros_like(x)
    mask = x < 0
    y[mask] = branches[0](x[mask])

    mask = x > 1
    y[mask] = branches[1](x[mask])
    return y

# x: torch tensor
# dhat: float, size of support
def log_barrier(x, dhat):
    y = torch.zeros_like(x)
    mask = x < dhat
    x_scaled = x[mask] / dhat
    y[mask] = -torch.log(x_scaled) * (x_scaled - 1)**2
    return y

def percentile(L, ratio=0.999):
    X = np.sort(L.flatten())
    id = int(len(X)*ratio)
    return X[id]

if __name__ == "__main__":

    #######################################
    # Example of autograd on cubic spline #
    #######################################

    x = torch.linspace(-3, 3, 100, requires_grad=True)

    # Compute the function value
    y = spline3(x)
    plt.plot(x.detach().numpy(), y.detach().numpy(), label="spline")

    z = torch.sum(y)
    grad_x = torch.autograd.grad(z, x, create_graph=True)[0]
    plt.plot(x.detach().numpy(), grad_x.detach().numpy(), label="1st deriv")

    h = torch.sum(grad_x)
    hess_x = torch.autograd.grad(h, x, create_graph=True)[0]
    plt.plot(x.detach().numpy(), hess_x.detach().numpy(), '--', label="2nd deriv")

    plt.legend()
    plt.savefig("tmp.png")

    # # Visualization of cubic spline

    # x = torch.linspace(-2.5, 2.5, 500, requires_grad=True)
    # y = spline3(x)
    # E = torch.sum(y)
    # deriv1 = torch.autograd.grad(E, x, create_graph=True)[0]
    # E = torch.sum(deriv1)
    # deriv2 = torch.autograd.grad(E, x, create_graph=True)[0]

    # fig = go.Figure(data=[
    #     go.Scatter(x=x.detach().numpy(), y=y.detach().numpy(), name="y(x)"),
    #     go.Scatter(x=x.detach().numpy(), y=deriv1.detach().numpy(), name="y'(x)"),
    #     go.Scatter(x=x.detach().numpy(), y=deriv2.detach().numpy(), name="y''(x)")
    #     ], layout=go.Layout(width=400, height=400))
    # fig.update_yaxes(exponentformat = 'E')
    # fig.show()