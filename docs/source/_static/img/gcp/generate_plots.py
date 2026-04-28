#!/usr/bin/env python3
"""Generate barriers.svg and heaviside.svg plots for GCP documentation."""

import matplotlib
import numpy as np

matplotlib.use("Agg")
import os

import matplotlib.pyplot as plt

# Publication-quality settings
plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "font.size": 10,
        "axes.labelsize": 12,
        "legend.fontsize": 9,
        "lines.linewidth": 1.8,
    }
)

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def B3(v):
    """Cubic B-spline basis function."""
    v = np.abs(v)
    result = np.zeros_like(v, dtype=float)
    mask1 = v < 1
    mask2 = (v >= 1) & (v < 2)
    result[mask1] = 2.0 / 3.0 - v[mask1] ** 2 + v[mask1] ** 3 / 2.0
    result[mask2] = (1.0 / 6.0) * (2.0 - v[mask2]) ** 3
    return result


def he(x, eps):
    """Spline envelope: he(x, eps) = B3(2x/eps) / B3(0)."""
    B3_0 = 2.0 / 3.0
    return B3(2.0 * x / eps) / B3_0


def ipc_barrier(d, dhat):
    """IPC log barrier: b(d, dhat) = -(d - dhat)^2 * ln(d / dhat) for 0 < d <= dhat."""
    result = np.zeros_like(d, dtype=float)
    mask = (d > 0) & (d <= dhat)
    result[mask] = -((d[mask] - dhat) ** 2) * np.log(d[mask] / dhat)
    return result


def heaviside_H(z):
    """
    Piecewise cubic Heaviside approximation:
      H(z) = 0                                 for z < -3
           = (1/6)(3+z)^3                      for -3 <= z < -2
           = (1/6)(3 - 9z - 9z^2 - 2z^3)      for -2 <= z < -1
           = 1 + (1/6)z^3                      for -1 <= z < 0
           = 1                                 for z >= 0
    """
    z = np.asarray(z, dtype=float)
    result = np.zeros_like(z)
    m1 = (z >= -3) & (z < -2)
    m2 = (z >= -2) & (z < -1)
    m3 = (z >= -1) & (z < 0)
    m4 = z >= 0
    result[m1] = (1.0 / 6.0) * (3.0 + z[m1]) ** 3
    result[m2] = (1.0 / 6.0) * (3.0 - 9.0 * z[m2] - 9.0 * z[m2] ** 2 - 2.0 * z[m2] ** 3)
    result[m3] = 1.0 + (1.0 / 6.0) * z[m3] ** 3
    result[m4] = 1.0
    return result


def H_ab(z, alpha, beta):
    """Scaled Heaviside: H_ab(z, alpha, beta) = H(3*(z - beta)/(alpha + beta))."""
    return heaviside_H(3.0 * (z - beta) / (alpha + beta))


# =============================================================================
# Plot 1: barriers.svg
# =============================================================================
def generate_barriers():
    fig, ax1 = plt.subplots(1, 1, figsize=(3, 2.5))

    # -- Left subplot: Barrier comparison --
    x = np.linspace(1e-5, 1.5, 2000)
    eps = 1.0
    dhat = 1.0

    # Curves
    h_vals = he(x, eps)
    p_vals = h_vals / x
    b_vals = ipc_barrier(x, dhat)

    ax1.plot(x, p_vals, color="#0072BD", label=r"$p_\epsilon$")
    ax1.plot(x, b_vals, color="#D95319", label=r"$p^{\mathrm{IPC}}$")
    ax1.plot(
        x, h_vals, color="black", linestyle="--", linewidth=1.4, label=r"$h_\epsilon$"
    )

    ax1.set_xlim([0, 1.5])
    ax1.set_ylim([-0.2, 3.5])
    ax1.set_xticks([0, 0.5, 1.0, 1.5])
    ax1.set_yticks([0, 1, 2, 3])
    ax1.set_xlabel(r"$z$")
    ax1.legend(loc="upper right", frameon=False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # -- Right subplot: Convergence as eps -> 0 --
    # x2 = np.linspace(1e-5, 1.5, 2000)

    # eps_vals = [1.0, 0.8, 0.5]
    # colors = ["#D95319", "#EDB120", "#7E2F8E"]
    # labels = [r"$\epsilon=1$", r"$\epsilon=0.8$", r"$\epsilon=0.5$"]

    # Vertical line at x=0 representing the discontinuous limit
    # ax2.plot([0, 0], [0, 3.5], color="#0072BD", linewidth=1.8, label="discontinuous")

    # for eps_i, col, lab in zip(eps_vals, colors, labels):
    #     h_i = he(x2, eps_i)
    #     p_i = h_i / x2
    #     ax2.plot(x2, p_i, color=col, label=lab)
    #     # Dashed vertical line at x=eps_i from bottom to y=1
    #     ax2.plot([eps_i, eps_i], [-0.2, 1.0], color=col, linestyle="--", linewidth=1.0)

    # ax2.set_xlim([0, 1.5])
    # ax2.set_ylim([-0.2, 3.5])
    # ax2.set_xticks([0, 0.5, 1.0, 1.5])
    # ax2.set_yticks([0, 1, 2, 3])
    # ax2.set_xlabel(r"$z$")
    # ax2.legend(loc="upper right", frameon=False)
    # ax2.spines["top"].set_visible(False)
    # ax2.spines["right"].set_visible(False)

    fig.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "barriers.svg")
    fig.savefig(
        out_path, format="svg", bbox_inches="tight", pad_inches=0, transparent=True
    )
    plt.close(fig)
    print(f"Saved {out_path}")


# =============================================================================
# Plot 2: heaviside.svg
# =============================================================================
def generate_heaviside():
    fig, ax = plt.subplots(figsize=(3, 2.5))

    alpha = 0.5
    beta = 0.0

    z = np.linspace(-0.6, 0.1, 2000)
    y = H_ab(z, alpha, beta)

    ax.plot(z, y, color="#0072BD", linewidth=2.0)

    ax.set_xlim([-0.6, 0.1])
    ax.set_ylim([-0.01, 1.1])
    ax.set_xticks([-0.5, 0])
    ax.set_xticklabels([r"$-\alpha$", r"$0$"])
    ax.set_yticks([0, 1])
    ax.set_xlabel(r"$z$")

    # Label the curve near the right end
    ax.annotate(
        r"$H^\alpha(z)$",
        xy=(0.02, 1.0),
        xytext=(0.04, 0.75),
        fontsize=11,
        color="#0072BD",
        arrowprops=dict(arrowstyle="->", color="#0072BD", lw=1.2),
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "heaviside.svg")
    fig.savefig(
        out_path, format="svg", bbox_inches="tight", pad_inches=0, transparent=True
    )
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    generate_barriers()
    generate_heaviside()
    print("Done.")
