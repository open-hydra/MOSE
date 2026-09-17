"""
Render the manufactured field (density, pressure, Mach number) as contour plots
for the V&V page.  The field is analytic, so this needs no solver output; it is
run offline and OUTPUT/mms-solution.svg is copied to docs/vv/images/.

Usage:  python plot_solution.py
"""
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

OM = 2.0 * np.pi
RHO0, ARHO = 1.0, 0.1
U0, AU = 40.0, 8.0
V0, AV = 30.0, 8.0
P0, AP = 8.0e3, 8.0e2
GAMMA, R = 1.4, 287.000001174

n = 200
x = np.linspace(0.0, 1.0, n)
y = np.linspace(0.0, 1.0, n)
X, Y = np.meshgrid(x, y)
rho = RHO0 + ARHO * np.sin(OM * X) * np.cos(OM * Y)
u = U0 + AU * np.sin(OM * X) * np.sin(OM * Y)
v = V0 + AV * np.cos(OM * X) * np.cos(OM * Y)
p = P0 + AP * np.cos(OM * X) * np.sin(OM * Y)
a = np.sqrt(GAMMA * p / rho)
mach = np.sqrt(u ** 2 + v ** 2) / a

mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

panels = [(r"density $\rho$", rho, "turbo"),
          (r"pressure $p$", p, "turbo"),
          (r"Mach number $M$", mach, "turbo")]
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
for ax, (title, field, cmap) in zip(axes, panels):
    cf = ax.pcolormesh(X, Y, field, cmap=cmap, shading="gouraud", rasterized=True)
    ax.contour(X, Y, field, levels=12, colors="k", linewidths=0.3, alpha=0.4)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=22)
    ax.set_xlabel("x", fontsize=20); ax.set_ylabel("y", fontsize=20)
    ax.tick_params(labelsize=20)
    cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=20)
plt.tight_layout()
os.makedirs("OUTPUT", exist_ok=True)
plt.savefig(os.path.join("OUTPUT", "mms-solution.svg"), bbox_inches="tight", transparent=True)
print("wrote OUTPUT/mms-solution.svg")
