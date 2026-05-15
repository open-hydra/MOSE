"""Plot Mach number from MOSE output (OUTPUT/field.tec)."""

import sys
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

output_dir = "OUTPUT"
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------------------
# Locate ORION Python package
# -------------------------------------------------------------
root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break

try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "ORION package not found. Run 'conda activate base' or "
        "install ORION in the current Python environment."
    )

# =============================================================
# Load MOSE solution (OUTPUT/field.tec)
# =============================================================
# Variable layout (x,y,z stripped by read_TEC):
#   0=rho  1=u   2=v   3=w   4=p   5=T   6=gamma  7=R
x_, y_, z_, var_, _ = read_TEC("OUTPUT/field.tec")

zone_xc, zone_yc, zone_mach = [], [], []

for zone in range(len(x_)):
    xn = np.asarray(x_[zone])
    yn = np.asarray(y_[zone])
    vn = np.asarray(var_[zone])

    # 2-D cell centres from nodal coordinates
    xc = 0.25 * (xn[:-1, :-1, 0] + xn[1:, :-1, 0] +
                 xn[:-1,  1:, 0] + xn[1:,  1:, 0])
    yc = 0.25 * (yn[:-1, :-1, 0] + yn[1:, :-1, 0] +
                 yn[:-1,  1:, 0] + yn[1:,  1:, 0])

    u = vn[1, :, :, 0]
    v = vn[2, :, :, 0]
    w = vn[3, :, :, 0]
    T = vn[5, :, :, 0]
    gamma = vn[6, :, :, 0]
    gas_constant = vn[7, :, :, 0]

    velocity_mag = np.sqrt(u ** 2 + v ** 2 + w ** 2)
    sound_speed = np.sqrt(np.maximum(gamma * gas_constant * T, 1.0e-30))
    mach = velocity_mag / sound_speed

    zone_xc.append(xc)
    zone_yc.append(yc)
    zone_mach.append(mach)

mach_min = min(np.min(m) for m in zone_mach)
mach_max = max(np.max(m) for m in zone_mach)
levels = np.linspace(mach_min, mach_max, 18)

# =============================================================
# Plot Mach number
# =============================================================
fig, ax = plt.subplots(1, 1, figsize=(9, 4.2))

for xc, yc, mach in zip(zone_xc, zone_yc, zone_mach):
    sc = ax.contourf(
        xc,
        yc,
        mach,
        levels=levels,
        cmap="turbo",
        antialiased=False,
    )

fig.colorbar(sc, ax=ax, label="Mach")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")

plt.tight_layout()

out_path = os.path.join(output_dir, "mach.svg")
plt.savefig(out_path, bbox_inches="tight", transparent=True)
print(f"  Figure saved to {out_path}")
