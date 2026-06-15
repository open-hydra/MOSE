"""
Shock Wave - Boundary Layer Interaction (SWBLI) validation case

Generates separate skin-friction comparison plots for the SA and SST
cases using the corresponding MOSE, SU2, Wind-US, and experiment data.

Usage:
    python verify.py           # saves figure to OUTPUT/verify_Cf.png
    python verify.py --plot    # shows interactive figure
"""

import argparse
import sys
from pathlib import Path
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import PowerNorm
import matplotlib as mpl
import warnings
warnings.filterwarnings("ignore")

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

# Increase font sizes
label_fontsize = 17
cbar_fontsize = 17

output_dir = "OUTPUT"

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

# -------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------
parser = argparse.ArgumentParser(description='Shock Wave - Boundary Layer Interaction (SWBLI) validation case')
parser.add_argument('--plot', action='store_true', help='Show interactive plots')
args = parser.parse_args()

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
    T = vn[6, :, :, 0]
    gamma = vn[7, :, :, 0]
    gas_constant = vn[8, :, :, 0]

    velocity_mag = np.sqrt(u ** 2 + v ** 2 + w ** 2)
    sound_speed = np.sqrt(np.maximum(gamma * gas_constant * T, 1.0e-30))
    mach = velocity_mag / sound_speed

    zone_xc.append(xc)
    zone_yc.append(yc)
    zone_mach.append(mach)

mach_min = 0.0
mach_max = 5.1

# =============================================================
# Plot Mach number
# =============================================================
fig, ax = plt.subplots(1, 1, figsize=(9, 4.2))

sc = None
all_xn, all_yn, all_mach = [], [], []

for xn, yn, vn in zip(x_, y_, var_):
    xn = np.asarray(xn)[:, :, 0]
    yn = np.asarray(yn)[:, :, 0]
    vn = np.asarray(vn)

    u = vn[1, :, :, 0]
    v = vn[2, :, :, 0]
    w = vn[3, :, :, 0]
    T = vn[6, :, :, 0]
    gamma        = vn[7, :, :, 0]
    gas_constant = vn[8, :, :, 0]

    velocity_mag = np.sqrt(u**2 + v**2 + w**2)
    sound_speed  = np.sqrt(np.maximum(gamma * gas_constant * T, 1e-30))
    mach = velocity_mag / sound_speed

    sc = ax.pcolormesh(
        xn, yn, mach,
        cmap="turbo",
        norm=PowerNorm(gamma=1.6, vmin=0, vmax=5),  # gamma < 1 clusters near 0
        shading="flat",
        edgecolors="none",
        linewidth=0,
        rasterized=True,
    )

    # save for contour lines
    xc = 0.5 * (xn[:-1, :] + xn[1:, :])
    xc = 0.5 * (xc[:, :-1] + xc[:, 1:])
    yc = 0.5 * (yn[:-1, :] + yn[1:, :])
    yc = 0.5 * (yc[:, :-1] + yc[:, 1:])
    all_xn.append(xc)
    all_yn.append(yc)
    all_mach.append(mach)

# ── contour lines ──────────────────────────────────────────
n_lines = 15
levels_c = np.linspace(mach_min, mach_max, n_lines)

for xc, yc, mach in zip(all_xn, all_yn, all_mach):
    ax.contour(
        xc, yc, mach,
        levels=levels_c,
        colors="black",
        linewidths=0.6,
        alpha=1.0,
    )

cbar = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.16, shrink=0.5, fraction=0.03)
cbar.set_label("Mach", fontsize=14)
cbar.ax.tick_params(labelsize=12)
ax.axis("off")
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "SWBLI-field.svg"), bbox_inches="tight", transparent=True, dpi=300)

# Extract freestream conditions from the last cell in the first zone

mu     = 1.1858685985e-5
U_inf   = var_[0][1][0,-1,0]
p_inf   = var_[0][4][0,-1,0]
rho_inf = var_[0][0][0,-1,0]
nu      = mu / rho_inf

print("Freestream conditions")
print(f"  U_inf   = {U_inf:.4f} m/s")
print(f"  rho_inf = {rho_inf:.4f} kg/m³")
print(f"  mu      = {mu:.2e} Pa·s")
print(f"  nu      = {nu:.4e} m²/s")
print(f"  Re/L    = {U_inf/nu:.1f} m⁻¹")

# -------------------------------------------------------------
# Read reference data (if available)
# -------------------------------------------------------------

def _load_xy_data(path):
    if not os.path.exists(path):
        return None, None
    try:
        ref_data = np.loadtxt(path)
    except ValueError:
        ref_data = np.loadtxt(path, delimiter=',')
    return ref_data[:, 0], ref_data[:, 1]

# -------------------------------------------------------------
# Helper: read a wall.tec-style file and return (x_cell, Cf)
# Cell-centered variables (indices 0-6):
#   0:y+  1:tauX  2:tauY  3:tauZ  4:pw  5:Tw  6:qw
# tauX is negative (wall drags fluid in -x), so take abs()
# -------------------------------------------------------------
def _read_wall_cf(path, rho_inf, U_inf):
    """Return (x_cell, Cf) arrays for x > 0 (plate region)."""
    xw_, yw_, zw_, varw_, _ = read_TEC(path)
    xw_n = xw_[0][:, 0, 0]
    xw_c = 0.5 * (xw_n[:-1] + xw_n[1:])
    tauX = (varw_[0][1][:, 0, 0])
    x_c  = xw_c
    Cf   = tauX / (0.5 * rho_inf * U_inf**2)
    return x_c, Cf

# SA / SST solutions
x_mose_sa, Cf_mose_sa = _read_wall_cf("reference/MOSE-SA.tec", rho_inf, U_inf)
x_mose_sst, Cf_mose_sst = _read_wall_cf("reference/MOSE-SST.tec", rho_inf, U_inf)

# SA reference data
x_ref_sa, Cf_ref_sa = _load_xy_data("reference/schulein.dat")
x_su2_sa, Cf_su2_sa = _load_xy_data("reference/SU2-SA.dat")
x_wind_sa, Cf_wind_sa = _load_xy_data("reference/wind-SA.dat")

# SST reference data
x_su2_sst, Cf_su2_sst = _load_xy_data("reference/SU2-SST.dat")
x_wind_sst, Cf_wind_sst = _load_xy_data("reference/wind-SST.dat")

# -------------------------------------------------------------
# Plots
# -------------------------------------------------------------

# --- Cf vs x ---
# Figure 1: Skin friction coefficient for SA
fig1 = plt.figure(figsize=(10, 6))
ax_cf = fig1.add_subplot(111)

if x_ref_sa is not None:
    ax_cf.plot(x_ref_sa, Cf_ref_sa, 'ko', ms=6, label='Schulein (exp)')

if x_su2_sa is not None:
    ax_cf.plot(x_su2_sa, Cf_su2_sa, 'r--', lw=3, label='SU2')

if x_wind_sa is not None:
    ax_cf.plot(x_wind_sa, Cf_wind_sa, 'b-.', lw=3, label='Wind-US')

ax_cf.plot(x_mose_sa, Cf_mose_sa, 'g', lw=4.0, label='MOSE')

ax_cf.set_xlim(0.32, 0.41)
ax_cf.set_ylim(-0.002, 0.0075)
ax_cf.set_xlabel(r'$x$  [m]', fontsize=label_fontsize)
ax_cf.set_ylabel(r'$C_f$', fontsize=label_fontsize)
ax_cf.tick_params(labelsize=label_fontsize - 2)
ax_cf.legend(loc='best', fontsize=label_fontsize - 2)
ax_cf.grid(True, alpha=0.3)

plt.tight_layout()

out_fig1 = Path(os.path.join(output_dir, "SWBLI-cf-sa.svg"))
fig1.savefig(out_fig1, dpi=150)

# Figure 2: Skin friction coefficient for SST
fig2 = plt.figure(figsize=(10, 6))
ax_cf = fig2.add_subplot(111)

if x_ref_sa is not None:
    ax_cf.plot(x_ref_sa, Cf_ref_sa, 'ko', ms=6, label='Schulein (exp)')

if x_su2_sst is not None:
    ax_cf.plot(x_su2_sst, Cf_su2_sst, 'r--', lw=3, label='SU2')

if x_wind_sst is not None:
    ax_cf.plot(x_wind_sst, Cf_wind_sst, 'b-.', lw=3, label='Wind-US')

ax_cf.plot(x_mose_sst, Cf_mose_sst, 'g', lw=4.0, label='MOSE')

ax_cf.set_xlim(0.32, 0.41)
ax_cf.set_ylim(-0.002, 0.0075)
ax_cf.set_xlabel(r'$x$  [m]', fontsize=label_fontsize)
ax_cf.set_ylabel(r'$C_f$', fontsize=label_fontsize)
ax_cf.tick_params(labelsize=label_fontsize - 2)
ax_cf.legend(loc='best', fontsize=label_fontsize - 2)
ax_cf.grid(True, alpha=0.3)

plt.tight_layout()

out_fig2 = Path(os.path.join(output_dir, "SWBLI-cf-sst.svg"))
fig2.savefig(out_fig2, dpi=150)

print(f"\nFigures saved to {out_fig1}")
print(f"Figures saved to {out_fig2}")
