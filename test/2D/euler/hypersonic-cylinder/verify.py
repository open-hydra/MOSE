"""
Hypersonic Cylinder test case

Loads:
  - MOSE output : OUTPUT/field.tec  (Tecplot BLOCK, 2 structured zones)

Usage:
    python verify.py           # prints error table, saves figure to OUTPUT/verify.png
    python verify.py --plot    # shows interactive figure instead
"""

import argparse
import sys
import os
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
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
parser = argparse.ArgumentParser(
    description="Hypersonic Cylinder — Riemann solver comparison"
)
parser.add_argument("--plot", action="store_true", help="Show interactive plots")
args = parser.parse_args()

# -------------------------------------------------------------
# Normal shock relations  (beta = 90 deg, bow shock at stagnation line)
# -------------------------------------------------------------
def normal_shock(M1, gamma):
    """Return M2, p2/p1, T2/T1, rho2/rho1 across a normal shock."""
    gp1  = gamma + 1.0
    gm1  = gamma - 1.0
    M1sq = M1 * M1

    M2sq = (M1sq + 2.0/gm1) / (2.0*gamma/gm1 * M1sq - 1.0)
    M2   = math.sqrt(M2sq)

    p2p1 = 1.0 + 2.0*gamma/gp1 * (M1sq - 1.0)
    r2r1 = gp1*M1sq / (gm1*M1sq + 2.0)
    t2t1 = p2p1 / r2r1

    return M2, p2p1, t2t1, r2r1


def total_pressure(p, M, gamma):
    return p * (1.0 + 0.5*(gamma - 1.0)*M*M) ** (gamma / (gamma - 1.0))


def total_temperature(T, M, gamma):
    return T * (1.0 + 0.5*(gamma - 1.0)*M*M)


# -------------------------------------------------------------
# Hard-coded inputs
# -------------------------------------------------------------
M1    = 8.1
gamma = 1.4
p1    = 101325.0
T1    = 288.14

# Analytical normal-shock solution (bow shock at stagnation line)
M2_ref, p2p1_ref, t2t1_ref, r2r1_ref = normal_shock(M1, gamma)
p2_ref = p2p1_ref * p1
T2_ref = t2t1_ref * T1

# Total quantities behind the shock (isentropic stagnation from post-shock state)
p0_ref = total_pressure(p2_ref, M2_ref, gamma)
T0_ref = total_temperature(T2_ref, M2_ref, gamma)   # equals T1*(1+(g-1)/2*M1^2)

# Probe location in the solver grid
I_PROBE = 80
J_PROBE = 0

# -------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------
parser = argparse.ArgumentParser(description='Normal shock (bow shock) verification')
parser.add_argument('--plot', action='store_true', help='Plot density contours comparison')
args = parser.parse_args()

# -------------------------------------------------------------
# Solver files
# Variable map (x,y,z stripped by read_TEC):
#   var[0]=rho  var[1]=u  var[2]=v  var[3]=w
#   var[4]=p    var[5]=T  var[6]=g  var[7]=R
# -------------------------------------------------------------
solvers = ['HLLC', 'HLLC+', 'HLLE', 'SLAU']
files   = [f'OUTPUT/field-{s}.tec' for s in solvers]

# -------------------------------------------------------------
# Load data and probe total quantities at I=I_PROBE, J=J_PROBE
# -------------------------------------------------------------
fields_data = {}   # store read data per solver for optional plotting

print(f"Normal-shock analytical reference  (M1={M1}, gamma={gamma})")
print(f"  M2    = {M2_ref:.6f}")
print(f"  p2/p1 = {p2p1_ref:.6f}")
print(f"  T2/T1 = {t2t1_ref:.6f}")
print(f"  r2/r1 = {r2r1_ref:.6f}")
print(f"  p0    = {p0_ref:.4f} Pa")
print(f"  T0    = {T0_ref:.4f} K")
print()
print(f"{'Solver':<8}  {'p0 [Pa]':>14}  {'err p0 [%]':>10}  "
      f"{'T0 [K]':>10}  {'err T0 [%]':>10}")
print("-" * 62)

for solver, fpath in zip(solvers, files):
    [x_, y_, z_, var_, vnames] = read_TEC(fpath)
    fields_data[solver] = {'x0': x_[0], 'y0': y_[0], 'p0': var_[0][4], 'rho0': var_[0][0]}

    # Extract primitive variables at probe cell
    p_pr = var_[0][4][I_PROBE, J_PROBE, 0]
    T_pr = var_[0][5][I_PROBE, J_PROBE, 0]
    u_pr = var_[0][1][I_PROBE, J_PROBE, 0]
    v_pr = var_[0][2][I_PROBE, J_PROBE, 0]
    w_pr = var_[0][3][I_PROBE, J_PROBE, 0]
    g_pr = var_[0][6][I_PROBE, J_PROBE, 0]   # local gamma
    R_pr = var_[0][7][I_PROBE, J_PROBE, 0]   # local gas constant

    speed = math.sqrt(u_pr**2 + v_pr**2 + w_pr**2)
    a_pr  = math.sqrt(g_pr * R_pr * T_pr)
    M_pr  = speed / a_pr

    p0_pr = total_pressure(p_pr, M_pr, g_pr)
    T0_pr = total_temperature(T_pr, M_pr, g_pr)

    err_p0 = 100.0 * (p0_pr - p0_ref) / p0_ref
    err_T0 = 100.0 * (T0_pr - T0_ref) / T0_ref

    print(f"{solver:<8}  {p0_pr:>14.4f}  {err_p0:>+10.4f}  "
          f"{T0_pr:>10.4f}  {err_T0:>+10.4f}")

print()

# -------------------------------------------------------------
# Optional plots
# -------------------------------------------------------------

# --- Figure 1: density contours ---
# Calculate aspect ratio from first solver to size the figure properly
d_first = fields_data[solvers[0]]
xn0_first = d_first['x0'][:, :, 0]; yn0_first = d_first['y0'][:, :, 0]
x_range = xn0_first.max() - xn0_first.min()
y_range = yn0_first.max() - yn0_first.min()
aspect_ratio = x_range / y_range

# Compute global vmin/vmax across all solvers for a shared colormap
vmin_global = min(
    (fields_data[s]['rho0'][:, :, 0] / fields_data[s]['rho0'][0, 0, 0]).min()
    for s in solvers
)
vmax_global = max(
    (fields_data[s]['rho0'][:, :, 0] / fields_data[s]['rho0'][0, 0, 0]).max()
    for s in solvers
)

# Create figure with size that preserves aspect ratio (1 row x N solvers)
subplot_width = 4
subplot_height = subplot_width / aspect_ratio
fig1, axs = plt.subplots(1, len(solvers), figsize=(subplot_width * len(solvers), subplot_height + 2.0))

# Increase font sizes
title_fontsize = 18
label_fontsize = 17
cbar_fontsize = 17

for ax, solver in zip(axs, solvers):
    d   = fields_data[solver]
    # Block 0
    xn0 = d['x0'][:, :, 0]; yn0 = d['y0'][:, :, 0]; vc0 = d['rho0'][:, :, 0]
    # Non-dimensional density
    vc0 = vc0 / d['rho0'][0, 0, 0]
    cf = ax.pcolormesh(xn0.T, yn0.T, vc0.T, shading='flat', cmap='OrRd', vmin=vmin_global, vmax=vmax_global, rasterized=True)
    xc0 = 0.5 * (0.5 * (xn0[:-1, :] + xn0[1:, :])[:, :-1] + 0.5 * (xn0[:-1, :] + xn0[1:, :])[:, 1:])
    yc0 = 0.5 * (0.5 * (yn0[:-1, :] + yn0[1:, :])[:, :-1] + 0.5 * (yn0[:-1, :] + yn0[1:, :])[:, 1:])
    ax.contour(xc0.T, yc0.T, vc0.T, levels=15, colors='black', linewidths=0.5, alpha=0.6)
    ax.set_title(f'{solver}', fontsize=title_fontsize, fontweight='bold')
    ax.set_xlabel('x', fontsize=label_fontsize)
    ax.set_ylabel('y', fontsize=label_fontsize)
    ax.tick_params(labelsize=label_fontsize - 2)
    ax.set_aspect('equal')

#cbar = fig1.colorbar(cf, ax=axs.tolist(), label='ρ', orientation='horizontal', pad=0.20, fraction=0.046, aspect=50)
#cbar.ax.tick_params(labelsize=cbar_fontsize - 2)
#cbar.set_label('ρ', fontsize=cbar_fontsize)
plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_path = os.path.join(output_dir, f"HC-fields.svg")
    plt.savefig(out_path, bbox_inches="tight", transparent=True, dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")

# --- Figure 2: density profile along J=0 vs radial angle ---
# Angle is 0 at y=2 (stagnation), ±90 deg at the extremes of the domain
Y_CENTER = 2.0
# Determine the half-span from the first solver
d0  = fields_data[solvers[0]]
yn0 = d0['y0'][:, J_PROBE, 0]
yc0 = 0.5 * (yn0[:-1] + yn0[1:])
Y_HALF = max(abs(yc0.max() - Y_CENTER), abs(yc0.min() - Y_CENTER))

fig2, ax2 = plt.subplots(figsize=(5, 4))
for solver in solvers:
    d    = fields_data[solver]
    yn   = d['y0'][:, J_PROBE, 0]
    y_j0 = 0.5 * (yn[:-1] + yn[1:])
    ang  = 90.0 * (y_j0 - Y_CENTER) / Y_HALF
    r_j0 = d['rho0'][:, J_PROBE, 0]
    ax2.plot(ang, r_j0, label=solver, lw=3)
ax2.set_xlabel('Radial angle  [deg]')
ax2.set_ylabel('ρ  [kg/m³]')
ax2.set_xlim(-90, 90)
ax2.set_xticks(range(-90, 91, 15))
ax2.legend()
ax2.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_path = os.path.join(output_dir, f"HC-slice.svg")
    plt.savefig(out_path, bbox_inches="tight", transparent=True, dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")
