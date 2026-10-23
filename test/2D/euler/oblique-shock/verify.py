"""
Oblique Shock test case

Loads:
  - MOSE outputs : OUTPUT/field-<solver>.tec  (Tecplot BLOCK, 2 structured zones), one per Riemann solver

Usage:
    python verify.py           # prints the error table, saves OUTPUT/OS-fields.svg and OUTPUT/OS-slice.svg,
                               # exits with 1 when a solver is outside the threshold
    python verify.py --plot    # shows the figures instead of saving them
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
    description="Oblique Shock — Riemann solver comparison"
)
parser.add_argument("--plot", action="store_true", help="Show interactive plots")
args = parser.parse_args()

# -------------------------------------------------------------
# theta-beta-Mach function
# -------------------------------------------------------------
def tbm(M, beta, gamma):
    c1 = M*M * math.sin(beta)**2 - 1.0
    c2 = M*M * (gamma + math.cos(2.0*beta)) + 2.0
    return 2.0 * (1.0 / math.tan(beta)) * (c1 / c2)


# -------------------------------------------------------------
# Oblique shock solver
# -------------------------------------------------------------
def sckobl(mach1, delta, gamma, tol):
    theta1 = 1.2 * delta
    delt1 = tbm(mach1, theta1, gamma)
    theta2 = 1.02 * theta1

    while True:
        delt2 = tbm(mach1, theta2, gamma)

        if abs((delt2 - math.tan(delta)) / math.tan(delta)) <= tol:
            break

        theta = theta2 - (math.tan(delta) - delt2) * (theta2 - theta1) / (delt1 - delt2)
        theta1 = theta2
        delt1 = delt2
        theta2 = theta

    theta = theta2

    # Post-shock perfect-gas relations
    rmn1 = mach1 * math.sin(theta)
    gp1 = gamma + 1.0
    gm1 = gamma - 1.0

    r2r1 = gp1 * rmn1**2 / (gm1 * rmn1**2 + 2.0)
    p2p1 = 1.0 + 2.0 * gamma * (rmn1**2 - 1.0) / gp1

    t1 = (rmn1**2 + 2.0 / gm1)
    t2 = (2.0 * gamma * rmn1**2 / gm1 - 1.0)
    rmn2 = math.sqrt(t1 / t2)

    t2t1 = p2p1 / r2r1
    mach2 = rmn2 / math.sin(theta - delta)

    return mach2, theta, p2p1, t2t1, r2r1


# -------------------------------------------------------------
# Hard-coded inputs
# -------------------------------------------------------------
M1    = 4.0
gamma = 1.4
tol   = 1e-6

# Measure ramp angle directly from mesh block-1 wall nodes (j=0)
_tmp = read_TEC('OUTPUT/field-HLLC.tec')
_xw  = _tmp[0][1][:, 0, 0]    # x of wall nodes along i-axis (block 1)
_yw  = _tmp[1][1][:, 0, 0]    # y of wall nodes
delta_deg = float(np.degrees(np.arctan2(_yw[-1] - _yw[0], _xw[-1] - _xw[0])))
del _tmp, _xw, _yw

delta = math.radians(delta_deg)

# Analytical solution
mach2_ref, beta, p2p1_ref, t2t1_ref, r2r1_ref = sckobl(M1, delta, gamma, tol)

print("Oblique Shock – Analytical Solution")
print(f"  M1    = {M1:.4f}   delta = {delta_deg:.1f} deg   gamma = {gamma}")
print(f"  beta  = {math.degrees(beta):.4f} deg")
print(f"  M2    = {mach2_ref:.6f}")
print(f"  p2/p1 = {p2p1_ref:.6f}")
print(f"  T2/T1 = {t2t1_ref:.6f}")
print(f"  r2/r1 = {r2r1_ref:.6f}")
print()

# -------------------------------------------------------------
# Solver files
# Variable map (x,y,z stripped by read_TEC):
#   var[0]=rho  var[1]=u  var[2]=v  var[3]=w
#   var[4]=p    var[5]=T  var[6]=g  var[7]=R
#
# Block layout:
#   block 0 – flat pre-shock section   (I-1=20, J-1=218, K-1=1)
#   block 1 – ramp post-shock section  (I-1=38, J-1=218, K-1=1)
#     axis 0 = i (streamwise along ramp, i=0 at wedge tip)
#     axis 1 = j (wall-normal, j=0 at wall, j→ far-field)
# -------------------------------------------------------------
solvers = ['HLLC', 'HLLC+Tramel', 'HLLE', 'SLAU']
files   = [f'OUTPUT/field-{s}.tec' for s in solvers]

# Pre-shock probe: block 0, interior cells away from interface
i_pre0 = slice(5, 15)
j_pre0 = slice(5, 15)

# Post-shock probe: block 1, i≥14 ensures shock is beyond j=80 for all rows;
# j=5..50 stays well below the shock and avoids the immediate wall layer.
i_post1 = slice(14, 35)
j_post1 = slice(5, 50)

# Error thresholds (%)
tol_p = 1.5
tol_T = 1.5
tol_r = 1.5
tol_M = 1.5

# -------------------------------------------------------------
# Load, compare, check
# -------------------------------------------------------------
fields_data = {}   # store read data per solver for optional plotting

hdr = (f"{'Solver':<12} {'p2/p1':>10} {'err%':>7}  "
       f"{'T2/T1':>10} {'err%':>7}  "
       f"{'r2/r1':>10} {'err%':>7}  "
       f"{'M2':>8} {'err%':>7}")
print(hdr)
print("-" * len(hdr))

any_fail = False

for solver, fpath in zip(solvers, files):
    [x_, y_, z_, var_, vnames] = read_TEC(fpath)

    # --- Block 0: pre-shock ---
    f0  = var_[0]
    p0  = f0[4][i_pre0, j_pre0, 0]
    rh0 = f0[0][i_pre0, j_pre0, 0]
    T0  = f0[5][i_pre0, j_pre0, 0]

    p1   = float(np.mean(p0))
    rho1 = float(np.mean(rh0))
    T1   = float(np.mean(T0))

    # --- Block 1: post-shock ---
    f1  = var_[1]
    p1b = f1[4][i_post1, j_post1, 0]
    rh1 = f1[0][i_post1, j_post1, 0]
    T1b = f1[5][i_post1, j_post1, 0]
    u2b = f1[1][i_post1, j_post1, 0]
    v2b = f1[2][i_post1, j_post1, 0]
    g2b = f1[6][i_post1, j_post1, 0]
    R2b = f1[7][i_post1, j_post1, 0]

    p2   = float(np.mean(p1b))
    rho2 = float(np.mean(rh1))
    T2   = float(np.mean(T1b))
    u2   = float(np.mean(u2b))
    v2   = float(np.mean(v2b))
    g2   = float(np.mean(g2b))
    R2   = float(np.mean(R2b))

    p2p1 = p2   / p1
    r2r1 = rho2 / rho1
    t2t1 = T2   / T1
    a2   = math.sqrt(g2 * R2 * T2)
    m2   = math.sqrt(u2**2 + v2**2) / a2

    e_p = (p2p1 - p2p1_ref) / p2p1_ref * 100.0
    e_T = (t2t1 - t2t1_ref) / t2t1_ref * 100.0
    e_r = (r2r1 - r2r1_ref) / r2r1_ref * 100.0
    e_m = (m2   - mach2_ref) / mach2_ref * 100.0

    print(f"{solver:<12} {p2p1:>10.6f} {e_p:>7.2f}%  "
          f"{t2t1:>10.6f} {e_T:>7.2f}%  "
          f"{r2r1:>10.6f} {e_r:>7.2f}%  "
          f"{m2:>8.4f} {e_m:>7.2f}%")

    if abs(e_p) > tol_p:
        print(f"  FAIL: {solver} p2/p1 error {e_p:.2f}% > threshold {tol_p}%")
        any_fail = True
    if abs(e_T) > tol_T:
        print(f"  FAIL: {solver} T2/T1 error {e_T:.2f}% > threshold {tol_T}%")
        any_fail = True
    if abs(e_r) > tol_r:
        print(f"  FAIL: {solver} r2/r1 error {e_r:.2f}% > threshold {tol_r}%")
        any_fail = True
    if abs(e_m) > tol_M:
        print(f"  FAIL: {solver} M2    error {e_m:.2f}% > threshold {tol_M}%")
        any_fail = True

    fields_data[solver] = {
        'x0': x_[0], 'y0': y_[0], 'p0': var_[0][4], 'rho0': var_[0][0],
        'x1': x_[1], 'y1': y_[1], 'p1': var_[1][4], 'rho1': var_[1][0],
        'p2p1': p2p1, 'r2r1': r2r1, 't2t1': t2t1, 'm2': m2,
        'rho_pre': rho1,
    }

print()

# -------------------------------------------------------------
# Optional plots
# -------------------------------------------------------------

# --- Figure 1: density of both blocks on their nodes (2×2 grid, one panel per solver) ---
fig1, axs = plt.subplots(2, 2, figsize=(14, 8))
axs = axs.ravel()

# Density scaled by the pre-shock density of each solver; one colour scale for all panels
vmin = min(min(d['rho0'][:, :, 0].min(), d['rho1'][:, :, 0].min()) / d['rho_pre'] for d in fields_data.values())
vmax = max(max(d['rho0'][:, :, 0].max(), d['rho1'][:, :, 0].max()) / d['rho_pre'] for d in fields_data.values())

# Increase font sizes
title_fontsize = 20
label_fontsize = 19
cbar_fontsize = 17

for ax, solver in zip(axs, solvers):
    d = fields_data[solver]
    # Each block on its own nodes: block 2 is the parallelogram whose lower side is the ramp
    for xb, yb, rb in ((d['x0'], d['y0'], d['rho0']), (d['x1'], d['y1'], d['rho1'])):
        cf = ax.pcolormesh(xb[:, :, 0], yb[:, :, 0], rb[:, :, 0] / d['rho_pre'], shading='flat', cmap='coolwarm',
                           vmin=vmin, vmax=vmax, rasterized=True)
    ax.set_title(f'{solver}', fontsize=title_fontsize, fontweight='bold')
    ax.set_xlabel('x', fontsize=label_fontsize)
    ax.set_ylabel('y', fontsize=label_fontsize)
    ax.tick_params(labelsize=label_fontsize - 2)
    ax.set_aspect('equal')

plt.tight_layout()
cbar = fig1.colorbar(cf, ax=axs.tolist(), orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label(r'$\rho/\rho_1$', fontsize=cbar_fontsize)
cbar.ax.tick_params(labelsize=cbar_fontsize - 2)

if args.plot:
    plt.show()
else:
    out_path = os.path.join(output_dir, f"OS-fields.svg")
    plt.savefig(out_path, bbox_inches="tight", transparent=True, dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")

# --- Figure 2: ρ profile along y at i=30 of block 1 ---
i_probe = 30
fig2, ax2 = plt.subplots(figsize=(8, 5))

for solver in solvers:
    d    = fields_data[solver]
    yn1  = d['y1']                          # nodes (39, 219, 2)
    rho1 = d['rho1']                        # cells (38, 218, 1)
    # Exact post-shock density: the pre-shock density times the exact ratio
    r2ref = d['rho_pre'] * r2r1_ref
    # Cell-centre y at i=i_probe: average the four surrounding nodes
    yc = 0.25 * (yn1[i_probe,   :-1, 0] + yn1[i_probe,   1:, 0] + yn1[i_probe+1, :-1, 0] + yn1[i_probe+1, 1:, 0])
    rho_prof = rho1[i_probe, :, 0]
    ax2.plot(yc, rho_prof/r2ref, label=solver, linewidth=2)

ax2.set_xlim(0.09, 0.12)
ax2.set_ylim(0.85, 1.15)
ax2.set_xlabel('y', fontsize=11)
ax2.set_ylabel('$\\rho/\\rho_2$', fontsize=11)
ax2.legend(loc='best', framealpha=0.9)
ax2.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()

if args.plot:
    plt.show()
else:
    out_path = os.path.join(output_dir, f"OS-slice.svg")
    plt.savefig(out_path, bbox_inches="tight", transparent=True, dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")

# Exit status of the check: 1 when a solver is outside the threshold
if any_fail:
    sys.exit(1)
