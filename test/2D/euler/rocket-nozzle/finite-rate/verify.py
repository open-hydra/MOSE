import argparse
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import warnings
warnings.filterwarnings("ignore")

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

output_dir = "../../../../docs/vv/images/"

# Try to make the script work even when conda is not activated
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
        "ORION package not found. Run 'conda activate base' or install ORION in the current Python environment."
    )

# -------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------
parser = argparse.ArgumentParser(description='Validation script')
parser.add_argument('--plot', action='store_true', help='Plot pressure contours and ratio comparison')
args = parser.parse_args()


[x_, y_, z_, var_, vnames] = read_TEC("OUTPUT/field.tec")
fields_data = {
    'x': x_[0],
    'y': y_[0],
    'p': var_[0][16],
    'rho': sum(var_[0][0:12]),
    'u': var_[0][13],
    'v': var_[0][14],
    'g': var_[0][18],
}

[x_, y_, z_, var_, vnames] = read_TEC("reference/field.tec")
ref_fields_data = {
    'x': x_[0],
    'y': y_[0],
    'p': var_[0][16],
    'rho': sum(var_[0][0:12]),
    'u': var_[0][13],
    'v': var_[0][14],
    'g': var_[0][18],
}

# Check that the reference and simulation grids match
for key in ['x', 'y']:
    if not np.allclose(fields_data[key], ref_fields_data[key]):
        raise ValueError(f"Grid mismatch for {key}: simulation and reference grids do not match.")
    
# Check that the reference and simulation fields match at the throat (x ≈ 0.0 m)
throat_idx = np.argmin(np.abs(fields_data['x'][:, 0, 0]))
for key in ['p', 'rho', 'u', 'v', 'g']:
    if not np.allclose(fields_data[key][throat_idx, 0, 0], ref_fields_data[key][throat_idx, 0, 0], rtol=1e-3):
        raise ValueError(f"Field mismatch for {key} at throat: simulation and reference values do not match within tolerance.")

# -------------------------------------------------------------
# Plots
# -------------------------------------------------------------

d = fields_data
xn = d['x'][:, :, 0]; yn = d['y'][:, :, 0]
a = np.sqrt(d['g'][:, :, 0] * d['p'][:, :, 0] / d['rho'][:, :, 0])
M = np.sqrt(d['u'][:, :, 0]**2 + d['v'][:, :, 0]**2) / a

# ----------------------------------------------------------
# Quasi-1D comparison: area-weighted averages vs CEA
# ----------------------------------------------------------
# CEA reference: frozen, Merlin 1D, RP-1/LOX O/F=2.35
# Source: CEA.out (Pin=108 bar, O/F=2.35)
# Stations: sub Ae/At=2&1.5, throat, sup Ae/At=2,4,8,16
P0_cea = 108.00    # BAR (stagnation)
T0_cea = 3644.12   # K   (stagnation)
#                     sub2.0    throat   sup2.0   sup4.0    sup8.0   sup16.0
AeAt_cea = np.array([ 2.0,        1.0,     2.0,     4.0,      8.0,     16.0])
M_cea    = np.array([ 0.312,    1.000,   2.066,   2.650,    3.184,    3.715])
P_cea    = np.array([101.84,   60.606,  12.676,   4.4664,   1.6840,   0.65259])
T_cea    = np.array([3606.55, 3289.48, 2477.34,  2037.39,  1686.54, 1393.69])
pp0_cea  = P_cea / P0_cea
TT0_cea  = T_cea / T0_cea

# Cell-face radii (axisymmetric: y = radius).
# Nodes: (I, J) = (121, 41); cells: (120, 40).
r_lo = 0.5 * (yn[:120, :40] + yn[1:, :40])   # inner radial face, shape (120, 40)
r_hi = 0.5 * (yn[:120, 1:]  + yn[1:, 1:])    # outer radial face, shape (120, 40)
w_ring = r_hi**2 - r_lo**2                    # proportional to ring area

def area_avg_1d(q2d):
    """Cross-sectional area-weighted average; q2d shape (Ncx, Ncy) -> (Ncx,)."""
    return np.sum(q2d * w_ring, axis=1) / np.sum(w_ring, axis=1)

# Cell-centre x coordinate
x_c = 0.5 * (xn[:120, 0] + xn[1:, 0])   # shape (120,)

# 2D field slices (k=0)
rho2d = sum(var_[0][k][:, :, 0] for k in range(13))  # all 13 species
p2d   = var_[0][16][:, :, 0]
T2d   = var_[0][17][:, :, 0]
g2d   = var_[0][18][:, :, 0]
u2d   = var_[0][13][:, :, 0]
v2d   = var_[0][14][:, :, 0]
a2d   = np.sqrt(g2d * p2d / rho2d)
M2d   = np.sqrt(u2d**2 + v2d**2) / a2d

# Cell-centre coordinates for Mach contour
x_cent = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
y_cent = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])

# Mach field plot
fig1, ax1 = plt.subplots(figsize=(10, 6))
contour = ax1.contourf(x_cent, y_cent, M2d, levels=50, cmap='turbo')
cax = inset_axes(ax1, width="40%", height="5%", loc='upper center')
cbar = plt.colorbar(contour, cax=cax, orientation='horizontal')
cbar.set_label("Mach number")
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_aspect('equal', adjustable='box')
plt.tight_layout()

# 1D averages and axis values
p_avg  = area_avg_1d(p2d)
T_avg  = area_avg_1d(T2d)
M_avg  = area_avg_1d(M2d)
p_axis = p2d[:, 0]
T_axis = T2d[:, 0]
M_axis = M2d[:, 0]

# Inlet stagnation recovered via isentropic relations
# (first cell is at Ae/At≈4.3, M≈0.14 — not at rest)
g0 = area_avg_1d(g2d)[0]
M0 = M_avg[0]
p0_sim = p_avg[0] * (1 + (g0 - 1) / 2 * M0**2)**(g0 / (g0 - 1))
T0_sim = T_avg[0] * (1 + (g0 - 1) / 2 * M0**2)

# CEA x-positions: locate each Ae/At station in the mesh geometry.
# Subsonic stations live in the converging section, supersonic in the diverging.
R_wall      = r_hi[:, -1]
i_throat    = int(np.argmin(R_wall))
R_t         = R_wall[i_throat]
AeAt_1d     = R_wall**2 / R_t**2            # local area ratio at each axial cell
sub_idx     = np.where(x_c < x_c[i_throat])[0]
sup_idx     = np.where(x_c >= x_c[i_throat])[0]

x_cea_list = []
for k, AeAt in enumerate(AeAt_cea):
    if M_cea[k] < 1.0:                      # subsonic → converging section
        region = sub_idx
    else:                                    # sonic/supersonic → diverging section
        region = sup_idx
    i = region[np.argmin(np.abs(AeAt_1d[region] - AeAt))]
    x_cea_list.append(x_c[i])
x_cea = np.array(x_cea_list)

fig3, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

# Mach
ax = axes[0]
ax.plot(x_c, M_avg,  label='MOSE 1D avg', linewidth=2)
ax.scatter(x_cea, M_cea, color='red', zorder=5, label='CEA', s=60)
ax.set_ylabel('Mach')
ax.legend()
ax.grid(True, alpha=0.3)

# p / p0
ax = axes[1]
ax.plot(x_c, p_avg  / p0_sim, label='MOSE 1D avg', linewidth=2)
ax.scatter(x_cea, pp0_cea, color='red', zorder=5, label='CEA', s=60)
ax.set_ylabel('p / p$_0$')
ax.legend()
ax.grid(True, alpha=0.3)

# T / T0
ax = axes[2]
ax.plot(x_c, T_avg  / T0_sim, label='MOSE 1D avg', linewidth=2)
ax.scatter(x_cea, TT0_cea, color='red', zorder=5, label='CEA', s=60)
ax.set_ylabel('T / T$_0$')
ax.set_xlabel('x  [m]')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()

if args.plot:
    plt.show()
else:
    os.makedirs(output_dir, exist_ok=True)

    out_path = os.path.join(output_dir, "Merlin-FZ-field.svg")
    fig1.savefig(out_path, bbox_inches="tight", transparent=True, dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")

    out_path = os.path.join(output_dir, "Merlin-FZ.svg")
    plt.savefig(out_path, bbox_inches="tight", transparent=True, dpi=100)
    fsize = os.path.getsize(out_path) / 1024
    print(f"  Figure saved to {out_path} ({fsize:.1f} KB)")
