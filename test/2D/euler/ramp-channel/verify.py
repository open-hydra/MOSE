import argparse
import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ORION import read_TEC
import warnings
warnings.filterwarnings("ignore")


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
M1    = 3.0
gamma = 1.4
tol   = 1e-6

# Measure ramp angle directly from mesh block-0 wall nodes (j=0)
_tmp = read_TEC('OUTPUT/field-HLLC.tec')
_xw  = _tmp[0][0][:, 0, 0]    # x of wall nodes along i-axis (block 0)
_yw  = _tmp[1][0][:, 0, 0]    # y of wall nodes
delta_deg = float(np.degrees(np.arctan2(_yw[50] - _yw[20], _xw[50] - _xw[20])))
del _tmp, _xw, _yw

delta = math.radians(delta_deg)

# Analytical solution
mach2_ref, beta, p2p1_ref, t2t1_ref, r2r1_ref = sckobl(M1, delta, gamma, tol)

# print("Oblique Shock – Analytical Solution")
# print(f"  M1    = {M1:.4f}   delta = {delta_deg:.1f} deg   gamma = {gamma}")
# print(f"  beta  = {math.degrees(beta):.4f} deg")
# print(f"  M2    = {mach2_ref:.6f}")
# print(f"  p2/p1 = {p2p1_ref:.6f}")
# print(f"  T2/T1 = {t2t1_ref:.6f}")
# print(f"  r2/r1 = {r2r1_ref:.6f}")
# print()

# -------------------------------------------------------------
# Command-line arguments
# -------------------------------------------------------------
parser = argparse.ArgumentParser(description='Oblique shock verification')
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
# Load, compare, check
# -------------------------------------------------------------
fields_data = {}   # store read data per solver for optional plotting

for solver, fpath in zip(solvers, files):
    [x_, y_, z_, var_, vnames] = read_TEC(fpath)
    fields_data[solver] = {'x0': x_[0], 'y0': y_[0], 'p0': var_[0][4], 'rho0': var_[0][0]}

# -------------------------------------------------------------
# Optional plots
# -------------------------------------------------------------
if args.plot:

    # --- Figure 1: density contours ---
    # Calculate aspect ratio from first solver to size the figure properly
    d_first = fields_data[solvers[0]]
    xn0_first = d_first['x0'][:, :, 0]; yn0_first = d_first['y0'][:, :, 0]
    x_range = xn0_first.max() - xn0_first.min()
    y_range = yn0_first.max() - yn0_first.min()
    aspect_ratio = x_range / y_range
    
    # Create figure with size that preserves aspect ratio (2x2 grid)
    # Each subplot scaled by aspect ratio
    subplot_width = 7
    subplot_height = subplot_width / aspect_ratio
    fig1, axs = plt.subplots(2, 2, figsize=(subplot_width * 2, subplot_height * 2))
    axs = axs.ravel()

    for ax, solver in zip(axs, solvers):
        d   = fields_data[solver]
        # Block 0
        xn0 = d['x0'][:, :, 0]; yn0 = d['y0'][:, :, 0]; vc0 = d['rho0'][:, :, 0]
        # Non-dimensional density
        vc0 = vc0 / d['rho0'][0, 0, 0]
        vmin = vc0.min()
        vmax = vc0.max()
        cf = ax.pcolormesh(xn0.T, yn0.T, vc0.T, shading='flat', cmap='turbo', vmin=vmin, vmax=vmax)
        fig1.colorbar(cf, ax=ax, label='ρ')
        ax.set_title(f'{solver}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')

    fig1.suptitle('Density Field')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4)
    plt.show()
