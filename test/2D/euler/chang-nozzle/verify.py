import argparse
import math
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

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
    'p': var_[0][4],
    'rho': var_[0][0],
    'u': var_[0][1],
    'v': var_[0][2],
    'g': var_[0][6],
}

# -------------------------------------------------------------
# Optional plots
# -------------------------------------------------------------
if args.plot:

    d = fields_data
    xn = d['x'][:, :, 0]; yn = d['y'][:, :, 0]
    a = np.sqrt(d['g'][:, :, 0] * d['p'][:, :, 0] / d['rho'][:, :, 0])
    M = np.sqrt(d['u'][:, :, 0]**2 + d['v'][:, :, 0]**2) / a

    fig2, ax2 = plt.subplots(figsize=(10, 5))
    # align grid coordinates with M shape
    x_axis = d['x'][: M.shape[0], 0, 0]
    x_wall = d['x'][: M.shape[0], -1, 0]
    ax2.plot(x_axis, M[:, 0], label='MOSE axis', linewidth=2)
    ax2.plot(x_wall, M[:, -1], label='MOSE wall', linewidth=2)

    # overlay OUTREF data
    ref1 = np.loadtxt('reference/ref1.dat')
    ref2 = np.loadtxt('reference/ref2.dat')
    ax2.plot(ref1[:, 0], ref1[:, 1], 'o', markersize=4, alpha=0.7, label='Ref. axis')
    ax2.plot(ref2[:, 0], ref2[:, 1], 's', markersize=4, alpha=0.7, label='Ref. wall')

    ax2.set_xlabel('x')
    ax2.set_ylabel('Mach')
    ax2.legend()
    ax2.grid(True)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
