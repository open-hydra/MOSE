"""
Figure of the four grid families of the planar-nozzle study on the coarse
level (level 0): O, H, ZZ and SB2, as make_case.py builds them, each of the
four blocks in its own colour.  No run is needed: the grids come from
make_case.grid().

Usage:
    python plot_grids.py            # writes OUTPUT/planar-nozzle-grids.svg
    python plot_grids.py --plot     # also shows the figure
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_case import grid, CELLS_F, LEVELS  # noqa: E402

parser = argparse.ArgumentParser(description="grid families of the planar nozzle, coarse level")
parser.add_argument("--plot", action="store_true", help="show the figure")
args = parser.parse_args()

import matplotlib as mpl  # noqa: E402
if not args.plot:
    mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

FAMILIES = ["O", "H", "ZZ", "SB2"]
nx = sum(LEVELS[0])
fig, axes = plt.subplots(2, 2, figsize=(15, 6.5), sharex=True, sharey=True)
for k, (ax, fam) in enumerate(zip(axes.ravel(), FAMILIES)):
    for b, N in enumerate(grid(fam, 0)):
        X, Y = N[:, :, 0, 0], N[:, :, 0, 1]
        col = "C%d" % b
        ax.plot(X, Y, "-", color=col, lw=0.5)
        ax.plot(X.T, Y.T, "-", color=col, lw=0.5)
    ax.set_aspect("equal")
    if k >= 2:
        ax.set_xlabel("x [m]", fontsize=20)
    if k % 2 == 0:
        ax.set_ylabel("y [m]", fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_title("%s, %d x %d cells" % (fam, nx, sum(CELLS_F[fam])), fontsize=18)
fig.legend([plt.Line2D([], [], color="C%d" % b) for b in range(4)], ["block %d" % (b + 1) for b in range(4)],
           loc="upper center", ncol=4, fontsize=18, bbox_to_anchor=(0.5, 1.04))
plt.tight_layout(rect=[0, 0, 1, 0.92])
os.makedirs("OUTPUT", exist_ok=True)
plt.savefig(os.path.join("OUTPUT", "planar-nozzle-grids.svg"), bbox_inches="tight", transparent=True)
if args.plot:
    plt.show()
