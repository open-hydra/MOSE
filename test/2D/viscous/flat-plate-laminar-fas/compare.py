#!/usr/bin/env python3
"""
Skin friction of the FAS cycle against that of the grid sequencing on the laminar flat plate.

Reads the wall shear stress tauX of RUN/fas/OUTPUT/wall.tec and RUN/seq/OUTPUT/wall.tec (the runs of
./MOSE.sh test), as verify.py of ../flat-plate-laminar reads it, and prints the RMS and the largest relative
difference of the two over the plate for x >= 0.05 m, the region that verify.py checks against Blasius.

Usage (from this folder, after ./MOSE.sh test):
    python3 compare.py
"""
import sys
from pathlib import Path

import numpy as np

# ORION Python package of the source tree
root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break
from ORION import read_TEC  # noqa: E402


def wall_tau(path):
    """Cell-centred x and |tauX| along the wall, for x >= 0.05 m."""
    xw, _, _, varw, _ = read_TEC(str(path))
    x, tau = xw[0], varw[0][1]
    # the index along which x varies, as verify.py picks it
    if x.shape[0] > 1 and (x.shape[1] == 1 or np.nanstd(x[:, 0, 0]) >= np.nanstd(x[0, :, 0])):
        xn, t = x[:, 0, 0], tau[:, 0, 0]
    else:
        xn, t = x[0, :, 0], tau[0, :, 0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    mask = xc >= 0.05
    return xc[mask], np.abs(t[mask])


x_fas, t_fas = wall_tau(root / "RUN" / "fas" / "OUTPUT" / "wall.tec")
x_seq, t_seq = wall_tau(root / "RUN" / "seq" / "OUTPUT" / "wall.tec")
if x_fas.shape != x_seq.shape or not np.allclose(x_fas, x_seq):
    sys.exit("compare.py: the two wall files have different cells")

rel = np.abs(t_fas - t_seq) / t_seq * 100.0
print(f"Skin friction, cycle against grid sequencing (x >= 0.05 m): RMS = {np.sqrt(np.mean(rel**2)):.2f} %, "
      f"max = {rel.max():.2f} %")
