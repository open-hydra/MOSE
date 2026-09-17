#!/usr/bin/env python3
"""
Build the initial condition for the 2D viscous Method-of-Manufactured-Solutions
verification (test N2).

Reads the GRIB-generated Tecplot mesh (nodal x,y,z) with the ORION library,
evaluates the MANUFACTURED primitive field (rho, u, v, w, p) at the cell centers
and writes them to ic.tec (self-contained: nodal mesh + cell data, the format
MOSE reads).

The manufactured field is a smooth, periodic solution on [0,1]^2 (om = 2*pi):

    rho = 1.0   + 0.1 *sin(om x) cos(om y)
    u   = 40.0  + 8.0 *sin(om x) sin(om y)
    v   = 30.0  + 8.0 *cos(om x) cos(om y)
    p   = 8.0e3 + 8.0e2*cos(om x) sin(om y)

It is made an EXACT *steady* solution of the compressible Navier-Stokes
equations by the analytic source term injected by the MOSE_mms driver
(src/app/mms.f90, generated + FD-checked by mms_gen.py).  These amplitudes MUST
match src/app/mms.f90, mms_gen.py and verify.py.  With mu = 10 the Reynolds
number is ~5, so the viscous terms are ~20 % of the convective ones and the
viscous operator (gradients, metrics) is genuinely exercised -- which the
inviscid isentropic vortex (N1) cannot do.  The low pressure level keeps the
sound speed low so the mean Mach number is ~0.47: high enough that the
compressible scheme does not suffer low-Mach accuracy degradation (which would
otherwise leave a grid-independent error in rho/u/v).

Usage:  python build_ic.py            # uses the mesh produced by 'GRIB meshgen'
"""

import os
import sys
import numpy as np

# --- portable ORION import (search upward for lib/ORION/src/python) ----------
from pathlib import Path
_here = Path(__file__).resolve().parent
for _p in [_here, *_here.parents]:
    _cand = _p / "lib" / "ORION" / "src" / "python"
    if _cand.exists():
        sys.path.insert(0, str(_cand))
        break
from ORION import ORION

# --- manufactured-field parameters (must match src/app/mms.f90 & verify.py) --
OM = 2.0 * np.pi
RHO0, ARHO = 1.0,   0.1
U0,   AU   = 40.0,  8.0
V0,   AV   = 30.0,  8.0
P0,   AP   = 8.0e3, 8.0e2

IC_FILE = "ic.tec"


def _find_mesh():
    for c in ("mesh.tec", os.path.join("MESH", "mesh.tec")):
        if os.path.exists(c):
            return c
    sys.exit("build_ic.py: mesh not found (run 'GRIB meshgen' first)")


def manufactured(x, y):
    sx, cx = np.sin(OM * x), np.cos(OM * x)
    sy, cy = np.sin(OM * y), np.cos(OM * y)
    rho = RHO0 + ARHO * sx * cy
    u   = U0   + AU   * sx * sy
    v   = V0   + AV   * cx * cy
    p   = P0   + AP   * cx * sy
    w   = np.zeros_like(u)
    return rho, u, v, w, p


def main():
    mesh = _find_mesh()
    xb, yb, zb, _vb, _names = ORION.read_TEC(mesh)

    xb_out, yb_out, zb_out, vb_out = [], [], [], []
    for b in range(len(xb)):
        xn, yn, zn = xb[b], yb[b], zb[b]

        def cell_center(a):
            return 0.125 * (a[:-1, :-1, :-1] + a[1:, :-1, :-1] +
                            a[:-1, 1:, :-1] + a[:-1, :-1, 1:] +
                            a[1:, 1:, :-1] + a[1:, :-1, 1:] +
                            a[:-1, 1:, 1:] + a[1:, 1:, 1:])

        xc = cell_center(xn)
        yc = cell_center(yn)
        rho, u, v, w, p = manufactured(xc, yc)

        vb_out.append([rho, u, v, w, p])
        xb_out.append(xn)
        yb_out.append(yn)
        zb_out.append(zn)

    var_names = ["x", "y", "z", "rho1", "u", "v", "w", "p"]
    ORION.write_TEC(IC_FILE, xb_out, yb_out, zb_out, vb_out, var_names)
    print(f"Wrote {IC_FILE} (manufactured viscous field, mu=10)")


if __name__ == "__main__":
    main()
