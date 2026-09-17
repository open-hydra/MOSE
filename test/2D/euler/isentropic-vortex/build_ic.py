#!/usr/bin/env python3
"""
Build the initial condition for the 2D isentropic Euler vortex (N1).

Reads the GRIB-generated Tecplot mesh (nodal x,y,z) with the ORION library,
evaluates the isentropic-vortex primitive variables (rho, u, v, w, p) at the
cell centers and writes them to ic.tec (self-contained: nodal mesh + cell
data, the format MOSE reads).

Isentropic vortex (Shu 1998 / Hu-Shu), a *smooth* (C-infinity) exact solution
of the compressible Euler equations superimposed on a uniform mean flow.  With
mean flow (u_inf, v_inf) the whole vortex is convected without change of shape;
after one domain pass it returns to its start, so the exact solution at that
time equals the initial condition.  Because it is smooth, a 2nd-order scheme
shows a clean 2nd-order error under grid refinement (that is the point of N1,
and the reason the kinked Gresho vortex cannot be used here).

    r^2 = (x-x0)^2 + (y-y0)^2
    du  = -(beta/2pi) (y-y0) exp((1-r^2)/2)
    dv  =  (beta/2pi) (x-x0) exp((1-r^2)/2)
    dT  = -(gamma-1) beta^2 / (8 gamma pi^2) exp(1-r^2)
    T   = 1 + dT,   rho = T^(1/(gamma-1)),   p = rho^gamma   (entropy = 1)
    u   = u_inf + du,   v = v_inf + dv,   w = 0

Domain [0,L]^2 with the vortex centered at (L/2, L/2); mean flow (1,1) so one
period is t = L / u_inf.

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

# --- parameters (must match verify.py) --------------------------------------
L      = 10.0          # square domain [0, L]^2
x0, y0 = 5.0, 5.0      # vortex center
gamma  = 1.4
beta   = 5.0           # vortex strength
uinf   = 1.0           # mean-flow velocity (u = v)
vinf   = 1.0

IC_FILE = "ic.tec"


def _find_mesh():
    for c in ("mesh.tec", "MESH/mesh.tec", os.path.join("MESH", "mesh.tec")):
        if os.path.exists(c):
            return c
    sys.exit("build_ic.py: mesh not found (run 'GRIB meshgen' first)")


def vortex(x, y):
    dx, dy = x - x0, y - y0
    r2 = dx * dx + dy * dy
    f = np.exp(0.5 * (1.0 - r2))
    du = -(beta / (2.0 * np.pi)) * dy * f
    dv = (beta / (2.0 * np.pi)) * dx * f
    dT = -((gamma - 1.0) * beta * beta) / (8.0 * gamma * np.pi * np.pi) * np.exp(1.0 - r2)
    T = 1.0 + dT
    rho = T ** (1.0 / (gamma - 1.0))
    p = rho ** gamma
    u = uinf + du
    v = vinf + dv
    w = np.zeros_like(u)
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
        rho, u, v, w, p = vortex(xc, yc)

        vb_out.append([rho, u, v, w, p])
        xb_out.append(xn)
        yb_out.append(yn)
        zb_out.append(zn)

    var_names = ["x", "y", "z", "rho1", "u", "v", "w", "p"]
    ORION.write_TEC(IC_FILE, xb_out, yb_out, zb_out, vb_out, var_names)
    print(f"Wrote {IC_FILE} (isentropic vortex, beta={beta}, mean flow=({uinf},{vinf}))")


if __name__ == "__main__":
    main()
