#!/usr/bin/env python3
"""
Build the initial conditions for the 2D low-Mach Gresho vortex.

Reads a Tecplot mesh (nodal x,y,z) with the ORION library, evaluates the
Gresho vortex primitive variables (rho, u, v, w, p) at the cell centers and
writes them to ic.tec (cell-centered).

Mach-scaled rotating Gresho vortex (Miczek, Roepke & Edelmann 2015). The
azimuthal velocity profile is fixed (peak |u| = 1); the target maximum Mach
number M sets the background pressure, so the same vortex can be run at any
Mach number by changing M only:

    rho = 1   everywhere

    u_phi(r) = 5 r                     0.0 <= r < 0.2
             = 2 - 5 r                 0.2 <= r < 0.4
             = 0                       r   >= 0.4

    p(r) = p0 + dp(r),   p0 = rho * u_max^2 / (gamma * M^2) = 1 / (gamma M^2)

    dp(r) = 12.5 r^2 + 2 + 4 ln(1/2)                  0.0 <= r < 0.2
          = 12.5 r^2 - 20 r + 6 + 4 ln(r/0.4)         0.2 <= r < 0.4
          = 0                                         r   >= 0.4

dp is the centrifugal-balance increment (dp/dr = rho u_phi^2 / r), zero in the
far field, so the relative pressure fluctuation dp/p0 ~ M^2.

The mesh spans [0,0.8] x [0,0.8], so the vortex is centered at (0.4, 0.4).
Usage:  python build_ic.py [Mach]      (default Mach = 0.35, the classic case)
"""

import sys
import numpy as np

sys.path.insert(0, '/Users/marcogrossi/Codici/MOSE/lib/ORION/src/python')
from ORION import ORION

# --- parameters --------------------------------------------------------------
MESH_FILE = 'mesh.tec'
IC_FILE   = 'INPUT/ic.tec'
xc0, yc0  = 0.4, 0.4          # vortex center
rho0      = 1.0               # constant density
gamma     = 1.4              # ratio of specific heats
mach      = 0.35             # target maximum Mach number (override via CLI arg)


def gresho_uphi(r):
    uphi = np.where(r < 0.2, 5.0 * r,
            np.where(r < 0.4, 2.0 - 5.0 * r, 0.0))
    return uphi


def gresho_pressure(r, p0):
    # Centrifugal-balance increment dp(r), zero for r >= 0.4
    lr = np.log(np.maximum(r, 1e-300))
    dp = np.where(r < 0.2,
                  12.5 * r**2 + 2.0 + 4.0 * np.log(0.5),
             np.where(r < 0.4,
                  12.5 * r**2 - 20.0 * r + 6.0 + 4.0 * (lr - np.log(0.4)),
                  0.0))
    return p0 + dp


def main():
    M = float(sys.argv[1]) if len(sys.argv) > 1 else mach
    p0 = rho0 * 1.0**2 / (gamma * M**2)   # u_max = 1
    print(f'Mach = {M},  background pressure p0 = {p0}')

    # Read the nodal mesh
    xb, yb, zb, vb_in, names = ORION.read_TEC(MESH_FILE)

    xb_out, yb_out, zb_out, vb_out = [], [], [], []

    for b in range(len(xb)):
        xn, yn, zn = xb[b], yb[b], zb[b]
        Ni, Nj, Nk = xn.shape          # node counts
        Nx, Ny, Nz = Ni - 1, Nj - 1, Nk - 1   # cell counts

        # Cell-center coordinates: average of the 8 surrounding nodes
        def cell_center(a):
            return 0.125 * (a[:-1, :-1, :-1] + a[1:, :-1, :-1] +
                            a[:-1, 1:, :-1] + a[:-1, :-1, 1:] +
                            a[1:, 1:, :-1] + a[1:, :-1, 1:] +
                            a[:-1, 1:, 1:] + a[1:, 1:, 1:])

        xc = cell_center(xn)
        yc = cell_center(yn)

        # Gresho vortex evaluated at cell centers
        dx = xc - xc0
        dy = yc - yc0
        r  = np.sqrt(dx**2 + dy**2)

        uphi = gresho_uphi(r)
        p    = gresho_pressure(r, p0)

        # Avoid division by zero at the center (velocity -> 0 there anyway)
        rsafe = np.where(r > 0.0, r, 1.0)
        u = -uphi * dy / rsafe
        v =  uphi * dx / rsafe
        w = np.zeros_like(u)
        rho = np.full_like(u, rho0)

        vb_out.append([rho, u, v, w, p])
        xb_out.append(xn)
        yb_out.append(yn)
        zb_out.append(zn)

    var_names = ['x', 'y', 'z', 'rho1', 'u', 'v', 'w', 'p']
    ORION.write_TEC(IC_FILE, xb_out, yb_out, zb_out, vb_out, var_names)
    print(f'Wrote {IC_FILE}')


if __name__ == '__main__':
    main()
