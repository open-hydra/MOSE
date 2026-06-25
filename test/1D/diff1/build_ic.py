#!/usr/bin/env python3
"""
Build the initial conditions for the 1D multicomponent diffusion test.

Reproduces section 4.2 ("Multicomponent diffusion test") of Forti et al.,
"Development and Validation of a Compressible Multicomponent Reactive Flow
Solver based on STREAmS-2.0" (EUCASS 2025).

Reads a Tecplot mesh (nodal x,y,z) with the ORION library, evaluates the
primitive variables (partial densities rho_k, u, v, w, p) at the cell centers
and writes them to ic.tec (cell-centered).

Setup (eq. 8-10 + Tab. 1):

    velocity   u = v = w = 0
    pressure   p = 101325 Pa (uniform)

    Y_k(x) = Y_{k,0} + (Y_{k,f} - Y_{k,0}) f(x),   k = O2, N2, H2O, CH4
    T(x)   = T_0     + (T_f     - T_0    ) f(x)
    f(x)   = 1 - exp( -(x - x0)^2 / d^2 )

with d = 2.5e-3 m and x0 = 25e-3 m. Subscript 0 is the oxidizer side
(f(x0) = 0 at the domain centre), subscript f the fuel side (f -> 1 far field).

Density follows from the ideal-gas EoS:  rho = p * W_mix / (Ru * T),
1/W_mix = sum_k Y_k / W_k.  Partial densities: rho_k = Y_k * rho.

The species order written to ic.tec MUST match INPUT/phase.txt:
    rho1 = CH4, rho2 = H2O, rho3 = N2, rho4 = O2.

Usage:  python build_ic.py
"""

import sys
import numpy as np

sys.path.insert(0, '/data10/grossi/tmp/MOSE/lib/ORION/src/python')
from ORION import ORION

# --- parameters --------------------------------------------------------------
MESH_FILE = 'mesh.tec'
IC_FILE   = 'INPUT/ic.tec'

x0   = 25.0e-3               # profile centre [m]
d    = 2.5e-3               # profile width  [m]
p0   = 101325.0            # uniform pressure [Pa]
Ru   = 8314.462618        # universal gas constant [J/(kmol K)]

# Species order = INPUT/phase.txt order: CH4, H2O, N2, O2
SPECIES = ['CH4', 'H2O', 'N2', 'O2']
W = {'CH4': 16.043, 'H2O': 18.015, 'N2': 28.014, 'O2': 31.998}   # [kg/kmol]

# Boundary states (Tab. 1).  '0' = oxidizer side (centre), 'f' = fuel side.
T_0, T_f = 1350.0, 320.0
Y_0 = {'CH4': 0.0,   'H2O': 0.1, 'N2': 0.758, 'O2': 0.142}   # oxidizer side
Y_f = {'CH4': 0.214, 'H2O': 0.0, 'N2': 0.591, 'O2': 0.195}   # fuel side


def profile(x):
    return 1.0 - np.exp(-(x - x0) ** 2 / d ** 2)


def main():
    # Read the nodal mesh
    xb, yb, zb, vb_in, names = ORION.read_TEC(MESH_FILE)

    xb_out, yb_out, zb_out, vb_out = [], [], [], []

    for b in range(len(xb)):
        xn, yn, zn = xb[b], yb[b], zb[b]

        # Cell-center coordinates: average of the 8 surrounding nodes
        def cell_center(a):
            return 0.125 * (a[:-1, :-1, :-1] + a[1:, :-1, :-1] +
                            a[:-1, 1:, :-1] + a[:-1, :-1, 1:] +
                            a[1:, 1:, :-1] + a[1:, :-1, 1:] +
                            a[:-1, 1:, 1:] + a[1:, 1:, 1:])

        xc = cell_center(xn)

        f = profile(xc)
        T = T_0 + (T_f - T_0) * f

        # Mass fractions and mean molar mass
        Yk = {s: Y_0[s] + (Y_f[s] - Y_0[s]) * f for s in SPECIES}
        inv_W = sum(Yk[s] / W[s] for s in SPECIES)
        W_mix = 1.0 / inv_W

        # Ideal-gas density and partial densities
        rho = p0 * W_mix / (Ru * T)
        rho_k = [Yk[s] * rho for s in SPECIES]

        u = np.zeros_like(xc)
        v = np.zeros_like(xc)
        w = np.zeros_like(xc)
        p = np.full_like(xc, p0)

        vb_out.append(rho_k + [u, v, w, p])
        xb_out.append(xn)
        yb_out.append(yn)
        zb_out.append(zn)

    var_names = (['x', 'y', 'z'] +
                 [f'rho{i + 1}' for i in range(len(SPECIES))] +
                 ['u', 'v', 'w', 'p'])
    ORION.write_TEC(IC_FILE, xb_out, yb_out, zb_out, vb_out, var_names)
    print(f'Wrote {IC_FILE}  (species order: {", ".join(SPECIES)})')


if __name__ == '__main__':
    main()
