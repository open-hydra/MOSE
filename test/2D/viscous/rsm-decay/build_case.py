#!/usr/bin/env python3
"""Generate the deck, initial condition and BC file for the Reynolds-stress decay case.

Geometry: an axisymmetric duct sector (a wedge of angle THETA about the x axis),
between the radii R_IN and R_OUT, far from any wall.

      face 4  (j = ny)   symmetry (slip)
    +---------------------------------+
    |                                 |
    |  inflow  ->  uniform stream  -> |  outflow     faces 5, 6 : axisymmetric
    |                                 |
    +---------------------------------+
      face 3  (j = 1)    symmetry (slip)

A uniform stream carries isotropic turbulence (R11 = R22 = R33, no shear stress)
in from the left. With no mean strain nothing produces it, and it decays along x
by the closed form that verify.py checks (see the module docstring there).

Usage:  python3 build_case.py <model>
"""

import math
import os
import sys

# --- models ------------------------------------------------------------------
# The decay is the same for all of them: the diffusion (generalized gradient or
# simple) is negligible at this Peclet number, and the length-scale correction of
# w2019 acts near walls only.
MODELS = ('SSGLRR', 'SSGLRR-SD', 'SSGLRR-w2019', 'SSGLRR-w2019-SD')

# --- discretisation ----------------------------------------------------------
NX, NY = 32, 4
LX          = 1.0                   # [m]
R_IN, R_OUT = 0.5, 1.5              # [m]
THETA       = math.radians(1.0)     # wedge angle

# --- stream ------------------------------------------------------------------
P0, T0, MACH = 1.0e5, 300.0, 0.2
GAMMA  = 1.4
R_AIR  = 8314.46261815324 / 28.970418    # gas constant of the single-species gas (FLINT Runiv)
RII    = 10.0       # normal Reynolds stresses at the inlet [m2/s2]: k = 1.5 RII, 0.3% of U^2
OMEGA  = 600.0      # specific dissipation rate at the inlet [1/s]
SHEAR  = 1.0e-8     # shear stresses at the inlet [m2/s2] (the value ATLAS writes)

# Boundary condition type codes (see src/lib/io/IO_BC.f90)
BC_AXISYMMETRIC, BC_SYMMETRY, BC_INLET_G_T, BC_OUTLET = 200, 300, 404, 406

DECK = """[MOSE-Parameters]
iter-threshold = 6000
newrun         = true

[MOSE-Physics]
equations  = navier-stokes
turbulence = {model}

[MOSE-Numerics]
cfl                   = 0.8
vnn                   = 0.3
time-accurate         = false
time-scheme           = RK3
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = SLAU

[MOSE-IO]
sol-diter      = 1000000000
shell-diter    = 2000
ini-diter      = 1000000000
"""


def stream():
    rho = P0 / (R_AIR * T0)
    u = MACH * math.sqrt(GAMMA * R_AIR * T0)
    return rho, u


def write_ic(path):
    """Uniform stream carrying the inlet turbulence."""
    ni, nj, nk = NX + 1, NY + 1, 2
    ncell = NX * NY
    rho, u = stream()

    x, y, z = [], [], []
    for k in range(nk):
        th = (k - 0.5) * THETA
        for j in range(nj):
            r = R_IN + (R_OUT - R_IN) * j / NY
            for i in range(ni):
                x.append(LX * i / NX); y.append(r * math.cos(th)); z.append(r * math.sin(th))

    # the turbulence bands of an initial condition are density-weighted
    turb = [rho * RII] * 3 + [rho * SHEAR] * 3 + [rho * OMEGA]
    with open(path, 'w') as f:
        f.write(" VARIABLES = x y z  rho1 u v w p ru'u' rv'v' rw'w' ru'v' ru'w' rv'w' rhoOmega\n")
        f.write(f' ZONE  T = Block1, I={ni}, J={nj}, K={nk}, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-3]=NODAL,[4-15]=CELLCENTERED)\n')
        for arr in (x, y, z):
            for v in arr:
                f.write(f' {v: .15E}\n')
        for value in [rho, u, 0.0, 0.0, P0] + turb:
            for _ in range(ncell):
                f.write(f' {value: .15E}\n')


def write_bc(path):
    """One record per boundary cell face: 'block i j k face type', plus a
    parameter line for the types that carry data."""
    rho, u = stream()
    per_face = {f: [] for f in range(1, 7)}

    def rec(face, i, j, k, bctype, params=None):
        per_face[face].append(f'{1:8d}{i:8d}{j:8d}{k:8d}{face:8d}{bctype:8d}')
        if params is not None:
            per_face[face].append(params)

    # face 1: inlet, T + mass flux; the turbulence values of an inlet are specific
    # T, g, alpha, beta, rel_fac, mass fractions(1:nsc), R11 R22 R33 R12 R13 R23 omega
    values = [T0, rho * u, None, None, 1.0, 1.0] + [RII] * 3 + [SHEAR] * 3 + [OMEGA]
    inlet = ''.join(f'{"normal":>15s},' if v is None else f'  {v:.16E},' for v in values)
    for j in range(1, NY + 1):
        rec(1, 1, j, 1, BC_INLET_G_T, inlet)

    # face 2: outlet at the back pressure (p_back, rel_fac)
    for j in range(1, NY + 1):
        rec(2, NX, j, 1, BC_OUTLET, f'  {P0:.16E},  {1.0:.16E},')

    # faces 3 / 4: slip walls
    for i in range(1, NX + 1):
        rec(3, i, 1, 1, BC_SYMMETRY)
    for i in range(1, NX + 1):
        rec(4, i, NY, 1, BC_SYMMETRY)

    # faces 5 / 6: the two sides of the wedge
    for f in (5, 6):
        for j in range(1, NY + 1):
            for i in range(1, NX + 1):
                rec(f, i, j, 1, BC_AXISYMMETRIC)

    lines = [ln for f in range(1, 7) for ln in per_face[f]]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in MODELS:
        sys.exit(f'usage: build_case.py [{"|".join(MODELS)}]')
    model = sys.argv[1]

    os.makedirs('INPUT', exist_ok=True)
    os.makedirs('OUTPUT', exist_ok=True)

    with open('input.ini', 'w') as f:
        f.write(DECK.format(model=model))
    write_ic('INPUT/ic.tec')
    write_bc('INPUT/bc.txt')

    # Record the *configuration* only: verify.py derives the reference solution itself.
    with open('INPUT/case.txt', 'w') as f:
        f.write(f'model {model}\n')
        f.write(f'Rii {RII!r}\n')
        f.write(f'omega {OMEGA!r}\n')

    print(f'model = {model}')
    print(f'inlet: Rii = {RII:g} m2/s2, omega = {OMEGA:g} 1/s, U = {stream()[1]:.3f} m/s')


if __name__ == '__main__':
    main()
