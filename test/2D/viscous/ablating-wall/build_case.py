#!/usr/bin/env python3
"""Generate the mesh, initial condition and BC file for the ablating-wall case.

Geometry: a thin vertical column of gas, one-dimensional in y.

      face 4  (j = ny)   outflow
    +-------------------+
    |                   |
    |  gas column       |   faces 1,2 : symmetry -> the solution is 1-D in y
    |                   |
    +-------------------+
      face 3  (j = 1)    ablating wall (BC 503 / 504)

The ablating wall blows gas into the column; it leaves through the outflow at
the top. The column is initialised isothermal at the wall temperature, which
makes the conductive heat flux at the wall vanish, so the surface energy
balance collapses to a closed form that verify.py checks (see the module
docstring there). The gas starts at rest: the velocity field must be generated
by the blowing itself.

Usage:  python3 build_case.py <config>
"""

import os
import sys

from gsi_models import (MELT, PYRO_ID, mdot_melting, mdot_pyrolysis,
                        wall_temperature_radiative_equilibrium)

# --- discretisation ----------------------------------------------------------
NX, NY = 4, 12
LX, LY = 2.0e-3, 6.0e-3         # [m]
LZ     = 5.0e-4                 # [m] (single cell, 2-D plane)

P0     = 101325.0               # ambient / outflow pressure [Pa]
R_AIR  = 8314.51 / 28.970418    # gas constant of the single-species gas

# Boundary condition type codes (see src/lib/io/IO_BC.f90)
BC_UNUSED, BC_SYMMETRY, BC_EXTRAPOLATION = 0, 300, 400
BC_MELTING, BC_PYROLYSIS, BC_SRM = 503, 504, 502

# Solid propellant grain (BC 502): Saint-Robert's law in mass-flux form,
#   mdot = a * (p/pRef)**n ,  a in kg/(m2 s).
# With pRef at the column pressure the reference mass flux is simply a.
SRM = dict(Taf=1500.0, a=0.2, n=0.35, rho_grain=1750.0, sf=1.0)

# --- configurations ----------------------------------------------------------
# Each entry fixes the wall model and the radiative flux driving the ablation.
# 'wall_face' selects which end of the column ablates: 3 is the low-j face,
# 4 the high-j one. The two are the same physical problem mirrored, but they
# exercise opposite branches of the face-orientation bookkeeping in the BC, so
# both are run.
CONFIGS = {
    'melting':     dict(bc=BC_MELTING,   Tw=800.0, qrad=2.0e5, wall_face=3),
    'htpb':        dict(bc=BC_PYROLYSIS, gsi_id=1, qrad=1.0e6, wall_face=3),
    'hdpe':        dict(bc=BC_PYROLYSIS, gsi_id=2, qrad=1.0e6, wall_face=3),
    'pp':          dict(bc=BC_PYROLYSIS, gsi_id=3, qrad=1.0e6, wall_face=3),
    'melting-top': dict(bc=BC_MELTING,   Tw=800.0, qrad=2.0e5, wall_face=4),
    'hdpe-top':    dict(bc=BC_PYROLYSIS, gsi_id=2, qrad=1.0e6, wall_face=4),
    # No radiation and a gas colder than the wall: the surface loses heat instead of
    # gaining it, the melting energy balance returns a non-positive mass flux, and
    # the boundary condition must fall back to an impermeable isothermal wall.
    'melting-inert': dict(bc=BC_MELTING, Tw=800.0, qrad=0.0, wall_face=3,
                          T_init=400.0, inert=True),
    # Burning propellant grain. The rate comes from the pressure, not from an
    # energy balance, and the wall sits at the flame temperature; the column is
    # still the same exact solution. This is the case that pins the enthalpy of
    # the injected gas: if the products were injected cold the column could not
    # stay at Taf.
    'srm':         dict(bc=BC_SRM, wall_face=3),
}


def expected_wall_state(cfg):
    """Wall temperature and blown mass flux of the isothermal-column solution."""
    if cfg['bc'] == BC_SRM:
        Tw = SRM['Taf']
        mdot = SRM['a']          # pRef is the column pressure, so (p/pRef)**n = 1
    elif cfg['bc'] == BC_MELTING:
        Tw = cfg['Tw']
        mdot = 0.0 if cfg.get('inert') else mdot_melting(cfg['qrad'], Tw=Tw)
    else:
        model = PYRO_ID[cfg['gsi_id']]
        Tw = wall_temperature_radiative_equilibrium(model, cfg['qrad'])
        mdot = mdot_pyrolysis(model, Tw)
    return Tw, mdot


def write_mesh(path):
    ni, nj, nk = NX + 1, NY + 1, 2
    xs = [LX * i / NX for i in range(ni)]
    ys = [LY * j / NY for j in range(nj)]
    zs = [0.0, LZ]

    x, y, z = [], [], []
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                x.append(xs[i]); y.append(ys[j]); z.append(zs[k])

    with open(path, 'w') as f:
        f.write(' VARIABLES ="x" "y" "z"\n')
        f.write(f' ZONE  T = Block1, I={ni}, J={nj}, K={nk}, '
                'DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL)\n')
        for arr in (x, y, z):
            for v in arr:
                f.write(f' {v: .15E}\n')


def write_ic(path, Tw):
    """Uniform state at rest, isothermal at the wall temperature."""
    ni, nj, nk = NX + 1, NY + 1, 2
    ncell = NX * NY * 1
    rho = P0 / (R_AIR * Tw)

    xs = [LX * i / NX for i in range(ni)]
    ys = [LY * j / NY for j in range(nj)]
    zs = [0.0, LZ]
    x, y, z = [], [], []
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                x.append(xs[i]); y.append(ys[j]); z.append(zs[k])

    with open(path, 'w') as f:
        f.write('VARIABLES = "x" "y" "z" "rho1" "u" "v" "w" "p"\n')
        f.write(f'ZONE T="Block1", I={ni}, J={nj}, K={nk}, DATAPACKING=BLOCK, '
                'VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED)\n')
        for arr in (x, y, z):
            for v in arr:
                f.write(f' {v: .15E}\n')
        for value in (rho, 0.0, 0.0, 0.0, P0):     # rho1, u, v, w, p
            for _ in range(ncell):
                f.write(f' {value: .15E}\n')


def write_bc(path, cfg):
    """One record per boundary cell face: 'block i j k face type', plus a
    parameter line for the types that carry data."""
    # Records must be grouped by face and in ascending face order: the solver
    # walks the boundary list with per-face counters, so an out-of-order file
    # silently mis-associates boundary cells with their faces.
    per_face = {f: [] for f in range(1, 7)}

    def rec(face, i, j, k, bctype, params=None):
        per_face[face].append(f'{1:8d}{i:8d}{j:8d}{k:8d}{face:8d}{bctype:8d}')
        if params is not None:
            per_face[face].append(params)

    # faces 1 / 2 : symmetry planes that keep the solution one-dimensional
    for j in range(1, NY + 1):
        rec(1, 1,  j, 1, BC_SYMMETRY)
    for j in range(1, NY + 1):
        rec(2, NX, j, 1, BC_SYMMETRY)

    # the ablating wall, with the outflow at the opposite end of the column
    wall_face = cfg['wall_face']
    wall_j    = 1 if wall_face == 3 else NY
    vent_face = 4 if wall_face == 3 else 3
    vent_j    = NY if wall_face == 3 else 1

    for i in range(1, NX + 1):
        # Full double precision: the verifier checks algebraic identities of the
        # wall model to round-off, which a truncated constant would break.
        if cfg['bc'] == BC_SRM:
            # Taf, a, n, pRef, rhoGrain, SF, mass fractions(1:nsc)
            p = (f"    {SRM['Taf']:.17E},    {SRM['a']:.17E},    {SRM['n']:.17E},"
                 f"    {P0:.17E},    {SRM['rho_grain']:.17E},    {SRM['sf']:.17E},"
                 f"    {1.0:.17E},")
        elif cfg['bc'] == BC_MELTING:
            # cp_wall, Tw, Ti_wall, dh_wall, qrad, eps_wall, mass fractions(1:nsc)
            p = (f"    {MELT['cp']:.17E},    {cfg['Tw']:.17E},    {MELT['Ti']:.17E},"
                 f"    {MELT['dh']:.17E},    {cfg['qrad']:.17E},    {0.0:.17E},"
                 f"    {1.0:.17E},")
        else:
            # pyrolysis model id, qrad, eps_wall, mass fractions(1:nsc)
            p = (f"    {cfg['gsi_id']:d},    {cfg['qrad']:.17E},    {0.0:.17E},"
                 f"    {1.0:.17E},")
        rec(wall_face, i, wall_j, 1, cfg['bc'], p)

    for i in range(1, NX + 1):
        rec(vent_face, i, vent_j, 1, BC_EXTRAPOLATION)

    # faces 5 / 6 : inactive third direction
    for f in (5, 6):
        for j in range(1, NY + 1):
            for i in range(1, NX + 1):
                rec(f, i, j, 1, BC_UNUSED)

    lines = [ln for f in range(1, 7) for ln in per_face[f]]
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in CONFIGS:
        sys.exit(f'usage: build_case.py [{"|".join(CONFIGS)}]')
    name = sys.argv[1]
    cfg = CONFIGS[name]

    os.makedirs('MESH', exist_ok=True)
    os.makedirs('INPUT', exist_ok=True)
    os.makedirs('OUTPUT', exist_ok=True)

    Tw, mdot = expected_wall_state(cfg)
    T_init = cfg.get('T_init', Tw)

    write_mesh('MESH/mesh.tec')
    write_ic('INPUT/ic.tec', T_init)
    write_bc('INPUT/bc.txt', cfg)

    # Record the *configuration* only. verify.py re-derives the reference solution
    # from gsi_models itself, so the check never degenerates into comparing the
    # generator with itself.
    rho = P0 / (R_AIR * Tw)
    with open('INPUT/case.txt', 'w') as f:
        f.write(f'config {name}\n')
        f.write(f'bc {cfg["bc"]}\n')
        f.write(f'qrad {cfg.get("qrad", 0.0)!r}\n')
        f.write(f'srm_Taf {SRM["Taf"]!r}\n')
        f.write(f'srm_a {SRM["a"]!r}\n')
        f.write(f'srm_n {SRM["n"]!r}\n')
        f.write(f'gsi_id {cfg.get("gsi_id", 0)}\n')
        f.write(f'Tw_prescribed {cfg.get("Tw", 0.0)!r}\n')
        f.write(f'wall_face {cfg["wall_face"]}\n')
        f.write(f'T_init {T_init!r}\n')
        f.write(f'inert {int(cfg.get("inert", False))}\n')
        f.write(f'p {P0!r}\n')
        f.write(f'R_gas {R_AIR!r}\n')

    print(f'{name}: Tw = {Tw:.4f} K, mdot = {mdot:.6f} kg/m2/s, '
          f'v_blow = {mdot / rho:.4f} m/s')


if __name__ == '__main__':
    main()
