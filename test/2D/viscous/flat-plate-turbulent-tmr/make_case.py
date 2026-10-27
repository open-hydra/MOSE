#!/usr/bin/env python3
"""
Write INPUT/bc.txt and INPUT/ic.tec of one run of the turbulent flat plate of
the NASA Turbulence Modeling Resource (TMR), for one model and one grid of the
TMR family (numpy is not needed: the files are written from the grid file and
the constants below).

Flow (TMR, 2D zero-pressure-gradient flat plate): M = 0.2, Re = 5 million per
unit length, plate from x = 0 to x = 2, inflow at x = -0.333, top boundary at
y = 1. Here T_inf = 555 K and p_inf = 1e5 Pa, with the constant viscosity of
test/common/Air-viscous, which gives Re = 5.0e6 per metre.

Grids: MESH/flatplate_clust2_<n>levelsdown_<NI>x<NJ>.p2dfmt, the files of the
TMR as published (one block, x and y of the nodes, double precision); each
coarser grid is every other node of the next finer one. The 2-D mesh of MOSE
has two planes of nodes, at z = -5.00000023748726e-4 and +5.00000023748726e-4:
the thickness that the mesh generator of ATLAS writes for a plane case
(0.0005 in single precision), kept so that the inputs are those of the runs
of the V&V page.

Boundaries, one record per boundary cell face (the six integers of a record of
bc.txt, as BCB of ATLAS writes them) and its property line, numbers that MOSE
reads list-directed (17 decimals for the inflow and the top boundary, six for
the outflow, the wall and the relaxation factor and mass fraction of the
inflow):
  face 1 (inflow, x = -0.333): type 401, total pressure and temperature of the
         TMR, p0/p_inf = 1.02828 and T0/T_inf = 1.008, flow normal to the face,
         and the free-stream turbulence of the TMR: SA nu~ = 3 nu_inf; SST and
         Wilcox 2006 k = 9e-9 a_inf^2 and omega = 1e-6 rho_inf a_inf^2/mu_inf;
  face 2 (outflow, x = 2) and face 4 (top, y = 1): type 406 at p_inf. MOSE has
         no far-field condition with a moving free stream (type 407 is a
         reservoir at rest), so the top boundary is an outflow at the
         free-stream pressure; where gas would enter, a type 406 face acts as
         a symmetry face;
  face 3 (y = 0): symmetry (type 300) for x < 0, adiabatic no-slip wall
         (type 301, q = 0) on the plate;
  faces 5 and 6: the two planes of the 2-D mesh (type 0).
The free-stream turbulence uses rho_inf = p_inf/(R T_inf) and
a_inf^2 = gamma R T_inf with R = 287 and gamma = 1.4, and mu_inf to eight
digits (1.1858686e-5 Pa s).

Initial field: the uniform state that ATLAS ICB writes for p0 = 1.02828e5 Pa,
T = 555 K and M = 0.2 with the gas of test/common/Air-viscous, with
nu~ = 1e-6 m^2/s (SA) or k = 2.25e-3 m^2/s^2 and omega = 6250 1/s (SST,
Wilcox 2006), density-weighted as MOSE reads them.

Usage: python3 make_case.py --model SA|SST|Wilcox2006 --grid 35x25|69x49|137x97 [--out DIR]
(DIR defaults to the current folder; DIR/INPUT/bc.txt and DIR/INPUT/ic.tec
are written, replacing files of the same name).
"""
import argparse
import math
import os

MODELS = ('SA', 'SST', 'Wilcox2006')
GRIDS = {'35x25': 4, '69x49': 3, '137x97': 2}   # levels below the finest grid of the TMR family (545x385)

T_INF, P_INF = 555.0, 1.0e5                     # free stream [K, Pa]
MU_INF = 1.1858686e-5                           # Pa s, the viscosity of test/common/Air-viscous to eight digits
R_GAS, GAMMA = 287.0, 1.4
P0_RATIO, T0_RATIO = 1.02828, 1.008             # TMR inflow
K_IC, W_IC = 2.25e-3, 6250.0                    # k [m^2/s^2] and omega [1/s] of the initial field (SST, Wilcox 2006)
Z_PLANE = '0.500000023748726E-003'              # half thickness of the 2-D mesh [m]

# uniform initial state of ATLAS ICB, 15 significant digits (Fortran E23.15E3)
IC_RHO = ' 0.627808435815472E+000'              # kg/m^3
IC_U = ' 0.944451633820311E+002'                # m/s
IC_ZERO = ' 0.000000000000000E+000'
IC_P = ' 0.999998972357490E+005'                # Pa
IC_RHO_NUT = ' 0.627808435815472E-006'          # rho nu~ [kg/(m s)], nu~ = 1e-6 m^2/s

# property lines of the outflow and of the wall, six decimals
BC_OUT = '    0.100000E+06,    0.100000E+01,'     # 406: p, relaxation factor
BC_WALL = '    0.000000E+00,    0.000000E+00,    0.000000E+00,'   # 301: adiabatic wall


def tmr_turbulence(model):
    """Free-stream turbulence of the TMR: [nu~] for SA, [k, omega] for SST and Wilcox 2006."""
    rho = P_INF / (R_GAS * T_INF)
    a2 = GAMMA * R_GAS * T_INF
    if model == 'SA':
        return [3.0 * MU_INF / rho]
    return [9.0e-9 * a2, 1.0e-6 * rho * a2 / MU_INF]


def fortran_e(v, digits=15, exp_digits=3):
    """A value as Fortran writes it with E(digits+8).(digits)E(exp_digits): ' 0.xxxE+eee' or '-0.xxxE+eee'."""
    if v == 0.0:
        return ' 0.' + '0' * digits + 'E+' + '0' * exp_digits
    m, e = ('%.*e' % (digits - 1, abs(v))).split('e')
    e = int(e) + 1
    return ('-' if v < 0 else ' ') + '0.' + m.replace('.', '') + 'E' + ('+' if e >= 0 else '-') + '%0*d' % (exp_digits, abs(e))


def ic_format(v):
    """A positive value of a turbulence band as the initial fields of these runs hold it, e.g. ' 0.141256898058481E-002'."""
    e = int(math.floor(math.log10(v))) + 1
    m = v / 10.0 ** e
    if round(m, 15) >= 1.0:
        m /= 10.0
        e += 1
    elif round(m, 15) < 0.1:
        m *= 10.0
        e -= 1
    return ' %.15fE%s%03d' % (m, '+' if e >= 0 else '-', abs(e))


def read_p2dfmt(path):
    """Nodes of a one-block formatted 2-D Plot3D file: (ni, nj, x, y), x and y lists with i fastest."""
    tok = open(path).read().split()
    if int(tok[0]) != 1:
        raise SystemExit('%s: one block expected' % path)
    ni, nj = int(tok[1]), int(tok[2])
    v = [float(t) for t in tok[3:]]
    if len(v) != 2 * ni * nj:
        raise SystemExit('%s: %d values, expected %d' % (path, len(v), 2 * ni * nj))
    return ni, nj, v[:ni * nj], v[ni * nj:]


def write_ic(path, model, ni, nj, x, y):
    turb = ['rhoNuTilde'] if model == 'SA' else ['rhoK', 'rhoOmega']
    ncell = (ni - 1) * (nj - 1)
    out = [' VARIABLES = x y z  rho1 u v w p ' + ' '.join(turb),
           ' ZONE  T = B1-IG, I=%d, J=%d, K=2, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL,[4-%d]=CELLCENTERED), '
           'SOLUTIONTIME=0.100000000000000E+002, STRANDID = 0' % (ni, nj, 8 + len(turb))]
    for c in (x, y):
        col = [fortran_e(v) for v in c]
        out += col + col                         # the two planes k = 1, 2
    out += ['-' + Z_PLANE] * (ni * nj) + [' ' + Z_PLANE] * (ni * nj)
    bands = [IC_RHO, IC_U, IC_ZERO, IC_ZERO, IC_P]
    if model == 'SA':
        bands.append(IC_RHO_NUT)
    else:
        rho = float(IC_RHO)
        bands += [ic_format(rho * K_IC), ic_format(rho * W_IC)]
    for b in bands:
        out += [b] * ncell
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')


def write_bc(path, model, ni, nj, x):
    nic, njc = ni - 1, nj - 1
    rec = lambda i, j, face, t: ''.join('%8d' % v for v in (1, i, j, 1, face, t))
    t0, p0 = '%.17E' % (T0_RATIO * T_INF), '%.17E' % (P0_RATIO * P_INF)
    inflow = '    ' + ', '.join([t0, p0, 'normal', 'normal', '0.100000E+01', '0.100000E+01'] +
                                ['%.17E' % v for v in tmr_turbulence(model)]) + ','
    top = '    %.17E, %.17E,' % (P_INF, 1.0)
    out = []
    for j in range(1, njc + 1):
        out += [rec(1, j, 1, 401), inflow]
    for j in range(1, njc + 1):
        out += [rec(nic, j, 2, 406), BC_OUT]
    for i in range(1, nic + 1):
        if 0.5 * (x[i - 1] + x[i]) < 0.0:          # face centre upstream of the plate (nodes of row j = 1)
            out.append(rec(i, 1, 3, 300))
        else:
            out += [rec(i, 1, 3, 301), BC_WALL]
    for i in range(1, nic + 1):
        out += [rec(i, njc, 4, 406), top]
    for face in (5, 6):
        for j in range(1, njc + 1):
            for i in range(1, nic + 1):
                out.append(rec(i, j, face, 0))
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--model', required=True, choices=MODELS)
    ap.add_argument('--grid', required=True, choices=list(GRIDS))
    ap.add_argument('--out', default='.')
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    mesh = os.path.join(here, 'MESH', 'flatplate_clust2_%dlevelsdown_%s.p2dfmt' % (GRIDS[a.grid], a.grid))
    ni, nj, x, y = read_p2dfmt(mesh)
    if '%dx%d' % (ni, nj) != a.grid:
        raise SystemExit('%s: %d x %d nodes' % (mesh, ni, nj))
    os.makedirs(os.path.join(a.out, 'INPUT'), exist_ok=True)
    write_ic(os.path.join(a.out, 'INPUT', 'ic.tec'), a.model, ni, nj, x, y)
    write_bc(os.path.join(a.out, 'INPUT', 'bc.txt'), a.model, ni, nj, x)
    print('make_case.py: %s on the %s grid: %s/INPUT/ic.tec and bc.txt' % (a.model, a.grid, a.out))


if __name__ == '__main__':
    main()
