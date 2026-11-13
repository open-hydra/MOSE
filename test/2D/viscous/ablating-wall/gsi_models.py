"""Closed-form gas-surface interaction models for the ablating-wall test case.

Single source of truth for both the case generator (build_case.py) and the
verifier (verify.py). The constants mirror src/lib/physics/Lib_GSI.f90 -- if
they are changed there, they must be changed here, and the test will say so by
failing.
"""

from math import exp, sqrt

RUNIV = 8314.46261815324  # universal gas constant [J/(kmol K)], as in FLINT
P_ATM = 101325.0

# --- melting (BC 503), paraffin-like -----------------------------------------
MELT = dict(cp=1946.0310193884, Ti=298.15, dh=1.698285789316e5)

# --- pyrolysis (BC 504) ------------------------------------------------------
HTPB = dict(cp=1632.0, dh=1.1e6, Ti=298.15, rho=960.0,
            A1=3.965, Ea1=55.8564e6, A2=11.04e-3, Ea2=20.54344e6, Ts=722.0)
HDPE = dict(cp=1255.2, dh=2.72e6, Ti=298.15, A=4588.8e3, Ea=251039.73e3)
PP   = dict(cp=1700.0, dh=2.4823e6, Ti=298.15, A=2.12e15, Ea=2.12e8)

PYRO_ID = {1: 'HTPB', 2: 'HDPE', 3: 'PP'}


def mdot_pyrolysis(model, Tw):
    """Blown mass flux [kg/(m2 s)] from the Arrhenius law of the given model."""
    if model == 'HTPB':
        c = HTPB
        A, Ea = (c['A1'], c['Ea1']) if Tw <= c['Ts'] else (c['A2'], c['Ea2'])
        return c['rho'] * A * exp(-Ea / (RUNIV * Tw))
    if model == 'HDPE':
        c = HDPE
        # note the factor 2 in the exponent denominator, as in GSI_HDPE
        return c['A'] * exp(-c['Ea'] / (2.0 * RUNIV * Tw))
    if model == 'PP':
        c = PP
        arrh = c['Ea'] / (RUNIV * Tw)
        rr = sqrt(4.60517 * (1.0 - c['Ti'] / Tw + c['dh'] / (c['cp'] * Tw))
                  - c['dh'] / (c['cp'] * Tw))
        return 910.0 * sqrt(c['A'] * exp(-arrh) / arrh * 0.2 / (910.0 * c['cp']) / rr)
    raise ValueError(model)


def q_pyrolysis(model, Tw):
    """Heat absorbed by the pyrolysing surface [W/m2]."""
    c = {'HTPB': HTPB, 'HDPE': HDPE, 'PP': PP}[model]
    return mdot_pyrolysis(model, Tw) * (c['dh'] + c['cp'] * (Tw - c['Ti']))


def wall_temperature_radiative_equilibrium(model, qrad, lo=300.0, hi=2500.0):
    """Wall temperature of a pyrolysing surface heated only by radiation.

    With the gas column isothermal at the wall temperature the conductive flux
    vanishes, so the surface energy balance collapses to the scalar equation

        q_pyro(Tw) = qrad

    which is solved here by bisection. q_pyro is monotonically increasing in Tw,
    so the root is unique.
    """
    f = lambda T: q_pyrolysis(model, T) - qrad
    if f(lo) > 0.0 or f(hi) < 0.0:
        raise ValueError(f'{model}: qrad={qrad:g} outside the bracketed range')
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def mdot_melting(qwall, cp=None, Tw=None, Ti=None, dh=None):
    """Melting mass flux: the heat reaching the surface divided by the heat
    needed to raise the solid to Tw and melt it (GSI_Melting)."""
    cp = MELT['cp'] if cp is None else cp
    Ti = MELT['Ti'] if Ti is None else Ti
    dh = MELT['dh'] if dh is None else dh
    return qwall / (dh + cp * (Tw - Ti))


# --- surface reactions of a carbon wall (BC 505) and the pyrolysing ablator (BC 506) --------------------------------
# GSI_Carbon_Bradley and GSI_Char_Pyrolysis of src/lib/physics/Lib_GSI.f90, on the gas of a phase.txt and a thermo.dat
# (test/common/CHO-N2-viscous), and the pyrolysis-gas table of FLINT_Lib_PyrolysisGas.
P_REF = 1.01325e5                 # Lib_GSI.f90 pref: the partial pressures of the rates in atmospheres
MW_C = 12.0107                    # Lib_GSI.f90 mwC
BRADLEY = ('O2', 'O', 'OH', 'H', 'H2', 'H2O', 'CO', 'CO2')


def read_phase(path):
    """Species names and molar masses of a phase.txt (line 1: the type of the phase), in the order of the file."""
    names, W = [], []
    with open(path) as f:
        next(f)
        for line in f:
            if line.strip():
                n, w = line.split()[:2]
                names.append(n)
                W.append(float(w))
    return names, W


class Thermo:
    """h_s(T) of a thermo.dat of ATLAS GPB (one zone per species, rows 'T cp h s dcp', 1 K step), read and interpolated as
    FLINT does (read_idealgas_thermo, f_tabT): a table that starts at 1 K gets a row at 0 K with h = 0; T is clamped to
    [Tmin, Tmax - 1], and h = h(i) + (h(i+1) - h(i)) (T - i), i = int(T)."""

    def __init__(self, path):
        self.h, cur, self.tmin = {}, None, None
        with open(path) as f:
            for line in f:
                if line.startswith('ZONE'):
                    cur = line.split('"')[1]
                    self.h[cur] = []
                elif cur is not None and line[:1].isdigit():
                    v = line.split()
                    if self.tmin is None:
                        self.tmin = int(float(v[0]))
                    self.h[cur].append(float(v[2]))
        self.tmax = self.tmin + len(next(iter(self.h.values()))) - 1
        if self.tmin == 1:
            for s in self.h:
                self.h[s].insert(0, 0.0)
            self.tmin = 0

    def __call__(self, s, T):
        tc = max(float(self.tmin), min(T, float(self.tmax - 1)))
        i = int(tc)
        j = min(i + 1, self.tmax)
        a, b = self.h[s][i - self.tmin], self.h[s][j - self.tmin]
        return a + (b - a) * (tc - i)


def bradley(rho_s, T, names, W):
    """GSI_Carbon_Bradley: the species sources omega (kg/(m2 s), on the order of names, positive into the gas: out of the
    solid) and the net mass flux mdot, the carbon removed, for the partial densities rho_s of the gas at the wall."""
    idx = {n: k for k, n in enumerate(names)}
    pi = [rho_s[k] * T * (RUNIV / W[k]) / P_REF for k in range(len(names))]
    RT, isq = RUNIV * T, 1.0 / sqrt(T)
    kH2O = 4.80e5 * exp(-288.00e6 / RT)
    kCO2 = 9.00e3 * exp(-285.00e6 / RT)
    kOH, kO = 3.61e2 * isq, 6.655e2 * isq
    k5 = 2.40e3 * exp(-125.60e6 / RT)
    k6 = 2.13e1 * exp(17.17e6 / RT)
    k7 = 5.35e-1 * exp(-63.64e6 / RT)
    k8 = 1.81e7 * exp(-406.10e6 / RT)
    p = {s: pi[idx[s]] for s in BRADLEY}
    x = 1.0 / (1.0 + k8 / (k7 * p['O2'])) if k7 * p['O2'] != 0.0 else 0.0
    d = dict(H2O=kH2O * sqrt(p['H2O']), CO2=kCO2 * sqrt(p['CO2']), OH=kOH * p['OH'], O=kO * p['O'],
             O2=x * k5 * p['O2'] / (1.0 + k6 * p['O2']) + k7 * p['O2'] * (1.0 - x))
    Wm = {s: W[idx[s]] for s in BRADLEY}
    om = [0.0] * len(names)
    om[idx['O2']] = -(Wm['O2'] / MW_C) * (0.5 * d['O2'])
    om[idx['H2O']] = -(Wm['H2O'] / MW_C) * d['H2O']
    om[idx['CO']] = (Wm['CO'] / MW_C) * (2.0 * d['CO2'] + d['H2O'] + d['O2'] + d['OH'] + d['O'])
    om[idx['CO2']] = -(Wm['CO2'] / MW_C) * d['CO2']
    om[idx['H2']] = (Wm['H2'] / MW_C) * d['H2O']
    om[idx['O']] = -(Wm['O'] / MW_C) * d['O']
    om[idx['H']] = (Wm['H'] / MW_C) * d['OH']
    om[idx['OH']] = -(Wm['OH'] / MW_C) * d['OH']
    return om, sum(om)


class PyrolysisTable:
    """One table of a pyrolysis.dat (format 1 of FLINT_Lib_PyrolysisGas, without the condensed column): y_pg(p, T) on the
    species order of names, bilinear in (p, T), clamped to the table, as pyrolysis_gas_composition."""

    def __init__(self, path, names, table=1):
        tok = [l.split('#')[0].split() for l in open(path)]
        tok = [t for t in tok if t]
        k = 2
        for _ in range(table - 1):
            while tok[k] != ['end', 'table']:
                k += 1
            k += 1
        nsp = int(tok[k + 4][1])
        self.species = [tok[k + 5 + i][0] for i in range(nsp)]
        k += 5 + nsp
        self.p = [float(v) for v in tok[k + 1]]
        self.T = [float(v) for v in tok[k + 3]]
        rows = tok[k + 5:k + 5 + len(self.p) * len(self.T)]
        self.y = [[[0.0] * len(names) for _ in self.T] for _ in self.p]
        for ip in range(len(self.p)):
            for iT in range(len(self.T)):
                for s, v in zip(self.species, rows[ip * len(self.T) + iT]):
                    self.y[ip][iT][names.index(s)] = float(v)

    @staticmethod
    def _cell(x, g):
        lo, hi = 0, len(g) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if x >= g[mid]:
                lo = mid
            else:
                hi = mid
        return lo

    def __call__(self, p, T):
        Tc = min(max(T, self.T[0]), self.T[-1])
        i = self._cell(Tc, self.T)
        wT = (Tc - self.T[i]) / (self.T[i + 1] - self.T[i])
        lin = lambda r: [(1.0 - wT) * a + wT * b for a, b in zip(r[i], r[i + 1])]
        if len(self.p) == 1:
            return lin(self.y[0])
        pc = min(max(p, self.p[0]), self.p[-1])
        j = self._cell(pc, self.p)
        wp = (pc - self.p[j]) / (self.p[j + 1] - self.p[j])
        return [(1.0 - wp) * a + wp * b for a, b in zip(lin(self.y[j]), lin(self.y[j + 1]))]


def ablator(rho_s, T, p, names, W, h, phi=0.0, h_v0=0.0, table=None):
    """GSI_Char_Pyrolysis: omega, mdot and the heat the surface absorbs, q = sum_s omega_s h_s(T) - mdot h_v0, with
    omega = omega_ch + phi mdot_ch y_pg(p, T) (table read only when phi > 0); with phi = 0 and h_v0 = 0, BC 505."""
    om, md = bradley(rho_s, T, names, W)
    if phi > 0.0:
        y = table(p, T)
        om = [o + phi * md * ys for o, ys in zip(om, y)]
        md = md + phi * md
    q = sum(o * h(s, T) for o, s in zip(om, names)) - md * h_v0
    return om, md, q


def wall_temperature_balance(f, T0, rel=0.3):
    """The root of f(T) = 0 in [T0 (1 - rel), T0 (1 + rel)] by bisection; ValueError if f does not change sign there."""
    lo, hi = T0 * (1.0 - rel), T0 * (1.0 + rel)
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0.0:
        raise ValueError(f'no root of the balance within {rel:g} of {T0:g} K')
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (f(mid) < 0.0) == (flo < 0.0):
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


# The pyrolysis-gas table of the 506 columns (build_case.py writes it as INPUT/pyrolysis.dat): one table, gas-only, the
# states of the shift CO + H2O = CO2 + H2 of 1 CO + 1 H2O + 0.5 H2 (moles) at the extent xi = 0.05 + 0.10 (iT - 1) +
# 0.05 (ip - 1), on 2 pressures x 3 temperatures: every row has the elements of that mixture (C 1, H 3, O 2 moles).
PYRO_SPECIES = ('CO', 'H2', 'H2O', 'CO2')
PYRO_P = (5.0e4, 2.0e5)
PYRO_T = (1000.0, 2000.0, 3000.0)


def write_pyrolysis_table(path, names, W, condensed=False):
    """INPUT/pyrolysis.dat (format 1 of FLINT_Lib_PyrolysisGas) with the molar masses W of the gas (names); condensed: the
    same table as treatment condensed, every row (1 - X_C) y with the condensed-carbon column X_C = 0.1 (refused by the
    wall 506)."""
    w = {s: W[names.index(s)] for s in PYRO_SPECIES}
    WC, WH, WO = 12.011, 1.008, 15.999            # the element weights of the molar masses of phase.txt
    M = WC + 3 * WH + 2 * WO
    e = lambda x: '%24.16E' % x
    xc = 0.1 if condensed else 0.0
    out = ['# pyrolysis-gas table of the ablating-wall columns (gsi_models.py)', 'format 1', 'tables 1', 'table shift',
           'treatment condensed' if condensed else 'treatment gas-only',
           f'elements 3 C {WC / M:.16E} H {3 * WH / M:.16E} O {2 * WO / M:.16E}', f'condensed {int(condensed)}',
           f'species {len(PYRO_SPECIES)}', *PYRO_SPECIES, f'pressures {len(PYRO_P)}', ' '.join(e(p) for p in PYRO_P),
           f'temperatures {len(PYRO_T)}', ' '.join(e(T) for T in PYRO_T), 'data']
    for ip in range(len(PYRO_P)):
        for iT in range(len(PYRO_T)):
            xi = 0.05 + 0.10 * iT + 0.05 * ip
            n = dict(CO=1 - xi, H2=0.5 + xi, H2O=1 - xi, CO2=xi)
            m = {s: n[s] * w[s] for s in PYRO_SPECIES}
            tot = sum(m.values())
            out.append(' '.join(e((1.0 - xc) * m[s] / tot) for s in PYRO_SPECIES) + (' ' + e(xc) if condensed else ''))
    out.append('end table')
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')
