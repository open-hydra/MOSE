"""Closed-form gas-surface interaction models for the ablating-wall test case.

Single source of truth for both the case generator (build_case.py) and the
verifier (verify.py). The constants mirror src/lib/physics/Lib_GSI.f90 -- if
they are changed there, they must be changed here, and the test will say so by
failing.
"""

from math import exp, sqrt

RUNIV = 8314.51          # universal gas constant [J/(kmol K)], as in FLINT
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
