"""Nozzle part of BoundaryConditions (run.sh): the face of a choked nozzle (420) below p_shock and from p_shock up to
p_sub, inputs and a reader of the outputs.

The pressure p of the boundary cell selects the state of the face of a nozzle 420 (docs/user/boundary-conditions.md):
below p_shock, the pressure behind a normal shock standing at the exit of the design flow, the supersonic design exit
state of its record, the isentropic state of the reservoir (T0, p0) with the mass flux g, on the tables of the gas; from
p_shock up to p_sub the state of mass flux g and total temperature T0 that the invariant of the cell reaches. The
records are those of real nozzles, as the boundary-condition builder of ATLAS writes them from an exit-to-throat area
ratio Ae_At: g = G*/Ae_At, G* the choked mass flux, the largest mass flux G(T) = rho(T) u(T) of the isentropic expansion
from (T0, p0); p_sub and p_sup the pressures of the subsonic and of the supersonic exit states, the roots of G = g. All
on the tables, with the arithmetic of the solver: h and cp interpolated linearly between the rows at 1 K, T clamped to
[Tmin, Tmax - 1]; p(T) = p0 exp(-int_T^T0 cp/(R T') dT') by the trapezoidal rule on the integer kelvins and the two
fractional ends; rho = p/(R T); u = sqrt(2 (h0 - h(T))). G* by a scan of 65 points and a golden-section search, the
roots by bisection on the two monotone branches of G, as ATLAS computes them.

The states of mass flux g and total enthalpy h0 lie on the line u(T) = sqrt(2 (h0 - h(T))), p(T) = g R T/u(T), whose
pressure grows with T: p_son is the pressure of its sonic state (u = a, a^2 = gamma R T with the gamma of the tables),
and p_shock that of its state on the subsonic side with the impulse p + g u of the design state, which a normal shock
keeps with g and h0; both by bisection. With the constant cp of the air of the tests, p_shock is the pressure behind a
normal shock at the Mach number M_e of the design state, p_e (1 + 2 gamma (M_e^2 - 1)/(gamma + 1)): the cases of the gas
of the tube on the line take this value, once the two agree within 1e-10 of it.

A 1-D tube along x: NX cells of DX, one cell of DX in y and in z; slip walls (300) on the j faces, null (0) k faces,
the nozzle on face 1 (i = 1) and an extrapolation face (400) at i = NX. The gas of the tube moves in. One forward Euler
step with first-order reconstruction: the change of the conserved variables of cell 1 times DX/DT, plus the flux of the
state of the cell through the interior face (the same state on both sides of it), is the flux through the nozzle face.
It must be that of the state the face takes, each of mass, momentum and energy within 1e-9 of its value:
  - below p_shock, the design state: g, p_e + g u_e and g h0 (h0 the total enthalpy of the reservoir; p_e and u_e of
    the root found here). HLLC gives the face the flux of the design state, on the outer side, when its estimate of
    the left wave speed, min(u_e - a_e, u~ - a~) with u~ and a~ the Roe averages of the two states, is not negative:
    so it is for every gas of the tube below (u~ - a~ is zero for the state behind the normal shock in a gas of
    constant cp);
  - from p_shock up to p_sub, for a gas of the tube on the line, the state of the gas of the tube itself.

Cases (CASES below), the gas linked by run.sh as INPUT/phase.txt and INPUT/thermo.dat:
  mach25, mach3, mach4  air of test/common/Air, T0 = 300 K, p0 = 2 bar, Ae_At that of the exit Mach number 2.5, 3 or
                        4 of a perfect gas with the gamma of this air; the gas of the tube the design state at 0.98 of
                        its pressure, at its temperature and velocity
  rocket45              the gas of the inlet 401 of test/2D/euler/rocket-nozzle/frozen (13 species, its mass fractions
                        normalized to sum 1), T0 = 3644.12 K, p0 = 108 bar, Ae_At = 40: an exit Mach number of 4.46; the
                        gas of the tube as for mach25
  sonic                 the air of mach25 with Ae_At = 1 and g = (1 + 0.5e-9) G*, the gas of the tube at 1.2 times the
                        velocity of the face state: the solver accepts g up to 1e-9 above its G* (the G* of ATLAS, to
                        the rounding), and the face takes the sonic state, the largest mass flux on the tables. Its
                        Mach number with the frozen speed of sound differs from 1 by a few parts in a million (the
                        trapezoidal isentrope), and on the subsonic side HLLC gives the face the mass and the energy
                        fluxes of the star state, which differ from those of the state by that fraction of the jump to
                        the gas of the tube; the momentum flux, the impulse p + g u, stationary at Mach 1, moves to
                        second order only: the momentum alone is checked
  overexp-sup           the air of mach25 with Ae_At = 4 (an exit Mach number of 2.94), the gas of the tube on the line
                        at 2 p_sup, supersonic: an over-expanded exit, the design state
  overexp-sub           the same with the gas of the tube on the line midway between p_son and p_shock, subsonic: the
                        design state
  shock-below           the same at (1 - 1e-6) p_shock: the design state
  shock-above           the same at (1 + 1e-6) p_shock: the state of the gas of the tube
  shock-in              the same at 1.02 p_shock, a normal shock in the divergent part: the state of the gas of the tube
shock-below and shock-above hold the p_shock of the solver within 1e-6 of the normal-shock pressure.

Usage (from the test directory):
  python3 -B design_state.py case NAME          INPUT/ic.tec and INPUT/bc.txt of case NAME
  python3 -B design_state.py check NAME FILE    the fluxes through the nozzle face in FILE (a field.tec) against those of
                                                the state the face takes; exit status 0 when each agrees within 1e-9
"""
import math
import os
import re
import sys

NX = 8
DX = 0.01
RUNIV = 8314.46261815324                      # universal gas constant of FLINT, J/(kmol K)
GAMMA_AIR = 1004.5 / (1004.5 - RUNIV / 28.970418)   # test/common/Air: cp = 1004.5 J/(kg K), molar mass 28.970418
ROCKET_Y = (0.445291, 0.236802, 0.473876e-04, 0.118430e-02, 0.580019e-04, 0.631957e-04, 0.103089e-01, 0.101750e-04,
            0.264504, 0.149121e-04, 0.353148e-02, 0.304376e-01, 0.774545e-02)   # the inlet of rocket-nozzle/frozen
NSCAN = 64
GOLDEN = 0.6180339887498949                   # golden ratio conjugate
P_CELL = 0.98                                 # pressure of the gas of the tube over that of the design state
TOL = 1.0e-9
TOL_SHOCK = 1.0e-10                           # p_shock on the line against the normal-shock pressure


def area_ratio(mach):
    """Exit-to-throat area ratio of the exit Mach number `mach` of a perfect gas with the gamma of the air."""
    g = GAMMA_AIR
    return ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * mach * mach)) ** ((g + 1.0) / (2.0 * (g - 1.0))) / mach


# name: (gas, T0, p0, Ae_At, g over G*/Ae_At, velocity of the gas of the tube over that of the face state, line), with
# line None for the gas of the tube at the temperature of the design state and 0.98 of its pressure, or (pressure,
# factor) for the gas of the tube on the line of g and h0 at factor times p_sup, p_shock or p_mid = (p_son + p_shock)/2
CASES = {
    "mach25": ("air", 300.0, 2.0e5, area_ratio(2.5), 1.0, 1.0, None),
    "mach3": ("air", 300.0, 2.0e5, area_ratio(3.0), 1.0, 1.0, None),
    "mach4": ("air", 300.0, 2.0e5, area_ratio(4.0), 1.0, 1.0, None),
    "rocket45": ("rocket", 3644.12, 1.08e7, 40.0, 1.0, 1.0, None),
    "sonic": ("air", 300.0, 2.0e5, 1.0, 1.0 + 0.5e-9, 1.2, None),
    "overexp-sup": ("air", 300.0, 2.0e5, 4.0, 1.0, 1.0, ("p_sup", 2.0)),
    "overexp-sub": ("air", 300.0, 2.0e5, 4.0, 1.0, 1.0, ("p_mid", 1.0)),
    "shock-below": ("air", 300.0, 2.0e5, 4.0, 1.0, 1.0, ("p_shock", 1.0 - 1.0e-6)),
    "shock-above": ("air", 300.0, 2.0e5, 4.0, 1.0, 1.0, ("p_shock", 1.0 + 1.0e-6)),
    "shock-in": ("air", 300.0, 2.0e5, 4.0, 1.0, 1.0, ("p_shock", 1.02)),
}


class Gas:
    """The tables of INPUT/phase.txt and INPUT/thermo.dat as the solver reads them, for the mass fractions y: rows at
    1 K from Tmin to Tmax, with a row at 0 K (cp of the row at 1 K, h = 0) when the file starts at 1 K; mixture values
    sum_s y_s v_s(T), each v_s interpolated linearly with T clamped to [Tmin, Tmax - 1]; R = sum_s y_s Runiv/W_s."""

    def __init__(self, y):
        names, weights = [], []
        with open("INPUT/phase.txt") as f:
            for line in f.read().splitlines()[1:]:
                words = line.split()
                if len(words) == 2:
                    names.append(words[0])
                    weights.append(float(words[1]))
        zones, cp, h, first = [], [], [], None
        with open("INPUT/thermo.dat") as f:
            for line in f:
                if line.startswith("ZONE"):
                    zones.append(line.split('"')[1])
                    cp.append([])
                    h.append([])
                elif zones and line.strip()[:1].isdigit():
                    words = line.split()
                    if first is None:
                        first = int(float(words[0]))
                    cp[-1].append(float(words[1]))
                    h[-1].append(float(words[2]))
        if zones != names or len(y) != len(names):
            sys.exit(f"the species of INPUT/thermo.dat ({zones}) are not those of INPUT/phase.txt ({names}) and the case")
        self.tmax = first + len(cp[0]) - 1
        self.tmin = first
        if first == 1:
            self.tmin = 0
            for s in range(len(names)):
                cp[s].insert(0, cp[s][0])
                h[s].insert(0, 0.0)
        self.y, self.cp_s, self.h_s = list(y), cp, h
        r = 0.0
        for ys, w in zip(self.y, weights):
            r = r + ys * (RUNIV / w)
        self.r = r
        self.r_s = [RUNIV / w for w in weights]
        # the mixture cp over R T at each integer T, the nodes of the isentrope: there the interpolation is the row
        self.f_node = [0.0] * (self.tmax - self.tmin + 1)
        for k in range(max(self.tmin, 1), self.tmax + 1):
            self.f_node[k - self.tmin] = self.mix(float(k), self.cp_s) / self.r / k

    def mix(self, t, table, y=None):
        tc = max(float(self.tmin), min(t, float(self.tmax - 1)))
        i = int(tc)
        d = tc - i
        k = i - self.tmin
        acc = 0.0
        for ys, col in zip(y or self.y, table):
            acc = acc + ys * (col[k] + (col[k + 1] - col[k]) * d)
        return acc

    def h(self, t, y=None):
        return self.mix(t, self.h_s, y)

    def gamma(self, t):
        cp = self.mix(t, self.cp_s)
        return cp / (cp - self.r)

    def isentrope(self, t, t0):
        """int_t^t0 cp/(R T) dT, t <= t0, by the trapezoidal rule on t, ceil(t), ..., floor(t0), t0."""
        t_prev, f_prev = t, self.mix(t, self.cp_s) / self.r / t
        integral = 0.0
        for k in range(math.ceil(t), math.floor(t0) + 1):
            f_k = self.f_node[k - self.tmin] if k <= self.tmax else self.mix(float(k), self.cp_s) / self.r / k
            integral = integral + 0.5 * (k - t_prev) * (f_prev + f_k)
            t_prev, f_prev = float(k), f_k
        f_end = self.mix(t0, self.cp_s) / self.r / t0
        return integral + 0.5 * (t0 - t_prev) * (f_prev + f_end)

    def pressure(self, t, t0, p0):
        return p0 / math.exp(self.isentrope(t, t0))

    def mass_flux(self, t, t0, p0, h0):
        rho = self.pressure(t, t0, p0) / (self.r * t)
        return rho * math.sqrt(max(2.0 * (h0 - self.h(t)), 0.0))

    def line_velocity(self, t, h0):
        """The velocity of the state at t of the line of total enthalpy h0."""
        return math.sqrt(max(2.0 * (h0 - self.h(t)), 0.0))

    def throat(self, t0, p0, h0):
        """Critical temperature and choked mass flux: a scan of 65 points of [Tlo, t0] and a golden-section search."""
        tlo = float(max(self.tmin, 1))
        kmax, gmax = 0, -1.0
        for k in range(NSCAN + 1):
            g = self.mass_flux(tlo + (t0 - tlo) * k / NSCAN, t0, p0, h0)
            if g > gmax:
                kmax, gmax = k, g
        if kmax in (0, NSCAN):
            sys.exit("the tables do not hold the sonic state of the reservoir")
        ta = tlo + (t0 - tlo) * (kmax - 1) / NSCAN
        tb = tlo + (t0 - tlo) * (kmax + 1) / NSCAN
        tc, td = tb - GOLDEN * (tb - ta), ta + GOLDEN * (tb - ta)
        gc, gd = self.mass_flux(tc, t0, p0, h0), self.mass_flux(td, t0, p0, h0)
        while tb - ta > 1.0e-10 * t0:
            if gc > gd:
                tb, td, gd = td, tc, gc
                tc = tb - GOLDEN * (tb - ta)
                gc = self.mass_flux(tc, t0, p0, h0)
            else:
                ta, tc, gc = tc, td, gd
                td = ta + GOLDEN * (tb - ta)
                gd = self.mass_flux(td, t0, p0, h0)
        t = 0.5 * (ta + tb)
        return t, self.mass_flux(t, t0, p0, h0)

    def root(self, t_pos, t_neg, g, t0, p0, h0):
        """The root of G(T) = g between t_pos (G > g) and t_neg (G < g), by bisection to the spacing of the doubles."""
        a, b = t_pos, t_neg
        for _ in range(200):
            if abs(b - a) <= 2.0 * sys.float_info.epsilon * t0:
                break
            m = 0.5 * (a + b)
            if self.mass_flux(m, t0, p0, h0) > g:
                a = m
            else:
                b = m
        return 0.5 * (a + b)


def bisect(f, lo, hi):
    """A root of f between lo (f < 0) and hi (f > 0), by bisection to the spacing of the doubles."""
    for _ in range(200):
        if hi - lo <= 4.0 * sys.float_info.epsilon * hi:
            break
        m = 0.5 * (lo + hi)
        if f(m) < 0.0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def nozzle(name):
    """The gas, the record (T0, p0, p_sub, p_sup, g, mass fractions), the design state, p_son and p_shock of the line of
    g and h0, and the state of the tube."""
    gas_name, t0, p0, ae_at, factor, ufac, line = CASES[name]
    if gas_name == "air":
        y = (1.0,)
    else:
        total = sum(ROCKET_Y)
        y = tuple(v / total for v in ROCKET_Y)
    gas = Gas(y)
    h0 = gas.h(t0)
    tlo = float(max(gas.tmin, 1))
    tstar, gstar = gas.throat(t0, p0, h0)
    g = factor * gstar / ae_at
    if ae_at == 1.0:                              # an exit at the throat: the two exit states are the sonic one
        t_sub = t_sup = tstar
    else:
        t_sub = gas.root(tstar, t0, g, t0, p0, h0)
        t_sup = gas.root(tstar, tlo, g, t0, p0, h0)
    p_sub, p_sup = gas.pressure(t_sub, t0, p0), gas.pressure(t_sup, t0, p0)
    rho_e = p_sup / (gas.r * t_sup)
    u_e = g / rho_e
    # the line of g and h0: its sonic state, and its state behind a normal shock at the design exit, whose impulse is
    # that of the design state: the root on [T*, T0] of f = g R T + g u^2 - I_e u = u (I - I_e), none if f(T*) >= 0
    t_son = bisect(lambda t: gas.h(t) + 0.5 * gas.gamma(t) * gas.r * t - h0, tlo, t0)
    p_son = g * gas.r * t_son / gas.line_velocity(t_son, h0)
    i_e = p_sup + g * u_e

    def impulse_gap(t):
        u = gas.line_velocity(t, h0)
        return g * gas.r * t + g * u * u - i_e * u

    if impulse_gap(tstar) >= 0.0:
        p_shock = p_sup
    else:
        t_2 = bisect(impulse_gap, tstar, t0)
        p_shock = g * gas.r * t_2 / gas.line_velocity(t_2, h0)
    if line is None:
        cell = dict(t=t_sup, p=P_CELL * p_sup, u=ufac * u_e)
    else:
        gam = gas.gamma(t_sup)
        mach_e = u_e / math.sqrt(gam * gas.r * t_sup)
        p_normal = p_sup * (1.0 + 2.0 * gam * (mach_e * mach_e - 1.0) / (gam + 1.0))
        if abs(p_shock - p_normal) > TOL_SHOCK * p_normal:
            sys.exit(f"{name}: p_shock on the line {p_shock:.15e} Pa differs from the normal-shock pressure "
                     f"{p_normal:.15e} Pa by more than {TOL_SHOCK:g} of it: the gas has no constant cp")
        where, fac = line
        p = fac * {"p_sup": p_sup, "p_shock": p_normal, "p_mid": 0.5 * (p_son + p_normal)}[where]
        t = bisect(lambda t: g * gas.r * t / gas.line_velocity(t, h0) - p, tlo, t0)
        cell = dict(t=t, p=p, u=gas.line_velocity(t, h0))
    cell["rho"] = cell["p"] / (gas.r * cell["t"])
    return dict(gas=gas, y=y, t0=t0, p0=p0, h0=h0, p_sub=p_sub, p_sup=p_sup, g=g, t_e=t_sup, p_e=p_sup, u_e=u_e,
                p_son=p_son, p_shock=p_shock, cell=cell)


def write_case(name):
    n = nozzle(name)
    c, y = n["cell"], n["y"]
    for p in ("INPUT/ic.tec", "INPUT/bc.txt"):
        if os.path.lexists(p):
            os.remove(p)
    with open("INPUT/ic.tec", "w") as f:
        rhos = " ".join(f"rho{s + 1}" for s in range(len(y)))
        f.write(f' VARIABLES ="x" "y" "z" {rhos} "u" "v" "w" "p"\n')
        f.write(f" ZONE  T = B1-IG, I={NX + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                f"VARLOCATION=([1-3]=NODAL,[4-{7 + len(y)}]=CELLCENTERED), SOLUTIONTIME=0.0\n")
        for axis in range(3):
            for k in range(2):
                for j in range(2):
                    for i in range(NX + 1):
                        f.write(f" {(DX * i, DX * j, DX * k)[axis]:.16E}\n")
        for v in [c["rho"] * ys for ys in y] + [c["u"], 0.0, 0.0, c["p"]]:
            f.write(f" {v:.16E}\n" * NX)
    with open("INPUT/bc.txt", "w") as f:
        f.write(f"       1       1       1       1       1     420\n")
        massf = ", ".join(f"{ys:.16E}" for ys in y)
        f.write(f"    {n['t0']:.16E}, {n['p0']:.16E}, {n['p_sub']:.16E}, {n['p_sup']:.16E}, {n['g']:.16E}, 1.0, "
                f"{massf},\n")
        f.write(f"       1     {NX:3d}       1       1       2     400\n")
        for face, t in ((3, 300), (4, 300), (5, 0), (6, 0)):
            for i in range(1, NX + 1):
                f.write(f"       1     {i:3d}       1       1       {face}     {t:3d}\n")
    print(f"    {name}: g = {n['g']:.10e} kg/(m2 s), p_sub = {n['p_sub']:.6e} Pa, p_sup = {n['p_sup']:.6e} Pa, "
          f"p_shock = {n['p_shock']:.6e} Pa, design exit T = {n['t_e']:.6f} K; the gas of the tube at "
          f"p = {c['p']:.6e} Pa, Mach {c['u'] / math.sqrt(n['gas'].gamma(c['t']) * n['gas'].r * c['t']):.4f}")


def read_field(path, nsp):
    """Solution time and the values of cell 1 of a one-zone MOSE ASCII field (BLOCK packing), names quoted or not."""
    with open(path) as f:
        lines = f.read().splitlines()
    names = [a or b for a, b in re.findall(r'"([^"]*)"|([^\s,"]+)', lines[0].split('=', 1)[1])]
    time = float(re.search(r'SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)', lines[1]).group(1).replace('D', 'E'))
    values = [float(v.replace('D', 'E')) for v in " ".join(l for l in lines[2:] if not l.startswith("DATASETAUXDATA")).split()]
    nn, first, pos = (NX + 1) * 2 * 2, {}, 0
    for k, var in enumerate(names):
        first[var] = values[pos]                    # a node value for the coordinates, that of cell 1 for the others
        pos += nn if k < 3 else NX
    rho_s = [first[f"rho({s + 1})"] for s in range(nsp)]
    return time, rho_s, first["u"], first["T"]


def conserved(gas, rho_s, u, t):
    """Density, momentum and total energy per unit volume of a state, its energy from the tables at its T."""
    rho = sum(rho_s)
    y = [r / rho for r in rho_s]
    r_mix = 0.0
    for ys, rs in zip(y, gas.r_s):
        r_mix = r_mix + ys * rs
    return rho, rho * u, rho * (gas.h(t, y) - r_mix * t) + 0.5 * rho * u * u


def check(name, path):
    n = nozzle(name)
    gas, c = n["gas"], n["cell"]
    dt, rho_s, u, t = read_field(path, len(n["y"]))
    new = conserved(gas, rho_s, u, t)
    old = conserved(gas, [c["rho"] * ys for ys in n["y"]], c["u"], c["t"])
    # the interior face takes the flux of the state of the cell, the same on both sides of it
    h_c = gas.h(c["t"])
    inner = (c["rho"] * c["u"], c["rho"] * c["u"] ** 2 + c["p"], c["rho"] * c["u"] * (h_c + 0.5 * c["u"] ** 2))
    face = [(a - b) * DX / dt + f for a, b, f in zip(new, old, inner)]
    if c["p"] < n["p_shock"]:                    # the design state
        want, state = (n["g"], n["p_e"] + n["g"] * n["u_e"], n["g"] * n["h0"]), "the design state"
    else:                                         # the state of the gas of the tube, on the line of g and h0
        want, state = inner, "the gas of the tube"
    err = [abs(a - b) / abs(b) for a, b in zip(face, want)]
    if CASES[name][3] == 1.0:                     # the sonic state: the momentum alone (module docstring)
        ok = err[1] <= TOL
        print(f"    {name}: through the nozzle face, momentum {face[1] / want[1]:.12f} (p* + g u*); difference "
              f"{err[1]:.1e}{'' if ok else '  <-- differs'}")
    else:
        ok = all(e <= TOL for e in err)
        print(f"    {name}: through the nozzle face, mass {face[0] / want[0]:.12f}, momentum {face[1] / want[1]:.12f}, "
              f"energy {face[2] / want[2]:.12f} of {state}; largest difference {max(err):.1e}"
              f"{'' if ok else '  <-- differs'}")
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2])
    elif what == "check":
        check(sys.argv[2], sys.argv[3])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
