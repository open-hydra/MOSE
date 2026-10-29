"""
Turbulent flat plate of the NASA Turbulence Modeling Resource (TMR): SA, SST
and Wilcox 2006 on the 35x25, 69x49 and 137x97 grids of the TMR family.

Quantities, from the wall output of the plate (face 3): c_f = |tau_x|/q_inf
at the face centres of the plate (x > 0), with q_inf = gamma p_inf M^2/2 =
2800 Pa; c_f(0.970084071) and c_f(1.90334) by linear interpolation between
the two nearest face centres; C_D = sum c_f dx / 2 over the plate (length 2).
Iterative error e_it of a quantity: largest minus smallest value over the last
three numbered solutions and the final one.

Three grids (f1 137x97, f2 69x49, f3 35x25, ratio 2), as in the uncertainty
summary of the TMR (Celik et al. 2008): e21 = f2 - f1, e32 = f3 - f2,
s = sign(e32 e21), observed order p = ln|e32/e21| / ln 2, limit
f_ext = (2^p f1 - f2)/(2^p - 1), GCI = 1.25 |e21/f1|/(2^p - 1),
uncertainty U = GCI |f1| + e_it (e_it the largest of the three grids); the
formal order 2 and the factor 3 replace p and 1.25 when e21 and e32 differ in
sign, when p < 0.5 or p > 4, when |e21| >= |e32|, or when |e21| or |e32| is
below 10 e_it.
Criteria of the V&V page:
  (d) e_it <= GCI |f1|/10, for c_f(0.970084071), c_f(1.90334) and C_D;
  (1) on each grid, |c_f - (c_CFL3D + c_FUN3D)/2| <= |c_CFL3D - c_FUN3D|/2
      + m c_f for c_f(0.970084071), against CFL3D and FUN3D on the same grid;
      required on the 137x97 grid, reported on the two coarser grids;
  (2) |f_ext - c_CFL3D,545x385| <= U + m c_f,1 for c_f(0.970084071);
with m = 0.5 % (SA) and 1 % (SST, Wilcox 2006), for the differences of the
case from the TMR (see the V&V page).

Short form (./MOSE.sh testfast): for each model, the restart from the
converged 35x25 solution (reference/field-35x25-<model>.tec) must
  1. have done the iterations of its input: its final solution at the global
     iteration count c0 + n (c0 that of the solution it starts from, n the
     iter-threshold of its input.ini), and n new rows in residual-history.dat;
  2. keep c_f(0.970084071), c_f(1.90334) and C_D within TOL_SHORT (relative)
     of the values of that solution (reference/levels.json).
Then it computes again, from the stored results of the full form (the nine
final wall outputs reference/wall-<model>-<grid>.tec and reference/levels.json),
the tables and the criteria of the V&V page, and writes its figures in
RUN/figures/: it fails if a stored wall output does not give its stored values
(to 1e-12, relative), or if (d), (1) on the 137x97 grid or (2) fails for a
model.
Full form (./MOSE.sh test): the same analysis on the nine runs in RUN/, the
figures, and RUN/reference/ with the files of reference/ that the runs give
(levels.json, the nine final wall outputs and the three converged 35x25
solutions), to be copied to reference/ when the solver changes the steady
state; it fails if a run is missing or if (d), (1) on the 137x97 grid or (2)
fails for a model.

Usage:
    python verify.py            # report, figures in RUN/figures/, exit 0/1

Exit 0 = pass, 1 = fail.
"""
import glob
import json
import math
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

MODELS = ("SA", "SST", "Wilcox2006")
GRIDS = ("35x25", "69x49", "137x97")                 # coarse to fine
N_TMR = {"35x25": 816, "69x49": 3264, "137x97": 13056}   # cells, as the TMR files count the grids
N_545 = 208896                                       # the 545x385 grid
Q_INF = 0.5 * 1.4 * 1.0e5 * 0.2 ** 2                 # gamma p_inf M^2 / 2 [Pa]
U_INF = 0.2 * math.sqrt(1.4 * 287.0 * 555.0)         # M a_inf [m/s]
STATIONS = (0.970084071, 1.90334)                    # c_f stations of the TMR [m]
KEYS = ("cf(0.970084071)", "cf(1.90334)", "C_D", "y+max")
UE_X, UE_Y = (0.97008, 1.90334), 0.1                 # stations and height of the TMR velocity profiles [m]
MARGIN = {"SA": 0.005, "SST": 0.01, "Wilcox2006": 0.01}
P_FORMAL = 2.0
TOL_SHORT = 1.0e-5
TMR_CONV = {"SA": "reference/tmr/SA/cf_convergence.dat", "SST": "reference/tmr/SST/cf_convergence_sstv.dat",
            "Wilcox2006": "reference/tmr/W06/cf_convergence_w06.dat"}
TMR_PLATE = {"SA": "reference/tmr/SA/cf_plate.dat", "SST": "reference/tmr/SST/cf_plate_sstv.dat",
             "Wilcox2006": "reference/tmr/W06/cf_plate_w06.dat"}
TMR_U = "reference/tmr/SA/flatplate_u.dat"
TMR_NAME = {"SA": "SA", "SST": "SST-Vm", "Wilcox2006": "Wilcox2006-klim-m"}

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break
try:
    from ORION import read_TEC
except ModuleNotFoundError:
    raise ModuleNotFoundError("ORION package not found. Run 'conda activate base' or install ORION.")


# ---------------------------------------------------------------------------- reading
def tec_time(path):
    """SOLUTIONTIME of the first zone and the DATASETAUXDATA records of a Tecplot ASCII file."""
    t, aux = None, {}
    with open(path) as fh:
        for line in fh:
            if t is None and line.lstrip().upper().startswith("ZONE"):
                m = re.search(r"SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)", line)
                t = float(m.group(1).replace("D", "E").replace("d", "e")) if m else None
            elif line.startswith("DATASETAUXDATA"):
                m = re.match(r'DATASETAUXDATA\s+(\S+)\s*=\s*"([^"]*)"', line)
                if m:
                    aux[m.group(1)] = m.group(2)
    return t, aux


def wall_quantities(path):
    """c_f at the stations, C_D and the largest y+ on the plate, from a wall output with one zone (face 3)."""
    xw, _, _, vw, names = read_TEC(path)
    if len(xw) != 1:
        raise SystemExit("verify.py: %s holds %d zones, one expected (the plate, face 3)" % (path, len(xw)))
    xn = xw[0][:, 0, 0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    dx = np.diff(xn)
    tau = np.abs(vw[0][names.index("tauX") - 3][:, 0, 0])
    m = xc > 0.0
    cf = tau[m] / Q_INF
    out = {KEYS[0]: float(np.interp(STATIONS[0], xc[m], cf)), KEYS[1]: float(np.interp(STATIONS[1], xc[m], cf))}
    out["C_D"] = float(np.sum(cf * dx[m]) / 2.0)
    out["y+max"] = float(np.max(vw[0][names.index("y+") - 3][m, 0, 0]))
    return out


def wall_cf(path):
    """x of the face centres of the plate and c_f there."""
    xw, _, _, vw, names = read_TEC(path)
    xn = xw[0][:, 0, 0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    m = xc > 0.0
    return xc[m], np.abs(vw[0][names.index("tauX") - 3][:, 0, 0])[m] / Q_INF


def outer_velocity(path):
    """u/U_inf at the stations UE_X, on the row of cells nearest y = UE_Y (interpolated along the row)."""
    xf, yf, _, vf, names = read_TEC(path)
    x, y, u = xf[0][:, :, 0], yf[0][:, :, 0], vf[0][names.index("u") - 3][:, :, 0]
    xc = 0.25 * (x[:-1, :-1] + x[1:, :-1] + x[:-1, 1:] + x[1:, 1:])
    yc = 0.25 * (y[:-1, :-1] + y[1:, :-1] + y[:-1, 1:] + y[1:, 1:])
    out = []
    for xs in UE_X:
        i = int(np.argmin(np.abs(xc[:, 0] - xs)))
        j = int(np.argmin(np.abs(yc[i, :] - UE_Y)))
        out.append(float(np.interp(xs, xc[:, j], u[:, j])) / U_INF)
    return out


def residuals(path):
    """Density residual of residual-history.dat: largest, last, fall in orders, range over the last 2000 rows."""
    r = np.loadtxt(path, usecols=(0, 2))
    last = r[-2000:, 1]
    return dict(rows=int(r.shape[0]), first_row=int(r[0, 0]), last_row=int(r[-1, 0]), max=float(r[:, 1].max()),
                last=float(r[-1, 1]), orders=float(math.log10(r[:, 1].max() / r[-1, 1])),
                last2000_min=float(last.min()), last2000_max=float(last.max()))


def tmr_convergence(path):
    """{code: {N: c_f(0.97)}} of a cf_convergence file of the TMR."""
    data, zone = {}, None
    for ln in open(path):
        s = ln.strip()
        if s.lower().startswith("zone"):
            zone = re.search(r't\s*=\s*"([^"]+)"', s, re.I).group(1)
            data[zone] = {}
        elif zone and s and s[0].isdigit():
            w = s.split()
            data[zone][int(float(w[0]))] = float(w[3])
    return data


def tmr_zones(path, cols=(0, 1)):
    """{zone title: array of the columns} of a multi-zone Tecplot ASCII data file of the TMR."""
    zones, label, rows = {}, None, []
    for ln in open(path):
        s = ln.strip()
        if not s or s.upper().startswith("VARIABLES"):
            continue
        if s.upper().startswith("ZONE"):
            if label is not None:
                zones[label] = np.array(rows)
            label, rows = re.search(r't\s*=\s*"([^"]+)"', s, re.I).group(1).strip(), []
        else:
            w = s.split()
            rows.append([float(w[c]) for c in cols])
    if label is not None:
        zones[label] = np.array(rows)
    return zones


def tmr_outer_velocity():
    """u/U_ref of CFL3D (545x385) at y = UE_Y at the two stations of FlatPlate/SA/flatplate_u.dat."""
    out = []
    for label, a in tmr_zones(TMR_U).items():
        a = a[np.argsort(a[:, 1])]
        out.append((float(label.split("=")[1]), float(np.interp(UE_Y, a[:, 1], a[:, 0]))))
    return [u for _, u in sorted(out)]


# ---------------------------------------------------------------------------- analysis
def run_series(d):
    """Values of the final wall output of a run of the full form, iterative error, counts and residuals."""
    fin = wall_quantities(os.path.join(d, "OUTPUT", "wall.tec"))
    nums = sorted(glob.glob(os.path.join(d, "OUTPUT", "wall[0-9]*.tec")),
                  key=lambda p: int(re.search(r"wall(\d+)\.tec$", p).group(1)))
    if len(nums) < 3:
        raise SystemExit("verify.py: %s holds %d numbered wall outputs, at least 3 needed" % (d, len(nums)))
    last = [wall_quantities(p) for p in nums[-3:]] + [fin]
    eit = {k: float(max(q[k] for q in last) - min(q[k] for q in last)) for k in fin}
    t, aux = tec_time(os.path.join(d, "OUTPUT", "field.tec"))
    return dict(final=fin, e_it=eit, cycle=int(aux.get("CYCLE", -1)), numbered_outputs=len(nums),
                u_e=outer_velocity(os.path.join(d, "OUTPUT", "field.tec")),
                residual=residuals(os.path.join(d, "OUTPUT", "residual-history.dat")))


def three_levels(f1, f2, f3, e1, pf=P_FORMAL):
    e21, e32 = f2 - f1, f3 - f2
    res = dict(f1=f1, f2=f2, f3=f3, e21=e21, e32=e32)
    s = math.copysign(1.0, e32 * e21) if e21 != 0 and e32 != 0 else 0.0
    res["s"] = s
    small = abs(e21) < 10 * e1 or abs(e32) < 10 * e1
    if abs(e21) < abs(e32) and e21 != 0:
        p = math.log(abs(e32 / e21)) / math.log(2.0)
        res["p"] = p
        bad = s < 0 or p < 0.5 or p > 2 * pf or small
    else:
        res["p"] = "div"
        bad = True
    pu = pf if bad else res["p"]
    fs = 3.0 if bad else 1.25
    res.update(p_used=pu, Fs=fs)
    res["f_ext"] = (2 ** pu * f1 - f2) / (2 ** pu - 1)
    res["GCI"] = fs * abs(e21 / f1) / (2 ** pu - 1)
    res["U"] = res["GCI"] * abs(f1) + e1
    res["d_ok"] = e1 <= res["GCI"] * abs(f1) / 10
    res["e_it"] = e1
    return res


def order3(c):
    """Observed order of a code of the TMR on the three grids of the case."""
    f3, f2, f1 = c[N_TMR["35x25"]], c[N_TMR["69x49"]], c[N_TMR["137x97"]]
    return math.log(abs((f3 - f2) / (f2 - f1))) / math.log(2.0)


def analyse(model, runs):
    """Tables and criteria of the V&V page for one model; runs[grid] as run_series gives it."""
    lines = ["%s (%s in the TMR data)" % (model, TMR_NAME[model])]
    for g in GRIDS:
        r = runs[g]
        lines.append("  %-6s iterations %d, %d numbered outputs: " % (g, r["cycle"], r["numbered_outputs"]) +
                     ", ".join("%s %.6e (e_it %.1e)" % (k, r["final"][k], r["e_it"][k]) for k in KEYS))
    res, ok = {}, True
    for q in KEYS[:3]:
        e1 = max(runs[g]["e_it"][q] for g in GRIDS)
        r = res[q] = three_levels(runs["137x97"]["final"][q], runs["69x49"]["final"][q], runs["35x25"]["final"][q], e1)
        lines.append("  %s: f3 %.6e, f2 %.6e, f1 %.6e; s %+d, p %s, F_s %.2f, f_ext %.6e, GCI %.3f %%, U %.3e, e_it %.1e, (d) %s"
                     % (q, r["f3"], r["f2"], r["f1"], r["s"], ("%.3f" % r["p"]) if r["p"] != "div" else "div", r["Fs"],
                        r["f_ext"], 100 * r["GCI"], r["U"], e1, "met" if r["d_ok"] else "NOT met"))
        ok &= r["d_ok"]
    T = tmr_convergence(TMR_CONV[model])
    c, f, m = T["CFL3D"], T["FUN3D"], MARGIN[model]
    crit1 = {}
    for g in GRIDS:
        cf, n = runs[g]["final"][KEYS[0]], N_TMR[g]
        mid, half = 0.5 * (c[n] + f[n]), 0.5 * abs(c[n] - f[n])
        dev = abs(cf - mid)
        crit1[g] = dict(cf=cf, cfl3d=c[n], fun3d=f[n], dev=dev, band=half + m * cf, ok=dev <= half + m * cf)
        lines.append("  (1) %-6s c_f %.6e; CFL3D %.6e, FUN3D %.6e; |c_f - mean| %.4e %s band %.4e: %s"
                     % (g, cf, c[n], f[n], dev, "<=" if crit1[g]["ok"] else ">", half + m * cf,
                        "inside" if crit1[g]["ok"] else "outside"))
    ok &= crit1["137x97"]["ok"]
    r = res[KEYS[0]]
    dlim = abs(r["f_ext"] - c[N_545])
    bound = r["U"] + m * r["f1"]
    ok2 = dlim <= bound
    crit2 = dict(f_ext=r["f_ext"], cfl3d=c[N_545], fun3d=f[N_545], delta=dlim, bound=bound, ratio=dlim / bound, ok=ok2,
                 difference=(r["f_ext"] - c[N_545]) / c[N_545])
    lines.append("  (2) f_ext %.6e against CFL3D 545x385 %.6e (FUN3D %.6e): |D| %.4e %s U + m c_f %.4e: %s; "
                 "|D|/(U + m c_f) %.3f; f_ext %+.3f %% from CFL3D"
                 % (r["f_ext"], c[N_545], f[N_545], dlim, "<=" if ok2 else ">", bound, "met" if ok2 else "NOT met",
                    crit2["ratio"], 100 * crit2["difference"]))
    ok &= ok2
    orders = dict(CFL3D=order3(c), FUN3D=order3(f))
    lines.append("  orders of CFL3D and FUN3D on the same three grids: %.3f, %.3f" % (orders["CFL3D"], orders["FUN3D"]))
    ue_tmr = tmr_outer_velocity()
    lines.append("  u/U_inf at y = %g, cells nearest it: " % UE_Y + "; ".join(
        "x = %g: %s, CFL3D 545x385 %.6f (SA)" % (x, ", ".join("%s %.6f" % (g, runs[g]["u_e"][k]) for g in GRIDS), ue_tmr[k])
        for k, x in enumerate(UE_X)))
    lines.append("  density residual, fall from its largest value to the last iteration: " + ", ".join(
        "%s %.1f orders" % (g, runs[g]["residual"]["orders"]) for g in GRIDS) +
        "; 35x25, last 2000 iterations: %.2e to %.2e" % (runs["35x25"]["residual"]["last2000_min"],
                                                         runs["35x25"]["residual"]["last2000_max"]))
    return ok, lines, dict(three_levels=res, criterion1=crit1, criterion2=crit2, orders_tmr=orders, u_e_tmr=ue_tmr)


def figure(model, levels, analysis, walls):
    """c_f(0.97) against h with the limit (left) and c_f(x) (right), MOSE on the three grids and the TMR data."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not available: no figure)")
        return
    plt.rcParams.update({"figure.facecolor": "none", "axes.facecolor": "none", "savefig.facecolor": "none",
                         "svg.fonttype": "none", "svg.hashsalt": "flat-plate-tmr"})
    fs = 15
    blue, orange, green = "#0072B2", "#E69F00", "#009E73"
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 6))
    T = tmr_convergence(TMR_CONV[model])
    for code, col, mk in (("CFL3D", orange, "s"), ("FUN3D", green, "^")):
        n = sorted(T[code])
        a.plot([1e3 * nn ** -0.5 for nn in n], [1e3 * T[code][nn] for nn in n], color=col, marker=mk, ms=8,
               linewidth=2.5, linestyle="--", label=code + " (" + TMR_NAME[model] + ")")
    h = [1e3 * N_TMR[g] ** -0.5 for g in GRIDS]
    a.plot(h, [1e3 * levels[g]["final"][KEYS[0]] for g in GRIDS], color=blue, marker="o", ms=10, linewidth=3.0,
           label="MOSE, 35x25, 69x49, 137x97")
    r = analysis["three_levels"][KEYS[0]]
    a.errorbar([0.0], [1e3 * r["f_ext"]], yerr=[1e3 * r["U"]], color=blue, marker="D", ms=9, capsize=6, linewidth=2.0,
               label="MOSE, limit $\\pm U$ (p = %.2f)" % r["p_used"])
    a.set_xlabel(r"$h = N^{-1/2}\ \times 10^3$", fontsize=fs)
    a.set_ylabel(r"$C_f(x = 0.97) \times 10^3$", fontsize=fs)
    a.set_xlim(-1.5, 37.5)
    a.tick_params(labelsize=fs - 2)
    a.grid(True, alpha=0.3)
    a.legend(loc="best", fontsize=fs - 3)
    zones = tmr_zones(TMR_PLATE[model])
    for g, ls in zip(GRIDS, (":", "--", "-")):
        x, cf = walls[g]
        b.plot(x, 1e3 * cf, color=blue, linestyle=ls, linewidth=3.0, label="MOSE, %s" % g)
    for code, col in (("CFL3D", orange), ("FUN3D", green)):
        z = [k for k in zones if k.upper().startswith(code)][0]
        b.plot(zones[z][:, 0], 1e3 * zones[z][:, 1], color=col, linewidth=2.0, linestyle="-.", label=code + ", 545x385")
    b.set_xlim(0.0, 2.0)
    b.set_ylim(2.0, 6.0)
    b.set_xlabel(r"$x$ [m]", fontsize=fs)
    b.set_ylabel(r"$C_f \times 10^3$", fontsize=fs)
    b.tick_params(labelsize=fs - 2)
    b.grid(True, alpha=0.3)
    b.legend(loc="best", fontsize=fs - 3)
    fig.tight_layout()
    os.makedirs("RUN/figures", exist_ok=True)
    out = "RUN/figures/flat-plate-tmr-%s.svg" % model.lower()
    fig.savefig(out, metadata={"Date": None})     # the same data give the same file
    plt.close(fig)
    print("  figure: %s" % out)


def page(levels, walls_of):
    """Tables, criteria and figures of the V&V page from the values of the nine runs; True if the criteria hold."""
    ok, out = True, {}
    for model in MODELS:
        okm, lines, an = analyse(model, levels[model])
        print("\n".join(lines))
        print("  criteria of the page: %s" % ("met" if okm else "NOT met"))
        figure(model, levels[model], an, {g: walls_of(model, g) for g in GRIDS})
        ok &= okm
        out[model] = an
    os.makedirs("RUN", exist_ok=True)
    with open("RUN/analysis.json", "w") as fh:
        json.dump(out, fh, indent=1)
    return ok


# ---------------------------------------------------------------------------- forms
def short_form(n_expected):
    ok = True
    ref = json.load(open("reference/levels.json"))["runs"]
    print("Short form: restarts from the converged 35x25 solutions (relative tolerance %.0e)" % TOL_SHORT)
    for model in MODELS:
        d = "RUN/%s-35x25-short" % model
        need = [os.path.join(d, f) for f in ("input.ini", "OUTPUT/field.tec", "OUTPUT/wall.tec", "OUTPUT/residual-history.dat")]
        missing = [f for f in need if not os.path.exists(f)]
        if missing:
            print("  %s: FAIL - %s missing (run './MOSE.sh testfast')" % (model, ", ".join(missing)))
            ok = False
            continue
        n = int(re.search(r"(?m)^\s*iter-threshold\s*=\s*(\d+)", open(os.path.join(d, "input.ini")).read()).group(1))
        _, aux0 = tec_time("reference/field-35x25-%s.tec" % model)
        c0 = int(aux0["CYCLE"])
        tf, auxf = tec_time(os.path.join(d, "OUTPUT", "field.tec"))
        tw, _ = tec_time(os.path.join(d, "OUTPUT", "wall.tec"))
        rows = np.loadtxt(os.path.join(d, "OUTPUT", "residual-history.dat"), usecols=(0, 2), ndmin=2)
        counted = (tf == -(c0 + n) and tw == -(c0 + n) and int(auxf.get("CYCLE", -1)) == c0 + n and
                   rows.shape[0] == n and np.array_equal(rows[:, 0], np.arange(c0 + 1, c0 + n + 1)))
        q = wall_quantities(os.path.join(d, "OUTPUT", "wall.tec"))
        r0 = ref[model]["35x25"]["final"]
        dev = {k: abs(q[k] - r0[k]) / abs(r0[k]) for k in KEYS[:3]}
        good = counted and max(dev.values()) <= TOL_SHORT
        print("  %-10s %d iterations from %d (%s); " % (model, n, c0, "counted" if counted else
                                                        "NOT the iterations of the input: SOLUTIONTIME %s, CYCLE %s, %d rows"
                                                        % (tf, auxf.get("CYCLE"), rows.shape[0])) +
              ", ".join("%s %.6e (%.1e)" % (k, q[k], dev[k]) for k in KEYS[:3]) +
              "; density residual %.2e to %.2e -> %s" % (rows[0, 1], rows[-1, 1], "PASS" if good else "FAIL"))
        if n != n_expected:
            print("  (%s ran %d iterations; MOSE.sh testfast runs %d)" % (model, n, n_expected))
        ok &= good
    print("\nThe V&V page from the stored results of the full form (reference/):")
    levels, same = {}, True
    for model in MODELS:
        levels[model] = ref[model]
        for g in GRIDS:
            q = wall_quantities("reference/wall-%s-%s.tec" % (model, g))
            if any(abs(q[k] - ref[model][g]["final"][k]) > 1e-12 * abs(ref[model][g]["final"][k]) for k in q):
                print("  %s %s: FAIL - reference/wall-%s-%s.tec does not give the values of reference/levels.json"
                      % (model, g, model, g))
                same = False
    if same:
        print("  the nine wall outputs give the values of reference/levels.json")
    okp = page(levels, lambda m, g: wall_cf("reference/wall-%s-%s.tec" % (m, g)))
    return ok and same and okp


def full_form():
    levels = {}
    for model in MODELS:
        levels[model] = {}
        for g in GRIDS:
            d = "RUN/%s-%s" % (model, g)
            if not os.path.exists(os.path.join(d, "OUTPUT", "wall.tec")):
                print("FAIL - %s/OUTPUT/wall.tec missing (run './MOSE.sh test')" % d)
                return False
            levels[model][g] = run_series(d)
    print("Full form: the nine runs of RUN/")
    ok = page(levels, lambda m, g: wall_cf("RUN/%s-%s/OUTPUT/wall.tec" % (m, g)))
    os.makedirs("RUN/reference", exist_ok=True)
    with open("RUN/reference/levels.json", "w") as fh:
        json.dump(dict(quantities=list(KEYS), stations=list(STATIONS), q_inf=Q_INF, u_inf=U_INF, runs=levels), fh, indent=1)
    for model in MODELS:
        for g in GRIDS:
            shutil.copyfile("RUN/%s-%s/OUTPUT/wall.tec" % (model, g), "RUN/reference/wall-%s-%s.tec" % (model, g))
        shutil.copyfile("RUN/%s-35x25/OUTPUT/field.tec" % model, "RUN/reference/field-35x25-%s.tec" % model)
    print("  RUN/reference/: levels.json, wall-<model>-<grid>.tec, field-35x25-<model>.tec")
    return ok


form = open("RUN/form.txt").read().split() if os.path.exists("RUN/form.txt") else ["unknown"]
if form[0] == "testfast":
    passed = short_form(int(form[1]))
elif form[0] == "test":
    passed = full_form()
else:
    print("verify.py: RUN/form.txt says %r; run './MOSE.sh testfast' or './MOSE.sh test' first" % " ".join(form))
    passed = False
print("  STATUS : %s" % ("PASS" if passed else "FAIL"))
sys.exit(0 if passed else 1)
