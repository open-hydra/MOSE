"""
2D manufactured solutions on skewed and non-smooth grids: order of accuracy.

The MOSE_mms driver (src/app/mms.f90) adds the analytic source that makes a
smooth field an exact steady solution of the compressible Navier-Stokes
equations, so the discretisation error of each grid is known everywhere.
MOSE.sh runs, on N = 10, 20 (short form, testfast) or N = 20, 40, 80 (full
form, test) cells per side:
  zigzag, periodic2d  every cell skewed by +-0.3 h, its neighbours the other way;
  iface,  periodic2d  a cell-size jump of 1.5 and kinked grid lines at y = 1/2;
  iface,  wall2d      the same grid with walls: no slip at y = 0 (type 302,
                      T = T_w(x)) and y = 1 (type 301, the exact heat flux).
On the zigzag family the MUSCL reconstruction of the convective fluxes, taken
along the grid lines, is not exact for linear fields, so the observed order
drops to about 1; on the iface family it stays close to 2 in L2.

The study passes if, for every case:
  1. the normalised L2 error of rho, u, v and p decreases on the finest pair;
  2. the observed L2 order on the finest pair is at least ORDER_MIN of the
     family: 70 % of the order the family is expected to show, 1 on zigzag
     (the MUSCL reconstruction along the grid lines is not exact for linear
     fields where every cell is skewed the other way of its neighbours) and 2
     on iface;
  3. with walls, the heat flux integrated over the wall y = 0 converges to the
     exact one with an order of at least ORDER_MIN of the family;
  4. on every grid, the totals sum(U V) of the conserved variables whose
     boundary flux is zero (all of them with periodic2d, the mass with
     wall2d) differ from those of the initial state, the exact field at the
     cell centres, by less than the smallest L2 error of that grid.  The
     discrete steady problem of these equations fixes their solution only up
     to a mean (the discrete source sums to zero); local time steps move the
     totals during the iteration, and this check makes sure that the mean so
     chosen does not dominate the measured error.

Figures, transparent and theme-aware like the other V&V figures:
  OUTPUT/mms-skewed-grids.svg   the coarsest grid of each family;
  OUTPUT/mms-skewed-<case>.svg  L2 and Linf errors against N (zigzag, iface,
                                iface-walls), as the figure of ../verify.py.

Usage:
    python verify.py            # report, write the SVGs, exit 0/1
    python verify.py --plot     # also show the figures

Exit 0 = pass, 1 = fail.
"""
import argparse
import glob
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_case import exact, dTdy, RGAS, CP, KCOND  # noqa: E402

ORDER_MIN = {"zigzag": 0.7, "iface": 1.4}   # minimum observed L2 order on the finest pair
GAMMA = CP / (CP - RGAS)

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

parser = argparse.ArgumentParser(description="manufactured solutions on skewed 2D grids")
parser.add_argument("--plot", action="store_true", help="show the figures")
args = parser.parse_args()

VARS = ["rho", "u", "v", "p"]
LABELS = {"rho": r"$\rho$", "u": r"$u$", "v": r"$v$", "p": r"$p$"}
CASES = [("zigzag", "periodic2d"), ("iface", "periodic2d"), ("iface", "wall2d")]
TITLES = {("zigzag", "periodic2d"): "zig-zag, periodic", ("iface", "periodic2d"): "kinked interface, periodic",
          ("iface", "wall2d"): "kinked interface, walls"}


def centres(x, y, z):
    """Cell centres (average of the 8 corner nodes) of one block."""
    def c(a):
        return 0.125 * (a[:-1, :-1, :-1] + a[1:, :-1, :-1] + a[:-1, 1:, :-1] + a[1:, 1:, :-1] +
                        a[:-1, :-1, 1:] + a[1:, :-1, 1:] + a[:-1, 1:, 1:] + a[1:, 1:, 1:])
    return c(x), c(y), c(z)


def volumes(x, y, z):
    """Volume of the cells of a quasi-2D block: area of the quadrilateral in the plane k = 0 (shoelace) times the depth."""
    xa, ya = x[:, :, 0], y[:, :, 0]
    x1, y1, x2, y2 = xa[:-1, :-1], ya[:-1, :-1], xa[1:, :-1], ya[1:, :-1]
    x3, y3, x4, y4 = xa[1:, 1:], ya[1:, 1:], xa[:-1, 1:], ya[:-1, 1:]
    area = 0.5 * np.abs((x1 * y2 - x2 * y1) + (x2 * y3 - x3 * y2) + (x3 * y4 - x4 * y3) + (x4 * y1 - x1 * y4))
    depth = np.mean(z[:, :, 1] - z[:, :, 0])
    return (area * depth)[:, :, None]


def conserved(rho, u, v, w, p):
    E = p / ((GAMMA - 1.0) * rho) + 0.5 * (u ** 2 + v ** 2 + w ** 2)
    return [rho, rho * u, rho * v, rho * E]


def read_case(path, field):
    """Errors, totals and grid of one field_<fam>-<field>-<N>.tec."""
    xb, yb, zb, vb, names = read_TEC(path)
    idx = {n: k - 3 for k, n in enumerate(names)}
    num = {v: [] for v in VARS}
    ex = {v: [] for v in VARS}
    tot_num = np.zeros(4)
    tot_ex = np.zeros(4)
    tot_abs = np.zeros(4)
    for x, y, z, var in zip(xb, yb, zb, vb):
        xc, yc, zc = centres(x, y, z)
        V = volumes(x, y, z)
        rho_e, u_e, v_e, w_e, p_e, _ = exact(field, xc, yc, zc)
        n = {"rho": var[idx["rho(1)"]], "u": var[idx["u"]], "v": var[idx["v"]], "p": var[idx["p"]]}
        e = {"rho": rho_e, "u": u_e, "v": v_e, "p": p_e}
        for q in VARS:
            num[q].append(np.asarray(n[q]).ravel())
            ex[q].append(np.broadcast_to(e[q], xc.shape).ravel())
        Un = conserved(n["rho"], n["u"], n["v"], var[idx["w"]], n["p"])
        Ue = conserved(rho_e, u_e, v_e, np.broadcast_to(w_e, xc.shape), p_e)
        for k in range(4):
            tot_num[k] += float(np.sum(Un[k] * V))
            tot_ex[k] += float(np.sum(Ue[k] * V))
            tot_abs[k] += float(np.sum(np.abs(Ue[k]) * V))
    if not all(np.all(np.isfinite(np.concatenate(num[q]))) for q in VARS):
        return None
    out = {}
    for q in VARS:
        a, b = np.concatenate(num[q]), np.concatenate(ex[q])
        out[q] = (math.sqrt(float(np.mean((a - b) ** 2))) / math.sqrt(float(np.mean(b ** 2))),
                  float(np.max(np.abs(a - b))) / float(np.max(np.abs(b))))
    out["totals"] = np.abs(tot_num - tot_ex) / tot_abs
    out["grid"] = [(x[:, :, 0], y[:, :, 0]) for x, y in zip(xb, yb)]
    out["cells"] = int(sum((x.shape[0] - 1) * (x.shape[1] - 1) for x in xb))
    return out


def zone_titles(path):
    return re.findall(r'ZONE\s+T\s*=\s*"?(\w+)"?', open(path).read())


def wall_flux(path, field):
    """Heat flux integrated over the wall y = 0, numerical and exact (|sum q A|; the sign of q follows the face)."""
    xb, yb, zb, vb, names = read_TEC(path)
    iq = names.index("qw") - 3
    num = ex = 0.0
    for x, y, z, var in zip(xb, yb, zb, vb):
        P = np.stack([x, y, z], axis=-1)
        flat = [d for d in range(3) if P.shape[d] == 1][0]
        Q = np.moveaxis(P, flat, 0)[0]
        c = 0.25 * (Q[:-1, :-1] + Q[1:, :-1] + Q[:-1, 1:] + Q[1:, 1:])
        if float(np.mean(c[..., 1])) > 0.5:
            continue
        a1 = 0.5 * (Q[1:, :-1] - Q[:-1, :-1] + Q[1:, 1:] - Q[:-1, 1:])
        a2 = 0.5 * (Q[:-1, 1:] - Q[:-1, :-1] + Q[1:, 1:] - Q[1:, :-1])
        A = np.linalg.norm(np.cross(a1, a2), axis=-1)
        q = np.asarray(var[iq]).reshape(A.shape)
        num += float(np.sum(q * A))
        ex += float(np.sum(KCOND * dTdy(field, c[..., 0], c[..., 1], c[..., 2]) * A))
    return abs(num), abs(ex)


form = open("OUTPUT/form.txt").read().strip() if os.path.exists("OUTPUT/form.txt") else "unknown"
form_name = {"testfast": "short form", "test": "full form"}.get(form, form)
results = {}
for fam, field in CASES:
    files = glob.glob("OUTPUT/field_%s-%s-*.tec" % (fam, field))
    Ns = sorted(int(re.search(r"-(\d+)\.tec$", f).group(1)) for f in files)
    if len(Ns) < 2:
        print("mms skewed: FAIL - %s %s: fewer than two grids (run './MOSE.sh testfast' or './MOSE.sh test')" % (fam, field))
        sys.exit(1)
    rows = []
    for N in Ns:
        r = read_case("OUTPUT/field_%s-%s-%d.tec" % (fam, field, N), field)
        if r is None:
            print("mms skewed: FAIL - non-finite field, %s %s N=%d" % (fam, field, N))
            sys.exit(1)
        if field.startswith("wall"):
            r["qw"] = wall_flux("OUTPUT/wall_%s-%s-%d.tec" % (fam, field, N), field)
        rows.append((N, r))
    results[(fam, field)] = rows

print("2D manufactured solutions on skewed grids - %s (normalised errors)" % form_name)
ok_all = True
for (fam, field), rows in results.items():
    Ns = np.array([N for N, _ in rows], float)
    cells = [r["cells"] for _, r in rows]
    print("  %s (%s): N = %s, cells = %s" % (TITLES[(fam, field)], field, ", ".join("%d" % N for N in Ns),
                                              ", ".join(str(c) for c in cells)))
    ok_case = True
    for q in VARS:
        e2 = np.array([r[q][0] for _, r in rows])
        ei = np.array([r[q][1] for _, r in rows])
        o2 = np.log(e2[:-1] / e2[1:]) / np.log(Ns[1:] / Ns[:-1])
        oi = np.log(ei[:-1] / ei[1:]) / np.log(Ns[1:] / Ns[:-1])
        print("    %-3s L2   : %s   orders %s" % (q, "  ".join("%.3e" % v for v in e2), " ".join("%.2f" % v for v in o2)))
        print("        Linf : %s   orders %s" % ("  ".join("%.3e" % v for v in ei), " ".join("%.2f" % v for v in oi)))
        ok_case &= bool(e2[-1] < e2[-2]) and bool(o2[-1] >= ORDER_MIN[fam])
    if field.startswith("wall"):
        eq = np.array([abs(a - b) / b for _, r in rows for a, b in [r["qw"]]])
        oq = np.log(eq[:-1] / eq[1:]) / np.log(Ns[1:] / Ns[:-1])
        print("    wall heat flux at y = 0, relative error: %s   orders %s" % ("  ".join("%.3e" % v for v in eq),
                                                                             " ".join("%.2f" % v for v in oq)))
        ok_case &= bool(oq[-1] >= ORDER_MIN[fam])
    keep = [0, 1, 2, 3] if field == "periodic2d" else [0]
    drift = [max(float(r["totals"][k]) for k in keep) for _, r in rows]
    floor = [min(r[q][0] for q in VARS) for _, r in rows]
    print("    totals with zero boundary flux, largest relative change from the initial state: %s   (smallest L2 error: %s)"
          % ("  ".join("%.1e" % d for d in drift), "  ".join("%.1e" % e for e in floor)))
    ok_case &= all(d < e for d, e in zip(drift, floor))
    print("    %s (L2 order on the finest pair >= %.1f for every variable)" % ("PASS" if ok_case else "FAIL", ORDER_MIN[fam]))
    ok_all &= ok_case
status = "PASS" if ok_all else "FAIL"
print("  STATUS : %s" % status)

# --- figures -------------------------------------------------------------------------------
import matplotlib as mpl  # noqa: E402
if not args.plot:
    mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})
os.makedirs("OUTPUT", exist_ok=True)

# coarsest grid of each family
fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
for ax, (fam, field) in zip(axes, [("zigzag", "periodic2d"), ("iface", "periodic2d")]):
    N, r = results[(fam, field)][0]
    for b, (X, Y) in enumerate(r["grid"]):
        col = "C0" if b == 0 else "C1"
        ax.plot(X, Y, "-", color=col, lw=1.0)
        ax.plot(X.T, Y.T, "-", color=col, lw=1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("x", fontsize=20)
    ax.set_ylabel("y", fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_title("%s, %d x %d cells" % (TITLES[(fam, field)].split(",")[0], N, N), fontsize=18)
fig.legend([plt.Line2D([], [], color="C0"), plt.Line2D([], [], color="C1")], ["block 1", "block 2"],
           loc="upper center", ncol=2, fontsize=18, bbox_to_anchor=(0.5, 1.04))
plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.savefig(os.path.join("OUTPUT", "mms-skewed-grids.svg"), bbox_inches="tight", transparent=True)

# L2 and Linf errors against N, one figure per case, as the figure of the Cartesian case (../verify.py)
from matplotlib.ticker import NullFormatter  # noqa: E402
SLUG = {("zigzag", "periodic2d"): "zigzag", ("iface", "periodic2d"): "iface", ("iface", "wall2d"): "iface-walls"}
for (fam, field), rows in results.items():
    Ns = np.array([N for N, _ in rows], float)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, norm, col in ((axes[0], "L_2", 0), (axes[1], "L_\\infty", 1)):
        data = {q: np.array([r[q][col] for _, r in rows]) for q in VARS}
        for q in VARS:
            ax.loglog(Ns, data[q], "o-", lw=2.5, ms=9, label=LABELS[q])
        ref = data["rho"][0] * (Ns[0] / Ns) ** 2
        ax.loglog(Ns, ref, "k--", lw=2.0, alpha=0.6, label="order 2")
        if fam == "zigzag":
            ax.loglog(Ns, data["rho"][0] * (Ns[0] / Ns), "k:", lw=2.0, alpha=0.6, label="order 1")
        ax.set_xlabel("N (cells per side)", fontsize=20)
        ax.set_ylabel(rf"normalised ${norm}$ error", fontsize=20)
        ax.grid(True, which="both", alpha=0.3)
        ax.tick_params(labelsize=20)
        ax.set_xticks(Ns)
        ax.set_xticklabels(["%d" % N for N in Ns])
        ax.xaxis.set_minor_formatter(NullFormatter())
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels), fontsize=20, bbox_to_anchor=(0.5, 1.05))
    fig.text(0.99, 0.0, "%s, %s; %s: %s cells" % (TITLES[(fam, field)], field, form_name,
                                                 ", ".join("%d x %d" % (N, N) for N in Ns)),
             ha="right", va="bottom", fontsize=14)
    plt.tight_layout(rect=[0, 0.04, 1, 0.93])
    plt.savefig(os.path.join("OUTPUT", "mms-skewed-%s.svg" % SLUG[(fam, field)]), bbox_inches="tight", transparent=True)

if args.plot:
    plt.show()

sys.exit(0 if status == "PASS" else 1)
