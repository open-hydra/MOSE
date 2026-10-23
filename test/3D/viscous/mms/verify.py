"""
3D manufactured solutions on skewed and non-smooth grids: order of accuracy.

The MOSE_mms driver (src/app/mms.f90) adds the analytic source that makes a
smooth field an exact steady solution of the compressible Navier-Stokes
equations, so the discretisation error of each grid is known everywhere.
MOSE.sh runs (short form, testfast: zigzag on N = 8, 12; full form, test:
zigzag on N = 16, 24, 36, iface with and without walls on N = 20, 30, 45):
  zigzag, periodic3d  every cell skewed by +-0.3 h in the three directions,
                      its neighbours the other way;
  iface,  periodic3d  a cell-size jump of 1.5 and kinked grid lines at y = 1/2;
  iface,  wall3d      the same grid with walls: no slip at y = 0 (type 302,
                      T = T_w(x, z)) and y = 1 (type 301, the exact heat flux).
The two blocks of each grid meet with permuted and reversed indices.  The
short form leaves out the iface grids: on N = 10, 15 the grid with walls
passes the checks below, but its convergence costs twenty to twenty-five
times the zigzag pair in CPU time; the full form runs it on N = 20, 30, 45.
On the zigzag family the MUSCL reconstruction of the convective fluxes, taken
along the grid lines, is not exact for linear fields, so the observed order
drops to about 1; on the iface family it stays close to 2 in L2.

The study passes if, for every case:
  1. the normalised L2 error of rho, u, v, w and p decreases on the finest
     pair;
  2. the observed L2 order on the finest pair is at least ORDER_MIN of the
     family: 70 % of the order the family is expected to show, 1 on zigzag
     (the MUSCL reconstruction along the grid lines is not exact for linear
     fields where every cell is skewed the other way of its neighbours) and 2
     on iface;
  3. with walls, the heat flux integrated over the wall y = 0 converges to the
     exact one with an order of at least ORDER_MIN of the family;
  4. on every grid, the totals sum(U V) of the conserved variables whose
     boundary flux is zero (all of them with periodic3d, the mass with
     wall3d) differ from those of the initial state, the exact field at the
     cell centres, by less than the smallest L2 error of that grid.  The
     discrete steady problem of these equations fixes their solution only up
     to a mean (the discrete source sums to zero); local time steps move the
     totals during the iteration, and this check makes sure that the mean so
     chosen does not dominate the measured error.

Figures, transparent and theme-aware like the other V&V figures:
  OUTPUT/mms3d-grids.svg        the bottom face of the coarsest grid of each
                                family, seen along z;
  OUTPUT/mms3d-<case>.svg       L2 and Linf errors against N (zigzag, iface,
                                iface-walls), as the figure of the 2D case.

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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "2D" / "viscous" / "mms" / "skewed"))
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

parser = argparse.ArgumentParser(description="manufactured solutions on skewed 3D grids")
parser.add_argument("--plot", action="store_true", help="show the figures")
args = parser.parse_args()

VARS = ["rho", "u", "v", "w", "p"]
LABELS = {"rho": r"$\rho$", "u": r"$u$", "v": r"$v$", "w": r"$w$", "p": r"$p$"}
CASES = [("zigzag", "periodic3d"), ("iface", "periodic3d"), ("iface", "wall3d")]
TITLES = {("zigzag", "periodic3d"): "zig-zag, periodic", ("iface", "periodic3d"): "kinked interface, periodic",
          ("iface", "wall3d"): "kinked interface, walls"}


def centres(x, y, z):
    """Cell centres (average of the 8 corner nodes) of one block."""
    def c(a):
        return 0.125 * (a[:-1, :-1, :-1] + a[1:, :-1, :-1] + a[:-1, 1:, :-1] + a[1:, 1:, :-1] +
                        a[:-1, :-1, 1:] + a[1:, :-1, 1:] + a[:-1, 1:, 1:] + a[1:, 1:, 1:])
    return c(x), c(y), c(z)


def volumes(x, y, z):
    """Cell volumes as the solver takes them (Compute_Metric_Tensor): the integral of det J of the trilinear cell,
    1/6 sum_d e_d . (S_d- + S_d+), e_d the mean edge of direction d and S_d-, S_d+ the area vectors (diagonal x diagonal)/2
    of the two faces of direction d, along +d (docs/theory/numerics.md, Cell volume); in absolute value."""
    X = np.stack([x, y, z], axis=-1)
    n1, n2, n3, n4 = X[:-1, :-1, :-1], X[:-1, :-1, 1:], X[:-1, 1:, :-1], X[:-1, 1:, 1:]
    n5, n6, n7, n8 = X[1:, :-1, :-1], X[1:, :-1, 1:], X[1:, 1:, :-1], X[1:, 1:, 1:]
    r2, r3, r4, r5, r6, r7, r8 = (n - n1 for n in (n2, n3, n4, n5, n6, n7, n8))

    def dot(a, b):
        return np.einsum('...i,...i->...', a, b)
    e1 = 0.25 * (r5 + (r6 - r2) + (r7 - r3) + (r8 - r4))
    e2 = 0.25 * (r3 + (r4 - r2) + (r7 - r5) + (r8 - r6))
    e3 = 0.25 * (r2 + (r4 - r3) + (r6 - r5) + (r8 - r7))
    v = (dot(e1, np.cross(r4, r2 - r3) + np.cross(r8 - r5, r6 - r7))
         + dot(e2, np.cross(r6, r5 - r2) + np.cross(r8 - r3, r7 - r4))
         + dot(e3, np.cross(r7, r3 - r5) + np.cross(r8 - r2, r4 - r6)))
    return np.abs(v) / 12.0


def conserved(rho, u, v, w, p):
    E = p / ((GAMMA - 1.0) * rho) + 0.5 * (u ** 2 + v ** 2 + w ** 2)
    return [rho, rho * u, rho * v, rho * w, rho * E]


def read_case(path, field):
    """Errors, totals and grid of one field_<fam>-<field>-<N>.tec."""
    xb, yb, zb, vb, names = read_TEC(path)
    idx = {n: k - 3 for k, n in enumerate(names)}
    num = {v: [] for v in VARS}
    ex = {v: [] for v in VARS}
    tot_num = np.zeros(5)
    tot_ex = np.zeros(5)
    tot_abs = np.zeros(5)
    for x, y, z, var in zip(xb, yb, zb, vb):
        xc, yc, zc = centres(x, y, z)
        V = volumes(x, y, z)
        rho_e, u_e, v_e, w_e, p_e, _ = exact(field, xc, yc, zc)
        n = {"rho": var[idx["rho(1)"]], "u": var[idx["u"]], "v": var[idx["v"]], "w": var[idx["w"]], "p": var[idx["p"]]}
        e = {"rho": rho_e, "u": u_e, "v": v_e, "w": w_e, "p": p_e}
        for q in VARS:
            num[q].append(np.asarray(n[q]).ravel())
            ex[q].append(np.broadcast_to(e[q], xc.shape).ravel())
        Un = conserved(n["rho"], n["u"], n["v"], n["w"], n["p"])
        Ue = conserved(rho_e, u_e, v_e, np.broadcast_to(w_e, xc.shape), p_e)
        for k in range(5):
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
    out["grid"] = bottom_faces(xb, yb, zb, path)
    out["cells"] = int(sum((x.shape[0] - 1) * (x.shape[1] - 1) * (x.shape[2] - 1) for x in xb))
    return out


def bottom_faces(xb, yb, zb, path):
    """(x, y) of the nodes of the bottom face of every block, the grid seen along z.  The bottom face of a block is the
    face that the period of the grids in z, (0, 0, 1), moves onto a face of a block: the faces that bc.txt pairs with the
    periodic records in z.  On zigzag it lies on the plane z = 0; on iface it is not plane, as the kink moves z by half the
    shift of x.  The blocks are turned, so it is not always k = 0.  A block without one such face stops the script."""
    faces = []
    for x, y, z in zip(xb, yb, zb):
        X = np.stack([x, y, z], axis=-1)
        faces.append((X[0], X[-1], X[:, 0], X[:, -1], X[:, :, 0], X[:, :, -1]))
    every = [G for six in faces for G in six]

    def same_nodes(A, B):
        """Whether the face B has the nodes of the face A, in one of the eight orders of the two index directions."""
        for C in (B, np.swapaxes(B, 0, 1)):
            for D in (C, C[::-1], C[:, ::-1], C[::-1, ::-1]):
                if D.shape == A.shape and float(np.max(np.abs(D - A))) < 1e-9:
                    return True
        return False

    out = []
    for b, six in enumerate(faces):
        found = [F for F in six if any(same_nodes(F + np.array([0.0, 0.0, 1.0]), G) for G in every)]
        if len(found) != 1:
            print("mms 3D: FAIL - %s, block %d: %d faces that the period (0, 0, 1) moves onto a face of a block, not one: "
                  "no bottom face to draw" % (path, b + 1, len(found)))
            sys.exit(1)
        out.append((found[0][..., 0], found[0][..., 1]))
    return out


def zone_titles(path):
    return re.findall(r'ZONE\s+T\s*=\s*"?(\w+)"?', open(path).read())


def wall_flux(path, field):
    """Heat flux integrated over the wall y = 0, numerical and exact, compared in absolute value: qw of the wall output
    is the heat flux into the gas on every face, -k dT/dy on this wall, and the exact flux is integrated as k dT/dy."""
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
    if not Ns:
        continue                      # the short form runs the zigzag family only
    if len(Ns) < 2:
        print("mms 3D: FAIL - %s %s: fewer than two grids (run './MOSE.sh testfast' or './MOSE.sh test')" % (fam, field))
        sys.exit(1)
    rows = []
    for N in Ns:
        r = read_case("OUTPUT/field_%s-%s-%d.tec" % (fam, field, N), field)
        if r is None:
            print("mms 3D: FAIL - non-finite field, %s %s N=%d" % (fam, field, N))
            sys.exit(1)
        if field.startswith("wall"):
            r["qw"] = wall_flux("OUTPUT/wall_%s-%s-%d.tec" % (fam, field, N), field)
        rows.append((N, r))
    results[(fam, field)] = rows

if not results:
    print("mms 3D: FAIL - no field_<case>.tec in OUTPUT (run './MOSE.sh testfast' or './MOSE.sh test')")
    sys.exit(1)
print("3D manufactured solutions on skewed grids - %s (normalised errors)" % form_name)
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
    keep = [0, 1, 2, 3, 4] if field == "periodic3d" else [0]
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

# bottom face of the coarsest grid of each family, seen along z
shown = [c for c in [("zigzag", "periodic3d"), ("iface", "periodic3d")] if c in results]
fig, axes = plt.subplots(1, len(shown), figsize=(5.5 * len(shown), 5.5), squeeze=False)
for ax, (fam, field) in zip(axes[0], shown):
    N, r = results[(fam, field)][0]
    for b, (X, Y) in enumerate(r["grid"]):
        col = "C0" if b == 0 else "C1"
        ax.plot(X, Y, "-", color=col, lw=1.0)
        ax.plot(X.T, Y.T, "-", color=col, lw=1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("x", fontsize=20)
    ax.set_ylabel("y", fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_title("%s, %d x %d x %d cells, bottom face" % (TITLES[(fam, field)].split(",")[0], N, N, N), fontsize=16)
fig.legend([plt.Line2D([], [], color="C0"), plt.Line2D([], [], color="C1")], ["block 1", "block 2"],
           loc="upper center", ncol=2, fontsize=18, bbox_to_anchor=(0.5, 1.04))
plt.tight_layout(rect=[0, 0, 1, 0.92])
plt.savefig(os.path.join("OUTPUT", "mms3d-grids.svg"), bbox_inches="tight", transparent=True)

# L2 and Linf errors against N, one figure per case, as the figure of the Cartesian case (../verify.py)
from matplotlib.ticker import NullFormatter  # noqa: E402
SLUG = {("zigzag", "periodic3d"): "zigzag", ("iface", "periodic3d"): "iface", ("iface", "wall3d"): "iface-walls"}
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
                                                 ", ".join("%d x %d x %d" % (N, N, N) for N in Ns)),
             ha="right", va="bottom", fontsize=14)
    plt.tight_layout(rect=[0, 0.04, 1, 0.93])
    plt.savefig(os.path.join("OUTPUT", "mms3d-%s.svg" % SLUG[(fam, field)]), bbox_inches="tight", transparent=True)

if args.plot:
    plt.show()

sys.exit(0 if status == "PASS" else 1)
