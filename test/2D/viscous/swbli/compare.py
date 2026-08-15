"""
SWBLI comparison figures for the V&V documentation. MOSE and OpenFOAM
(rhoCentralFoam) are run on the same 237,440-cell grid with the same freestream,
thermodynamics and Sutherland transport. Figures are written straight into
docs/vv/images as transparent, theme-aware SVG.

Per model (--model SA|SST):
    SWBLI-cf-<model>-openfoam.svg     verification: MOSE vs OpenFOAM, skin friction
    SWBLI-pw-<model>-openfoam.svg     verification: MOSE vs OpenFOAM, wall pressure
    SWBLI-cf-<model>-validation.svg   validation: MOSE, OpenFOAM, Wind-US, experiment
    SWBLI-cf-sst-omegaBC-mose.svg     (SST) MOSE, practical vs asymptotic omega BC

Inputs (per <MODEL>/):
    OUTPUT/wall.tec                        MOSE final (Sutherland) wall
    OPENFOAM/bottomWall.xy                 OpenFOAM final (Sutherland) wall
    reference/{wind,SU2}-<MODEL>.dat, reference/schulein.dat   reference data
    reference/MOSE-SST{,-asymptotic}.tec   MOSE constant-mu both-BC runs (omega-BC fig)

The OpenFOAM wall is the raw bottom-wall patch sample:
    postProcess -func wallShearStress -latestTime
    postProcess -func "patchSurface(patch=bottomWall, fields=(wallShearStress p), \
        interpolate=false, surfaceFormat=raw)" -latestTime
    cp postProcessing/patchSurface*/*/patch.xy  OPENFOAM/bottomWall.xy

Usage:
    python compare.py --model SA
    python compare.py --model SST [--plot]
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# Same theme-aware SVG setup as the verify.py scripts: transparent canvas and
# real <text> elements, so the docs CSS can recolor them for the page theme.
mpl.rcParams.update({
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "svg.fonttype": "none",
})

root = Path(__file__).resolve().parent
for parent in [root, *root.parents]:
    candidate = parent / "lib" / "ORION" / "src" / "python"
    if candidate.exists():
        sys.path.insert(0, str(candidate))
        break
from ORION import read_TEC   # noqa: E402

# Figures are the V&V documentation figures: write them straight into the docs.
OUT = root.parents[3] / "docs" / "vv" / "images"
LABEL_FS = 17

# Freestream (<MODEL>/input.ini): M=5, p=4000 Pa, T=68.3 K, air
R, GAMMA = 287.0, 1.4
T_INF, P_INF, MACH = 68.3, 4000.0, 5.0
RHO_INF = P_INF / (R * T_INF)
U_INF = MACH * np.sqrt(GAMMA * R * T_INF)
Q_INF = 0.5 * RHO_INF * U_INF ** 2

# Categorical trio (validated CVD-safe in light & dark, dataviz palette slots 1-3):
# MOSE = orange (hero), OpenFOAM = aqua, Wind-US = blue.
MOSE_C, OF_C, WIND_C = "#eb6834", "#1baf7a", "#2a78d6"
OF_D_C = "#4a3aa7"   # violet: OpenFOAM asymptotic in the legacy omega-BC figure

# Interaction region shown in the iso-line figure
XLIM, YLIM = (0.15, 0.45), (0.0, 0.09)

# For verification the two codes must solve the *same* model, so the OpenFOAM SST
# run uses beta1 = 0.0075 in omegaWallFunction, i.e. omega_w = 800*nu/y^2 — the
# Menter practical form MOSE applies (Lib_SST.f90:260).  bottomWall-defaultBC.xy
# is the same case with OpenFOAM's default (asymptotic, 80*nu/y^2) condition; it
# is a validation variant and is only used in the omega-BC study figure.
OF_LABEL = {"SA": "OpenFOAM (rhoCentralFoam)",
            "SST": r"OpenFOAM (Menter practical $\omega$ BC)"}


# ----------------------------------------------------------------------
# OpenFOAM readers (binary case files, no OpenFOAM installation needed)
# ----------------------------------------------------------------------
def _next_list(buf, pos, itemsize, ncomp=1, dtype="<f8"):
    """Read the next '<count>\\n(<binary blob>)' list starting at/after pos."""
    m = re.compile(rb"\n(\d+)\s*\n\(").search(buf, pos)
    if m is None:
        raise ValueError("no binary list found")
    n = int(m.group(1))
    start = m.end()
    nbytes = n * ncomp * itemsize
    a = np.frombuffer(buf, dtype=dtype, count=n * ncomp, offset=start)
    end = start + nbytes
    return (a.reshape(n, ncomp) if ncomp > 1 else a), end


def of_cell_centres(case):
    """Cell centres of an OpenFOAM polyMesh, from points/faces/owner/neighbour.

    For a hex, every vertex is shared by 3 of its faces, so summing the points
    of all the faces of a cell and dividing by the count gives the vertex
    centroid exactly.
    """
    pm = case / "constant" / "polyMesh"

    pts, _ = _next_list((pm / "points").read_bytes(), 0, 8, ncomp=3)

    buf = (pm / "faces").read_bytes()
    off, pos = _next_list(buf, 0, 4, dtype="<i4")          # compact-list offsets
    fpts, _ = _next_list(buf, pos, 4, dtype="<i4")         # face point indices

    owner, _ = _next_list((pm / "owner").read_bytes(), 0, 4, dtype="<i4")
    neigh, _ = _next_list((pm / "neighbour").read_bytes(), 0, 4, dtype="<i4")

    face_sum = np.add.reduceat(pts[fpts], off[:-1], axis=0)
    face_cnt = np.diff(off).astype(float)

    ncells = int(max(owner.max(), neigh.max())) + 1
    csum = np.zeros((ncells, 3))
    ccnt = np.zeros(ncells)
    np.add.at(csum, owner, face_sum)
    np.add.at(ccnt, owner, face_cnt)
    np.add.at(csum, neigh, face_sum[: neigh.size])
    np.add.at(ccnt, neigh, face_cnt[: neigh.size])
    return csum / ccnt[:, None]


def of_field(path):
    """internalField of a binary volScalarField."""
    buf = path.read_bytes()
    m = re.compile(rb"internalField\s+nonuniform\s+List<scalar>").search(buf)
    if m is None:
        m = re.compile(rb"internalField\s+uniform\s+([-\d.eE+]+)").search(buf)
        raise ValueError(f"{path.name}: uniform field ({float(m.group(1))})")
    a, _ = _next_list(buf, m.end(), 8)
    return a


def of_latest_time(case):
    times = [d for d in case.iterdir()
             if d.is_dir() and re.fullmatch(r"\d+", d.name) and d.name != "0"]
    if not times:
        raise FileNotFoundError(f"no time directory in {case} — is the run copied over?")
    return max(times, key=lambda d: int(d.name))


def of_wall(path):
    """(x, Cf, p) on the bottom wall from a raw patchSurface sample."""
    d = np.loadtxt(path)
    o = np.argsort(d[:, 0])
    # wallShearStress is the traction on the fluid, hence negative under an
    # attached BL: flip the sign to match the MOSE tauX convention.
    return d[o, 0], -d[o, 3] / Q_INF, d[o, 6]


def ref_xy(path):
    """(x, Cf) from a literature reference file, or (None, None) if absent."""
    if not path.exists():
        return None, None
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1]


# ----------------------------------------------------------------------
# MOSE readers
# ----------------------------------------------------------------------
def mose_wall(path):
    """(x, Cf, p) from wall.tec. Vars are ['y+','tauX','tauY','tauZ','pw',...]."""
    xw, _, _, vw, _ = read_TEC(str(path))
    xn = xw[0][:, 0, 0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    return xc, vw[0][1][:, 0, 0] / Q_INF, vw[0][4][:, 0, 0]


def mose_mach(path):
    """[(x, y, Mach)] per block from field.tec.

    Coordinates are nodal, variables are cell-centred: average the four nodes
    of each cell so the contour grid matches the data.
    """
    xb, yb, _, vb, names = read_TEC(str(path))
    iu, iv, iw = (names.index(n) - 3 for n in ("u", "v", "w"))
    iT, ig, iR = (names.index(n) - 3 for n in ("T", "g", "R"))
    out = []
    for x, y, v in zip(xb, yb, vb):
        u, w2, w3 = v[iu][:, :, 0], v[iv][:, :, 0], v[iw][:, :, 0]
        T, g, Rg = v[iT][:, :, 0], v[ig][:, :, 0], v[iR][:, :, 0]
        M = np.sqrt(u ** 2 + w2 ** 2 + w3 ** 2) / np.sqrt(g * Rg * T)
        xn, yn = x[:, :, 0], y[:, :, 0]
        xc = 0.25 * (xn[:-1, :-1] + xn[1:, :-1] + xn[:-1, 1:] + xn[1:, 1:])
        yc = 0.25 * (yn[:-1, :-1] + yn[1:, :-1] + yn[:-1, 1:] + yn[1:, 1:])
        out.append((xc, yc, M))
    return out


# ----------------------------------------------------------------------
def sep_reatt(x, cf, lo=0.30, hi=0.42):
    m = (x > lo) & (x < hi)
    xs, cs = np.asarray(x)[m], np.asarray(cf)[m]
    n = np.where(cs < 0)[0]
    return (xs[n[0]], xs[n[-1]]) if n.size else (None, None)


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, fontsize=LABEL_FS)
    ax.set_ylabel(ylabel, fontsize=LABEL_FS)
    ax.tick_params(labelsize=LABEL_FS - 2)
    ax.legend(loc="best", fontsize=LABEL_FS - 2)
    ax.grid(True, alpha=0.3)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / name
    fig.savefig(p, dpi=150, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote docs/vv/images/{name}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=("SA", "SST"), default="SA")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    M = args.model
    m = M.lower()
    case = root / M
    ofc = case / "OPENFOAM"

    print(f"Model: {M}")
    print(f"Freestream: U={U_INF:.2f} m/s  rho={RHO_INF:.4f} kg/m3  q={Q_INF:.0f} Pa")

    # Final MOSE and OpenFOAM results (Sutherland), plus the Wind-US reference.
    xm, cfm, pwm = mose_wall(case / "OUTPUT" / "wall.tec")
    xo, cfo, pwo = of_wall(ofc / "bottomWall.xy")
    xw, cfw = ref_xy(case / "reference" / f"wind-{M}.dat")

    print("\nSeparation / reattachment (Cf = 0):")
    three = {"MOSE": (xm, cfm), "OpenFOAM": (xo, cfo)}
    if xw is not None:
        three["Wind-US"] = (xw, cfw)
    for name, (x, cf) in three.items():
        s, r = sep_reatt(x, cf)
        L = f"{1e3 * (r - s):.1f} mm" if s else "attached"
        print(f"  {name:9s}: x_sep = {s}  x_reatt = {r}  (L = {L})")

    print("\nMean |difference| vs MOSE:")
    for lo, hi, lbl in ((0.10, 0.30, "upstream   "),
                        (0.35, 0.40, "post-reatt."),
                        (0.45, 0.52, "downstream ")):
        xs = np.linspace(lo, hi, 300)
        cm, pm = np.interp(xs, xm, cfm), np.interp(xs, xm, pwm)
        ec = 100 * np.mean(np.abs(np.interp(xs, xo, cfo) - cm) / np.abs(cm))
        ep = 100 * np.mean(np.abs(np.interp(xs, xo, pwo) - pm) / pm)
        print(f"  {lbl} ({lo:.2f}-{hi:.2f}):  Cf {ec:5.1f}%   p_w {ep:5.2f}%")

    # --- skin friction: verification (MOSE vs OpenFOAM, same grid) ----
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(xm, cfm, MOSE_C, lw=4.0, label="MOSE")
    ax.plot(xo, cfo, OF_C, lw=2.5, label="OpenFOAM")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlim(0.32, 0.41)
    ax.set_ylim(-0.002, 0.0075)
    style(ax, r"$x$  [m]", r"$C_f$")
    save(fig, f"SWBLI-cf-{m}-openfoam.svg")

    # --- wall pressure ------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(xm, pwm / P_INF, MOSE_C, lw=4.0, label="MOSE")
    ax.plot(xo, pwo / P_INF, OF_C, lw=2.5, label="OpenFOAM")
    ax.set_xlim(0.05, 0.52)
    style(ax, r"$x$  [m]", r"$p_w / p_\infty$")
    save(fig, f"SWBLI-pw-{m}-openfoam.svg")

    # --- validation: MOSE and OpenFOAM vs the reference data ----------
    # For SST both codes use the asymptotic omega wall condition (the correct
    # form on this wall-resolved mesh); for SA there is a single run each.
    # Labels stay short (the asymptotic-BC note goes in the docs caption).
    xmv, cfmv, xov, cfov = xm, cfm, xo, cfo   # base (Sutherland) MOSE + OpenFOAM

    fig, ax = plt.subplots(figsize=(10, 6))
    xe, cfe = ref_xy(case / "reference" / "schulein.dat")
    xw, cfw = ref_xy(case / "reference" / f"wind-{M}.dat")
    if xe is not None:
        ax.plot(xe, cfe, "ko", ms=6, label="Schulein (exp)")
    if xw is not None:
        ax.plot(xw, cfw, WIND_C, ls="-.", lw=2.5, label="Wind-US")
    ax.plot(xmv, cfmv, MOSE_C, lw=4.0, label="MOSE")
    ax.plot(xov, cfov, OF_C, lw=2.0, label="OpenFOAM")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlim(0.32, 0.41)
    ax.set_ylim(-0.002, 0.0075)
    style(ax, r"$x$  [m]", r"$C_f$")
    save(fig, f"SWBLI-cf-{m}-validation.svg")

    # --- omega wall-BC effect in MOSE (SST): same grid, only the BC changes -
    # MOSE run with both wall conditions (the `omega-wall-bc` option). Same
    # colour (MOSE) for both curves; line style carries the BC.
    mose_prac = case / "reference" / "MOSE-SST.tec"             # practical 800 nu/y^2
    mose_asy2 = case / "reference" / "MOSE-SST-asymptotic.tec"  # asymptotic 80 nu/y^2
    if M == "SST" and mose_prac.exists() and mose_asy2.exists():
        xp, cfp, _ = mose_wall(mose_prac)
        xa, cfa, _ = mose_wall(mose_asy2)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(xp, cfp, MOSE_C, ls="--", lw=2.5,
                label=r"MOSE, practical $\omega$ BC  ($800\,\nu/y^2$)")
        ax.plot(xa, cfa, MOSE_C, lw=3.5,
                label=r"MOSE, asymptotic $\omega$ BC  ($80\,\nu/y^2$)")
        ax.axhline(0.0, color="0.6", lw=0.8)
        ax.set_xlim(0.32, 0.45)
        ax.set_ylim(-0.002, 0.0075)
        style(ax, r"$x$  [m]", r"$C_f$")
        save(fig, "SWBLI-cf-sst-omegaBC-mose.svg")

    if args.plot:
        plt.show()


if __name__ == "__main__":
    main()
