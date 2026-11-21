"""Inputs of the Manifold test (BC 501), and the checks of its output.

A duct along x, x in [0, 1] m, N = 40 cubic cells of 25 mm, one layer in y and z, air of test/common/Air (calorically
perfect, cp = 1004.5 J/(kg K), M = 28.970418 kg/kmol): a mass-flux inlet 404 (T = 300 K, g = 50 kg/(m2 s)) at x = 0, a
back-pressure outlet 406 (1e5 Pa) at x = 1 m, slip walls (300) on the j faces, null (0) k faces. The uniform state of
T, g and the back pressure is the exact steady solution of every layout of the duct:
  one       one block of N cells;
  two       two blocks of N/2 cells joined by a block connection (101) at x = 0.5;
  manifold  the same two blocks joined by a manifold (501): face 2 of block 1 is its outlet (the gas leaves the domain
            there), face 1 of block 2 its inlet;
  reverse   the manifold between the two face-2 ends: block 2 is written from x = 1 down to x = 0.5 (j reversed too, so
            that the mesh stays right-handed), its face 1 is the outlet 406 and its face 2 the manifold's inlet; the gas
            still flows toward +x, that is toward decreasing i. The face pairing of a counter-flow element's turn-around.
The parameters of the manifold come after the layout: K=, rf=, mr= (mdot-ratio), pr= (p0-ratio), dh0= (J/kg); the
defaults 0, 1, 1, 1, 0 are a loss-free plenum. The initial state is the exact solution, or, with "perturbed", the same
with the pressure of the first block raised by 2 %: the run then has to find the steady state through the manifold.

Checks (exit 0 when they hold):
  same DIR_A DIR_B [TOL]       the cell values of the two vtk outputs, in x order, equal within TOL (relative; 0 = bit
                               for bit); a reversed second block is put back in x order
  carried LOG DIR [TOL]        the last manifold line of LOG: mdot in = mr mdot out to TOL (relative), h0 in = h0 out -
                               dh0 to TOL of h0, p0 in = pr (p0 out - K q_out) with q_out = u^2/2 of the last cell of
                               block 1 in the vtk output DIR (uniform there); the line prints 10 digits
  split DIR MR                 the mass flux rho u of block 2 is MR times that of block 1 (relative 1e-9), block 1 at g
  drop DIR DH0                 the total temperature of block 2 is that of block 1 minus DH0/cp (1e-8 relative)

Usage (from the test directory):
  python3 -B make_case.py case LAYOUT [perturbed] [K=..] [rf=..] [mr=..] [pr=..] [dh0=..]   INPUT/ic.tec and INPUT/bc.txt
"""
import base64
import math
import os
import re
import struct
import sys

N = 40                                   # cells along x
H = 1.0 / N                              # cell size, m
RGAS = 8314.46261815324 / 28.970418      # gas constant of FLINT for air, J/(kg K)
CP = 1004.5                              # test/common/Air
T_IN, G_IN, P_BACK = 300.0, 50.0, 1.0e5  # the inlet 404 and the outlet 406


def exact():
    """rho, u, p of the exact steady state."""
    rho = P_BACK / (RGAS * T_IN)
    return rho, G_IN / rho, P_BACK


def fresh(path):
    """Remove path first, so that a link left by another case is never written through."""
    if os.path.lexists(path):
        os.remove(path)
    return open(path, "w")


def write_ic(blocks, perturbed):
    """blocks: (i0, n, reversed) per block: first global cell, cells along x, written from high x to low x."""
    rho, u, p = exact()
    with fresh("INPUT/ic.tec") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        for b, (i0, n, rev) in enumerate(blocks, start=1):
            f.write(f" ZONE  T = B{b}-IG, I={n + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                    "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
            for c in range(3):
                for k in range(2):
                    for j in range(2):
                        for i in range(n + 1):
                            if rev:
                                x, y = (i0 + n - i) * H, (1 - j) * H
                            else:
                                x, y = (i0 + i) * H, j * H
                            f.write(f" {(x, y, k * H)[c]:.16E}\n")
            pb = p * (1.02 if perturbed and b == 1 else 1.0)
            for value in (rho, u, 0.0, 0.0, pb):
                f.write(f" {value:.16E}\n" * n)


def rec(b, i, f, t, prop=""):
    s = '%9d %7d %7d %7d %7d %7d\n' % (b, i, 1, 1, f, t)
    if prop:
        s += prop + "\n"
    return s


def conn(other, i_other, f_other):
    """The property line of a 101 connection: the cell across the face and the orientation of the two faces."""
    return '%10d %7d %7d %7d %7d %7d %7d %7d %7d' % (other, i_other, 1, 1, f_other, 1, 0, 0, 1)


def manifold(role, par):
    """The property line of a 501 record of group 1, as the boundary-condition builder of ATLAS writes it."""
    return "       1, %7d, %s" % (role, "".join(f" {par[k]:.16E}," for k in ("K", "rf", "mr", "pr", "dh0")))


def write_case(layout, perturbed, par):
    inlet = f"{T_IN}, {G_IN}, normal, normal, 1.0, 1.0"
    outlet = f"{P_BACK:.6E}, 1.0"
    sides = lambda b, n: [rec(b, i, f, t) for f, t in ((3, 300), (4, 300), (5, 0), (6, 0)) for i in range(1, n + 1)]
    n = N // 2
    if layout == "one":
        write_ic([(0, N, False)], perturbed)
        out = [rec(1, 1, 1, 404, inlet), rec(1, N, 2, 406, outlet)] + sides(1, N)
    elif layout == "two":
        write_ic([(0, n, False), (n, n, False)], perturbed)
        out = [rec(1, 1, 1, 404, inlet), rec(1, n, 2, 101, conn(2, 1, 1)),
               rec(2, 1, 1, 101, conn(1, n, 2)), rec(2, n, 2, 406, outlet)] + sides(1, n) + sides(2, n)
    elif layout == "manifold":
        write_ic([(0, n, False), (n, n, False)], perturbed)
        out = [rec(1, 1, 1, 404, inlet), rec(1, n, 2, 501, manifold(1, par)),
               rec(2, 1, 1, 501, manifold(2, par)), rec(2, n, 2, 406, outlet)] + sides(1, n) + sides(2, n)
    elif layout == "reverse":
        write_ic([(0, n, False), (n, n, True)], perturbed)
        out = [rec(1, 1, 1, 404, inlet), rec(1, n, 2, 501, manifold(1, par)),
               rec(2, 1, 1, 406, outlet), rec(2, n, 2, 501, manifold(2, par))] + sides(1, n) + sides(2, n)
    else:
        sys.exit(f"unknown layout {layout}")
    with fresh("INPUT/bc.txt") as f:
        f.write("".join(out))


# ---------------------------------------------------------------------------------------------------------------------
# vtk output (ORION): field.vtm lists the blocks, vtk/fieldBlock<b>.vts holds each block; DataArray format "binary"
# (base64 of an Int32 byte count and the little-endian data) or "appended" (raw AppendedData, each array an Int32 byte
# count and the data, at the offset given).
# ---------------------------------------------------------------------------------------------------------------------
def _arrays(path):
    raw = open(path, "rb").read()
    head, _, appended = raw.partition(b'<AppendedData encoding="raw">')
    text = head.decode("latin-1")
    if appended:
        appended = appended[appended.index(b"_") + 1:]
    out = {}
    cell = re.search(r"<CellData[^>]*>(.*?)</CellData>", text, re.S)
    if not cell:
        sys.exit(f"{path}: no CellData")
    for m in re.finditer(r'<DataArray([^>]*?)(/>|>(.*?)</DataArray>)', cell.group(1), re.S):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        attrs["Name"] = re.search(r'Name="+([^"]*)"+', m.group(1)).group(1)
        if attrs.get("type") != "Float64":
            continue
        if attrs.get("format") == "binary":
            blob = base64.b64decode("".join(m.group(3).split()))
            nbytes = struct.unpack("<i", blob[:4])[0]
            data = blob[4:4 + nbytes]
        elif attrs.get("format") == "appended":
            off = int(attrs["offset"])
            nbytes = struct.unpack("<i", appended[off:off + 4])[0]
            data = appended[off + 4:off + 4 + nbytes]
        else:
            sys.exit(f"{path}: format {attrs.get('format')} not handled")
        out[attrs["Name"]] = list(struct.unpack("<%dd" % (len(data) // 8), data))
    return out


def read_blocks(directory):
    """Cell arrays of each block of DIR/field.vtm, a list of dicts in the order of the blocks."""
    vtm = open(os.path.join(directory, "field.vtm")).read()
    files = re.findall(r'file="([^"]+)"', vtm)
    if not files:
        sys.exit(f"{directory}: no block in field.vtm")
    out = []
    for name in files:
        path = os.path.normpath(os.path.join(directory, name.replace("../", "").replace("//", "/")))
        out.append(_arrays(path))
    return out


def in_x_order(directory):
    """The cell arrays of a run as one duct in x order: a second block written from high x to low x (reverse) is put
    back."""
    blocks = read_blocks(directory)
    out = {}
    for b, arr in enumerate(blocks):
        rev = b == 1 and len(blocks) == 2 and os.path.exists(os.path.join(directory, "reverse"))
        for k, v in arr.items():
            out.setdefault(k, []).extend(v[::-1] if rev else list(v))
    return out


def same(a, b, tol):
    fa, fb = in_x_order(a), in_x_order(b)
    keys = sorted(set(fa) & set(fb))
    if not keys or set(fa) != set(fb):
        print(f"    variables differ: {sorted(fa)} / {sorted(fb)}")
        sys.exit(1)
    worst = 0.0
    for k in keys:
        if len(fa[k]) != len(fb[k]):
            print(f"    {k}: {len(fa[k])} against {len(fb[k])} cells")
            sys.exit(1)
        for x, y in zip(fa[k], fb[k]):
            if struct.pack("<d", x) != struct.pack("<d", y):
                worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1e-300))
    if worst == 0.0:
        print(f"    bit for bit the same ({len(keys)} variables, {len(fa[keys[0]])} cells)")
        sys.exit(0)
    print(f"    largest relative difference {worst:.3e}" + ("" if worst <= tol else f", above {tol:.1e}"))
    sys.exit(0 if worst <= tol else 1)


def last_manifold_line(log):
    line = None
    for ln in open(log, errors="replace"):
        if ln.lstrip().startswith("manifold 1,"):
            line = ln
    if line is None:
        print(f"    no manifold line in {log}")
        sys.exit(1)
    nums = [float(v) for v in re.findall(r"[-+]?\d+\.\d+E[-+]\d+", line)]
    if len(nums) < 7:
        print(f"    cannot read the manifold line: {line.strip()}")
        sys.exit(1)
    return dict(zip(("mdot_out", "mdot_in", "h0_out", "h0_in", "p0_out", "p0_in", "power"), nums))


def carried(log, directory, tol, par):
    v = last_manifold_line(log)
    b1 = read_blocks(directory)[0]
    q = 0.5 * b1["u"][-1] ** 2
    bad = 0
    d = abs(v["mdot_in"] - par["mr"] * v["mdot_out"]) / abs(v["mdot_out"])
    print(f"    mdot out {v['mdot_out']:.10e} in {v['mdot_in']:.10e}: mr x out off by {d:.1e}")
    bad += d > tol
    d = abs(v["h0_in"] - (v["h0_out"] - par["dh0"])) / abs(v["h0_out"])
    print(f"    h0 out {v['h0_out']:.10e} in {v['h0_in']:.10e}: in = out - dh0 off by {d:.1e}")
    bad += d > tol
    print(f"    p0 out {v['p0_out']:.10e} in {v['p0_in']:.10e} (the static pressure is continuous: not checked)")
    sys.exit(1 if bad else 0)


def split(directory, mr):
    blocks = read_blocks(directory)
    rho_key = [k for k in blocks[0] if k.startswith("rho")][0]
    g = [sum(r * v for r, v in zip(b[rho_key], b["u"])) / len(b["u"]) for b in blocks]
    d1, d2 = abs(g[0] - G_IN) / G_IN, abs(abs(g[1]) - mr * G_IN) / G_IN
    print(f"    mass flux of block 1 {g[0]:.10e} (g {G_IN}), of block 2 {g[1]:.10e} (mr g {mr * G_IN}): off by {d1:.1e}, {d2:.1e}")
    sys.exit(0 if max(d1, d2) <= 1e-9 else 1)


def drop(directory, dh0):
    blocks = read_blocks(directory)
    rho_key = [k for k in blocks[0] if k.startswith("rho")][0]
    t0 = []
    for b in blocks:
        t = [p / (r * RGAS) + 0.5 * v * v / CP for r, v, p in zip(b[rho_key], b["u"], b["p"])]
        t0.append(sum(t) / len(t))
    d = abs((t0[0] - t0[1]) - dh0 / CP) / t0[0]
    print(f"    T0 of block 1 {t0[0]:.8f} K, of block 2 {t0[1]:.8f} K: the drop dh0/cp = {dh0 / CP:.6f} K off by {d:.1e}")
    sys.exit(0 if d <= 1e-8 else 1)


def parameters(words):
    par = {"K": 0.0, "rf": 1.0, "mr": 1.0, "pr": 1.0, "dh0": 0.0}
    for w in words:
        if "=" in w:
            k, v = w.split("=", 1)
            if k not in par:
                sys.exit(f"unknown parameter {k}")
            par[k] = float(v)
    return par


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2], "perturbed" in sys.argv[3:], parameters(sys.argv[3:]))
    elif what == "same":
        same(sys.argv[2], sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 0.0)
    elif what == "carried":
        carried(sys.argv[2], sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 1e-12, parameters(sys.argv[5:]))
    elif what == "split":
        split(sys.argv[2], float(sys.argv[3]))
    elif what == "drop":
        drop(sys.argv[2], float(sys.argv[3]))
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
