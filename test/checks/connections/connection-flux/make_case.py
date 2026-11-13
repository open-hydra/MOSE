"""Inputs of the ConnectionFlux tests, and the checks of their output.

A periodic 1-D domain along x, x in [0, 1] m, N = 40 cubic cells of 25 mm, written either as one block periodic on
itself (connection 201 between its faces 1 and 2) or as two blocks of 20 cells (block cut 101 at x = 0.5, periodic
connection 201 between the outer faces). Faces 3 to 6 are null (type 0). Air of test/common/Air (calorically perfect,
cp = 1004.5 J/(kg K), M = 28.970418 kg/kmol). Two states:
  sod   at rest, Sod's states: rho = 1 kg/m3, p = 1e5 Pa for x < 0.5, rho = 0.125 kg/m3, p = 1e4 Pa for x > 0.5. A shock
        leaves each discontinuity, at the block cut and at the periodic connection, so the shock sensors switch at the
        connection cells;
  wave  low Mach number (u = 10 m/s, p = 1e5 Pa), with an entropy wave rho = 1.2 (1 + 0.1 sin(2 pi x/L)) and a pressure
        pulse of 1e-3 p over x in [0.35, 0.65]: for the preconditioned solver (HLLC-PC), whose reference velocity Ur
        matters at the connection cells.

A face shared by two blocks is computed twice, once by each block. With one face flux, identical on the two sides and
equal to that of an interior face (the shock sensor of the cell on the low-index side of the face, the reference
velocity of the preconditioner of both cells):
  - the two-block run gives the one-block run bit for bit (each cell of this 1-D case sums two face fluxes, so the order of
    the sums does not matter);
  - mass, momentum and energy are conserved to round-off in time-accurate runs (a periodic domain has no boundary flux).

Usage (from the test directory):
  python3 -B make_case.py case one|two sod|wave INPUT/ic.tec and INPUT/bc.txt
  python3 -B make_case.py same DIR_A DIR_B       exit 0 when the cell values of the vtk outputs in DIR_A and DIR_B, in x
                                                 order, are bit for bit the same (prints the largest difference otherwise)
  python3 -B make_case.py conserved DIR          exit 0 when the totals of mass, momentum and energy of the final field in
                                                 DIR equal those of the initial sod state to 1e-13 (relative)
"""
import base64
import math
import os
import re
import struct
import sys

N = 40                          # cells along x
H = 1.0 / N                     # cell size, m
RGAS = 8314.46261815324 / 28.970418      # gas constant of FLINT for air, J/(kg K)
CP = 1004.5                     # test/common/Air
LEFT = (1.0, 1.0e5)             # sod: rho, p for x < 0.5
RIGHT = (0.125, 1.0e4)          # sod: rho, p for x > 0.5


def state(x, kind):
    """rho, u, p at x."""
    if kind == "sod":
        rho, p = LEFT if x < 0.5 else RIGHT
        return rho, 0.0, p
    rho = 1.2 * (1.0 + 0.1 * math.sin(2.0 * math.pi * x))
    p = 1.0e5 * (1.0 + (1.0e-3 if 0.35 < x < 0.65 else 0.0))
    return rho, 10.0, p


def fresh(path):
    """Remove path first, so that a link left by another case is never written through."""
    if os.path.lexists(path):
        os.remove(path)
    return open(path, "w")


def write_ic(blocks, kind):
    """blocks: list of (i0, n): first global cell and number of cells along x of each block."""
    with fresh("INPUT/ic.tec") as f:
        f.write(' VARIABLES ="x" "y" "z"  "rho1" "u" "v" "w" "p"\n')
        for b, (i0, n) in enumerate(blocks, start=1):
            f.write(f" ZONE  T = B{b}-IG, I={n + 1}, J=2, K=2, DATAPACKING=BLOCK, "
                    "VARLOCATION=([1-3]=NODAL,[4-8]=CELLCENTERED), SOLUTIONTIME=0.0\n")
            for c in range(3):
                for k in range(2):
                    for j in range(2):
                        for i in range(n + 1):
                            f.write(f" {((i0 + i) * H, j * H, k * H)[c]:.16E}\n")
            cells = [state((i0 + i + 0.5) * H, kind) for i in range(n)]
            for col in ([c[0] for c in cells], [c[1] for c in cells], [0.0] * n, [0.0] * n, [c[2] for c in cells]):
                f.write("".join(f" {v:.16E}\n" for v in col))


def rec(b, i, f, t, prop=None):
    s = '%9d %7d %7d %7d %7d %7d\n' % (b, i, 1, 1, f, t)
    if prop is not None:
        s += '%10d %7d %7d %7d %7d %7d %7d %7d %7d\n' % prop
    return s


def write_case(layout, kind):
    if kind not in ("sod", "wave"):
        sys.exit(f"unknown state {kind}")
    if layout == "one":
        write_ic([(0, N)], kind)
        out = [rec(1, 1, 1, 201, (1, N, 1, 1, 2, 1, 0, 0, 1)), rec(1, N, 2, 201, (1, 1, 1, 1, 1, 1, 0, 0, 1))]
        out += [rec(1, i, f, 0) for f in (3, 4, 5, 6) for i in range(1, N + 1)]
    elif layout == "two":
        n = N // 2
        write_ic([(0, n), (n, n)], kind)
        out = []
        for b, other in ((1, 2), (2, 1)):
            # face 1 of block 1 and face 2 of block 2: periodic connection; face 2 of block 1 and face 1 of block 2: cut
            out.append(rec(b, 1, 1, 201 if b == 1 else 101, (other, n, 1, 1, 2, 1, 0, 0, 1)))
            out.append(rec(b, n, 2, 101 if b == 1 else 201, (other, 1, 1, 1, 1, 1, 0, 0, 1)))
            out += [rec(b, i, f, 0) for f in (3, 4, 5, 6) for i in range(1, n + 1)]
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
        # the name may come doubly quoted (Name=""rho(1)""): read it on its own
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


def read_field(directory):
    """Cell arrays of all blocks of DIR/field.vtm, concatenated in the order of the blocks (x order here)."""
    vtm = open(os.path.join(directory, "field.vtm")).read()
    files = re.findall(r'file="([^"]+)"', vtm)
    if not files:
        sys.exit(f"{directory}: no block in field.vtm")
    out = {}
    for name in files:
        path = os.path.normpath(os.path.join(directory, name.replace("../", "").replace("//", "/")))
        for k, v in _arrays(path).items():
            out.setdefault(k, []).extend(v)
    return out


def same(a, b):
    fa, fb = read_field(a), read_field(b)
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
    print(f"    largest relative difference {worst:.3e}")
    sys.exit(1)


def totals(rho, u, p):
    gm1 = RGAS / (CP - RGAS)
    mass = sum(rho) * H ** 3
    mom = sum(r * v for r, v in zip(rho, u)) * H ** 3
    ener = sum(q / gm1 + 0.5 * r * v * v for r, v, q in zip(rho, u, p)) * H ** 3
    return mass, mom, ener


def conserved(directory):
    f = read_field(directory)
    rho_key = [k for k in f if k.startswith("rho")][0]
    m1, q1, e1 = totals(f[rho_key], f["u"], f["p"])
    init = [state((i + 0.5) * H, "sod") for i in range(N)]
    m0, _, e0 = totals([c[0] for c in init], [c[1] for c in init], [c[2] for c in init])
    q_scale = m0 * (1.0e5 / 1.0) ** 0.5           # mass times a velocity scale sqrt(p/rho) of the left state
    dm, dq, de = abs(m1 - m0) / m0, abs(q1) / q_scale, abs(e1 - e0) / e0
    ok = max(dm, dq, de) < 1.0e-13
    print(f"    drift of mass {dm:.1e}, momentum {dq:.1e}, energy {de:.1e}")
    sys.exit(0 if ok else 1)


def main():
    what = sys.argv[1]
    if what == "case":
        write_case(sys.argv[2], sys.argv[3])
    elif what == "same":
        same(sys.argv[2], sys.argv[3])
    elif what == "conserved":
        conserved(sys.argv[2])
    else:
        sys.exit(f"unknown request {what}")


if __name__ == "__main__":
    main()
