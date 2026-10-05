"""Wall outputs in vtk raw: every wallB<b>F<f>.vts has the extent of face f of block b and its values.

For each OUTPUT/vtk/wallB<b>F<f>.vts the extent of block b is the WholeExtent of
fieldBlock<b>.vts (0 NI 0 NJ 0 NK), in whatever vtk format the field is written; the wall
file must have the extent of face f, the one of the block with the direction normal to
the face collapsed to its first index (0 0 0 NJ 0 NK for faces 1 and 2, 0 NI 0 0 0 NK for
faces 3 and 4, 0 NI 0 NJ 0 0 for faces 5 and 6), as WholeExtent and as Piece Extent; its
Points array must hold three values for each of the points of that extent, and each
CellData array one value for each of its cells (a collapsed direction counts one layer of
cells, as in the cell-centered ordered data of the Tecplot format and in the writers of
ORION). Every array of a wall file must be raw appended data, whose values are counted
from its size header.

Usage: python check_wall_vts.py <OUTPUT/vtk directory> <label>
Exit 0 = pass, 1 = fail.
"""
import re
import struct
import sys
from pathlib import Path

vtk, label = Path(sys.argv[1]), sys.argv[2]
SIZE = {"Float64": 8, "Float32": 4, "Int32": 4, "UInt32": 4, "Int64": 8, "UInt64": 8}


def whole_extent(path):
    data = path.read_bytes()
    cut = data.find(b"<AppendedData")
    head = (data[:cut] if cut >= 0 else data).decode("latin-1")
    return [int(x) for x in re.search(r'WholeExtent="([^"]*)"', head).group(1).split()]


def parse(path):
    data = path.read_bytes()
    cut = data.find(b"<AppendedData")
    head = (data[:cut] if cut >= 0 else data).decode("latin-1")
    whole = [int(x) for x in re.search(r'WholeExtent="([^"]*)"', head).group(1).split()]
    piece = [int(x) for x in re.search(r'<Piece[^>]*Extent="([^"]*)"', head).group(1).split()]
    header = re.search(r'header_type="([^"]*)"', head)
    hsize = SIZE[header.group(1)] if header else 4
    start = data.find(b"_", cut) + 1 if cut >= 0 else -1
    arrays = {}
    for section in ("Points", "CellData", "PointData"):
        m = re.search(r"<%s>(.*?)</%s>" % (section, section), head, re.S)
        if not m:
            continue
        for a in re.finditer(r"<DataArray([^>]*)/?>", m.group(1)):
            attrs = dict(re.findall(r'(\w+)="([^"]*)"', a.group(1)))
            if "offset" not in attrs:
                arrays[(section, attrs.get("Name", "?"))] = (None, 0)
                continue
            offset = int(attrs["offset"])
            nbytes = struct.unpack("<Q" if hsize == 8 else "<I", data[start + offset:start + offset + hsize])[0]
            arrays[(section, attrs.get("Name", "?"))] = (nbytes // SIZE[attrs["type"]], int(attrs.get("NumberOfComponents", "1")))
    return whole, piece, arrays


fail, nwall = [], 0
for wall in sorted(vtk.glob("wallB*F*.vts")):
    m = re.fullmatch(r"wallB(\d+)F(\d)\.vts", wall.name)
    b, f = int(m.group(1)), int(m.group(2))
    field = vtk / ("fieldBlock%d.vts" % b)
    if not field.is_file():
        fail.append(f"{wall.name}: {field.name} not found")
        continue
    nwall += 1
    block = whole_extent(field)
    want = list(block)
    d = (f - 1) // 2                      # 0: i, 1: j, 2: k
    want[2 * d + 1] = want[2 * d]
    whole, piece, arrays = parse(wall)
    npts = 1
    ncell = 1
    for e in range(3):
        n = want[2 * e + 1] - want[2 * e]
        npts *= n + 1
        ncell *= max(n, 1)
    if whole != want or piece != want:
        fail.append(f"{wall.name}: extent {whole} (piece {piece}), {want} expected")
    if not any(s == "Points" for s, _ in arrays):
        fail.append(f"{wall.name}: no Points array")
    for (section, name), (nval, ncomp) in arrays.items():
        if nval is None:
            fail.append(f"{wall.name}: {section} {name}: not raw appended data")
            continue
        expected = (npts if section in ("Points", "PointData") else ncell) * ncomp
        if nval != expected:
            fail.append(f"{wall.name}: {section} {name}: {nval} values, {expected} expected")
    if not any(s == "CellData" for s, _ in arrays):
        fail.append(f"{wall.name}: no CellData")

for line in fail:
    print("  " + line)
if nwall == 0:
    print(f"  {label}: no wall file in {vtk}")
    sys.exit(1)
print(f"  {label}: {nwall} wall files, {len(fail)} faults")
sys.exit(1 if fail else 0)
