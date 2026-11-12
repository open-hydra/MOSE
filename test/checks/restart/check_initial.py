"""Initial condition from a solution written by MOSE (Restart test, part 3).

Compares the cell values of the final solution of two runs: the restart (newrun = false) from a solution, and a new run
(newrun = true) from the same solution copied into INPUT/ as the initial condition. The two runs start from the same values
and make the same steps (constant CFL number, no time-dependent boundary values), so their cell values must be the same,
bit for bit; only the time and the counters stored with the files differ.

- vtk (ascii, base64 or raw appended): every CellData array of every block file listed by OUTPUT/field.vtm, decoded to its
  bytes (raw, base64) or to its tokens (ascii);
- Tecplot ASCII: the numbers of OUTPUT/field.tec after the header of each zone (VARIABLES, ZONE and auxiliary data lines are
  not compared: they hold the solution time).

Usage: python3 -B check_initial.py FORMAT <OUTPUT of the restart> <OUTPUT of the new run> [LABEL]
       (FORMAT: vtk or tecplot; LABEL: the name printed, FORMAT by default)
Exit 0 = the cell values are the same, 1 = they differ or a file is missing.
"""
import base64
import re
import sys
from pathlib import Path


def vtk_cells(vtm):
    """{(block entry, array name): bytes or tokens} of the CellData arrays of the block files of a .vtm."""
    out = {}
    for entry in re.findall(r'file="([^"]+)"', vtm.read_text()):
        data = (vtm.parent / entry).read_bytes()
        cut = data.find(b"<AppendedData")
        head = data[:cut] if cut >= 0 else data
        cells = head[head.index(b"<CellData"):head.index(b"</CellData>")]
        for tag in re.finditer(rb'<DataArray[^>]*?(/?)>', cells):
            attrs = tag.group(0)
            name = re.search(rb'Name="([^"]+)"', attrs).group(1).decode()
            fmt = re.search(rb'format="([^"]+)"', attrs).group(1).decode()
            if fmt == "appended":
                start = data.index(b"_", cut) + 1
                off = start + int(re.search(rb'offset="(\d+)"', attrs).group(1))
                nbytes = int.from_bytes(data[off:off + 4], "little")
                value = data[off + 4:off + 4 + nbytes]
            else:
                body = cells[tag.end():cells.index(b"</DataArray>", tag.end())]
                value = base64.b64decode(b"".join(body.split())) if fmt == "binary" else tuple(body.split())
            out[(entry, name)] = value
    return out


def tecplot_numbers(tec):
    """The lines of numbers of a Tecplot ASCII file, without the header lines of its zones."""
    keep = []
    for line in tec.read_text().splitlines():
        word = line.strip().split("=")[0].strip().upper()
        if word.startswith(("VARIABLES", "ZONE", "TITLE", "DATASETAUXDATA", "AUXDATA", "DT=")):
            continue
        keep.append(line.split())
    return keep


fmt, restart, new = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
label = sys.argv[4] if len(sys.argv) > 4 else fmt
name = "field.vtm" if fmt == "vtk" else "field.tec"
missing = [str(p / name) for p in (restart, new) if not (p / name).is_file()]
if missing:
    print(f"  {label}: FAIL — missing {', '.join(missing)}")
    sys.exit(1)
if fmt == "vtk":
    a, b = vtk_cells(restart / name), vtk_cells(new / name)
    same = a.keys() == b.keys() and all(a[k] == b[k] for k in a)
    what = f"{len(a)} cell arrays"
else:
    a, b = tecplot_numbers(restart / name), tecplot_numbers(new / name)
    same = a == b
    what = f"{sum(len(x) for x in a)} values"
if same:
    print(f"  {label}: ok — the new run from the solution has the cell values of the restart from it ({what}, bit for bit)")
    sys.exit(0)
print(f"  {label}: FAIL — the new run from the solution and the restart from it have different cell values")
sys.exit(1)
