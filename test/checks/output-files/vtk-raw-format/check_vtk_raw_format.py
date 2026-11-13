"""vtk raw outputs of MOSE: every DataArray of every block file is raw appended data.

For each OUTPUT/<name>.vtm the block files it lists are looked for relative to
the .vtm and in OUTPUT/vtk (where MOSE writes them); the XML header of each
(everything before <AppendedData>) must name format="appended" in every
DataArray, and the file must hold the raw appended data.

Usage: python check_vtk_raw_format.py <OUTPUT directory> <label>
Exit 0 = pass, 1 = fail.
"""
import re
import sys
from pathlib import Path

out, label = Path(sys.argv[1]), sys.argv[2]
fail = []
nfiles = 0
vtms = sorted(out.glob("*.vtm"))
if not vtms:
    fail.append("no .vtm in " + str(out))
for vtm in vtms:
    for f in re.findall(r'file="([^"]+)"', vtm.read_text()):
        vts = vtm.parent / f
        if not vts.is_file():
            vts = out / "vtk" / Path(f).name
        if not vts.is_file():
            fail.append(f"{vtm.name}: block file {f} not found")
            continue
        data = vts.read_bytes()
        cut = data.find(b"<AppendedData")
        head = (data[:cut] if cut >= 0 else data).decode("latin-1")
        formats = re.findall(r'<DataArray[^>]*format="([^"]+)"', head)
        nfiles += 1
        if not formats:
            fail.append(f"{vts.name}: no DataArray")
        elif any(x != "appended" for x in formats):
            fail.append(f"{vts.name}: formats {sorted(set(formats))}, not all appended")
        elif cut < 0 or b'encoding="raw"' not in data[cut:cut + 64]:
            fail.append(f"{vts.name}: no raw appended data")

for line in fail:
    print("  " + line)
print(f"  {label}: {nfiles} block files, {len(fail)} not raw")
sys.exit(1 if fail else 0)
