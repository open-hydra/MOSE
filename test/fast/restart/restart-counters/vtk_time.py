"""Field data of a MOSE solution: of a .vtm (TIME, or a counter such as MOSE_OUTPUT_NUMBER), read from its first block
file, in any vtk format (ascii, base64 or raw appended); and the solution time of a Tecplot ASCII .tec file, the
SOLUTIONTIME of its first zone header."""
import base64
import re
import struct
from pathlib import Path


def field_of(vtm, name):
    vtm = Path(vtm)
    entry = re.search(r'file="([^"]+)"', vtm.read_text()).group(1)
    vts = vtm.parent / entry
    if not vts.is_file():
        vts = vtm.parent / "vtk" / Path(entry).name
    data = vts.read_bytes()
    cut = data.find(b"<AppendedData")
    head = data[:cut] if cut >= 0 else data
    tag = re.search(rb'<DataArray[^>]*Name="' + name.encode() + rb'"[^>]*?(/?)>', head)
    attrs = tag.group(0)
    fmt = re.search(rb'format="([^"]+)"', attrs).group(1).decode()
    if fmt == "appended":
        start = data.index(b"_", cut) + 1
        off = start + int(re.search(rb'offset="(\d+)"', attrs).group(1))
        return struct.unpack("<d", data[off + 4:off + 12])[0]
    body = head[tag.end():head.index(b"</DataArray>", tag.end())]
    if fmt == "ascii":
        return float(body.split()[0])
    return struct.unpack("<d", base64.b64decode(b"".join(body.split()))[4:12])[0]


def time_of(vtm):
    vtm = Path(vtm)
    if vtm.suffix == ".tec":
        zone = next(l for l in vtm.read_text().splitlines() if l.lstrip().upper().startswith("ZONE"))
        return float(re.search(r"SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)", zone).group(1).replace("D", "E"))
    return field_of(vtm, "TIME")
