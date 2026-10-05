"""TIME field data of a MOSE .vtm: read from its first block file, in any vtk format (ascii, base64 or raw appended)."""
import base64
import re
import struct
from pathlib import Path


def time_of(vtm):
    vtm = Path(vtm)
    entry = re.search(r'file="([^"]+)"', vtm.read_text()).group(1)
    vts = vtm.parent / entry
    if not vts.is_file():
        vts = vtm.parent / "vtk" / Path(entry).name
    data = vts.read_bytes()
    cut = data.find(b"<AppendedData")
    head = data[:cut] if cut >= 0 else data
    tag = re.search(rb'<DataArray[^>]*Name="TIME"[^>]*?(/?)>', head)
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
