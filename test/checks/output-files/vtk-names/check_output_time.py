"""Solution time of the field, wall and diagnostic outputs of a steady run, in vtk and in Tecplot.

A steady run stamps each of its outputs with minus its global iteration count, the same number in every format, so that
a reader tells a steady solution from a time-accurate one by the sign of its time:
  vtk      the field data TIME of every block file of field.vtm, wall.vtm and diagnostic.vtm is minus the count, which
           the field solution also holds as CYCLE;
  tecplot  every zone header of field.tec, wall.tec and diagnostic.tec gives SOLUTIONTIME = minus the count and no
           STRANDID: the zones of a steady solution are static, they belong to no strand of a time series.

Usage: python check_output_time.py <OUTPUT directory> vtk|tecplot <iteration count>
Exit 0 = pass, 1 = fail.
"""
import base64
import re
import struct
import sys
from pathlib import Path

out, fmt, count = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])


def field_data(vts):
    """{name: first value} of the field data of a .vts file, in any vtk format (ascii, base64 or raw appended)."""
    data = vts.read_bytes()
    cut = data.find(b"<AppendedData")
    head = data[:cut] if cut >= 0 else data
    fd = re.search(rb"<FieldData>(.*?)</FieldData>", head, re.S)
    values = {}
    if not fd:
        return values
    for tag in re.finditer(rb"<DataArray[^>]*?>", fd.group(1)):
        attrs = tag.group(0)
        name = re.search(rb'Name="([^"]+)"', attrs).group(1).decode()
        kind = re.search(rb'type="([^"]+)"', attrs).group(1).decode()
        form = re.search(rb'format="([^"]+)"', attrs).group(1).decode()
        code, size = {"Float64": ("<d", 8), "Float32": ("<f", 4), "Int32": ("<i", 4), "Int64": ("<q", 8)}[kind]
        if form == "appended":
            start = data.index(b"_", cut) + 1
            off = start + int(re.search(rb'offset="(\d+)"', attrs).group(1))
            values[name] = struct.unpack(code, data[off + 4:off + 4 + size])[0]
        else:
            body = fd.group(1)[tag.end():fd.group(1).index(b"</DataArray>", tag.end())]
            if form == "ascii":
                values[name] = float(body.split()[0])
            else:
                values[name] = struct.unpack(code, base64.b64decode(b"".join(body.split()))[4:4 + size])[0]
    return values


fail = []
for name in ("field", "wall", "diagnostic"):
    if fmt == "vtk":
        vtm = out / f"{name}.vtm"
        blocks = [vtm.parent / f for f in re.findall(r'file="([^"]+)"', vtm.read_text())] if vtm.is_file() else []
        if not blocks:
            fail.append(f"{vtm}: missing or lists no block")
            continue
        stamps = [field_data(b) for b in blocks]
        times = sorted({s.get("TIME") for s in stamps}, key=str)
        cycles = sorted({s.get("CYCLE") for s in stamps}, key=str) if name == "field" else [count]
        good = times == [-count] and cycles == [count]
        print(f"  {name}.vtm: TIME {times} of {len(blocks)} block files" +
              (f", CYCLE {cycles}" if name == "field" else "") + f"  [{'ok' if good else 'FAIL'}]")
    else:
        tec = out / f"{name}.tec"
        zones = [l for l in tec.read_text().splitlines() if l.lstrip().upper().startswith("ZONE")] if tec.is_file() else []
        if not zones:
            fail.append(f"{tec}: missing or holds no zone")
            continue
        times = sorted({float(m.group(1).replace("D", "E")) if m else None
                        for m in (re.search(r"SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)", z) for z in zones)}, key=str)
        strands = [z for z in zones if re.search(r"STRANDID", z, re.I)]
        good = times == [-count] and not strands
        print(f"  {name}.tec: SOLUTIONTIME {times} in {len(zones)} zone headers, {len(strands)} with a STRANDID"
              f"  [{'ok' if good else 'FAIL'}]")
    if not good:
        fail.append(f"{name}: expected the time {-count} (minus the iteration count){'' if fmt == 'vtk' else ' and no STRANDID'}")

if fail:
    for line in fail:
        print("  " + line)
    print(f"VtkNames, solution time ({fmt}): FAIL")
    sys.exit(1)
print(f"VtkNames, solution time ({fmt}): ok — field, wall and diagnostic stamped with minus the iteration count, {-count}")
sys.exit(0)
