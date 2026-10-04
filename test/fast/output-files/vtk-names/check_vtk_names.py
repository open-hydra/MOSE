"""VTK outputs of MOSE: well-formed XML and the right variable names.

For each multi-block file OUTPUT/<name>.vtm (field, wall, diagnostic) the .vtm
must parse as XML, every block file it lists must exist relative to it or to
OUTPUT/vtk (where MOSE writes them), the XML part of every .vts (everything
before the raw appended data, if any) must parse, and the names of the cell
data arrays must be the variables MOSE writes for this case, without quotes.

Usage: python check_vtk_names.py <OUTPUT directory>
Exit 0 = pass, 1 = fail.
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

out = Path(sys.argv[1])
EXPECTED = {
    # single-species air, Navier-Stokes, sol-variables = thermo (the default)
    "field": ["rho(1)", "u", "v", "w", "p", "T", "g", "R"],
    # wall-variables = mech thermal (the default)
    "wall": ["y+", "tauX", "tauY", "tauZ", "pw", "Tw", "qw"],
    # steady run: residuals, local time step, shock sensor
    "diagnostic": ["rho", "rhou", "rhov", "rhow", "rhoe", "dt", "beta"],
}


def xml_part(path):
    """The XML of a .vts: up to the raw appended data, closed again if it was cut there."""
    data = path.read_bytes()
    cut = data.find(b"<AppendedData")
    if cut < 0:
        return data
    return data[:cut] + b"</VTKFile>"


fail = []
for name, expected in EXPECTED.items():
    vtm = out / f"{name}.vtm"
    if not vtm.is_file():
        fail.append(f"{vtm}: missing")
        continue
    try:
        files = [d.get("file") for d in ET.parse(vtm).getroot().iter("DataSet")]
    except ET.ParseError as e:
        fail.append(f"{vtm}: not well-formed XML ({e})")
        continue
    if not files:
        fail.append(f"{vtm}: lists no block")
    for f in files:
        vts = vtm.parent / f
        if not vts.is_file():
            vts = out / "vtk" / Path(f).name
        if not vts.is_file():
            fail.append(f"{vtm}: block file {f} not found")
            continue
        try:
            root = ET.fromstring(xml_part(vts))
        except ET.ParseError as e:
            fail.append(f"{vts.name}: not well-formed XML ({e})")
            continue
        names = [a.get("Name") for a in root.iter("DataArray")
                 if a.get("Name") not in ("Points", "TIME")]
        if names != expected:
            fail.append(f"{vts.name}: variable names {names}, expected {expected}")

if fail:
    for line in fail:
        print("  " + line)
    print("VtkNames: FAIL")
    sys.exit(1)
print("VtkNames: ok — field, wall and diagnostic VTK files are well-formed XML with the expected names")
sys.exit(0)
