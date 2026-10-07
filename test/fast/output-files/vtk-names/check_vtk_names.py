"""VTK outputs of MOSE: well-formed XML and the right variable names.

For each multi-block file OUTPUT/<name>.vtm (field, wall, diagnostic) the .vtm
must parse as XML, every block file it lists must exist relative to it or to
OUTPUT/vtk (where MOSE writes them), the XML part of every .vts (everything
before the raw appended data, if any) must parse, and the names of the cell
data arrays must be the variables MOSE writes for this case, without quotes.
Field data (TIME, CYCLE and the counters a restart reads back) are global
values, not variables, and are not counted.

The three writers of MOSE (field, wall, diagnostic) each decide the quotes of
the names, the paths and the vtk format: they must agree. Every .vtm lists its
block files as vtk/<file>.vts, relative to its own folder (no leading or double
slash), and that path must resolve; every data array of the three outputs has
the format of sol-format: appended (vtk raw) or ascii (vtk ascii).

With the second argument sa-transport: the run with the Spalart-Allmaras model and sol-variables = thermo
transport, whose field adds the turbulence variable, the laminar viscosity and conductivity and the eddy viscosity, and
whose diagnostic adds the residual of the turbulence variable.

Usage: python check_vtk_names.py <OUTPUT directory> [sa-transport] [ascii]
(ascii: the run with sol-format = vtk ascii; vtk raw otherwise)
Exit 0 = pass, 1 = fail.
"""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

out = Path(sys.argv[1])
flags = set(sys.argv[2:])
FORMAT = "ascii" if "ascii" in flags else "appended"
EXPECTED = {
    # single-species air, Navier-Stokes, sol-variables = thermo (the default)
    "field": ["rho(1)", "u", "v", "w", "p", "T", "g", "R"],
    # wall-variables = mech thermal (the default)
    "wall": ["y+", "tauX", "tauY", "tauZ", "pw", "Tw", "qw"],
    # steady run: residuals, local time step, shock sensor
    "diagnostic": ["rho", "rhou", "rhov", "rhow", "rhoe", "dt", "beta"],
}
if "sa-transport" in flags:
    EXPECTED["field"] = ["rho(1)", "u", "v", "w", "p", "rhoNuTilde", "T", "g", "R", "mil", "kl", "mut"]
    EXPECTED["diagnostic"] = ["rho", "rhou", "rhov", "rhow", "rhoe", "rhoNuTilde", "dt", "beta"]


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
        if not re.fullmatch(r"vtk/[^/]+\.vts", f or ""):
            fail.append(f"{vtm}: block file listed as {f!r}, not vtk/<file>.vts relative to the .vtm")
        vts = vtm.parent / f
        if not vts.is_file():
            fail.append(f"{vtm}: block file {f} not found relative to the .vtm")
            continue
        try:
            root = ET.fromstring(xml_part(vts))
        except ET.ParseError as e:
            fail.append(f"{vts.name}: not well-formed XML ({e})")
            continue
        field_data = {id(a) for fd in root.iter("FieldData") for a in fd.iter("DataArray")}
        names = [a.get("Name") for a in root.iter("DataArray")
                 if id(a) not in field_data and a.get("Name") != "Points"]
        if names != expected:
            fail.append(f"{vts.name}: variable names {names}, expected {expected}")
        formats = sorted({a.get("format") for a in root.iter("DataArray")})
        if formats != [FORMAT]:
            fail.append(f"{vts.name}: data array formats {formats}, expected {FORMAT} for every array")

if fail:
    for line in fail:
        print("  " + line)
    print("VtkNames: FAIL")
    sys.exit(1)
print("VtkNames: ok — field, wall and diagnostic VTK files are well-formed XML with the expected names")
sys.exit(0)
