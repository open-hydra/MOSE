#!/usr/bin/env python3
"""
Convert a structured SU2 mesh (2D quad elements) to structured Tecplot ASCII format.
The mesh is assumed to be NI x NJ structured, stored as unstructured in SU2.
The element connectivity (type 9 quads) is used to detect NI and NJ automatically.
The output is a K=2 extruded zone (Z=0 and Z=1) to match Tecplot 3D-structured convention.

Usage:
    python3 su2_to_tec.py <input.su2> [output.tec]
"""

import sys
import math


def parse_su2(filename):
    """Parse SU2 file and return (nelem, npoin, elements, points)."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    idx = 0

    # NDIME
    while not lines[idx].strip().startswith('NDIME'):
        idx += 1
    ndime = int(lines[idx].split('=')[1])
    idx += 1

    # NELEM
    while not lines[idx].strip().startswith('NELEM'):
        idx += 1
    nelem = int(lines[idx].split('=')[1])
    idx += 1

    elements = []
    for _ in range(nelem):
        parts = lines[idx].split()
        # type, n0, n1, n2, n3, elem_id
        elem_type = int(parts[0])
        if elem_type == 9:  # quad
            nodes = [int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])]
        else:
            raise ValueError(f"Expected quad elements (type 9), got type {elem_type}")
        elements.append(nodes)
        idx += 1

    # NPOIN
    while not lines[idx].strip().startswith('NPOIN'):
        idx += 1
    npoin = int(lines[idx].split('=')[1])
    idx += 1

    points = [None] * npoin
    for _ in range(npoin):
        parts = lines[idx].split()
        x = float(parts[0])
        y = float(parts[1])
        pt_id = int(parts[2])
        points[pt_id] = (x, y)
        idx += 1

    return nelem, npoin, elements, points


def detect_structured_dims(nelem, npoin):
    """Infer NI x NJ from nelem = (NI-1)*(NJ-1), npoin = NI*NJ."""
    # Try integer factorizations
    for ni in range(2, npoin):
        if npoin % ni == 0:
            nj = npoin // ni
            if (ni - 1) * (nj - 1) == nelem:
                return ni, nj
    raise ValueError("Cannot detect structured NI x NJ dimensions.")


def write_tecplot(filename, points, ni, nj, title="SU2 mesh"):
    """Write structured Tecplot BLOCK format with K=2 extrusion."""
    nk = 2
    ntotal = ni * nj

    # points are ordered: node_id = i + j*ni  (i=0..ni-1 fast, j=0..nj-1 slow)
    # This matches Tecplot BLOCK I-fastest ordering exactly.

    x_vals = [points[k][0] for k in range(ntotal)]
    y_vals = [points[k][1] for k in range(ntotal)]

    with open(filename, 'w') as f:
        f.write(f' TITLE     = "{title}"\n')
        f.write(' VARIABLES = "X", "Y", "Z"\n')
        f.write(f' ZONE T="BLOCCO 1"\n')
        f.write(f' I={ni:3d}, J={nj:3d}, K={nk:3d}, ZONETYPE=Ordered\n')
        f.write(' DATAPACKING=BLOCK\n')

        # X block: k=1 layer then k=2 layer (same x)
        for k in range(nk):
            for val in x_vals:
                f.write(f'{val}\n')

        # Y block: k=1 layer then k=2 layer (same y)
        for k in range(nk):
            for val in y_vals:
                f.write(f'{val}\n')

        # Z block: k=1 -> Z=0, k=2 -> Z=1
        for z in [0.0, 1.0]:
            for _ in range(ntotal):
                f.write(f'{z}\n')

    print(f"Written: {filename}")
    print(f"  Grid: I={ni}, J={nj}, K={nk}")
    print(f"  Total points: {ni * nj * nk}")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.su2> [output.tec]")
        sys.exit(1)

    su2_file = sys.argv[1]
    tec_file = sys.argv[2] if len(sys.argv) > 2 else su2_file.replace('.su2', '.tec')

    print(f"Reading: {su2_file}")
    nelem, npoin, elements, points = parse_su2(su2_file)
    print(f"  NELEM={nelem}, NPOIN={npoin}")

    ni, nj = detect_structured_dims(nelem, npoin)
    print(f"  Detected structured grid: NI={ni}, NJ={nj}")

    # Verify that the first element matches expected node ordering i + j*ni
    # Element (0,0) should have nodes: 0, 1, ni+1, ni
    expected = [0, 1, ni + 1, ni]
    if elements[0] != expected:
        # Try transposed: node_id = j + i*nj
        expected_t = [0, nj, nj + 1, 1]
        if elements[0] == expected_t:
            print(f"  Node ordering: j + i*{nj} (j-fast). Reordering to i-fast for Tecplot.")
            # Reorder: build i-fast array from j-fast storage
            reordered = [None] * npoin
            for i in range(ni):
                for j in range(nj):
                    reordered[i + j * ni] = points[j + i * nj]
            points = reordered
        else:
            print(f"  Warning: first element nodes {elements[0]} don't match expected {expected}.")
            print(f"  Proceeding with natural node ordering (check output carefully).")

    write_tecplot(tec_file, points, ni, nj)


if __name__ == '__main__':
    main()
