# Boundary conditions

Boundary conditions are described in the case's boundary-condition file, normally `INPUT/bc.txt`. In multigrid cases, the file names for the coarser levels use the `bc<n>.txt` suffix. Multi-phase cases add the phase prefix.

ATLAS can generate these files from the mesh and boundary-condition definitions. You can also prepare them manually; the format is record-based and depends on the face type.

## Before editing `bc.txt`

- Check that the file belongs to the mesh you are running. A boundary record refers to a particular block and face.
- Use the face numbering and record format documented in the [complete boundary-condition reference](../user/boundary-conditions.md).
- Keep the boundary data and mesh consistent. A file made for a different block ordering or mesh split may refer to the wrong face.

## Two boundary types worth understanding

**Extrapolation (type `400`) is not a prescribed inflow.** It takes the flux from the state in the boundary cell. If gas flows into the domain, it enters with that same cell state. Use a specified inlet or the ambient boundary when the incoming state is known or the surroundings must be able to drive inflow.

**Type `200` is geometric as well as a boundary setting.** With one cell in the `k` direction it defines an axisymmetric wedge; with multiple cells in `k` it can close a rotationally periodic sector. The mesh geometry must match the selected interpretation.

## Full reference

The [boundary-condition reference](../user/boundary-conditions.md) contains the face numbering, record syntax, all boundary types, required parameters, and examples. Use it when writing or diagnosing individual records.
