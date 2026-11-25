# Manifold boundary (`501`)

A manifold represents a lumped plenum that connects one or more outlet faces to one or more inlet faces. Gas leaving the outlet side is mixed in the plenum and returned through the inlet side.

The data line is:

```text
gid, role, K, rel_fac, mdot_ratio, p0_ratio, dh0
```

- `gid` identifies the manifold group.
- `role = 1` marks an outlet face; `role = 2` marks an inlet face.
- All records in a group share the same parameter values.

At each stage, MOSE gathers the outgoing mass flow, total enthalpy, dynamic pressure, species, and turbulence quantities. It uses these to construct the inlet state and the outlet pressure for the connected group.

With defaults (`K = 0`, `mdot_ratio = 1`, `p0_ratio = 1`, `dh0 = 0`), the plenum is loss-free and perfectly mixed: mass, total enthalpy, and total pressure are passed through. `K` models a loss on outlet dynamic pressure; `mdot_ratio`, `p0_ratio`, and `dh0` can represent a simplified network element, bleed, pressure ratio, or specific work.

Every group must contain at least one inlet and one outlet. MOSE rejects invalid role values, negative loss coefficient, non-positive ratios, a relaxation factor outside `(0, 1]`, or inconsistent parameters within a group.

The run prints one manifold summary per group at `shell-diter` intervals. It reports the inlet/outlet mass flow, total enthalpy, total pressure, and the work term.

For the exact averaging definitions, flow-reversal behaviour, and diagnostic-line format, see [BC types in the complete reference](../boundary-conditions.md#bc-types).
