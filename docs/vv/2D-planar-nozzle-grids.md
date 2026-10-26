# Planar Nozzle on Grid Families

A turbulent flow through a planar converging–diverging nozzle, computed with the
Spalart–Allmaras and SST models on grid families that differ in the direction of
their grid lines and in their smoothness. The geometry and the flow are
mirror-symmetric about the centre line, while the four blocks of every grid have
different index orientations: a discretisation that does not depend on the
orientation of the indices keeps the solution mirror-symmetric at every
iteration. The full form runs the families on three levels of refinement and
writes the field of every run; this page gives the result of the short form.

## Problem setup

Geometry, in units of the throat half-height $R_t = 0.01$ m: the half-height is

$$
h(x) = \sqrt{1 + (A_i - 1)\sin^2\!\Big(\frac{\pi x}{2L_c}\Big)},\quad -L_c \le x \le 0,
\qquad
h(x) = \sqrt{1 + (A_e - 1)\sin^2\!\Big(\frac{\pi x}{2L_d}\Big)},\quad 0 \le x \le L_d,
$$

with $A_i = 2.25$, $A_e = 1.8$, $L_d = 4$ and $L_c = L_d\sqrt{(A_i-1)/(A_e-1)} = 5$
(continuous second derivative at the throat): inlet and exit half-heights 1.5 and
1.342, zero wall slope at inlet and exit, isentropic exit Mach number 1.70; one
cell in depth.

| Boundary | Type | Condition |
|---|---|---|
| Inlet, $x = -L_c$ | `407` | ambient inlet at total conditions, $p_0 = 10^5$ Pa, $T_0 = 600$ K, flow normal to the face, freestream turbulence of the model |
| Outlet, $x = L_d$ | `400` | supersonic, extrapolated |
| Walls | `302` | isothermal, $T_w = 300$ K |
| $k$ faces | `300` | symmetry |

The initial field is the quasi-1D isentropic flow (subsonic in the convergent,
supersonic in the divergent).

## Grid families

With $F = y/h$ from $-1$ (lower wall) to $+1$ (upper wall); every level doubles
the cells of the coarse one in both directions:

- **O**: grid lines normal to the walls and to the centre line;
- **H**: vertical grid lines;
- **ZZ**: the O grid with every other row moved in $x$ and every other column in
  $F$ ($\pm0.15$ of the local spacing) and a period-6 triangle in $F$
  ($\pm0.10$): non-smooth at every level, and not symmetric;
- **SB2**: the ZZ grid plus random terms (grid lines that meet the walls at angles
  up to $\pm50$ degrees, kinks, noise), drawn anew on every level, so that no two
  levels belong to one family.

Next to the walls the first two layers of the coarse grid are $10^{-4}\,h$ thick,
then grow by 1.2 per coarse cell up to a cap; O and H have 48 cells from each wall
to the centre line on the coarse grid, ZZ and SB2 44 below and 56 above. The grid
is split into four blocks (convergent and divergent, lower and upper half) with
different index orientations. `make_case.py` in the test folder builds the mesh,
the initial field and the boundary records; `plot_grids.py` draws the four families on
the coarse level, each block in its own colour:

<figure>
  {% include "vv/images/planar-nozzle-grids.svg" %}
</figure>

## Setup

| Parameter | Value |
|---|---|
| Gas | air: perfect gas with $\gamma = 1.4$ and $R = 287$ J/(kg·K) (`test/common/Air-viscous`, the files that ATLAS GPB writes at ATLAS commit `ebeaf4f` for an ideal gas `Air` with `gamma = 1.4`, `R = 287.0`, `mil = 1.1858685985e-5` and `Pr = 0.72`), viscosity of Sutherland's law $\mu = 1.458\times10^{-6}\,T^{3/2}/(T + 110.4)$ Pa·s (`test/common/Air-sutherland/transport.dat`), $\mathrm{Pr} = 0.72$, $\mathrm{Pr}_t = 0.9$ |
| Turbulence | Spalart–Allmaras; SST (`omega-wall-bc = practical`) |
| Time scheme | RK3, local time stepping (steady) |
| CFL / VNN | 0.8 / 0.3 |
| Space reconstruction | MUSCL, van Leer |
| Riemann solver | HLLC |

## Checks

- **Short form** (`./MOSE.sh testfast`, run by CTest): 300 iterations on the
  coarse level from the quasi-1D field.
  1. On O and H, with SA and SST, the largest mirror asymmetry of every variable
     ($\rho$, $u$, $p$, $T$ and the turbulence variables even in $y$, $v$ odd),
     relative to its largest magnitude, is below $10^{-10}$.
  2. On ZZ, with SA, the field is finite, with positive density and pressure.
- **Full form** (`./MOSE.sh -p <threads> test`): O, H, ZZ and SB2 on the three
  levels, SA and SST, up to the residual threshold (at most 200 000 iterations);
  it writes `OUTPUT/field_<family>-L<level>-<model>.tec` and the wall files of each
  run, and `verify.py` applies the two checks above to the solutions of level 0.

## Results

- **Short form** (`PlanarNozzleGrids`, validation tier of CTest): passes both
  checks.

## Running the case

```bash
cd test/2D/viscous/planar-nozzle
./MOSE.sh testfast        # or: ./MOSE.sh -p <threads> test (full form, long)
python verify.py
python plot_grids.py      # the figure of the grid families (no run needed)
```
