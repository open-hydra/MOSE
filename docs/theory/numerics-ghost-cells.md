# Ghost Cells at Walls and Extrapolation Faces

This page details how the two layers of ghost cells behind viscous walls and extrapolation faces are filled. The summary is in
[Spatial Discretization](numerics.md#ghost-cells-of-walls-and-extrapolation-faces).

## Ghost Cells of Walls and Extrapolation Faces

Behind a viscous wall (adiabatic or heat-flux, isothermal, gas–surface interaction, and the walls
coupled to another solver) and an extrapolation face the two layers of ghost cells are filled once per stage, before the
fluxes, and the viscous fluxes, the MUSCL reconstruction and the source terms
of the turbulence models read the same values. Extrapolated are the quantities
whose differences the viscous fluxes and the sources take, the mass fractions
$Y_s$, the temperature $T$, the velocity $\mathbf{u}$ and the specific
turbulence and scalar variables $\phi = (\rho\phi)/\rho$, and the pressure
$p$; the ghost state follows from them, $\rho = p/(R(Y)\,T)$,
$\rho_s = Y_s\rho$, $(\rho\phi) = \rho\phi$. With $q_1$, $q_2$ the values in the
wall cell and in the next one across the face:

| Class | Quantities | First ghost cell | Position of its value |
|---|---|---|---|
| wall value | velocity at the coupled, adiabatic or heat-flux and isothermal walls (zero, or the wall velocity of a wall at rest in a rotating frame); the turbulence variables at the walls ($k$, $R_{ij}$, $\tilde\nu$ on a smooth wall: zero; $\omega$, $\tilde\nu$ on a rough wall: the wall value of the model) | $2q_w - q_1$ | $2\mathbf{x}_f - \mathbf{x}_1$ |
| wall temperature | $T$ at isothermal walls | $\max(2T_w - T_1,\ T_w/2)$ | $2\mathbf{x}_f - \mathbf{x}_1$ |
| extrapolated | $Y_s$, $p$; $T$ at coupled, adiabatic or heat-flux and gas–surface walls and at extrapolation faces; the velocity at extrapolation faces and gas–surface walls; turbulence and scalars at extrapolation faces, scalars at the walls | $q_1 + \psi(q_1 - q_2)$ | $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$ |

$\mathbf{x}_f$ is the centre of the face (the mean of its four nodes) and
$s = l_0/l_1$ the ratio of the distances of the MUSCL reconstruction at the
first interior face, $l_0 = (\Delta l_g + \Delta l_1)/2 = \Delta l_1$ (the ghost
cell has the edge across the face of the first cell, see *Gradient Computation
for Diffusive Fluxes*), $l_1 = (\Delta l_1 + \Delta l_2)/2$:
$s = 2\Delta l_1/(\Delta l_1 + \Delta l_2)$, 1 on uniform cells.

$\psi \le s$ is
the largest value that keeps each species, the pressure, the temperature and
each scalar within $[q_1/2,\ 2q_1]$; the species are then normalised to sum
one; the velocity is not limited ($\psi = s$). The wall values of the turbulence
variables are those of the wall fluxes, from the same wall state (isothermal: $T_w$;
adiabatic or heat-flux: the wall temperature estimated from the heat flux;
gas–surface walls, whose surface model sets the wall state at every call, and the
coupled walls, whose viscous flux the other solver gives: the
state of the wall cell, as an adiabatic wall).

With one cell across the face the
extrapolated quantities are copied (position $\mathbf{x}_1$). In an inviscid
run the walls carry no wall values and every quantity is extrapolated, as
behind an extrapolation face. The second layer takes $2q_w - q_2$,
$\max(2T_w - T_2,\ T_w/2)$ and $q_1 + 2\psi_2(q_1 - q_2)$; no flux reads it.

Properties, for every mesh and in 3-D, with the hypotheses written:

- **Positivity.** The ghost state has $\rho_s \ge 0$, $p > 0$, $T > 0$ and
  so $\rho > 0$; the extrapolation amplifies $q_1$ at most by $1 + s$. The
  temperature and the mass fractions are what the viscous fluxes take
  differences of, so they must be those of a gas. The extrapolation of every
  primitive as $3q_1 - 3q_2 + q_3$, which the other types keep (below), can be
  negative (density that rises steeply away from a hot wall) and amplifies by 3
  and 6 on the two layers.
- **Wall value.** For the first two classes without the floor, the value
  interpolated at the face centre, $(G_1 + q_1)/2$, is the wall value, for any
  profile (zero velocity, zero $k$, the model value of $\omega$). The floor
  $T_w/2$ acts only where $T_1 > 1.5\,T_w$, a first cell much hotter than a
  cooled wall that the mesh does not resolve; there the wall flux still takes
  $T_w$, and the ghost cell stays on the cold side of the wall.
- **Linear fields.** The rules are affine combinations of the values whose
  coefficients sum to one, so for a linear field the ghost value is the field at
  the position of the table: for the extrapolated class for every linear field,
  for the other two for the linear fields that take the wall value at the face
  centres (zero velocity on the wall, $T(\mathbf{x}_f) = T_w$ of each face).
  With these positions in the face metrics (stencil) the face gradients next
  to walls and extrapolation faces are exact for those fields, whatever the density:
  the quantities differenced are extrapolated themselves. So is the MUSCL
  reconstruction at the first interior face along a straight grid line: both
  positions are $l_0 = \Delta l_1$ from $\mathbf{x}_1$, the distance it assumes.
  The density gradient of the source terms takes the position of the
  temperature. At adiabatic or heat-flux walls, extrapolation faces and gas–surface walls the pressure, the composition
  and the temperature all sit at $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$,
  so the density of the ghost cell is that of the field there, exact for linear
  $p$, $Y_s$ and $T$ where no limit acts, with no further hypothesis. At isothermal walls
  the temperature sits at $2\mathbf{x}_f - \mathbf{x}_1$, $\Delta l_1$ from
  $\mathbf{x}_1$ along the edge across the face of the wall cell, and the
  pressure and the composition at $\mathbf{x}_1 + s(\mathbf{x}_1 - \mathbf{x}_2)$,
  along the line from $\mathbf{x}_2$ to $\mathbf{x}_1$, $\Delta l_1$ from
  $\mathbf{x}_1$ on a straight grid line: there the two positions are the same
  point and the density is exact; elsewhere they differ by about $\Delta l_1$
  times the angle $\theta$ between the two directions, and the density is exact
  where the pressure and the composition do not vary across the wall cell (a
  boundary layer), its relative error otherwise of the order of
  $\Delta l_1\,\theta\,|\nabla p|/p$, and the same for the composition.
- **Fluxes and sources.** One set of ghost cells: the wall flux routines write
  no ghost values, and the source terms take the positions of the table.

Where a limit acts the ghost value is no longer the value of the field at its
position, and the gradients next to the face are first order there: it happens
only where a quantity changes by more than a factor two between the first two
cells, or where the first cell is 1.5 times hotter than an isothermal wall; the run reports how often
either acted.

Declared limits:

- gas–surface interaction walls: the velocity and the temperature are extrapolated (the surface
  may blow, its state is set at every call of its routine) and the wall values
  of the turbulence variables are those of an adiabatic wall at the state of
  the cell: exact for $k$, $R_{ij}$ and $\tilde\nu$ (zero), while that of
  $\omega$ differs, at uniform pressure, from the value at the state of the
  surface by the factor $\mu(T_1)T_1/(\mu(T_s)T_s)$ (about $(T_1/T_s)^{1.7}$
  with $\mu \propto T^{0.7}$): of order one but bounded, it enters only the
  source terms (the gradient of $\omega$ in the cross-diffusion of SST,
  $\nabla k\cdot\nabla\omega$) and the limited MUSCL reconstruction of the
  first face, while the wall flux keeps the state of the surface; the wall flux
  of the pyrolysis and surface-reaction walls while they ablate takes the value of a wall that blows,
  $\omega_\text{wall} = 800\,\nu_w/y_c^2$ ($y_c$ the wall distance of the centre
  of the wall cell) for SST and Wilcox 2006 with either
  [wall condition of $\omega$](turbulence.md#wall-boundary-conditions), so that with the asymptotic one it is $800/80$ (SST) or
  $800/84.7$ (Wilcox 2006) times that of the ghost cell on top of this factor;
- the coupled walls take, as the gas–surface walls, the wall values of the
  turbulence variables of an adiabatic wall at the state of the cell, with the
  same factor for $\omega$, $T_s$ the temperature of the coupled surface; their
  viscous flux comes from the other solver;
- the temperature of an adiabatic or heat-flux wall is extrapolated (consistent, exact for every
  linear field) rather than set from the heat flux;
- the solid propellant grain as a wall that blows its products in (as an inflow it is a mass-flux inlet)
  and the other boundaries keep the extrapolation
  $3q_1 - 3q_2 + q_3$ of the primitives; in a RANS run the turbulence variables
  of the grain take their wall values instead, $2\phi_w - \phi_l$ at
  $2\mathbf{x}_f - \mathbf{x}_1$ as behind a wall, with the wall values of the
  flux of the grain (the blowing wall while it burns, the wall at the flame
  temperature when it does not).

## References
