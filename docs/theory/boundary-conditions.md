# Boundary Conditions

This page gives the state that the inlet, outlet and thin-direction boundaries impose, and why. The data lines, the fields
and the rules of use are in the [User Guide](../user/boundary-conditions.md). The gas–surface boundaries are in
[Gas–Surface Interaction](gsi.md).

Notation: $p$, $a$ and $w$ are the pressure, the speed of sound and the velocity along the outward normal of the gas of the
boundary cell; $\gamma$ is the ratio of specific heats.

## Outlet (406)

Where the gas of the boundary cell moves into the domain, the face is a slip wall.

Elsewhere the face takes the exit pressure

$$
p_e = \texttt{rel\_fac}\;p_\text{back} + (1 - \texttt{rel\_fac})\,p .
$$

Where the normal Mach number $M_n$ of the cell is below 1 and $p_e > 0$, the face takes $p_e$, raised to

$$
p\,\Bigl[\tfrac{2}{\gamma+1}\Bigl(1 + \tfrac{\gamma-1}{2}M_n^2\Bigr)\Bigr]^{\gamma/(\gamma-1)}
$$

where that is higher. The other quantities of the face are:

- the density of the isentropic expansion of the gas of the cell to that pressure;
- the composition, the turbulence variables per unit mass and the tangential velocity of the cell;
- the normal velocity, from the invariant that reaches the face from inside.

At a supersonic exit, and where $p_e = 0$, every quantity of the face is that of the cell.

## Inlets at total conditions (401, 402, 407)

### When the gas enters

The gas enters when $p_0$ is above

$$
p_\text{stop} = p\,\Bigl(1 + \tfrac{\gamma-1}{2}\,\tfrac{w}{a}\Bigr)^{2\gamma/(\gamma-1)},
$$

the pressure that the gas of the boundary cell reaches when the outgoing acoustic wave brings it to rest:

- a cell at rest takes the inflow as soon as $p_0 > p$;
- a cell moving out takes it when $p_0$ is high enough to turn the flow back;
- a cell moving in at $2a/(\gamma-1)$ or faster, which no pressure brings to rest ($p_\text{stop} = 0$), takes it for every
  $p_0 > 0$.

### The face state

The face state follows from $T_0$, $p_0$, the flow direction and the invariant $a + \tfrac{\gamma-1}{2}\,w$ that reaches the
face from inside, whatever its sign: the wave that links the cell to the face is a rarefaction or a shock, as for slower
gas. The mass flux is continuous where the invariant changes sign.

### Choking

When that state would be supersonic across the face, the invariant can no longer reach it and the face is choked: normal
Mach number 1, with the total temperature and pressure of the reservoir. With the direction normal this is the critical
state $a_0\sqrt{2/(\gamma+1)}$, the largest mass flux the reservoir can give.

The face is choked too when no inflow state exists: when the gas of the cell moves in so fast that, even expanded to
vacuum, it leaves the face faster than the gas of the reservoir can follow. With the direction normal, that is when
$a + \tfrac{\gamma-1}{2}\,w < -a_0\sqrt{(\gamma-1)/2}$.

### When the gas does not enter

- `401` and `402` make the face a slip wall: they let no gas out, and reflect a wave that reaches them with a stagnation
  pressure above $p_0$.
- `407` is then an outlet at $p_\text{back}$. With the gas of the cell at rest or moving out, it takes the outlet
  treatment of `406`. With the gas moving in, the face takes the state that the outgoing wave gives at the pressure
  $\min(p_\text{back}, p_\text{stop})$, where its velocity points out (at $p_\text{stop}$ the gas is at rest at the face).
- A `407` face whose gas moves in at $2a/(\gamma-1)$ or faster is an outlet only with $p_0 \le 0$ or $T_0 \le 0$, a
  reservoir that feeds nothing, and then has no flux.

`rel_fac` has no effect on these three types: when the gas enters, the face state exists with $p_0$ as given.

## Choked nozzle (420)

The face is the exit of a nozzle fed from the reservoir at $T_0$ and $p_0$. The pressure $p$ of the boundary cell selects
its state:

| Pressure of the cell | State of the face |
|----------------------|-------------------|
| $p_\text{sub} \le p < p_0$ | The inflow at the total conditions $T_0$ and $p_0$. |
| $p_\text{shock} \le p < p_\text{sub}$ | The inflow with the mass flux $g$, the total temperature $T_0$ and the invariant that reaches the face from inside: a normal shock stands in the divergent part and the exit is subsonic. |
| $p$ below both | The supersonic design exit state: the isentropic state of the reservoir with the mass flux $g$, and the flux of the Riemann solver between that state and the cell. |

The last row is the exit of a real nozzle whose back pressure is below $p_\text{shock}$: over-expanded down to $p_\text{sup}$,
under-expanded below it.

$p_\text{shock}$ is the pressure behind a normal shock standing at the exit of the design flow: the pressure of the state
of the mass flux $g$, the total enthalpy of the reservoir and the impulse $p + g\,u$ of the design exit state, on the
subsonic side. MOSE computes the design exit state and $p_\text{shock}$ once at set-up from $T_0$, $p_0$, $g$ and the gas,
on the tables of the gas, so $p_\text{sup}$ selects no branch.

## Directions of one cell

A direction of one cell (the $k$ direction of a 2-D case, $j$ and $k$ of a 1-D case) has two faces across it. Let
$u_\perp$ be the velocity across, $\Delta$ the thickness of the cell across and $a$ the speed of sound.

With symmetry (`300`) on both faces, the two mirror planes one cell apart give $u_\perp$ the mode

$$
\frac{du_\perp}{dt} = -\frac{2a}{\Delta}\,u_\perp ,
$$

which the explicit schemes cannot follow at the time step of the cells, about $\text{cfl}\,\Delta x/a$: round-off in
$u_\perp$ would grow until the run fails. So in a planar case MOSE does not integrate $u_\perp$ and sets it to zero, when
every face of the direction is a plane normal to its axis.

Where the faces are not normal to the axis, symmetry asks the gas to move along them, across the cell, and setting
$u_\perp$ to zero would change the flow (the mass flow, the pressure): there $u_\perp$ is integrated, and the mode can stop
the run on a thin cell. With extrapolation (`400`) on a face, the face lets the gas through and $u_\perp$ is integrated.

| Faces across | $u_\perp$ | Mode | Stable step, Heun / RK3 | The run fails for cells thinner across than, Heun / RK3 |
|---|---|---|---|---|
| `300` + `300`, planes normal to the axis | set to zero, not integrated | — | — | — |
| `300` + `300`, planes inclined to the axis | integrated | $-2\,a\,u_\perp/\Delta$ | $\Delta/a$ / $1.26\,\Delta/a$ | $\text{cfl}\,\Delta x$ / $0.8\,\text{cfl}\,\Delta x$ |
| `300` + `400` | integrated | $-a\,u_\perp/\Delta$ | $2\,\Delta/a$ / $2.51\,\Delta/a$ | $\text{cfl}\,\Delta x/2$ / $0.4\,\text{cfl}\,\Delta x$ |
| `400` + `400`, parallel planes | integrated | none: the fluxes of the two faces cancel | — | — |

Cells as wide across as along $x$ are within all these limits up to cfl 1.
