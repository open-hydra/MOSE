# Reynolds-stress model: decay of isotropic turbulence

Verification case for the SSG/LRR-omega Reynolds-stress model (`turbulence =
SSGLRR` and its variants), against an analytic solution. The flat plate of the
NASA TMR (`../flat-plate-turbulent-tmr`, `FlatPlateTMR`) runs the model too,
against CFL3D and TAU.

## Setup

A uniform stream (Mach 0.2) through a sector of an axisymmetric duct, far from
any wall:

```
   face 4 (j = ny)   symmetry (slip)
 +---------------------------------+
 |                                 |
 |  inflow  ->  uniform stream  -> |  outflow      faces 5, 6 : axisymmetric
 |                                 |
 +---------------------------------+
   face 3 (j = 1)    symmetry (slip)
```

The inlet carries isotropic turbulence, `R11 = R22 = R33 = 10 m2/s2`
(`k` = 0.3% of `U^2`) and `omega = 600 1/s`.

## Why it has an analytic solution

With no mean strain there is no production, and with no anisotropy the
pressure-strain model vanishes: the model reduces to

```
U dk/dx     = - C_mu k omega
U domega/dx = - beta omega^2
```

with the solution

```
omega(x) = omega_0 / (1 + beta omega_0 x/U)
k(x)     = k_0 (1 + beta omega_0 x/U)^(-C_mu/beta)
```

Far from walls `F1 = 0`, so `beta` is that of the SSG set (0.0828). Over the
duct `k` falls to 56% of its inlet value; the turbulent diffusion along `x` is
three orders of magnitude below the convection.

## What is checked

| check | what it exercises |
|---|---|
| `k(x)`, `omega(x)` follow the closed form | convection and source terms of the seven equations, inlet values (specific, times the inlet density) |
| `R11 = R22 = R33`, no shear stress | the stress equations stay consistent with each other |
| no radial velocity, `k` uniform across the duct | the azimuthal Reynolds stress on the two side faces of the wedge, which must balance the radial one as the pressure does; the mirror image of the stresses on the slip walls |
| `H + k + R11` constant along `x`, `dH/dk = -5/3` | the energy the mean flow exchanges with the turbulence, and the Reynolds stress on the inflow and outflow faces |

`H` is the total enthalpy of the mean flow. The Reynolds stresses do work on
the mean flow as a pressure does, so the conserved sum carries `R11` next to
`k`. Without the energy exchange `dH/dk` would be -2/3; without the Reynolds
stress on the open faces the pressure jumps by `rho*R11` in the two boundary
cells.

The two cells next to the inlet and the one next to the outlet are first-order
(their ghost cells are zero-gradient copies): the closed form is checked to
0.2% on the inner cells and to 2% on all of them.

## Running

```sh
./run.sh                           # the four variants, as CTest runs it
MODELS=SSGLRR-SD ./run.sh          # just one

python3 build_case.py SSGLRR       # generate deck / IC / BC by hand
./MOSE.sh solve
python3 verify.py
```

## Variants

`SSGLRR`, `SSGLRR-SD` (simple diffusion), `SSGLRR-w2019` (length-scale
correction) and `SSGLRR-w2019-SD` (both) have the same solution here: the
diffusion is negligible and the correction acts near walls only. Running the
four checks that each is
selectable and goes through its own diffusive-flux and source routines.
