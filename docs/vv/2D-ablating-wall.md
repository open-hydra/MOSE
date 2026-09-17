# 2D Ablating Wall – Blown Column

Analytic verification of the gas–surface interaction wall boundary conditions — melting, pyrolysis and solid propellant grain combustion — plus the fallback that closes an inert surface as an impermeable wall.

A thin gas column is blown from one end by an ablating wall and vented at the
other.  Initialised isothermal at the wall temperature and **at rest**, this
configuration is an exact steady solution of the governing equations, and the
wall state it settles on is available in closed form.

---

## Problem setup

```
   face 4 (j = ny)   outflow (extrapolation)
 +---------------------+
 |                     |
 |     gas column      |   faces 1, 2 : symmetry  ->  the solution is 1-D in y
 |                     |
 +---------------------+
   face 3 (j = 1)    ablating wall
```

The wall blows gas into the column; it leaves through the vent at the opposite
end.  The gas is a single-species calorically perfect air
($c_p = 1004.5$ J/(kg·K)), so the enthalpy the injected mass carries is exactly
$c_p T_w$ and the analysis below is free of table interpolation error.

| Parameter | Value |
|---|---|
| Domain | 2 × 6 mm, 4 × 12 cells |
| Pressure | 101325 Pa |
| Initial state | isothermal at $T_w$, at rest |
| Gas | single-species air, $R$ = 287.0 J/(kg·K) |

## Numerical setup

| Parameter | Value |
|---|---|
| Equations | Navier–Stokes (laminar) |
| Time scheme | RK2, local time stepping |
| CFL / VNN | 0.5 / 0.5 |
| Space reconstruction | MUSCL, Van Leer limiter |
| Riemann solver | HLLC |
| Iterations | 4000 |

---

## Analytic solution

With the column isothermal at the wall temperature the wall-normal temperature
gradient vanishes, so $q_\mathrm{cond} = 0$ and the surface energy balance
collapses to a single scalar equation.

**Melting.**  The wall temperature is prescribed, so the mass flux
follows directly:

$$
\dot m = \frac{q_\mathrm{rad}}{\Delta h_m + c_s\,(T_w - T_i)}
$$

**Pyrolysis.**  The mass flux obeys an Arrhenius law in $T_w$, and the
wall temperature is the root of

$$
q_\mathrm{abl}(T_w) \;=\; \dot m(T_w)\,\bigl[\Delta h_p + c_s (T_w - T_i)\bigr] \;=\; q_\mathrm{rad}
$$

$q_\mathrm{abl}$ increases monotonically with $T_w$, so the root is unique.

**The field.**  With $p$ and $T$ uniform the density is uniform, so continuity
fixes the blowing velocity at $v = \dot m / \rho$ everywhere in the column.  The
isothermal blown column is therefore an exact steady solution: the test checks
both that the solver *reproduces* it and that it *preserves* it, the latter
being sensitive to any sign error in the surface bookkeeping.

---

## Results

Wall quantities are reported against the closed form:

| Case | Model | Face | $T_w$ [K] | error | $\dot m$ [kg/m²s] | error |
|---|---|:---:|---|---|---|---|
| `melting` | melting, $q_\mathrm{rad}$ = 2 × 10⁵ | 3 | 800.000 | 0 (prescribed) | 0.1744523 | 1.1 × 10⁻⁶ |
| `htpb` | HTPB, $q_\mathrm{rad}$ = 10⁶ | 3 | 816.41786 | 2.3 × 10⁻⁷ | 0.5139235 | 7.1 × 10⁻⁷ |
| `hdpe` | HDPE, $q_\mathrm{rad}$ = 10⁶ | 3 | 910.05829 | 1.8 × 10⁻⁷ | 0.2866909 | 2.9 × 10⁻⁶ |
| `pp` | PP, $q_\mathrm{rad}$ = 10⁶ | 3 | 824.07440 | 1.7 × 10⁻⁷ | 0.2961751 | 2.7 × 10⁻⁶ |
| `melting-top` | melting, mirrored | 4 | 800.000 | 0 (prescribed) | 0.1744523 | 1.1 × 10⁻⁶ |
| `hdpe-top` | HDPE, mirrored | 4 | 910.05829 | 1.8 × 10⁻⁷ | 0.2866909 | 2.9 × 10⁻⁶ |

The field reproduces the exact solution to the same order: the column stays
isothermal to $3$–$9 \times 10^{-4}$ in $|T - T_w|/T_w$, and the mass flux
$\rho v$ matches $\dot m$ to $3$–$9 \times 10^{-4}$ through **every** cell, the
wall-adjacent one included.

### Solid propellant grain

The `srm` configuration replaces the ablating wall with a burning propellant grain
(BC 502).  The rate comes from Saint-Robert's law rather than from an energy
balance, and the wall is held at the flame temperature, but the blown column is the
same exact solution.  BC 502 is not a viscous wall, so it writes no `wall.tec`;
everything is checked on the field, which is the stronger test here.

| Quantity | Reference | Result |
|---|---|---|
| Column temperature | $T_{\mathrm{af}}$ = 1500 K | 4.8 × 10⁻⁴ relative |
| Mass flux $\rho v$ | $a\,(p/p_{\mathrm{ref}})^n$ = 0.2001 kg/m²s | 4.7 × 10⁻⁴ relative |
| Blowing speed | $\dot m/\rho$ = 0.8502 m/s | 1.7 × 10⁻³ relative |

The temperature check is what pins the enthalpy of the injected products: if the
grain injected its combustion products cold, the column would settle near 485 K
instead of 1500 K, and the blowing speed would fall by a factor of four.

### Inert surface

The `melting-inert` configuration removes the radiative input and starts the gas
at 400 K against an 800 K melting wall.  The surface now *loses* heat to the gas,
the melting energy balance returns a non-positive mass flux, and the boundary
condition must fall back to an impermeable wall:

| Quantity | Result |
|---|---|
| Injected mass flux | exactly 0 |
| Wall temperature | 800 K (the prescribed melting point) |
| Wall heat flux | −1.17 × 10⁴ W/m² — the wall heats the gas |
| Flow field | conduction layer developing, only weak expansion venting |

---

## Running the test

```bash
cd test/2D/viscous/ablating-wall
./run.sh                       # all configurations, as CTest runs them
CONFIGS=hdpe ./run.sh          # a single configuration

python3 build_case.py hdpe     # generate mesh / IC / BC by hand
./MOSE.sh solve
python3 verify.py
```

The case is registered in the fast tier as `AblatingWall` (~8 s):

```bash
ctest -R AblatingWall --output-on-failure
```

Mesh, initial condition and boundary-condition file are generated by
`build_case.py`, so nothing but the scripts is stored in the repository.