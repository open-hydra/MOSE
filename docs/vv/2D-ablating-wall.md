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
| Iterations | 4000 (`melting-inert`: 8000) |

The run is a pseudo-time march (local time stepping) that is stopped before full convergence. The
inert-surface check (no blowing, $|\rho v| < 5\times10^{-2}$ kg m⁻² s⁻¹) looks at the start-up flow
while it decays, so `melting-inert` runs 8000 iterations, in its own directory: the largest
$|\rho v|$ is 5.450e-2 after 4000 iterations and 4.455e-2 after 8000, when
the largest gas temperature is 530.93 K. The other configurations pass their checks with 4000.

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

Among its checks, `verify.py` requires that, in every blown configuration,
$\dot m$ and the $T_w$ of the pyrolysis walls are within $10^{-3}$ (relative) of
the closed form, the melting wall is at its prescribed $T_w$, the wall state is
uniform along the wall within $10^{-8}$, and the column stays isothermal, the
pressure uniform and $|\rho v|$ equal to $\dot m$ in every cell, within
$5 \times 10^{-3}$, and that `srm` passes the same field checks against
Saint-Robert's law.

Wall quantities after the 4000 iterations of the case, against the closed form:

| Case | Model | Face | $T_w$ [K] | error | $\dot m$ [kg/m²s] | error |
|---|---|:---:|---|---|---|---|
| `melting` | melting, $q_\mathrm{rad}$ = 2 × 10⁵ | 3 | 800.000 | 0 (prescribed) | 0.1744567 | 2.4 × 10⁻⁵ |
| `htpb` | HTPB, $q_\mathrm{rad}$ = 10⁶ | 3 | 816.42127 | 7.0 × 10⁻⁷ | 0.5139212 | 2.1 × 10⁻⁶ |
| `hdpe` | HDPE, $q_\mathrm{rad}$ = 10⁶ | 3 | 910.06371 | 1.9 × 10⁻⁷ | 0.2866921 | 3.1 × 10⁻⁶ |
| `pp` | PP, $q_\mathrm{rad}$ = 10⁶ | 3 | 824.07923 | 1.8 × 10⁻⁷ | 0.2961761 | 3.0 × 10⁻⁶ |
| `melting-top` | melting, mirrored | 4 | 800.000 | 0 (prescribed) | 0.1744567 | 2.4 × 10⁻⁵ |
| `hdpe-top` | HDPE, mirrored | 4 | 910.06371 | 1.9 × 10⁻⁷ | 0.2866921 | 3.1 × 10⁻⁶ |

The field reproduces the exact solution: the column stays
isothermal to $3$–$9 \times 10^{-4}$ in $|T - T_w|/T_w$, and the mass flux
$\rho v$ matches $\dot m$ to $3$–$9 \times 10^{-4}$ through **every** cell, the
wall-adjacent one included.

### Solid propellant grain

The `srm` configurations replace the ablating wall with a burning propellant grain
(BC 502).  The rate comes from Saint-Robert's law rather than from an energy
balance, and the wall is held at the flame temperature, but the blown column is the
same exact solution.  The grain runs in both modes of its record: a wall that blows
the products in (mode 1: `srm`, `srm-top`) and an inflow, the mass-flux inlet 404 at
the total temperature $T_{\mathrm{af}}$ (mode 2: `srm-inflow`, `srm-inflow-top`),
each on the low and on the high end of the column (faces 3 and 4: the two orientations
of the normal give the same numbers).  The gas starts at rest, which an inflow grain
must let the products into.  BC 502 is neither a viscous wall nor writes one, so it
writes no `wall.tec`; everything is checked on the field, which is the stronger test
here.

| Quantity | Reference | Mode 1 | Mode 2 |
|---|---|---|---|
| Column temperature | $T_{\mathrm{af}}$ = 1500 K | 4.6 × 10⁻⁴ relative | 4.0 × 10⁻⁴ relative |
| Mass flux $\rho v$ | $a\,(p/p_{\mathrm{ref}})^n$ = 0.2001 kg/m²s | 4.7 × 10⁻⁴ relative | 4.0 × 10⁻⁴ relative |
| Blowing speed | $\dot m/\rho$ = 0.8502 m/s | 1.6 × 10⁻³ relative | 1.4 × 10⁻³ relative |

The temperature check is what pins the enthalpy of the injected products: if the
grain injected its combustion products cold, the column would settle near 485 K
instead of 1500 K, and the blowing speed would fall by a factor of four.

### Surface reactions of a carbon wall (505)

The `carbon` configurations run the surface reactions of `GSI_Carbon_Bradley` on the gas `test/common/CHO-N2-viscous`
(O2, O, OH, H, H2, H2O, CO, CO2 and N2), in `gsi/<config>/`, 1000 iterations:

| Case | Gas (mass fractions) | $T_0$ [K] | $q_\mathrm{rad}$ [W/m²] | $\dot m(T_0)$ [kg/m²s] | Face |
|---|---|---|---|---|:---:|
| `carbon` | N2 0.8, H2O 0.1, CO2 0.1 | 1800 | 11445.13 (endothermic) | 8.354923 × 10⁻⁴ | 3 |
| `carbon-top` | the same | 1800 | 11445.13 | 8.354923 × 10⁻⁴ | 4 |
| `carbon-o2` | N2 0.98, O2 0.02 | 1500 | −5488.90 (exothermic) | 7.166338 × 10⁻⁴ | 3 |

$q_\mathrm{rad}$ is the heat the reactions absorb at the initial state, so the balance starts with its root at $T_0$.
The reactions consume the oxidisers of the wall cell and the blown products displace them: unlike the columns above,
this one has **no steady state**. `verify.py` checks the relations of the wall model at the state the solver reached,
with the wall temperature and pressure of `wall.tec` and the partial densities of the wall cell: the rates of an
independent copy of the model (`gsi_models.py`) give the reported $\dot m$ within $10^{-11}$ (the 15 printed digits,
times the sensitivity $E/\mathcal R T_w \le 24$ of the rates), the surface energy balance
$\sum_s \omega_s h_s(T_w) = q_\mathrm{rad} + q_\mathrm{cond}$ holds within $10^{-11}$, and its root in $T_w$ is the
reported one; the blowing is directed away from the wall and the pressure uniform. The values at the end of the run are
in the output of `verify.py` (the drift from $T_0$ and $\dot m(T_0)$ is printed there).

### Pyrolysing ablator (506)

The `ablator` configurations run 506 on the column of `carbon`, in `gsi/<config>/`. `build_case.py` writes the table
`INPUT/pyrolysis.dat` (`gsi_models.write_pyrolysis_table`: the water-gas shift of 1 CO + 1 H2O + 0.5 H2 at six states on
2 pressures × 3 temperatures, every row with the elements of that mixture):

| Case | $\phi$ | $h_{v,0}$ [J/kg] | $q_\mathrm{rad}$ [W/m²] | Face | Check |
|---|---|---|---|:---:|---|
| `ablator` | 0.25 | −8.44 × 10⁵ | 11345.04 | 3 | the relations of the wall model at the state reached, $\dot m = (1+\phi)\dot m_{ch}$ with the table read at $(p_w, T_w)$ |
| `ablator-top` | 0.25 | −8.44 × 10⁵ | 11345.04 | 4 | the same |
| `ablator-phi0` | 0 | 0 | 11445.13 | 3 | `field.tec` and `wall.tec` byte for byte those of `carbon` (one thread) |
| `ablator-hv0` | 0 | −8.44 × 10⁵ | 12150.29 | 3 | the relations; its wall file differs from that of `carbon` (the control of the comparison above) |

$q_\mathrm{rad}$ moves by $-\dot m_{ch}(T_0)\,h_{v,0} = 705.16$ W/m² from `carbon` to `ablator-hv0`: the closed form of the
virgin-material term at the initial state. Three more configurations must stop at set-up with their message:
`ablator-nofile` (no table file), `ablator-pyro2` (`pyro_id = 2`, one table), `ablator-cond` (a table with a
condensed-carbon column). The element balance of the blown gas and the pressure of the lookup are checked by the unit
check `GSIAblator`.

### Inert surface

The `melting-inert` configuration removes the radiative input and starts the gas
at 400 K against an 800 K melting wall.  The surface now *loses* heat to the gas,
the melting energy balance returns a non-positive mass flux, and the boundary
condition must fall back to an impermeable wall:

| Quantity | Result |
|---|---|
| Injected mass flux | exactly 0 |
| Wall temperature | 800 K (the prescribed melting point) |
| Wall heat flux | +1.78 × 10⁴ W/m² into the gas (after the 8000 iterations of `melting-inert`) — the wall heats the gas |
| Flow field | conduction layer developing, only weak expansion venting |

### Default Schmidt number

The `default-schmidt` configuration runs the `melting` column again, in its own directory, with no `schmidt` in
`input.ini`.  The default, `schmidt = 0`, selects mixture-averaged multicomponent diffusion, which needs the binary
diffusion table `INPUT/diffusion.dat` only when the gas has more than one species (see
[thermodynamics](../theory/thermo.md)): a single species does not diffuse, so the run needs no table and gives
the `melting` column, and every check of `verify.py` prints the same numbers as for `melting` in the table above.

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

The case is registered in the fast tier as `AblatingWall` (about 11 s of CPU time at one thread; see
`test/TESTING-PLAN.md`, section 4):

```bash
ctest -R AblatingWall --output-on-failure
```

Mesh, initial condition and boundary-condition file are generated by
`build_case.py`, so nothing but the scripts is stored in the repository.