# HyShot II Scramjet Combustor

The results are discussed in the documentation, [HyShot II Scramjet Combustor](../../../docs/vv/3D-hyshotII.md); this file holds the setup of the case and how to run it.

## Initial condition

ATLAS ICB builds the initial field, `INPUT/ic.szplt`, from the `[ICB-*]` sections of `input.ini` (see [What ships, and what you build](#what-ships-and-what-you-build)). Every block takes the state of `[ICB-Block*]`, except blocks 6 to 9, the four blocks of the porthole that run down the pipe: their own sections split them along $z$ at the wall plane, and a cell whose centre lies below it takes the state of `[hydroIC]`, any other cell the state of `[airIC]`, equal to that of `[ICB-Block*]`.

| Region | Gas | Velocity | $p$ | $T$ | $\tilde\nu$ |
|---|---|---|---|---|---|
| Every cell above the wall plane | air, $Y_{\mathrm{N_2}} = 0.7686$, $Y_{\mathrm{O_2}} = 0.2314$ | $M = 2.49$ along $x$ | 130.2 kPa | 1377 K | $10^{-6}$ m$^2$/s |
| The pipe, below the wall plane | pure H$_2$ | $w = 1204.2$ m/s, into the duct | 252.8 kPa | 248.4 K | $10^{-6}$ m$^2$/s |

The hydrogen moves at its speed of sound with the mass flux of the `403` inflow of the pipe: with the tables of the case, 1204.2 m/s is the speed of sound of hydrogen at 248.4 K and $\rho w = 297.15$ kg/(m$^2$ s) is the mass flux $g = 297.134$ kg/(m$^2$ s) of that inflow, each to within 0.01 %, and the stagnation temperature of this state is 299.5 K, against $T_0 = 300$ K of the inflow. The air is uniform: it does not follow the profile of the inlet, and next to the 300 K walls it starts at 1377 K and $M = 2.49$ too.

## Numerical setup

### Grid

Three independent, overlapping structured regions - a stretched upstream box, the homogeneous chamber background, and a refined O-grid patch around the porthole with its pipe - coupled by **chimera overlap**, with no conformal interface between the regions. The overlap connectivity (donor/overset, type 102) is computed by the ATLAS pre-processor from the geometry, so the regions can be re-meshed independently.

| Parameter | Value |
|---|---|
| Blocks | 11 |
| Cells | 4,563,968 |
| Boundary faces | 414,752 |
| of which chimera / connection / symmetry / wall / inlet-outlet | 28,768 / 59,712 / 224,224 / 95,008 / 7,040 |
| Chamber spacing $\Delta x = \Delta y$ | 0.25 mm |
| First wall cell | 1 $\mu$m, growth 1.2, 28 cells (mirrored at the ceiling) |
| Wall resolution, bottom wall | $y^+$ median 1.3, 95th percentile 1.9, max 2.1 |

### Physics and numerics

| Parameter | Value |
|---|---|
| Equations | Navier-Stokes |
| Equation of state | ideal gas, thermally perfect (variable $c_p(T)$) |
| Turbulence model | Spalart-Allmaras |
| Species diffusion | constant laminar Schmidt number, $Sc = 0.7$ |
| Turbulent Schmidt number | $Sc_t = 0.7$ |
| Chemistry | finite-rate, Gerlinger-9 (9 species, 19 reactions) |
| Chemistry-flow coupling | Strang operator splitting |
| ODE solver | `H-radau5` |
| Riemann solver | HLLC (Batten) |
| Space reconstruction | MUSCL, Van Leer limiter |
| Time scheme | RK2, conservative variables, steady (local time stepping) |
| Implicit residual smoothing | $\beta = 0.3$ |
| CFL / VNN | 0.8 / 0.5 |
| Multigrid | 3 levels, 100,000 (L3) + 20,000 (L2) + 20,000 (L1) iterations |

### Chemistry tolerances

The ODE error control is set explicitly and tightly:

```ini
[MOSE-Chemistry]
ode-absolute-tol-species     = 1e-12 (x9)
ode-relative-tol-species     = 1e-8  (x9)
ode-absolute-tol-temperature = 1e-4
ode-relative-tol-temperature = 1e-8
ode-max-steps                = 100000
```

This is not a convergence detail. MOSE integrates the chemistry on **partial densities** $\rho_s$ [kg/m$^3$], and the registered default absolute tolerance sits one to two decades *above* the entire pre-ignition radical pool. At the default the induction zone is integrated as if it were zero, heat release is displaced downstream, and the case shows no OH anywhere near the porthole - while the flow field (separation size, reversed-flow fraction, peak backflow velocity) is identical to a fraction of a percent. Tightening the tolerance alone raises peak $Y_{\mathrm{OH}}$ in the porthole region by two to three orders of magnitude.

Getting ignition to occur **ahead** of the jet, as the literature reports, additionally requires the profiled inflow: the uniform inlet delivers near-wall gas about 145 K colder, which weakens the upstream separation and roughly quarters the amount of hot flammable gas there. Both changes are necessary; neither is sufficient.

## Running the case

### What ships, and what you build

Everything the case needs is committed **except the initial condition**. At 4.6 M cells an initial field is about 375 MB even as a single-precision `.szplt`, and 2.2 GB as Tecplot ASCII, so it is not something a repository can carry.

| Committed | Build yourself |
|---|---|
| `MESH/mesh.szplt` - the overset grid, and `MESH/hyshot_mesh.py` that generates it | `INPUT/ic.szplt` - the initial field |
| `INPUT/phase.txt`, `thermo.dat`, `transport.dat`, `chemistry-*` - the 9-species Gerlinger mixture | |
| `INPUT/bc.txt`, `bc2.txt`, `bc3.txt` - boundary conditions and chimera connectivity, one file per multigrid level | |
| `input.ini` - the solver setup documented above | |

So the first step is to build the initial condition with the ATLAS pre-processor, from the shipped mesh and phase data. ICB looks for the mesh (`mesh.tec`, then `mesh.p3d`, then `mesh.szplt`), the phase file and its thermodynamic table **in the working directory**, and writes into `fromATLAStoSolver/`:

```bash
cd test/3D/hyshotII
ln -s MESH/mesh.szplt .          # ICB searches the cwd, not MESH/
ln -s INPUT/phase.txt .          # and needs the phase file beside it
ln -s INPUT/thermo.dat .         # and its thermodynamic table
ulimit -s unlimited
ATLAS ICB                        # -> fromATLAStoSolver/ic.szplt
mv fromATLAStoSolver/ic.szplt INPUT/
```

Then solve:

```bash
./MOSE.sh compile              # or copy an existing bin/MOSE
./MOSE.sh -p 60 solve          # 60 OpenMP threads (the option goes before the task); ~30 h for the full MG schedule
```

The boundary-condition files are committed, so ATLAS BCB does not need to run. They already encode everything in the [boundary conditions](../../../docs/vv/3D-hyshotII.md#boundary-conditions) table - the profiled supersonic inlet, the sonic hydrogen porthole, the 300 K isothermal walls and the chimera overlap - which is why the physical settings are documented on this page rather than being readable from `input.ini`.

The wall-pressure figure is drawn from `reference/wall.tec`, a MOSE wall solution of this case, which `wall_pressure.py` reads in its own folder:

```bash
cd test/3D/hyshotII/reference
python wall_pressure.py                  # -> docs/vv/images/HyShot-pw.svg
python wall_pressure.py --preset paper   # -> wall_pressure.pdf / .png
```

`wall_pressure.py` stitches the between-injector line across whatever bottom-wall zones `wall.tec` contains, accepting only rows that genuinely sample $y = W$, so it survives a re-mesh without editing. It also re-audits the thesis-to-paper coordinate and pressure-normalisation transform on every run, and prints the worst disagreement between the two digitisations of the same measurements.
