# SWBLI — Schülein shock-wave / boundary-layer interaction

Mach-5 turbulent boundary layer on an isothermal flat plate, hit by an oblique
shock from a generator on the upper wall. The interaction separates the boundary
layer and it reattaches downstream; the metric is the wall skin friction through
separation and reattachment.

The case is run with two turbulence models (SA, SST), each in its own folder, and
each MOSE run is paired with an **OpenFOAM (`rhoCentralFoam`) run on the same grid**
for code-to-code verification.

```
swbli/
├── compare.py          your run vs the references; --docs: the V&V figures
├── plot_sutherland.py  constant-μ vs Sutherland MOSE Cf (analysis helper)
├── SA/                 Spalart–Allmaras
│   ├── input.ini       MOSE case (ATLAS + GPB + GRIB sections)
│   ├── MOSE.sh         ./MOSE.sh solve
│   ├── INPUT/          bc.txt + multigrid levels, ic.szplt, transport tables
│   ├── MESH/script.geo gmsh mesh definition (GRIB builds the mesh from it)
│   ├── reference/      reference wall solutions: MOSE-SA.tec, OF-SA.xy,
│   │                   wind-SA.dat, SU2-SA.dat, schulein.dat
│   ├── OUTPUT/         your MOSE run (untracked)
│   └── OPENFOAM/       rhoCentralFoam companion case (run inputs)
└── SST/                k-omega SST — same layout, plus
    ├── reference/      three MOSE and four OpenFOAM SST wall solutions (below)
    └── OPENFOAM/
        ├── variant.sh      switch the case between the four OpenFOAM SST variants
        ├── kOmegaSSTMOSE/  OpenFOAM model with MOSE's default production (μₜS²)
        └── system/omegaFaceBC  MOSE's ω wall condition for OpenFOAM
```

Tracked: everything needed to re-run, and the reference wall solutions in
`reference/`. Untracked: your runs — MOSE `OUTPUT/`, and on the OpenFOAM side the
mesh (`mesh.msh`, `constant/polyMesh/`), time directories, `postProcessing/` and
your own `bottomWall.xy` sample.

## Flow conditions

| Quantity | Value |
|---|---|
| Mach | 5.0 |
| p∞ | 4000 Pa |
| T∞ | 68.3 K |
| U∞ | 828.29 m/s |
| ρ∞ | 0.2041 kg/m³ |
| μ | Sutherland (As = 1.458e-6, Ts = 110.4) → μ∞ = 4.605e-6 Pa·s |
| unit Reynolds | 36.7 × 10⁶ /m (= NASA NPARC/Wind-US reference) |
| Pr / Prt | 0.69 (Eucken) / 0.85 |
| Walls | isothermal, T = 300 K |
| Wall resolution | y⁺ ≈ 0.3 (wall-resolved, no wall functions) |

Dynamic pressure q∞ = ½ρ∞U∞² = 70008 Pa.

## Running

**MOSE**

```sh
cd SA            # or SST
./MOSE.sh solve  # -> OUTPUT/wall.tec, OUTPUT/field.tec
```

**OpenFOAM**

```sh
cd SA/OPENFOAM
gmsh -3 mesh.geo -o mesh.msh          # build the mesh
./Allrun.pre                          # gmshToFoam + patch types + checkMesh
# then run rhoCentralFoam with LTS (see run.slurm)
```

The SST case compiles its ω wall condition (`system/omegaFaceBC`) on the first
time step, so the compute nodes need the OpenFOAM build environment (`wmake`).
`gmsh` is not needed there: build `mesh.msh` wherever gmsh is available and copy it.

The committed SST case is the `asymptotic compressible` variant. To run another
one, switch it before submitting:

```sh
cd SST/OPENFOAM
./variant.sh wallfunction incompressible   # <asymptotic|wallfunction> <compressible|incompressible>
```

The `incompressible` variants use `kOmegaSSTMOSE`, which `variant.sh` builds into
`$FOAM_USER_LIBBIN` the first time (`wmake libso kOmegaSSTMOSE`); the compute
nodes must see that directory. `./variant.sh asymptotic compressible` restores the
committed case exactly.

Sample the OpenFOAM wall from the **raw patch faces** (not an interpolated line,
which smooths away the cell-to-cell noise):

```sh
postProcess -func wallShearStress -latestTime
postProcess -func "patchSurface(patch=bottomWall, fields=(wallShearStress p), \
    interpolate=false, surfaceFormat=raw)" -latestTime
cp postProcessing/patchSurface*/*/patch.xy  bottomWall.xy
```

## Comparing your run with the references

```sh
python compare.py --model SA     # or SST
```

reads your `OUTPUT/wall.tec` and compares it with the MOSE and OpenFOAM reference
solutions **for the same setup**: for SST, `omega-wall-bc` and `sst-production`
are read from `SST/input.ini` and select `reference/MOSE-SST-<bc>-<production>.tec`
and `reference/OF-SST-<bc>-<production>.xy` (where they exist). It prints
separation, reattachment and bubble length, and the mean C_f and p_w differences
upstream, after reattachment and downstream, and writes `OUTPUT/compare-cf.svg`
and `OUTPUT/compare-pw.svg`. A run of the committed `input.ini` should reproduce
its MOSE reference to within convergence.

```sh
python compare.py --model SST --docs   # regenerate docs/vv/images from reference/ only
python compare.py --model SST --fields # + BL and shock metrics of the V&V page
```

`--docs` needs no run: the V&V figures are drawn from the tracked reference files.
`--fields` needs your `OUTPUT/field.tec` and an OpenFOAM run (latest time directory
with `p`, `U`, `rho`, and `constant/polyMesh`) in `OPENFOAM/`. OpenFOAM data are
read directly from the case files, so no OpenFOAM install is needed to
post-process.

## Reference solutions

All references are converged steady states on the 237,440-cell grid: separation
and reattachment fixed, and between saves the wall C_f changing by < 0.1% (MOSE)
or by 1–2% cell-to-cell without trend (OpenFOAM).

**MOSE** — `reference/MOSE-<model>[-<omega-wall-bc>-<sst-production>].tec`, the
`OUTPUT/wall.tec` of `input.ini` as committed, or with the one line given changed:

| File | Setup | Iterations (fine grid) |
|---|---|---|
| `SA/reference/MOSE-SA.tec` | `SA/input.ini` | 100k |
| `SST/reference/MOSE-SST-asymptotic-compressible.tec` | `SST/input.ini` (base case) | 200k |
| `SST/reference/MOSE-SST-asymptotic-incompressible.tec` | `sst-production = incompressible` | 200k |
| `SST/reference/MOSE-SST-practical-compressible.tec` | `omega-wall-bc = practical` | 200k |

**OpenFOAM** — `reference/OF-<model>[-<wall>-<production>].xy`, raw bottom-wall
samples of the companion case. For SST the two parts of the name are those of
`OPENFOAM/variant.sh`: `asymptotic` = MOSE's face condition (`system/omegaFaceBC`),
`wallfunction` = OpenFOAM's `omegaWallFunction` (same value fixed in the
wall-adjacent cell); `compressible` = `kOmegaSST`, `incompressible` =
`kOmegaSSTMOSE` (P = μₜS²). Same name, same model as the MOSE file:

| File | `variant.sh` | How it was run |
|---|---|---|
| `SA/reference/OF-SA.xy` | — (`SA/OPENFOAM`) | 300k from the initial condition |
| `SST/reference/OF-SST-asymptotic-compressible.xy` | `asymptotic compressible` (committed) | 300k from the initial condition |
| `SST/reference/OF-SST-wallfunction-compressible.xy` | `wallfunction compressible` | 300k from the initial condition |
| `SST/reference/OF-SST-wallfunction-incompressible.xy` | `wallfunction incompressible` | restart of `wallfunction compressible` at 300k, to 450k |
| `SST/reference/OF-SST-asymptotic-incompressible.xy` | `asymptotic incompressible` | restart of `wallfunction compressible` at 300k, to 450k |

The verification of the V&V page pairs `MOSE-SST-asymptotic-compressible` with
`OF-SST-asymptotic-compressible`; its attribution table uses the four OpenFOAM
SST files.

**Updating a reference** (after a change that is meant to move it): copy the
converged `OUTPUT/wall.tec` over the MOSE file, or the OpenFOAM sample (see Running)
over the `OF-` file, then `python compare.py --model <M> --docs`.

## The OpenFOAM companion case

Targets **OpenFOAM-10 (openfoam.org)**, whose `rhoCentralFoam` supports local time
stepping. org-10 dictionary names are used (`constant/physicalProperties`,
`momentumTransport`, `thermophysicalTransport`).

**Same grid.** `OPENFOAM/mesh.geo` is `MESH/script.geo` extruded one hexahedral
layer in z (front/back `empty`), with the top wall split at x = 0 into `symTop`
(symmetry) and `topWall` (isothermal wall) to reproduce the MOSE `[up]` multipatch.
It yields **237,440 hexes — exactly the MOSE quad count**, the check that the two
codes see the same mesh. The node counts in `mesh.geo` are duplicated from
`MESH/script.geo` and **must be kept in sync with it**.

**Same transport.** Both codes use Sutherland's law with standard-air `As`/`Ts`;
OpenFOAM's `sutherland` transport takes the Eucken conductivity (Pr ≈ 0.69), and
the MOSE transport table is built to match.

**Same freestream.** MOSE's turbulence inputs are *density-weighted*
(`mit` = ρν̃, `kappa` = ρk, `omega` = ρω) for the initial condition and, since
the inlet fix in `Lib_BC_Fluxes_Inflow.f90`, for the inlet boundary as well, so
the OpenFOAM values are converted:

| | MOSE `input.ini` | OpenFOAM `0/` |
|---|---|---|
| SA | `mit = 1e-6` | `nuTilda = 4.9005e-06` |
| SST | `kappa = 5.040e-5`, `omega = 248.1` | `k = 2.4699e-04`, `omega = 1215.8` |

The SST levels are the NASA TMR freestream (k = 9e-9·a∞², ω = 1e-6·ρ∞a∞²/μ∞).
Before that fix MOSE multiplied the inlet values by ρ, and the SST `INPUT/bc*.txt`
still carried an older ω, so MOSE actually ran k∞ 4.9x and ω∞ 12.6x below
OpenFOAM (SA: ν̃∞ 4.9x below). It turned out not to change the boundary layer,
but the two codes now see the same numbers.

**Same SST model.** Two implementation choices are aligned explicitly:

- *Production.* MOSE runs `sst-production = compressible`, i.e. P = τᵢⱼ ∂uᵢ/∂xⱼ
  with the dilatation terms (−⅔μₜ(∇·u)², −⅔ρk∇·u, and −⅔γρω∇·u in the ω
  equation) — the form of OpenFOAM's `kOmegaSST`. The default incompressible form
  (μₜS², NASA TMR "SST-2003m") gives a ~16% longer separation bubble here; for
  OpenFOAM runs of that form, `kOmegaSSTMOSE` is `kOmegaSST` with exactly the
  three dilatation terms removed.
- *ω wall condition.* MOSE (`omega-wall-bc = asymptotic`) sets ω = 80ν_w/y_c² on
  the wall **face** and solves the wall-adjacent cell; OpenFOAM's
  `omegaWallFunction` would fix the wall-adjacent **cell** to the same number,
  which moves separation by ~1 mm. The OpenFOAM case therefore uses
  `system/omegaFaceBC`, a `codedFixedValue` that applies MOSE's face condition
  (compiled at run time into `dynamicCode/`).

**Numerics.** MOSE uses HLLC + MUSCL/van Leer + RK3; `rhoCentralFoam` has no HLLC
and uses the Kurganov central-upwind flux with van Leer reconstruction. Steady
state is reached by LTS (`localEuler`, `maxCo 0.2`, `rDeltaTSmoothingCoeff 0.1`,
300k iterations) — read "The LTS trap" before loosening those.

## The LTS trap

An early OpenFOAM run (`maxCo 0.4`, `rDeltaTSmoothingCoeff 0.02`) *looked* plausible
with flat residuals but was **not converged**: downstream of shock impingement the
wall fields oscillated 13–20% cell-to-cell. Differencing two saved times exposed it —
wall pressure there still moved 26% between iterations while the upstream half was
steady to 0.05%. `maxCo 0.2` + `rDeltaTSmoothingCoeff 0.1` removes it entirely.

**Flat residuals are not convergence here.** Verify by differencing successive saved
times (hence `purgeWrite 0`, `writeInterval 25000`).

With the settings above the OpenFOAM solution is converged: separation and
reattachment are fixed from ~75k iterations on, and the wall fields change by
≲1–2% between saves with no trend. The staircase still visible in the OpenFOAM
C_f after reattachment (steps 3.5–4.5 mm apart) is **steady** — identical at
every saved time — and comes from a fan of weak stationary waves next to the wall
that the central-upwind flux leaves undamped; MOSE's HLLC shows the same feature
at about a third of the amplitude.

## Regenerating `INPUT/`

`INPUT/bc*.txt` (one file per multigrid level, ATLAS BCB) and `INPUT/ic.szplt`
(ATLAS ICB) are committed, so the case runs as is. To regenerate them, build the
mesh with GRIB from `MESH/script.geo`, then run `ATLAS BCB` / `ATLAS ICB` and copy
`fromATLAStoSolver/*` into `INPUT/`. Check the result: `bc.txt` must hold 237,440
interior-cell records (and blocks of 849 × {73, 113, 97} nodes). As of October 2026
the Python GRIB builds block 2 with 537 instead of 849 nodes in i for this
geometry, and an ifx + OpenMP build of BCB segfaults (ASSOCIATE opened inside an
OpenMP region); a serial BCB on the 237,440-cell mesh reproduces the committed
files exactly.

## Full analysis

The verification (MOSE vs OpenFOAM), validation (vs Wind-US, SU2, experiment) and
the ω wall-condition and production-form studies are written up, with figures, in
the V&V page: `docs/vv/2D-swbli.md`.
