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
├── compare.py          comparison figures for the V&V docs (--model SA|SST)
├── plot_sutherland.py  constant-μ vs Sutherland MOSE Cf (analysis helper)
├── SA/                 Spalart–Allmaras
│   ├── input.ini       MOSE case (ATLAS + GPB + GRIB sections)
│   ├── MOSE.sh         ./MOSE.sh solve
│   ├── INPUT/          bc.txt + multigrid levels, ic.tec, transport tables
│   ├── MESH/script.geo gmsh mesh definition (GRIB builds the mesh from it)
│   ├── reference/      MOSE-SA.tec, SU2-SA.dat, wind-SA.dat, schulein.dat
│   ├── OUTPUT/         MOSE result: wall.tec (untracked)
│   └── OPENFOAM/       rhoCentralFoam companion case: run inputs + bottomWall.xy
└── SST/                k-omega SST — same layout (+ MOSE-SST-asymptotic.tec)
```

Results and regenerable data are untracked: MOSE `OUTPUT/`, and the OpenFOAM mesh
(`mesh.msh`, `constant/polyMesh/`), time directories and `postProcessing/`. What
is kept is the final wall data (`OUTPUT/wall.tec`, `OPENFOAM/bottomWall.xy`) and
everything needed to re-run.

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

Sample the OpenFOAM wall from the **raw patch faces** (not an interpolated line,
which smooths away the cell-to-cell noise):

```sh
postProcess -func wallShearStress -latestTime
postProcess -func "patchSurface(patch=bottomWall, fields=(wallShearStress p), \
    interpolate=false, surfaceFormat=raw)" -latestTime
cp postProcessing/patchSurface*/*/patch.xy  bottomWall.xy
```

**Comparison figures** (written straight into `docs/vv/images/`):

```sh
python compare.py --model SA
python compare.py --model SST
```

It reads `OUTPUT/wall.tec` (MOSE) and `OPENFOAM/bottomWall.xy` (OpenFOAM) — the
OpenFOAM wall is read straight from the raw sample, so no OpenFOAM install is
needed to post-process.

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
(`mit` = ρν̃, `kappa` = ρk, `omega` = ρω), so the OpenFOAM values are converted:

| | MOSE `input.ini` | OpenFOAM `0/` |
|---|---|---|
| SA | `mit = 1e-6` | `nuTilda = 4.9005e-06` |
| SST | `kappa = 5.040e-5`, `omega = 248.1` | `k = 2.4699e-04`, `omega = 1215.8` |

The SST levels are the NASA TMR freestream (k = 9e-9·a∞², ω = 1e-6·ρ∞a∞²/μ∞).

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

## Full analysis

The verification (MOSE vs OpenFOAM), validation (vs Wind-US, SU2, experiment) and
the ω wall-BC study are written up, with figures, in the V&V page:
`docs/vv/2D-swbli.md`.
