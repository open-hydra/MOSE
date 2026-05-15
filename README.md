<p align="center">
  <h1 align="center">MOSE</h1>
  <p align="center"><b>Multipurpose Open Solver for ideal-gas Equations</b></p>
</p>

<p align="center">
  <a href="https://github.com/open-hydra/MOSE/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
  <a href="https://open-hydra.github.io/MOSE/"><img src="https://img.shields.io/badge/docs-online-brightgreen.svg" alt="Documentation"></a>
  <img src="https://img.shields.io/badge/language-Fortran-734f96.svg" alt="Language: Fortran">
</p>

---

MOSE is an open-source compressible flow solver for the Euler and Navier–Stokes equations on multi-block structured grids. Written in modern Fortran, it targets a wide range of inviscid and viscous compressible flow problems — from classical shock-tube benchmarks to turbulent reacting flows with finite-rate chemistry and soot formation.

## Features

- **Compressible Euler & Navier–Stokes** — steady-state and time-accurate computations on multi-block structured grids.
- **High-order numerics** — MUSCL reconstruction with multiple flux limiters and Riemann solvers.
- **Turbulence modeling** — RANS closures: Spalart–Allmaras (with rotation/curvature corrections), Menter SST, Wilcox k-ω 2006, SSGLRR, QCR2000.
- **Finite-rate chemistry & soot** — multi-species reacting flows with finite-rate kinetics or chemical equilibrium (NASA CEA). Soot formation models included.
- **Parallel execution** — shared-memory parallelism via OpenMP; MPI support for distributed-memory runs.
- **Flexible I/O** — solution output in Tecplot (ASCII and binary) and VTK formats. Point probes for time-history recording. Restart capability.

## Quick Start

### Prerequisites

| Requirement | Details |
|---|---|
| **CMake** | ≥ 3.23 |
| **Fortran compiler** | GNU (`gfortran`) or Intel/oneAPI (`ifort` / `ifx`) |
| **C/C++ compiler** | Required only for optional components (Cantera, TecIO, SUNDIALS) |

### Build

```bash
git clone --recurse-submodules https://github.com/open-hydra/MOSE.git
cd MOSE

# Build with GNU compilers and OpenMP
./install.sh build --compiler=gnu --use-openmp

# — or with Intel compilers and full feature set —
./install.sh build --compiler=intel --use-openmp --use-mpi --use-tecio --use-sundials --use-cantera
```

The executable is placed in `bin/MOSE`.

See the [Installation Guide](https://open-hydra.github.io/MOSE/getting-started/installation/) for all build options, CMake presets, and troubleshooting.

### Run the Sod Shock Tube

```bash
cd test/1D/Sod79
./MOSE.sh test
```

See the [Quick Start](https://open-hydra.github.io/MOSE/getting-started/quick-start/) for a full walkthrough.


## Dependencies

MOSE is built on top of several companion libraries, included as Git submodules:

| Library | Role |
|---|---|
| [FLINT](https://github.com/MarcoGrossi92/FLINT) | Thermodynamic and chemical kinetics database |
| [ORION](https://github.com/MarcoGrossi92/ORION) | Multi-format I/O (Tecplot, VTK, Plot3D) |
| [OSlo](https://github.com/MarcoGrossi92/OSlo) | ODE solver library (stiff chemistry integration) |
| [FiNeR](https://github.com/szaghi/FiNeR) | INI configuration file parser |
| [ExactPack](https://github.com/lanl/ExactPack) | Exact solutions for verification |

Optional external libraries: **OpenMP**, **MPI**, **SUNDIALS**, **Cantera**, **TecIO**.

## Project Structure

```
MOSE/
├── src/
│   ├── app/           # Main application
│   ├── lib/           # Solver library sources
│   └── test/          # Unit tests
├── lib/               # Git submodule dependencies
├── test/              # Verification & validation cases
│   ├── 1D/            #   1-D test problems
│   ├── 2D/            #   2-D test problems
│   └── ...
├── docs/              # MkDocs documentation source
├── cmake/             # CMake modules
├── install.sh         # Build helper script
└── CMakeLists.txt
```

## Documentation

Full documentation is available at **[open-hydra.github.io/MOSE](https://open-hydra.github.io/MOSE/)**, covering:

- Installation & quick start
- User guide & input file reference
- Verification & validation results
- Theory guide (governing equations, numerical methods, turbulence models)

## License

MOSE is free and open-source software released under the [GNU General Public License v3.0](LICENSE).
