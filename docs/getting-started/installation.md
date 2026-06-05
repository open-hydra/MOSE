# Installation

This document describes how to obtain and build **MOSE**. The instructions describe the current `install.sh` script, CMake configuration and Git submodule layout.

!!! note
    MOSE has a dual nature: it is both a library and an executable. The installation process described here produces both the static library `libMOSE.a` and the main executable `bin/MOSE`. If you are only interested in using MOSE as a library, you can link against `libMOSE.a` without caring about the executable.

## Prerequisites

Before attempting to build MOSE make sure your system provides the following external tools and compilers:

- **CMake** – 3.23 or newer.
- **Fortran compiler** – either the GNU toolchain (`gfortran`) or Intel/oneAPI (`ifort`/`ifx`) are supported.
- **C++ compiler** – required only if you wish to enable optional components that depend on C++ (Cantera, TecIO).
- **C compiler** – needed when building SUNDIALS (optional).
- **OpenMP** – needed for optional parallelization support.

### Git submodules

MOSE depends on several repositories that are included as Git submodules.

| Path                               | Repository URL                                        | Purpose                                  |
|------------------------------------|-------------------------------------------------------|------------------------------------------|
| `lib/third_party/ExactPack`        | `https://github.com/lanl/ExactPack.git`              | Exact solutions database for reference    |
| `lib/ORION`                        | `https://github.com/MarcoGrossi92/ORION.git`         | I/O routines (TecIO, VTK, Plot3D, etc.)  |
| `lib/third_party/FiNeR`            | `https://github.com/szaghi/FiNeR.git`                | INI file parser                      |
| `build/lib/OSLO` (installed by FLINT) | `https://github.com/MarcoGrossi92/OSLO.git`       | ODE solvers library               |
| `lib/FLINT`                        | `https://github.com/MarcoGrossi92/FLINT.git`         | Thermodynamic database utilities          |

## Build methods

First clone the repository with submodules:

```bash
git clone https://github.com/MarcoGrossi92/MOSE.git
cd MOSE
# initialise submodules
git submodule update --init --recursive
```

To fully install MOSE, you may either use the bundled install script or invoke CMake manually. The script is convenient and is the preferred route for most users.

### Build with `install.sh` (recommended)

The script exposes three commands: `build`, `compile`, and `update`.  It also maintains a `CMakePresets.json` file that records the configuration used for
the most recent `build` invocation.

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global options**

* `-v`, `--verbose` – enable verbose logging.

**`build` command**

Performs a clean configure+build cycle.  Example usage:

```bash
# minimal GNU build with OpenMP enabled
./install.sh build --compiler=gnu --use-openmp

# full configuration with Intel compilers and all optional features
./install.sh build --compiler=intel \
                  --use-openmp --use-mpi \
                  --use-tecio --use-sundials --use-cantera
```

Options accepted by `build`

* `--compiler=<gnu|intel>` – select the compiler family (default: `gnu`).
* `--use-openmp` – enable OpenMP parallelization.
* `--use-mpi` – enable MPI parallelization.
* `--use-tecio` – enable TecIO support (requires a C++ compiler).
* `--use-sundials` – build and link the bundled SUNDIALS library (requires a C compiler).
* `--use-cantera` – build and link Cantera (requires a C++ compiler; if enabled, SUNDIALS should also be enabled).
* `--include-orion=PATH` – use an external ORION tree instead of the submodule.
* `--include-flint=PATH` – same for FLINT.
* `--include-oslo=PATH` – same for OSLO.
* `--include-finer=PATH` – same for FiNeR.

The script sets up environment variables for the chosen compilers and then invokes CMake. After a successful build, a `CMakePresets.json` file is written in the source root so that subsequent compilations can reuse the configuration.

<!-- The release path, build artifacts etc.--> 

**`compile` command**

Re‑runs CMake using the previously generated preset and rebuilds the project without clearing the build directory. This is useful during development when only the source has changed. Example usage:

```bash
./install.sh compile
```

**`update` command**

Synchronises the Git submodules. By default it checks out the commit recorded in `.gitmodules`; passing `--remote` will fetch the latest commit from each remote branch.

```bash
./install.sh update            # sync to recorded commit
./install.sh update --remote   # update to newest remote commit
```

### Build with CMake

If you prefer fine‑grained control, perform the configuration yourself. This is essentially what `install.sh` does under the hood.

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \  # or ifx
    -DUSE_OPENMP=ON \                    # or OFF
    -DUSE_TECIO=ON \                     # optional (needs C++ compiler)
    -DUSE_SUNDIALS=ON \                  # optional (needs C compiler)
    -DUSE_CANTERA=ON \                   # optional (needs C++ compiler)
    -DORION_PATH=/path/to/ORION \        # optional
    -DFLINT_PATH=/path/to/FLINT \        # optional
    -DOSLO_PATH=/path/to/OSLO \          # optional
    -DFINER_PATH=/path/to/FiNeR          # optional
cmake --build . --parallel
```

The resulting artifacts are placed in `build/` by default. The static library is `lib/libMOSE.a` and the main executable is `bin/MOSE` (inside the build
directory unless you set `CMAKE_INSTALL_PREFIX`).

## CMake presets

The file `CMakePresets.json` produced by the install script records all of the cache variables that were used during configuration. You can build the project later simply with

```bash
cmake --preset default
cmake --build build
```

or using the `compile` command of the install script as described above.

## Optional components

### TecIO

MOSE can be built with support for TecIO, a library for writing Tecplot files. This is an optional feature shipped with ORION. If enabled, MOSE will be able to read and write output in Tecplot binary format, which can be useful for post‑processing and visualization of large datasets. Enabling TecIO requires a working C++ compiler.

### SUNDIALS

SUNDIALS is a suite of solvers for ODEs. The library is included as a submodule in OSLo and will be built automatically if you select the `--use-sundials` option during installation. Enabling SUNDIALS requires a working C compiler.

### Cantera

MOSE’s core thermodynamic library does **not** depend directly on Cantera, but only on FLINT. Same warnings and observations valid for FLINT are applicable to MOSE: Cantera should be used for small projects when advanced chemistry not already included in FLINT is required, but for large‑scale simulations with parallel support the FLINT fortran interface should be preferred. See the [FLINT documentation](https://marcogrossi92.github.io/FLINT/) for details.

!!! note
    If Cantera is enabled, a working C++ compiler is mandatory and `--use-sundials` should also be selected because Cantera itself relies on the SUNDIALS library.

## Library linking (advanced)

To use MOSE from an external Fortran program you can compile as follows:

```bash
gfortran -I/path/to/MOSE/include \
         -L/path/to/MOSE/lib \
         -lMOSE \
         your_program.f90 -o your_program
```

or, from a CMake project:

```cmake
find_package(MOSE REQUIRED)
add_executable(myapp main.f90)
target_link_libraries(myapp MOSE::MOSE)
```

Installation prefix and other details may be customised via standard CMake
variables such as `CMAKE_INSTALL_PREFIX` and by running `cmake --install`.

## Next steps

* **[Quick Start Tutorial](quick-start.md)** – build and verify the installation.

---