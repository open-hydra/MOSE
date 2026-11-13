# Installation

This document describes how to obtain and build **MOSE**. The instructions describe the current `install.sh` script, CMake configuration and Git submodule layout.

!!! note
    MOSE has a dual nature: it is both a library and an executable. The installation process described here produces both the static library `libMOSEL.a` and the main executable `bin/MOSE`. If you are only interested in using MOSE as a library, you can link against `libMOSEL.a` without caring about the executable.

## Prerequisites

Before attempting to build MOSE make sure your system provides the following external tools and compilers:

- **CMake** – 3.23 or newer.
- **Fortran compiler** – the GNU toolchain (`gfortran`) or Intel oneAPI (`ifort`/`ifx`). MOSE is tested with gfortran 12 and 13 and with ifx 2024.1. ifx 2023.0 is not supported: it stops with an internal compiler error on StringiFor, a library that FiNeR needs, and with Intel MPI 2021.8 its runs on several MPI ranks with several OpenMP threads each do not repeat the bits of the run on one rank.
- **C and C++ compilers** – required, whatever the options: the CMake project of MOSE enables C and C++ (`project` in `CMakeLists.txt`).
- **OpenMP** – needed for optional parallelization support.

**ifx 2024.1 and IPO.** ifx 2024.1 (tested: 2024.1.0, build 20240308) compiles a loop of the ghost metrics of MOSE wrongly at `-O2` and `-O3` without interprocedural optimization: a minimal program in standard Fortran gives wrong values at `-O2` and `-O3`, also at the default level with no `-O` option, and the right ones at `-O0`, `-O1` and with `-ipo`. MOSE built so stops at set-up on the check of the ghost metrics, or runs with wrong ghost metrics. CMake therefore refuses ifx 2024.1 at `-O2` or above without IPO, with a message: a RELEASE build with `-DUSE_IPO=OFF` (`--no-ipo`), or whose toolchain cannot do IPO, and the TESTING build. The RELEASE build with IPO (the default) and the DEBUG build (`-O0`) are accepted. Other versions of ifx have not been tested.

### ORION, OSLO and FLINT

MOSE requires these three libraries by a range of releases, as FLINT, ATLAS and ICE do (`cmake/SemVerDeps.cmake`). The configure checks the release out in `lib/`.

| Path         | Repository URL                                | Releases           | Purpose                                 |
|--------------|-----------------------------------------------|--------------------|-----------------------------------------|
| `lib/ORION`  | `https://github.com/MarcoGrossi92/ORION.git` | 2.1 or a later 2.x | I/O routines (TecIO, VTK, Plot3D, etc.) |
| `lib/OSLO`   | `https://github.com/MarcoGrossi92/OSLO.git`  | 1.1 or a later 1.x | ODE solvers library, required by FLINT  |
| `lib/FLINT`  | `https://github.com/MarcoGrossi92/FLINT.git` | 2.0 or a later 2.x | Thermodynamic database utilities        |

The first configure needs the network: it lists the tags of each library (`git ls-remote`, https) and clones the newest release that satisfies the range. `deps.lock`, in the source root and not tracked, keeps the releases selected: the next configures use them without the network (`-DFETCHCONTENT_FULLY_DISCONNECTED=ON` works then). `-DSEMVER_UPDATE=ALL` (`./install.sh build --update-deps`) selects the newest releases again; `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<path>` (`--include-orion|flint|oslo=<path>`) builds a local checkout instead; a local checkout outside the range stops the configure (`-DMOSE_PIN_CHECK=WARNING` builds it on purpose). When a project embeds MOSE (Hydra), that project chooses the releases, and MOSE checks them against its ranges. See [Versioning](../development/versioning.md).

### FiNeR

FiNeR is a Git submodule, pinned at one release.

| Path                               | Repository URL                                        | Release  | Purpose                                  |
|------------------------------------|-------------------------------------------------------|----------|------------------------------------------|
| `lib/FiNeR`                        | `https://github.com/szaghi/FiNeR.git`                 | v2.3.0   | INI file parser                          |

The configure stops when `lib/FiNeR` is at another commit than the one MOSE records. FiNeR's own libraries (PENF, FACE, FLAP, BeFoR64 and StringiFor) are not tracked by MOSE: the configure clones them into `lib/FiNeR/src/third_party/`, as FiNeR's author does (`-DMOSE_FINER_DEPS=FREE`, the default: the default branch of each). If an upstream change breaks MOSE, `-DMOSE_FINER_DEPS=BACKUP` checks out the releases MOSE is tested with. `-DFETCHCONTENT_SOURCE_DIR_FINER=<path>` (`--include-finer=<path>`) builds an external FiNeR clone instead; it must be v2.3.0 and have its libraries in place.

[ExactPack](https://github.com/lanl/ExactPack) is a pip dependency of the test
suite only, installed from LANL's GitHub because it is not published on PyPI:

```bash
pip install "ExactPack @ git+https://github.com/lanl/ExactPack.git@1.7.11"
```

## Build methods

First clone the repository with submodules:

```bash
git clone --recurse-submodules https://github.com/open-hydra/MOSE.git
cd MOSE
```

!!! note "An existing clone"
    After pulling new commits into an existing clone, run `git submodule sync && git submodule update --init`: CMake initializes a submodule that is missing, but does not update one that is checked out at an older commit. If `lib/FiNeR` is not at the commit MOSE records, CMake stops and gives this command. FiNeR moved from `lib/third_party/FiNeR` to `lib/FiNeR`: the sync is what makes git know the new path, and the folder `lib/third_party` can be deleted.

    `lib/ORION` and `lib/FLINT` were submodules before they were taken by release: git leaves the two folders in place after that pull. The configure moves each one to its release when it is a detached checkout without changes, and otherwise stops and names the folder.

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
./install.sh build --compilers=gnu --use-openmp

# full configuration with Intel compilers and all optional features
./install.sh build --compilers=intel \
                  --use-openmp --use-mpi \
                  --use-tecio --use-sundials --use-cantera
```

Options accepted by `build`

* `--compilers=<gnu|intel>` – select the compiler family: `gnu` sets `gfortran`, `gcc` and `g++`, `intel` sets `ifx`, `icx` and `icpx`; without it the script sets no compiler and CMake chooses.
* `--use-openmp` – enable OpenMP parallelization.
* `--use-mpi` – enable MPI parallelization.
* `--use-tecio` – enable TecIO support (requires a C++ compiler).
* `--use-sundials` – build and link the bundled SUNDIALS library (requires a C compiler).
* `--use-cantera` – build and link Cantera (requires a C++ compiler; if enabled, SUNDIALS should also be enabled).
* `--no-ipo` – disable interprocedural (link-time) optimization in the RELEASE build.
* `--include-orion=PATH` – build a local ORION checkout instead of a release (`-DFETCHCONTENT_SOURCE_DIR_ORION`).
* `--include-flint=PATH` – same for FLINT.
* `--include-oslo=PATH` – same for OSLO.
* `--include-finer=PATH` – build an external FiNeR clone instead of the submodule `lib/FiNeR` (`-DFETCHCONTENT_SOURCE_DIR_FINER`).
* `--update-deps` – select again the newest releases of ORION, OSLO and FLINT (`-DSEMVER_UPDATE=ALL`, this build only).

The script sets up environment variables for the chosen compilers and then invokes CMake. After a successful build, a `CMakePresets.json` file is written in the source root so that subsequent compilations can reuse the configuration.

<!-- The release path, build artifacts etc.--> 

**`compile` command**

Re‑runs CMake using the previously generated preset and rebuilds the project without clearing the build directory. This is useful during development when only the source has changed. Example usage:

```bash
./install.sh compile
```

**`update` command**

Synchronises the Git submodule `lib/FiNeR`. By default it checks out the commit MOSE records; passing `--remote` fetches the latest commit of its remote branch, which the configure then refuses, as FiNeR is pinned.

```bash
./install.sh update            # sync to recorded commit
./install.sh update --remote   # update to newest remote commit
```

### Build with CMake

If you prefer fine‑grained control, perform the configuration yourself. This is essentially what `install.sh` does under the hood.

```bash
mkdir build && cd build
# Fortran compiler: gfortran or ifx. USE_OPENMP: ON or OFF. Optional: USE_TECIO and
# USE_CANTERA (need a C++ compiler), USE_SUNDIALS (needs a C compiler), FETCHCONTENT_SOURCE_DIR_<NAME> (ORION, OSLO, FLINT, FINER).
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DUSE_OPENMP=ON \
    -DUSE_TECIO=ON \
    -DUSE_SUNDIALS=ON \
    -DUSE_CANTERA=ON \
    -DFETCHCONTENT_SOURCE_DIR_ORION=/path/to/ORION \
    -DFETCHCONTENT_SOURCE_DIR_FLINT=/path/to/FLINT \
    -DFETCHCONTENT_SOURCE_DIR_OSLO=/path/to/OSLO \
    -DFETCHCONTENT_SOURCE_DIR_FINER=/path/to/FiNeR
cmake --build . --parallel
```

The build directory is `build/` in these examples. The static library is `build/src/lib/libMOSEL.a`; the executables (`MOSE`, `DocGen`) go to `bin/` of the source tree,
whatever the build directory; the programs of the test suite, the checks of `test/unit` and the
manufactured-solution driver `MOSE_mms`, are built in `build/test/`.

In a new build directory, an option that the command line does not set takes its default. With no `-D` option, `cmake ..`
configures a RELEASE build with `USE_IPO` ON, `USE_OPENMP`, `USE_MPI` and `USE_CANTERA` OFF, and `USE_TECIO` and `USE_SUNDIALS`
ON, the defaults of ORION and OSLO, which declare these two options (see [TecIO](#tecio) and [SUNDIALS](#sundials)).
`install.sh build` sets each of these options, and turns TecIO and SUNDIALS on only with `--use-tecio` and `--use-sundials`;
with CMake, `-DUSE_TECIO=OFF -DUSE_SUNDIALS=OFF` leave them out.

## CMake presets

The file `CMakePresets.json` produced by the install script records all of the cache variables that were used during configuration. You can build the project later simply with

```bash
cmake --preset default
cmake --build build
```

or using the `compile` command of the install script as described above.

## Optional components

### MPI

`-DUSE_MPI=ON` (`--use-mpi`) builds MOSE with MPI through the Fortran wrapper of the MPI library that CMake finds (`mpif90`, `mpiifx`). With GNU Fortran 10 or later the build compiles `Mod_MPI.f90` and `Mod_Multirate.f90` with `-fallow-argument-mismatch`: they call the same MPI routine with buffers of different types and ranks.

**gfortran with Intel MPI.** The `mpif90` of Intel MPI 2021.12 puts on the include path the `mpi` module built for ifx, which gfortran cannot read ("Reading module .../include/mpi/mpi.mod ... Unexpected EOF"). Give CMake the modules built for gfortran:

```bash
cmake .. -DCMAKE_Fortran_COMPILER=gfortran -DUSE_MPI=ON \
    "-DMPI_Fortran_COMPILER_INCLUDE_DIRS=$I_MPI_ROOT/include/mpi/gfortran/11.1.0;$I_MPI_ROOT/include"
```

The recipe was tested with gfortran 13.2.1: where `gfortran` is an older version, name the newer compiler, for example
`-DCMAKE_Fortran_COMPILER=gfortran-13`.

### TecIO

MOSE can be built with support for TecIO, a library for writing Tecplot files. This is an optional feature shipped with ORION. If enabled, MOSE will be able to read and write output in Tecplot binary format, which can be useful for post‑processing and visualization of large datasets. Enabling TecIO requires a working C++ compiler.
`USE_TECIO` is an option of ORION (`lib/ORION/CMakeLists.txt`), ON by default: a configuration of a new build directory that
does not set it builds TecIO from its sources in `lib/ORION/lib/TecIO` when that folder holds no install of TecIO for the
compilers of the configuration and its choice of MPI.
`install.sh build` enables it only with `--use-tecio`.

### SUNDIALS

SUNDIALS is a suite of solvers for ODEs. The library is included as a submodule in OSLO (`lib/OSLO/lib/sundials`).
`USE_SUNDIALS` is an option of OSLO (`lib/OSLO/CMakeLists.txt`), ON by default: a configuration of a new build
directory that does not set it builds SUNDIALS from that submodule when it finds in OSLO no install of SUNDIALS, or one built
with other Fortran flags (the OpenMP flag of `USE_OPENMP`), and stops with
"Failed to initialize git submodules." when the submodule is missing and OSLO's `git submodule update --init lib/sundials` cannot
fetch it. `install.sh build` enables it only with `--use-sundials`. Enabling SUNDIALS requires a working C compiler.

### Cantera

MOSE’s core thermodynamic library does **not** depend directly on Cantera, but only on FLINT. Same warnings and observations valid for FLINT are applicable to MOSE: Cantera should be used for small projects when advanced chemistry not already included in FLINT is required, but for large‑scale simulations with parallel support the FLINT fortran interface should be preferred. See the [FLINT documentation](https://marcogrossi92.github.io/FLINT/) for details.

!!! note
    If Cantera is enabled, a working C++ compiler is mandatory and `--use-sundials` should also be selected because Cantera itself relies on the SUNDIALS library.

## Library linking (advanced)

The build makes the static library `libMOSEL.a` (CMake target `MOSEL`, alias
`MOSEL::MOSEL`) in `build/src/lib/`, and writes the Fortran modules of MOSE and of
its libraries to `build/modules/`. MOSE has no install rules and no CMake package
file: `cmake --install` installs none of its files, and `find_package(MOSE)` does
not find it.

## Next steps

* **[Quick Start Tutorial](quick-start.md)** – build and verify the installation.

---