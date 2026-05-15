# Contributing

Thank you for your interest in contributing to MOSE.  This page
describes the development workflow, coding conventions, and review
process.

---

## Development Workflow

```mermaid
gitGraph
    commit id: "main"
    branch feature/my-change
    checkout feature/my-change
    commit id: "implement"
    commit id: "test"
    checkout main
    merge feature/my-change id: "merge"
    commit id: "release"
```

1. **Fork or branch** — create a feature branch from `main`.
2. **Implement** — make your changes in small, focused commits.
3. **Test** — run the relevant test suite (see [Testing](testing.md)).
4. **Push** — push your branch and open a pull request.
5. **Review** — address feedback; maintainers will merge when ready.

---

## Setting Up a Development Environment

```bash
# 1. Clone with submodules
git clone --recurse-submodules https://github.com/open-hydra/MOSE.git
cd MOSE

# 2. Build in debug mode
./install.sh build --compilers=gnu

# 3. (Optional) Update submodules to latest
./install.sh update --remote

# 4. Iterative recompilation
./install.sh compile
```

!!! tip "CMakePresets.json"

    After the first `install.sh build`, a `CMakePresets.json` is
    generated with your compiler paths and library locations.
    Subsequent builds only need `./install.sh compile` (or
    `cmake --build build`).

---

## Coding Conventions

### Fortran Style

| Rule | Detail |
|------|--------|
| **Standard** | Fortran 2008 (`-std=f2008`) |
| **Line length** | Hard limit of **132 characters**; use `&` continuation |
| **Indentation** | 2 spaces; no tabs |
| **Module naming** | `MOSE_<Name>` for public modules |
| **Variable naming** | `snake_case` for locals; `obj_` prefix for global singletons |
| **Implicit typing** | Always use `implicit none` |
| **Intent** | Declare `intent(in)`, `intent(out)`, or `intent(inout)` for all arguments |
| **Precision** | Use `iso_fortran_env` kinds: `int32`, `real64` |

### File Organisation

- One module per file (exceptions: small helper modules).
- File name matches module name minus the `MOSE_` prefix
  (e.g. `Lib_SST.f90` → module `MOSE_Lib_SST`).
- Computational routines go in `Lib_*.f90`; types and pointers in
  `Mod_*.f90`.

### Commit Messages

Use imperative mood, present tense:

```
Add SLAU2 Riemann solver

Implement the Kitamura & Shima (2013) pressure-flux variant of SLAU.
Includes unit test for normal-shock pressure recovery.
```

---

## Adding a New Riemann Solver

As a concrete example of extending MOSE, here is how to add a new
Riemann solver:

```mermaid
flowchart TD
    A["1. Create<br/>Lib_Riemann_New.f90<br/>in numerics/fluxes/"] --> B["2. Implement<br/>riemann_New(...)"]
    B --> C["3. Register in<br/>Mod_Riemann.f90"]
    C --> D["4. Add case to<br/>Assign_Setup.f90"]
    D --> E["5. Add test case<br/>in test/basic/Riemann/"]
    E --> F["6. Document in<br/>docs/theory/riemann-solvers.md"]

    style A fill:#1a237e,stroke:#7986cb,color:#e8eaf6
    style F fill:#1b5e20,stroke:#a5d6a7,color:#fff
```

1. Create `src/lib/numerics/fluxes/Lib_Riemann_New.f90` with a
   subroutine matching the `Riemann` procedure-pointer interface.
2. Implement the solver; the interface receives left/right primitive
   states, the face normal, and returns the numerical flux.
3. In `Mod_Riemann.f90`, add the new solver name to the selection
   logic.
4. In `config/Assign_Setup.f90`, bind the procedure pointer
   `Riemann => riemann_New` for the new keyword.
5. Add a test case that exercises the solver on a known Riemann
   problem.
6. Document the solver's theory and properties.

---

## Adding a New Turbulence Model

```mermaid
flowchart TD
    A["1. Create<br/>Lib_MyModel.f90<br/>in physics/turbulence/"] --> B["2. Implement<br/>Eddy_Viscosity<br/>Diffusive_Flux<br/>Source_Terms<br/>Wall_Values"]
    B --> C["3. Register in<br/>Mod_RANS.f90"]
    C --> D["4. Add case to<br/>Assign_Setup.f90"]
    D --> E["5. Validate on<br/>flat-plate / channel"]
    E --> F["6. Document in<br/>docs/theory/turbulence.md"]

    style A fill:#4a148c,stroke:#ce93d8,color:#f3e5f5
    style F fill:#1b5e20,stroke:#a5d6a7,color:#fff
```

The procedure-pointer interface in `Mod_RANS` requires four routines:

| Pointer | Signature purpose |
|---------|-------------------|
| `Eddy_Viscosity` | Compute $\mu_t$ from model variables |
| `RANS_Diffusive_Flux` | Turbulent diffusion for model equations |
| `RANS_Source_Terms` | Production / destruction source |
| `RANS_Set_Wall_Values` | Wall boundary conditions |

---

## Layer Boundaries

Respect the separation between layers:

```mermaid
graph TB
    A["driver/"] -->|calls| B["numerics/ · physics/"]
    B -->|calls| C["base/ · io/"]
    A -->|calls| C
    D["config/"] -->|configures| B
    D -->|reads| C

    A x-.-x|never| D
    B x-.-x|never| A

    style A fill:#263238,stroke:#90a4ae,color:#eceff1
    style B fill:#004d40,stroke:#80cbc4,color:#e0f2f1
    style C fill:#3e2723,stroke:#a1887f,color:#efebe9
    style D fill:#1a237e,stroke:#7986cb,color:#e8eaf6
```

- **driver/** calls numerics/physics and io but never config.
- **numerics/** and **physics/** never call upward into driver.
- **config/** configures numerics/physics but does not call driver.
- **base/** is a leaf layer with no upward dependencies.

---

## Pull Request Checklist

Before submitting a PR, please verify:

- [ ] Code compiles with `-std=f2008 -Wall -Wextra` without warnings
- [ ] No line exceeds 132 characters
- [ ] All new subroutines have `implicit none` and argument intents
- [ ] Existing tests pass (`test/test.sh`)
- [ ] New functionality has at least one test case
- [ ] Documentation updated if user-facing behaviour changed
- [ ] Commit messages follow the convention above
