# Acknowledgements

MOSE is built upon several projects.

## FiNeR

**Fortran INI ParseR and generator**

- **Repository:** [github.com/szaghi/FiNeR](https://github.com/szaghi/FiNeR)
- **License:** GPL v3.0

FiNeR is a pure Fortran 2003+ OOP library for reading and writing INI configuration files.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

ORION provides built-in functions to read files in different formats. MOSE uses ORION to read the initial condition,
which holds the mesh, and to write the solution, wall and diagnostic files (Tecplot and VTK).

## FLINT

**Fortran Library for INtegrated Thermochemistry**

- **Repository:** [github.com/MarcoGrossi92/FLINT](https://github.com/MarcoGrossi92/FLINT)
- **License:** GPL v3.0

FLINT is a lightweight, modular Fortran library for thermochemistry in reactive flow simulations.

## OSLO

**ODE SoLver toolbOx**

- **Repository:** [github.com/MarcoGrossi92/OSLO](https://github.com/MarcoGrossi92/OSLO)
- **License:** GPL v3.0

OSLO is a Fortran toolbox for systems of ordinary differential equations; its integrators are those that `ode-solver` selects for the chemistry.

## Documentation Tools

### MkDocs

**Static Site Generator**

- **Website:** [mkdocs.org](https://www.mkdocs.org)
- **License:** BSD-2-Clause
- **Contribution:** Documentation site generation

MkDocs transforms MOSE's documentation into a beautiful, searchable website.

### Material for MkDocs

**Modern Documentation Theme**

- **Website:** [squidfunk.github.io/mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
- **License:** MIT
- **Contribution:** Professional documentation theme with advanced features

Material for MkDocs provides the sleek, modern interface for MOSE's documentation.

## Institutional Support

While MOSE is an independent open-source project, it has benefited from:

- Academic research environments
- Access to various CFD datasets for validation

## License Compliance

MOSE respects all licenses of dependencies and foundations:

- **GPL v3.0** - FLINT, ORION, OSLO, FiNeR

See the [License](license.md) page for MOSE's full license text.

---