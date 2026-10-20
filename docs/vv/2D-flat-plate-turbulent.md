# 2D Flat Plate – Turbulent Boundary Layer

Turbulent flow over a flat plate at Mach 0.2, a compressible flow at low Mach number. This configuration compares MOSE's turbulent boundary layer modeling with well-established compressible Navier–Stokes solvers from the NASA suite (CFL3D, FUN3D), on grids of the NASA zero-pressure-gradient flat plate case but with its own boundary conditions (see [Differences from the TMR case](#differences-from-the-tmr-case)). Multiple turbulence closure models are tested: Spalart–Allmaras (SA), Shear Stress Transport (SST), and Wilcox (2006) k-ω.

**Reference**: [NASA Turbulence Modeling Resource](https://www.nasa.gov/nasa-turbulence-modeling-resource/)

---

## Problem setup

Uniform flow at Mach 0.2 enters from the left over a flat plate that starts at $x = 0$. The boundary layer rapidly transitions from laminar to turbulent flow, developing a fully turbulent boundary layer downstream. A symmetry condition is imposed upstream of the leading edge (pre-stagnation region), and a no-slip adiabatic wall condition is applied over the plate.

The skin friction coefficient is the primary validation metric:

$$C_f = \frac{\tau_w}{0.5 \rho_\infty U_\infty^2}$$

where $\tau_w = \mu \frac{\partial u}{\partial y}\big|_{\text{wall}}$ is the wall shear stress.

**Freestream conditions**

| Parameter | Value |
|---|---|
| Mach number | 0.2 |
| Temperature | 555 K |
| Pressure | 1.0 × 10⁵ Pa (total pressure 1.02828 × 10⁵ Pa) |
| Viscosity | 1.1859 × 10⁻⁵ Pa·s, constant |
| Re/m | 5.0 × 10⁶ m⁻¹ |

Air is a perfect gas with $\gamma = 1.4$ and $R = 287$ J/(kg·K) (`test/common/Air-viscous`), with a constant laminar Prandtl number of 0.72 and MOSE's default turbulent Prandtl number (`Prt` = 0.85).

## Numerical setup

All three cases march in pseudo-time (RK3 with local time stepping) for the number of iterations below.

| Parameter | SA | SST | Wilcox 2006 |
|---|---|---|---|
| Grid (nodes) | 137 × 97 | 69 × 49 | 69 × 49 |
| Iterations | 200 000 | 100 000 | 100 000 |
| CFL / VNN | 0.8 / 0.3 | 0.8 / 0.3 | 1.0 / 1.0 |
| IRS $\beta$ | 0.5 | 0.5 | 0.5 |
| Space reconstruction | MUSCL, Van Leer | MUSCL, Van Leer | MUSCL, Van Leer |
| Riemann solver | AUSM+ | SLAU | HLLC |
| Wall value of $\omega$ | — | asymptotic | asymptotic |

## Grid structure

The meshes are grids of the TMR family `flatplate_clust2`, node for node: the 137 × 97 grid (`flatplate_clust2_2levelsdown_137x97`) for SA and the 69 × 49 grid (`flatplate_clust2_3levelsdown_69x49`) for SST and Wilcox 2006, two and three levels below the finest TMR grid of 545 × 385 nodes. The nodes cluster near the leading edge and near the wall.

Boundary conditions:

- **Inlet** ($x = -0.33333$, type `404`): mass flux 59.2934 kg m⁻² s⁻¹ and static temperature 555 K, flow normal to the face; inflow turbulence $\tilde\nu = 10^{-6}$ m² s⁻¹ (SA), $k = 2.25\times10^{-3}$ m² s⁻² and $\omega = 6250$ s⁻¹ (SST, Wilcox 2006)
- **Adiabatic wall** (x ≥ 0, type `301`): No-slip, no-heat-flux
- **Symmetry** (x < 0, type `300`): Symmetry boundary condition upstream
- **Top** ($y = 1$, type `300`): symmetry
- **Outlet** ($x = 2$, type `406`): back pressure 1.0 × 10⁵ Pa

## Turbulence models

1. **Spalart–Allmaras (SA)**: One-equation model; computationally efficient for wall-bounded flows
2. **Shear Stress Transport (SST)**: Two-equation model; blends k-ω near-wall and k-ε in the freestream
3. **Wilcox k-ω**: Traditional two-equation model; excellent wall-layer resolution

## Differences from the TMR case

The case runs at the Mach and Reynolds numbers of the TMR zero-pressure-gradient flat plate ($M = 0.2$, $Re = 5\times10^6$ per unit length) on two of its grids, but it is not the TMR setup:

| | This case | TMR zero-pressure-gradient flat plate |
|---|---|---|
| Inflow | mass flux and static temperature (`404`) | total pressure and temperature, $p_t/p_\mathrm{ref} = 1.02828$, $T_t/T_\mathrm{ref} = 1.008$, one quantity from the interior |
| Top boundary | symmetry (`300`): a channel of height 1 | far-field Riemann condition |
| Gas | 555 K; constant viscosity and conductivity ($\mathrm{Pr} = 0.72$); $\mathrm{Pr}_t = 0.85$ | $T_\mathrm{ref}$ = 540 °R (300 K); Sutherland's law, $\mathrm{Pr} = 0.72$ and $\mathrm{Pr}_t = 0.9$, as the TMR states for its 2D Finite Flat Plate case, source of the SA reference curve |
| Inflow turbulence, SA | $\tilde\nu_\infty = 10^{-6}$ m² s⁻¹ $= 0.053\,\nu_\infty$ | $\tilde\nu_\infty \geq 3\,\nu_\infty$ ($5.67\times10^{-5}$ m² s⁻¹ here) |
| Inflow turbulence, SST and Wilcox 2006 | $k_\infty = 2.25\times10^{-3}$ m² s⁻², $\omega_\infty = 6250$ s⁻¹ | $k_\infty = 9\times10^{-9}\,a_\infty^2$, $\omega_\infty = 10^{-6}\,\rho_\infty a_\infty^2/\mu_\infty$ ($2.007\times10^{-3}$ m² s⁻² and 11 806 s⁻¹ here) |
| Grids | 137 × 97 (SA), 69 × 49 (SST, Wilcox 2006) | five nested grids, 35 × 25 to 545 × 385 |

The reference curves of `verify.py` (`reference/cf.dat`) are TMR data:

- **SA**: CFL3D on the finest grid (level 1 of 8) of the TMR *2D Finite Flat Plate* numerics case, whose grid extends upstream and downstream of the plate (`all_cf_global.dat`);
- **SST**: FUN3D with SST-Vm (the name the TMR now gives to these results, first reported as SST-V) on the 545 × 385 grid cut into triangles (`cf_plate_fun3d_sstv_tri.dat`);
- **Wilcox 2006**: CFL3D with Wilcox2006-klim-m on the 545 × 385 grid (`cf_plate_w06.dat`, which also holds FUN3D).

## Results and verification

For each turbulence model, MOSE computes the wall shear stress distribution along the flat plate. The skin friction coefficient is extracted and compared pointwise against reference solutions.

`verify.py` measures the relative difference of $c_f$ from the reference curve at the cell centres of the plate from $x = 0.02$ m to the last one ($x = 1.979$ m on the 137 × 97 grid, $1.957$ m on the 69 × 49 grid). With MOSE at commit `10b59003b72c` (the cases as they are, the differences above included):

| Model | RMS difference | Largest difference | Tolerance on the RMS |
|---|---|---|---|
| SA | 0.77 % | 3.65 % | 5 % |
| SST | 0.95 % | 3.15 % | 8 % |
| Wilcox 2006 | 0.50 % | 2.08 % | 8 % |

<figure>
  {% include "vv/images/TP-cf-sa.svg" %}
  Spalart–Allmaras
</figure>

<figure>
  {% include "vv/images/TP-cf-sst.svg" %}
  Shear Stress Transport (SST)
</figure>

<figure>
  {% include "vv/images/TP-cf-wilcox.svg" %}
   Wilcox k-ω
</figure>

!!! note "Figure and code version"
    The figures come from July 2026 (commit `bd7ce104`), and the code of 13 October 2026 (commit
    `b29130b`) gives them again, within 0.2 % of the range of the axes with SST and 0.4 % with Wilcox
    2006. MOSE at commit `10b59003b72c` moves the SST curve by 1.8 % and the Wilcox 2006 curve by 1.7 % of
    that range.

    - SST: the curve moves at three commits, by 0.4 %, 1.1 % and 3.0 % of the range of the axes, steps
      that partly cancel: the viscous face metric, from the mean of the metrics of the two cells to the
      inverse of the mean of their inverses (commit `208f6f14`); the point-implicit factor of the
      turbulence equations, which now divides their whole increment (commit `2cc7bc91`; before it the
      steady state depended on the local time step); the face metric from the differences of the cell
      centres (commit `000d28a12`, measured as the switch of `diffusive-metric` at `10b59003b72c`). A bisection over the commits between `b29130b` and `10b59003b72c`
      finds no other step above 0.1 %.
    - Wilcox 2006: with the inverse-mean face metric of commit `208f6f14`
      (`diffusive-metric = inverse-mean`) the curve is within 0.5 % of the figure, so most of the change
      comes from the face metric of commit `000d28a12`.

    These numbers compare the skin friction after the 100 000 iterations of the inputs; the change of the
    time step rule alone (commit `667fc3c0`) moves the SST curve by 0.01 %.
