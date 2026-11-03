# 2D Flat Plate – Turbulent Boundary Layer

Turbulent flow over a flat plate at Mach 0.2, a compressible flow at low Mach number. This configuration compares MOSE's turbulent boundary layer modeling with well-established compressible Navier–Stokes solvers from the NASA suite (CFL3D, FUN3D), on grids of the NASA zero-pressure-gradient flat plate case but with its own boundary conditions (see [Differences from the TMR case](#differences-from-the-tmr-case)). Multiple turbulence closure models are tested: Spalart–Allmaras (SA), Shear Stress Transport (SST), Wilcox (2006) k-ω, and the SSG/LRR-ω Reynolds-stress model.

**Reference**: [NASA Turbulence Modeling Resource](https://www.nasa.gov/nasa-turbulence-modeling-resource/)

---

## Problem setup

Uniform flow at Mach 0.2 enters from the left over a flat plate that starts at $x = 0$. A symmetry condition is imposed upstream of the leading edge (pre-stagnation region), and a no-slip adiabatic wall condition is applied over the plate.

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
- **SST**: FUN3D with SST-Vm (the name the TMR gives to these results, first reported as SST-V) on the 545 × 385 grid cut into triangles (`cf_plate_fun3d_sstv_tri.dat`);
- **Wilcox 2006**: CFL3D with Wilcox2006-klim-m on the 545 × 385 grid (`cf_plate_w06.dat`, which also holds FUN3D).
- **SSG/LRR-ω**: CFL3D and DLR's TAU with SSGLRR-RSM-w2012 on the 545 × 385 grid, from the TMR page of the model's flat plate (`SSGLRR/reference/cf.dat`), and the profile $u^+(y^+)$ of CFL3D at $x = 0.97$ (`SSGLRR/reference/uplus.dat`).

## Results and verification

For each turbulence model, MOSE computes the wall shear stress distribution along the flat plate. The skin friction coefficient is extracted and compared pointwise against reference solutions.

`verify.py` measures the relative difference of $c_f$ from the reference curve at the cell centres of the plate from $x = 0.02$ m to the last one ($x = 1.979$ m on the 137 × 97 grid, $1.957$ m on the 69 × 49 grid). The case passes when the RMS difference is below the tolerance of its model. For the cases as they are (the differences above included):

| Model | RMS difference | Largest difference | Tolerance on the RMS |
|---|---|---|---|
| SA | 0.77 % | 3.65 % | 5 % |
| SST | 0.95 % | 3.15 % | 8 % |
| Wilcox 2006 | 0.50 % | 2.08 % | 8 % |

On the grids of the case no laminar stretch is seen. From the first cell on the plate ($x = 2$ mm on the 137 × 97 grid, 4 mm on the 69 × 49 grid), where the skin friction is 1.7 to 2.3 times the laminar (Blasius) value, it stays above that value along the whole plate, and it comes within 2 % of the fully turbulent reference curve from $x = 0.044$ m (SA) and $x = 0.029$ m (SST, Wilcox 2006).

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

### SSG/LRR-ω Reynolds-stress model

`turbulence = SSGLRR` (TMR's SSGLRR-RSM-w2012, see the [theory page](../theory/turbulence.md#ssglrr-reynolds-stress-model)), run on the 69 × 49 grid of the TMR series with the TMR wall condition for $\omega$ and compared with the CFL3D and TAU results of the [TMR case page](https://tmbwg.github.io/turbmodels/flatplate_ssglrrrsm.html).

<figure>
  {% include "vv/images/TP-cf-ssglrr.svg" %}
  SSG/LRR-ω: skin friction of MOSE on the 137 × 97 grid (see below) and of CFL3D and TAU on the finest TMR grid (545 × 385)
</figure>

On the grid of the test, 69 × 49, $C_f$ at $x = 0.97$ is 0.002612: 1.6% below CFL3D and 0.7% below TAU on the same grid, which differ by 0.9% from each other. Over the plate it differs from the finest-grid CFL3D curve by 3.4% RMS.

These differences are an effect of the grid. The case was also run on the next coarser and the next finer grid of the TMR series, 35 × 25 and 137 × 97, with the same `input.ini` (300,000 iterations), and `INPUT/ic.tec` and `INPUT/bc.txt` written by ATLAS ICB and BCB on the grid files of `test/2D/viscous/flat-plate-turbulent-tmr/MESH`. On the 69 × 49 grid this gives the files of the case, byte for byte. These two runs are not part of the test.

| Grid | MOSE | CFL3D | TAU | MOSE vs CFL3D | MOSE vs TAU | $C_f$ RMS vs CFL3D 545 × 385 | $u^+$ RMS |
|---|---|---|---|---|---|---|---|
| 35 × 25 | 0.002460 | 0.002626 | 0.002555 | −6.3% | −3.7% | 9.1% | 7.5% |
| 69 × 49 | 0.002612 | 0.002654 | 0.002631 | −1.6% | −0.7% | 3.4% | 3.1% |
| 137 × 97 | 0.002670 | 0.002674 | 0.002667 | −0.1% | +0.1% | 1.0% | 1.2% |
| 545 × 385 | | 0.002690 | 0.002690 | | | | |

The gap between MOSE and the TMR codes closes 4 to 12 times at each refinement. On the 137 × 97 grid MOSE lies between CFL3D and TAU, within 0.13% of both. Richardson extrapolation of the three MOSE values gives 0.002706 (observed order 1.4), 0.6% above the 545 × 385 value. The same extrapolation of the CFL3D and TAU values on these three grids overshoots their own 545 × 385 value by 0.9% and 0.5%: the 35 × 25 grid is outside the asymptotic range for all three codes. On the 137 × 97 grid $C_f$ at $x = 0.97$ changes by $2\times10^{-6}$ over the last 60,000 iterations.

<figure>
  {% include "vv/images/TP-uplus-ssglrr.svg" %}
  SSG/LRR-ω: velocity profile in wall units in the cell column nearest x = 0.97 of the 137 × 97 grid (x = 0.951) against CFL3D at x = 0.97 (finest grid), and the Reynolds stresses of the model
</figure>

The velocity profile follows the CFL3D one within 1.2% RMS over $2 < y^+ < 1000$ (3.1% on the 69 × 49 grid). In the log layer the model returns the anisotropy it is built to resolve, $\overline{u'u'} : \overline{v'v'} : \overline{w'w'} = 0.91 : 0.44 : 0.66$ of $k$ with $-\overline{u'v'}/k = 0.30$, where an eddy-viscosity model gives 2/3 to each normal stress.

The variants of the model, on the 69 × 49 grid:

| `turbulence` | `omega-wall-bc` | $C_f$ at $x = 0.97$ | Difference |
|---|---|---|---|
| `SSGLRR` | `practical` (TMR) | 0.002612 | - |
| `SSGLRR-SD` | `practical` | 0.002629 | +0.7% |
| `SSGLRR-w2019` | `practical` | 0.002612 | none (the length-scale correction is inactive in an attached boundary layer) |
| `SSGLRR` | `asymptotic` | 0.002698 | +3.3% |

The solution does not depend on the `point-implicit` treatment of the source terms: with it off the wall shear stress changes by less than $10^{-6}$ anywhere on the plate. The density residual does not fall below $7\times10^{-7}$: over the last 10,000 of the 100,000 iterations it stays between $7.3\times10^{-7}$ and $7.8\times10^{-7}$. The two runs reach the same wall shear stress, within $10^{-6}$, by different iteration paths, so the residual that is left does not show in it.
