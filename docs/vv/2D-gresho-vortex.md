# Gresho Vortex

A stationary, rotating vortex held in exact centrifugal balance. Because the
Gresho vortex is a **steady solution of the Euler equations**, the exact answer
is simply the initial field preserved unchanged in time — any departure of the
computed solution is *pure numerical error*. Running the same vortex at
decreasing Mach numbers turns it into a stringent probe of a Riemann solver's
**low-Mach accuracy**: standard upwind fluxes carry an acoustic dissipation that
scales as $\mathcal O(1/M)$ and spin the vortex down, whereas all-speed schemes
must keep it alive. This case therefore validates the low-Mach treatment of the
[Riemann solvers](../theory/riemann-solvers.md) and the role of the shared
cutoff Mach $M_{co}$.

**Reference**: Gresho & Chan (1990); Mach-scaled formulation of Miczek, Röpke &
Edelmann, *Astron. Astrophys.* 576 (2015) A50; low-Mach Roe analysis of Rieper,
*J. Comput. Phys.* 230 (2011).

---

## Problem setup

A vortex of unit density is centred in a unit square. The azimuthal velocity is a
fixed triangular profile (peak $|u_\phi| = 1$), and the pressure carries exactly
the radial increment needed for centrifugal balance, $\mathrm{d}p/\mathrm{d}r =
\rho\,u_\phi^2/r$:

$$
u_\phi(r) =
\begin{cases}
5r & 0 \le r < 0.2 \\
2 - 5r & 0.2 \le r < 0.4 \\
0 & r \ge 0.4
\end{cases}
\qquad
\rho = 1 .
$$

The background pressure sets the Mach number. With the peak velocity fixed at
unity,

$$
p_0 = \frac{\rho\,u_{\max}^2}{\gamma M^2} = \frac{1}{\gamma M^2},
$$

so the **same vortex** is run at any target Mach number $M$ by changing $p_0$
only. The balance increment $\mathrm{d}p$ is independent of $M$, hence the
relative pressure fluctuation scales as $\mathrm{d}p/p_0 \sim M^2$ — the
hallmark of the low-Mach regime that the scheme must resolve without being
swamped by acoustic dissipation.

Three reference Mach numbers are tested: $M_{\rm ref} = 0.1,\ 0.01,\ 0.001$.

## Numerical setup

| Parameter | Value |
|---|---|
| Time scheme | RK3 (time-accurate) |
| CFL | 0.3 |
| Space reconstruction | MUSCL |
| Flux limiter | minmod |
| Low-Mach cutoff | $M_{co} = 0.005$ ($\approx 5\,M_{\rm peak}$) |
| **Riemann solvers tested** | **HLLC (baseline), LMRoe, AUSM+M** |

The cutoff $M_{co}$ is the shared low-Mach dissipation floor described in
[Low-Mach Cutoff Mach $M_{co}$](../theory/riemann-solvers.md#low-mach-cutoff-mach-m_co).
**HLLC** is the *control*: it carries no low-Mach treatment and shows what the
baseline upwind dissipation does to the vortex. **LMRoe** and **AUSM+M** are the
low-Mach-corrected schemes.

## Grid structure

A uniform $100 \times 100$ Cartesian grid on the unit square $[0,1]^2$, vortex
centred at $(0.5, 0.5)$. All four boundaries are **periodic**, so the only error
sources are the spatial flux and the time integration — there is no boundary
contamination.

## Results and verification

### Mach-number field: solver vs Mach matrix

The normalised Mach field $M/M_{\rm ref}$ is shown for the three solvers (rows)
across the three reference Mach numbers (columns). A perfectly preserved vortex
reproduces the same ring pattern at every Mach number.

<figure>
  {% include "vv/images/gresho_mach_matrix.svg" %}
</figure>

- **HLLC (baseline)** preserves the vortex acceptably at $M_{\rm ref}=0.1$ and
  $0.01$ but **collapses at $M_{\rm ref}=0.001$**: the $\mathcal O(1/M)$ acoustic
  dissipation drains the rotation, halving the peak velocity. This is exactly the
  low-Mach failure the all-speed schemes are designed to remove.
- **LMRoe** and **AUSM+M** keep the vortex intact across the **entire** Mach
  range — the ring is essentially Mach-independent, confirming the correct
  low-Mach scaling of their dissipation.

### Azimuthal velocity profiles

The tangential velocity along the horizontal centreline ($y = 0.5$), compared to
the exact triangular profile:

<figure>
  {% include "vv/images/gresho_profiles_by_mach.svg" %}
</figure>

At $M_{\rm ref}=0.001$ the HLLC profile is strongly clipped (peak well below
$1$), while LMRoe and AUSM+M track the exact peak closely. At higher Mach all
three are reasonable, the low-Mach schemes always being less dissipative.

### Quantitative error metrics

For each run: peak azimuthal velocity decay (exact $|u_\phi|_{\max} = 1$),
kinetic-energy loss relative to the analytic vortex, and the $L_2$ error of
$u_\phi(r<0.4)$.

| Solver | $M_{\rm ref}$ | Peak $u_\phi$ decay | KE lost | $L_2(u_\phi)$ |
|---|---|---|---|---|
| HLLC (baseline) | 0.1   | 14.6 % | 9.4 % | 4.90e-2 |
| HLLC (baseline) | 0.01  | 12.7 % | 5.2 % | 3.41e-2 |
| HLLC (baseline) | 0.001 | **52.3 %** | **54.0 %** | **3.21e-1** |
| LMRoe | 0.1   | 8.3 % | 1.9 % | 1.75e-2 |
| LMRoe | 0.01  | 5.2 % | 0.4 % | 9.52e-3 |
| LMRoe | 0.001 | 11.2 % | 5.6 % | 3.83e-2 |
| AUSM+M | 0.1   | 9.0 % | 2.4 % | 2.13e-2 |
| AUSM+M | 0.01  | 5.5 % | 0.3 % | 1.01e-2 |
| AUSM+M | 0.001 | 5.4 % | 0.9 % | 2.83e-2 |

The contrast at $M_{\rm ref}=10^{-3}$ is decisive: HLLC loses **54 %** of the
kinetic energy, whereas LMRoe and AUSM+M lose only **5.6 %** and **0.9 %**. Both
low-Mach schemes also keep density variations below $\sim 0.1\%$ and reproduce
the expected $\mathrm{d}p/p_0 \sim M^2$ pressure scaling, confirming a clean
incompressible-limit behaviour.

### Effect of the cutoff Mach number$

The cutoff Mach has a **two-sided** optimum, illustrated by sweeping $M_{co}$ for
LMRoe at the hardest case $M_{\rm ref}=10^{-3}$:

<figure>
  {% include "vv/images/gresho_lmroe_sweep.svg" %}
</figure>

| $M_{co}$ | Peak $u_\phi$ decay | KE change | $L_2(u_\phi)$ | Behaviour |
|---|---|---|---|---|
| 0.500 | 54.2 % | −55.2 % | 3.23e-1 | **Too high** — over-dissipation; vortex spun down |
| 0.005 | 11.2 % | −5.6 % | 3.83e-2 | **Optimal** — vortex preserved |
| 0.001 | −17.1 % | **+109 %** | 4.73e-1 | **Too low** — odd–even instability; energy pumped in |

- **Too high** ($M_{co}=0.5$): the floor reinstates a large acoustic
  dissipation $\sim M_{co}/M_{\rm local}$, draining the vortex almost as badly as
  plain HLLC. In the limit $M_{co}\to 1$ the low-Mach correction switches off and
  the baseline scheme is recovered — stable but inaccurate.
- **Too low** ($M_{co}=0.001$): the acoustic coupling is too weak to suppress
  the checkerboard/odd–even mode. Kinetic energy *grows* by 109 % (peak velocity
  exceeds the exact value) and the radial velocity is dominated by grid-scale
  noise — the signature of the spurious mode pumping energy into the solution.
- **Optimal** ($M_{co}=0.005 \approx 5\,M_{\rm peak}$): sits in the window,
  above the stability threshold and only mildly over-dissipative.

The asymmetry is the practical take-away: *too high* is stable but inaccurate,
*too low* is outright unstable. Set $M_{co}$ to a few times the flow's peak Mach,
and err slightly high when uncertain.

## Conclusions

- Plain HLLC exhibits the textbook $\mathcal O(1/M)$ low-Mach failure, losing
  half the vortex kinetic energy at $M=10^{-3}$.
- **LMRoe and AUSM+M preserve the vortex across the full $M = 0.1\!-\!0.001$
  range** once the shared cutoff $M_{co}$ is set appropriately, validating their
  all-speed dissipation scaling.
- The $M_{co}$ sweep confirms the predicted stability/accuracy window: a single,
  Mach-scaled cutoff (here $M_{co}=0.005$) makes all low-Mach solvers behave
  consistently.
