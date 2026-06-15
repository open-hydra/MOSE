# Numerical Dissipation and Low-Mach Flows

This page explains **how the numerical dissipation of a Riemann solver behaves
as the Mach number drops**, why most "standard" solvers lose accuracy in the
incompressible limit, and what that means in practice when you pick a solver and
the settings coupled to it (CFL, reconstruction, steady vs time-accurate).

It is a companion to the [Riemann Solvers](riemann-solvers.md) reference: that
page describes *what each solver computes*; this page describes *how much
artificial viscosity it adds and when that hurts you*.

!!! note "Scope"
    Everything below concerns **smooth, low-speed flow** ($M \lesssim 0.1$),
    e.g. natural convection, combustor cold-flow, vortex transport, aeroacoustic
    near-fields. For shock-dominated flows the opposite concern (enough
    dissipation for stability) dominates.

!!! warning "Low dissipation costs CFL"
    The solvers that are most accurate at low Mach add the least artificial
    viscosity — which is also what stabilises the explicit march. **SLAU/SLAU2**
    in particular (~400× less momentum dissipation than HLLC/Roe) typically need
    **CFL ≤ 0.2**, where HLLC/LMRoe tolerate ~0.5. If a low-dissipation solver
    diverges, lower the CFL before changing solver. See
    [Numerical Dissipation and Low-Mach Flows](low-mach-dissipation.md) for the
    full analysis and the dissipation/CFL/time-integration coupling.

---

## The low-Mach accuracy problem

As $M \to 0$ the compressible Euler equations converge to the incompressible
equations. A *good* scheme should reproduce that limit on a fixed mesh. Most
upwind schemes do **not**: their artificial viscosity on the velocity field
grows like $1/M$, so the discrete solution drifts away from the incompressible
limit as the Mach number falls.

A discrete asymptotic analysis (Guillard & Viozat 1999; Rieper 2011) identifies
the culprit precisely. Write the interface dissipation of a Roe-type flux and
expand in $M$. The momentum equation picks up a term proportional to the **jump
in the normal velocity across the face**, $\Delta v_n$, scaled by the acoustic
speed:

$$
\text{(momentum dissipation)} \;\sim\; \tfrac{1}{2}\,\rho\,a\,\Delta v_n .
$$

Since $\Delta v_n = \mathcal{O}(M a)$ while the *physical* momentum-flux
variation is $\mathcal{O}(M^2 a^2)$, this dissipation is $\mathcal{O}(1/M)$ too
large. The consequence is a spurious pressure field of the wrong order
($\mathcal{O}(M)$ instead of $\mathcal{O}(M^2)$) and rapid, unphysical decay of
vortical structures.

!!! warning "This is not cured by mesh refinement"
    The artificial kinematic viscosity is $\nu_{\text{num}} \sim a\,\Delta x$.
    A vortex of size $L$ resolved by $N = L/\Delta x$ cells survives one
    turnover only if $N \gtrsim 1/M$. At $M = 10^{-3}$ that is ~1000 cells
    across the vortex — impractical. The fix must come from the *scheme*, not
    the mesh.

---

## Key result: HLLC is **not** a low-Mach solver

A common misconception is that HLLC, being a "better" (contact-resolving)
solver than Roe, behaves better at low Mach. **It does not.** A direct
computation of the two interface fluxes shows they share the *same* offending
term.

For a subsonic interface in the low-Mach limit, the HLLC star pressure is

$$
p^\ast \;=\; \overline{p} \;-\; \tfrac{1}{2}\,\rho\,a\,\Delta v_n ,
$$

which is exactly the $\tfrac{1}{2}\rho a\,\Delta v_n$ term Rieper flags in the
Roe scheme. Evaluating both fluxes numerically on an identical low-Mach state
($M = 10^{-3}$, uniform density) confirms the equivalence to round-off:

| normal-momentum dissipation (flux − central) | value |
|---|---|
| analytic prediction $-\tfrac{1}{2}\rho a\,\Delta v_n$ | $-3.50000\times10^{-5}$ |
| plain Roe | $-3.50002\times10^{-5}$ |
| HLLC (Batten) | $-3.50010\times10^{-5}$ |
| **ratio HLLC / Roe** | **1.00002** |

**Interpretation.** Plain HLLC and plain Roe add identical leading-order
low-Mach dissipation. HLLC often *looks* better in practice only because it
fails **gracefully** — it over-diffuses into a smooth, symmetric (but decaying)
solution — whereas an unfixed or mis-tuned low-Mach scheme can fail
**catastrophically** (odd–even / checkerboard decoupling). Looking intact is not
the same as being accurate.

The genuine low-Mach cures change the dissipation itself:

- **LMRoe** scales the normal-velocity jump by $\min(M,1)$, directly removing
  the $\rho a\,\Delta v_n$ term while leaving the eigen-structure untouched.
- **Preconditioned solvers** (HLLC-PC) rescale the whole acoustic dissipation to
  the convective speed.
- **SLAU / SLAU2** multiply the pressure-diffusion by $\chi = (1-\tilde M)^2$,
  which vanishes as $M \to 0$.

See the [Riemann Solvers](riemann-solvers.md) page for the formulation of each.

---

## How the cures rescale the dissipation

Measuring the full dissipation vector $|A|(\mathbf U_R-\mathbf U_L)$ on a
low-Mach vortex state ($M=10^{-3}$) shows what each fix does, component by
component, relative to plain Roe:

| component | plain Roe | LMRoe ($\phi=\min(M,1)$) |
|---|---|---|
| density / mass | 1.00 | ≈ 1.0 (unchanged) |
| **normal momentum** | 1.00 | **≈ $M$ (× $10^{-3}$)** |
| tangential momentum | 1.00 | ≈ 1.0 (unchanged) |
| energy | 1.00 | ≈ 1.0 (unchanged) |

LMRoe surgically removes the over-dissipation **only** on the normal-velocity
(acoustic) component, exactly as the asymptotic analysis prescribes. Density,
shear, and energy dissipation are left at their physically correct level.

!!! tip "A floor is required"
    Pushing the scaling all the way to zero at stagnation points removes the
    acoustic pressure–velocity coupling and re-introduces checkerboard modes.
    LMRoe therefore keeps a small floor, $\phi = \max(\phi_{\min},\min(M,1))$
    with $\phi_{\min}\approx 0.05$ (as in SU2's L2Roe/LMRoe). Lower the floor for
    less diffusion (and more risk of noise), raise it for more robustness.

---

## Case study: Gresho vortex at $M = 10^{-3}$

The Gresho vortex is a **steady** solution (centrifugal balance), so the exact
answer is the initial vortex preserved unchanged; any departure is pure
numerical error. Run on an 80×80 mesh, MUSCL + limiter, to a pseudo-steady
state, the solvers rank exactly as the dissipation analysis predicts:

| solver | peak $|u_\phi|$ (exact 1.0) | kinetic energy vs exact | character |
|---|---|---|---|
| **LMRoe** | **1.03** (preserved) | **+1.6 %** | amplitude preserved; mildly under-dissipative (slight overshoot/scatter) |
| **HLLC-PC** | 0.92 | **−0.6 %** | cleanest convergence + best KE retention; rounds the velocity peak |
| **HLLC** | 0.91 | −2.9 % | smooth, symmetric, **but spinning down** |
| **SLAU / SLAU2** | — | — | most accurate dissipation, but **diverges unless CFL is lowered** (see below) |

Both low-Mach schemes (LMRoe, HLLC-PC) clearly beat plain HLLC: HLLC loses ~9 %
of its peak velocity and ~3 % of kinetic energy to the $\rho a\,\Delta v_n$
viscosity, while LMRoe and HLLC-PC hold the vortex. The visual "intactness" of
plain HLLC is smooth over-diffusion, not accuracy.

---

## Coupled choices — read this before changing the solver

Choosing a low-dissipation solver is never an isolated decision. The amount of
artificial viscosity a scheme adds is also what **stabilises the explicit time
march and suppresses odd–even modes**, so lowering it forces compensating
changes elsewhere.

### Dissipation level ↔ CFL number

**Less dissipation ⇒ smaller stable CFL.** The schemes that are most accurate at
low Mach are also the least forgiving:

| solver | relative dissipation | practical CFL (explicit, this case) |
|---|---|---|
| HLLE / LLF | very high | large, very robust |
| HLLC, HLLC-PC, LMRoe | moderate | ~0.5 |
| **SLAU / SLAU2** | **lowest (by design)** | **~0.05** (cfl 0.2 diverges slowly; see note) |

SLAU's momentum dissipation was measured at **~400× smaller** than HLLC/Roe on
the same state — that is its entire selling point at low Mach, and also why it
is the first to go unstable. On the steady Gresho case SLAU2 needed **CFL ≈ 0.05**
for clean convergence; at CFL 0.2 it is *slowly* unstable — the residual only
turns and grows after ~500–600 iterations, so a short run looks fine and a long
one diverges. If a low-dissipation solver diverges, **lower the CFL first, and
always confirm stability over a long run**, not a few hundred steps.

### Reconstruction ↔ dissipation

Higher-order MUSCL reconstruction shrinks the face jumps $\Delta v_n$, which *reduces* the low-Mach error for every solver (Thornber et al. 2008) — but it also reduces the stabilising dissipation. A low-Mach scheme on a smooth flow at 2nd order can look fine even with HLLC; the differences sharpen as you refine or run longer.

### Steady vs time-accurate

- **Time-accurate low-Mach:** use **LMRoe** or **SLAU2**. These are pure flux
  modifications — they do not touch the time integration and remain
  time-accurate at the standard CFL.
- **Steady / pseudo-time low-Mach:** **HLLC-PC** and other preconditioned paths
  are available and converge cleanly. Note that **preconditioning is a
  steady-state technique in MOSE**: it pairs the modified flux with a
  preconditioned residual (`integration-variables = prec`) and is *not*
  available for time-accurate runs. Selecting a preconditioned flux with
  `integration-variables = cons` couples a slow (preconditioned) dissipation to
  a fast (acoustic) time march — an inconsistent combination that is unstable.

### Accuracy ↔ robustness

The low-Mach fix trades one failure mode for another:

- **Too much** dissipation (plain HLLC): the vortex survives visually but
  **decays** (loss of amplitude / kinetic energy).
- **Too little** dissipation (un-floored LMRoe, SLAU at high CFL): the amplitude
  is preserved but **odd–even/checkerboard noise** appears.

The floor on the LMRoe factor and the CFL on SLAU are the knobs that balance
these two. Start conservative (floor $\approx 0.05$, CFL $\approx 0.2$) and relax
only after confirming a clean solution.

---

## Quick guidance for low-Mach simulations

| Situation | Recommended | Why |
|---|---|---|
| Time-accurate, $M \lesssim 0.1$ | **LMRoe** or **SLAU2** | flux-only fix, time-accurate, robust |
| Steady / pseudo-time, $M \lesssim 0.1$ | **HLLC-PC** or **LMRoe** | cleanest convergence; HLLC-PC best KE retention |
| All-speed (low + transonic + shocks) | **SLAU2** | single solver across the whole range |
| Best vortex/amplitude preservation | **LMRoe** | removes exactly the offending term |
| A solver diverges at low Mach | lower **CFL** first, then raise the LMRoe floor | low dissipation needs a smaller step |
| "It looks fine with HLLC" | verify quantitatively | HLLC over-diffuses; check peak velocity / KE, not just contours |

---

## References

1. F. Rieper, "A low-Mach number fix for Roe's approximate Riemann solver,"
   *J. Comput. Phys.*, 230, 2011.
2. H. Guillard, C. Viozat, "On the behaviour of upwind schemes in the low Mach
   number limit," *Comput. Fluids*, 28, 1999.
3. B. Thornber, A. Mosedale, D. Drikakis, D. Youngs, R. Williams, "An improved
   reconstruction method for compressible flows with low Mach number features,"
   *J. Comput. Phys.*, 227, 2008.
4. E. Shima, K. Kitamura, "Parameter-free simple low-dissipation AUSM-family
   scheme for all speeds," *AIAA J.*, 49(8), 2011.
5. J. M. Weiss, W. A. Smith, "Preconditioning applied to variable and constant
   density flows," *AIAA J.*, 33(11), 1995.
