# Riemann Solvers

The numerical convective flux at each cell interface is obtained by
solving a **Riemann problem** between the reconstructed left ($L$) and right
($R$) states.  MOSE provides several solvers grouped into five families,
each with distinct characteristics regarding accuracy, computational cost, and robustness.

Solvers are selected at run time via the input file and can optionally be combined with
the [shock detector](numerics.md#shock-detection)

**Why multiple solvers?**

- **Accuracy:** Different solvers capture different flow features (shocks, contacts, shear layers)
- **Robustness:** Dissipativeness trades accuracy for stability
- **Cost:** Approximate solvers are much faster than the exact solver
- **Regime:** Low-Mach flows, hypersonic flows, and boundary layers have different requirements

This page describes the mathematical foundations and practical considerations for choosing
an appropriate solver for your simulation.

---

## Notation

| Symbol | Definition |
|--------|-----------|
| $v_n$ | Normal velocity $\mathbf{v}\!\cdot\!\hat{\mathbf{n}}$ |
| $a$ | Speed of sound |
| $M$ | Mach number $v_n / a$ |
| $H_0$ | Total specific enthalpy |
| $\dot{m}$ | Interface mass flux per unit area |

Subscripts $L$ and $R$ denote left and right states; subscript $\frac{1}{2}$
denotes the interface value.

---

## AUSM Family

**Physical principle:**
AUSM solvers exploit the **splitting of flux into convective (mass-driven) and pressure-driven components**.
Instead of decomposing along eigenvectors, AUSM splits the interface velocity into upwinded parts
and the interface pressure via smoothed functions, avoiding expensive eigenvalue calculations.

**Mathematical foundation:**  
- **Convective part:** Advects mass and species based on directional velocity  
- **Pressure part:** Local pressure from left and right, weighted smoothly across incompressible--hypersonic Mach range  
- **Advantage:** Naturally handles low-Mach flows by decoupling velocity and pressure solutions

The Advection Upstream Splitting Method (AUSM) family splits the
convective flux into a **convective** (mass-flux) part and a **pressure**
part, treated independently.  This avoids a full eigenvalue decomposition
and naturally handles low-Mach flows.

### AUSM+ (Liou, 2006)

Mach-based splitting extending standard AUSM (+) with improved low-Mach behavior.

**Key features:**

1. **Pressure-dependent Mach splitting** (for $|M| < 1$):

   $$ \mathcal{M}^{\pm}(M) =
     \pm\tfrac{1}{4}(M \pm 1)^2 + \beta\,(M^2 - 1)^2,\quad\beta = \tfrac{1}{8} $$

2. **Pressure splitting** (for $|M| < 1$):
   $$ \mathcal{P}^{\pm}(M) =
     \tfrac{1}{4}(M \pm 1)^2(2 \mp M) + \alpha\,M\,(M^2 - 1)^2, \quad \alpha = \tfrac{3}{16} $$

3. **Interface pressure:** Combination of left and right contributions,
   $$ p_{\frac{1}{2}} = \mathcal{P}^+(M_L)\,p_L + \mathcal{P}^-(M_R)\,p_R $$

**Advantages:** Conservative, handles subsonic and supersonic regions, simple.  
**Disadvantages:** Can exhibit odd–even decoupling at very low Mach.

### AUSM+-up

Extends AUSM+ with:

- **Pressure diffusion** in the mass flux — scales with $1/M^2$ and
  eliminates odd–even decoupling at low Mach.
- **Velocity diffusion** in the interface pressure — improves robustness
  across the subsonic/transonic transition.

Recommended for flows with mixed subsonic and supersonic regions.

### AUSM+-up2

Further improvement for hypersonic flows with refined Mach-dependent
dissipation scaling.  Preserves the accuracy of AUSM+-up in the
low-Mach limit while adding stability at strong-shock crossings.

---

## HLL Family

**Physical principle:**
HLL solvers approximate the exact Riemann solution by **bracketing the entire wave fan between two bounding waves**.
Instead of resolving each acoustic/contact wave individually (as in Roe), HLL computes a single
intermediate state across the entire fan, leading to unconditional stability.

**Key variants:**  
- **HLLE:** Two-wave solver; very robust but highly dissipative  
- **HLLEM:** Adds anti-diffusion correction to capture contacts sharper  
- **HLLC:** Three waves (contact explicitly resolved); sharp contacts, may oscillate at shocks  
- **HLLE++, HLLC+:** Tramel variants with improved eigenvalue handling and shock detection  

The Harten–Lax–van Leer (HLL) solvers form a family of approximate Riemann solvers
that compute the flux by integrating over a simplified wave structure. They avoid
the expensive eigenvalue decomposition of the Roe solver while maintaining strong
stability properties. The key idea is to **replace the exact Riemann fan with two or
three waves** (left shock/rarefaction, contact/shear, and right shock/rarefaction)
and integrate the conserved variables across these waves.

### Wave Speed Estimates (Batten)

All HLL variants solve a limited Riemann problem by estimating **two bounding wave speeds**
$S_L$ (leftmost) and $S_R$ (rightmost) that bracket the entire wave structure.
The Batten et al. (1997) estimate is the industry standard:

$$
S_L = \min\!\bigl(0,\;v_{n,L} - a_L,\;\tilde{v}_n - \tilde{a}\bigr),\qquad
S_R = \max\!\bigl(0,\;v_{n,R} + a_R,\;\tilde{v}_n + \tilde{a}\bigr)
$$

where $\tilde{v}_n$ and $\tilde{a}$ are **Roe-averaged** values. The Roe average with
sqrt-weighting $w = \sqrt{\rho_L/\rho_R}$ is:

$$\tilde{v}_n = \frac{w\,v_{n,L} + v_{n,R}}{w + 1}, \quad \tilde{a} = a(\tilde{\rho}, \tilde{T}, \tilde{Y}_s)$$

**Purpose of Batten's bounds:**
- $\max(0, \ldots)$ ensures $S_R \geq 0$ (rightmost wave travels right or stalls)  
- $\min(0, \ldots)$ ensures $S_L \leq 0$ (leftmost wave travels left or stalls)  
- The Roe averages provide entropy-consistent estimates that capture expansion fans correctly  

### HLLE

The simplest two-wave variant: assumes the entire Riemann fan collapses to a single
intermediate state. The flux is computed as a conservative average:

$$
\mathbf{F}_{\frac{1}{2}} =
\frac{S_R\,\mathbf{F}_L - S_L\,\mathbf{F}_R
+ S_L\,S_R\,(\mathbf{U}_R - \mathbf{U}_L)}{S_R - S_L}
$$

**Characteristics:**
- **Unconditionally stable** (mathematically proven; often called entropy-stable)  
- Highly dissipative, especially on **contact discontinuities** and **shear layers**  
- Diffuses small-amplitude acoustic waves and material interfaces  
- Recommended as a fallback when more sophisticated solvers fail  

**Use cases:** Extremely strong shocks, near-vacuum flows, severe transients.

### HLLEM (Einfeldt's Modification)

Improves HLLE's contact-wave resolution via an **anti-diffusion correction** inspired by
the exact Roe solver. When $S_L < 0 < S_R$ (subsonic flow):

1. Decompose left/right state jumps into eigencomponents (density, shear, acoustic):
   $$ \mathbf{U}_R - \mathbf{U}_L = \sum_i \alpha_i\,\mathbf{r}_i $$
   where $\mathbf{r}_i$ are Roe eigenvectors.

2. Compute anti-diffusion weight:
   $$ \delta = \frac{\tilde{a}}{\tilde{a} + |\tfrac{1}{2}(S_L + S_R)|} $$
   (stronger in subsonic regions, weaker when mean wave speed is large)

3. Subtract correction from HLLE flux:
   $$ \mathbf{F}_{\frac{1}{2}}^{\text{HLLEM}} = \mathbf{F}_{\frac{1}{2}}^{\text{HLLE}} - C\,\delta\,\sum_i |\lambda_i|\,\alpha_i\,\mathbf{r}_i $$
   where $C = S_L S_R / (S_R - S_L)$ is a normalization factor.

**Advantages:** Significantly sharper contact resolution than HLLE.  
**Disadvantages:** Weaker at strong shocks than HLLC; intermediate cost.

### HLLC (Batten, 1997)

**Three-wave solver:** Explicitly resolves the **contact discontinuity** to overcome HLLE's 
excessive diffusion of material interfaces.

**Contact-wave speed:** Computed from pressure and momentum balance:

$$
S^\ast =
\frac{p_R - p_L + \rho_L\,v_{n,L}(S_L - v_{n,L})
                 - \rho_R\,v_{n,R}(S_R - v_{n,R})}
     {\rho_L(S_L - v_{n,L}) - \rho_R(S_R - v_{n,R})}
$$

**Intermediate pressure:** Pressure continuity constraint gives

$$
p^\ast =
\rho_L\,(v_{n,L} - S_L)(v_{n,L} - S^\ast) + p_L
$$

**Advantages:**  
- Exact contact-wave resolution for smooth flows  
- Sharp material interface tracking  
- Moderate cost increase over HLLE

**Disadvantages:**  
- Can show oscillations at strong shocks (carbuncle)  
- Requires stabilization for robustness  

**Recommendation:** Use HLLC+ instead for adaptive shock-aware fallback.

### HLLE++ (Tramel, 2009)

Tramel's improved variant of HLLE for better shear-layer resolution.

**Key modification:**

Instead of the standard Batten wave-speed bound using $\max(\tilde{v}_n, \tilde{a})$,
HLLE++ uses $\tilde{v}_n$ directly as the lower eigenvalue bound.
This change preserves **Eulerian shear waves** (pure tangential velocity jumps)
that would otherwise be artificially diffused by standard HLLE.

**Eigenvalue scaling:**

HLLE++ blends Roe eigenvalues with HLLE eigenvalues via a parameter $\beta$:

$$\lambda_i^{++} = \beta\,\lambda_i^{\text{Roe}} + (1-\beta)\,\lambda_i^{\text{HLLE}}$$

For the normal eigenvalue: $ \lambda_1^{\text{Roe}} = |\tilde{v}_n|\quad \text{(not } \max(|\tilde{v}_n|, \tilde{a})\text{)} $

For acoustic eigenvalues, Harten–Hyman entropy correction is applied.

**Advantages:**  
- Sharp contact and shear resolution at moderate additional cost

**Disadvantages:**  
- Less robust at strong normal shocks than HLLE

### HLLC+ (Tramel, 2009) with Shock Detection

Hybrid solver that **adaptively blends HLLC and HLLE** using a shock-detection parameter.

**Hybrid blending formula:**

Given a shock-detection indicator $\beta \in [0,1]$ computed from pressure gradients
(see [Shock Detection](numerics.md#shock-detection)):

$$\mathbf{F}_{\frac{1}{2}} = \beta\,\mathbf{F}_{\text{HLLC}} + (1-\beta)\,\mathbf{F}_{\text{HLLE}}$$

| Regime | $\beta$ | Solver | Characteristic |
|:------:|:-------:|:------:|----------------|
| Smooth flow (no shock) | $\approx 1.0$ | HLLC | Sharp contact waves, lower dissipation |
| Weak shock | $0.4 \text{--} 1.0$ | HLLC dominant | Balanced |
| Strong shock / discontinuity | $\approx 0.0$ | HLLE | Maximum dissipation, unconditionally stable |

**Use cases:**  
- Flows with shock-boundary layer interactions  
- Unsteady shock-contact interactions  
- Carbuncle-prone geometries (e.g., blunt bodies) where HLLC exhibits instabilities  

**Recommended:** Excellent all-around choice for compressible flows with mixed subsonic/supersonic regions and shocks.

### Rotated HLLC / HLLE

Alternative hybrid approach using a **frame rotation** aligned with the local
velocity-difference direction.

**Concept:**

Instead of a single scalar $\beta$, this variant rotates the Riemann problem
into a frame where the velocity vector $\mathbf{v}_R - \mathbf{v}_L$ is normal.
The rotated problem is then solved with a blend of HLLC and HLLE, and the result
is rotated back.

**Advantages over scalar-$\beta$ blending:**  

- Better handling of **oblique shocks** and **shear layers** at non-normal angles  
- Natural detection of flow-aligned discontinuities  
- Smooth transition between solvers in multi-dimensional problems  

**Trade-off:** Additional rotation overhead; primarily beneficial for flows with
complex shock orientations.

**Comparison:**
While HLLC+ uses shock detection in the original frame, the rotated variant
senses shock proximity through the velocity-difference orientation, providing
a complementary robustness mechanism for angled discontinuities.

---

## Godunov (Exact Solver)

The exact Riemann solver computes the entropy-satisfying weak solution to the
Riemann problem by iterating on the contact velocity and pressure simultaneously.

**Algorithm**:

1. Rotate left/right states into normal–tangential frame.
2. Initial guess from Riemann invariants:
   $\displaystyle R_L^+ = v_{n,L} + \frac{2\,a_L}{\gamma - 1},\qquad
   R_R^- = v_{n,R} - \frac{2\,a_R}{\gamma - 1}$
3. Newton–Raphson on $v_n^\ast$ (contact velocity) until the
   pressure jump $|p_L^\ast - p_R^\ast|$ vanishes.
   For each iteration the 1-D shock/rarefaction relations of each wave
   family are evaluated.
4. Sample the full wave structure at the interface ($x/t = 0$).
5. Compute fluxes from the interface state.

Most accurate solver available; higher computational cost
(Newton–Raphson with up to 1000 iterations).

---

## Lax–Friedrichs Family

Simple, globally stable first-order solvers based on Lax–Friedrichs averaging.
They serve as robust fallback options and as dissipative building blocks in
more sophisticated schemes.

### Local Lax–Friedrichs (LLF) / Rusanov

Average flux with global dissipation proportional to the maximum eigenvalue magnitude:

$$
\mathbf{F}_{\frac{1}{2}} =
\tfrac{1}{2}\bigl(\mathbf{F}_L + \mathbf{F}_R\bigr)
- \tfrac{1}{2}\,\Lambda_{\max}\,(\mathbf{U}_R - \mathbf{U}_L)
$$

where the spectral radius (worst-case wave speed) is

$$
\Lambda_{\max} = \max\!\bigl(|v_{n,L} - a_L|,\;|v_{n,L} + a_L|,\;
                              |v_{n,R} - a_R|,\;|v_{n,R} + a_R|\bigr)
$$

**Characteristics:**  
- **Unconditionally stable** and robust  
- Simplest implementation (no eigenvalue decomposition)  
- Extreme dissipation; smears all features (shocks, contacts, acoustics)  
- Useful as a debugging baseline and emergency fallback  

**Cost:** Very low.

---

## SLAU Family

**Physical principle:**
SLAU (Simple Low-dissipation AUSM) combines aspect-ratio-preserving Mach splitting with compressibility-dependent dissipation. The key innovation is a compressibility parameter $\chi = (1 - \tilde{M})^2$ that suppresses artificial viscosity near $M \to 0$ while maintaining shock-capturing at high Mach.

**Mathematical structure:**  
- **Mass flux:** Upwind-weighted convection based on Mach-dependent splitting  
- **Pressure reconstruction:** Blends AUSM+ splitting with compressibility-aware dissipation  
- **Low-Mach limit:** $\chi \to 1$ reduces dissipation to near-incompressible accuracy  
- **High-Mach limit:** $\chi \to 0$ recovers good shock resolution

The SLAU solvers are designed for all-speed flows, from incompressible subsonic through hypersonic regimes. They provide superior accuracy at low Mach while maintaining robustness at shocks, a difficult balance to achieve.

**Unique advantage:** Single solver suitable for entire Mach range without solver switching.

### Common Low-Mach Stabilization Parameters

Both SLAU variants use Mach-dependent scaling to suppress dissipation in low-speed flows:

**Density-weighted velocity scale:**
$\bar{V}_n = \frac{\rho_L\,|v_{n,L}| + \rho_R\,|v_{n,R}|}{\rho_L + \rho_R}$

**Reference Mach number** (capped at unity):
$\tilde{M} = \min\!\left(1,\;\sqrt{\frac{|\mathbf{v}|^2}{2\,a_F^2}}\right)$
where $a_F = \sqrt{a_L\,a_R}$ is the face sound speed.

**Compressibility parameter** (vanishes at $\tilde{M} \to 0$):
$\chi = (1 - \tilde{M})^2$

### SLAU (Shima & Kitamura, 2011)

Simple low-dissipation AUSM with Mach-dependent scaling.

**Mass flux:** Weighted average of left/right contributions with velocity correction:  

$$ \dot{m}_L = \rho_L\,v_{n,L},\quad \dot{m}_R = \rho_R\,v_{n,R} $$

$$ \dot{m} = \frac{1}{2}(\dot{m}_L + \dot{m}_R) + \frac{1}{2}(\dot{m}_L - \dot{m}_R)\,\text{sgn}(v_{n,L} + v_{n,R}) - \frac{\chi}{2a_F}(p_R - p_L) $$

where the third term provides compressibility-dependent pressure correction.

**Interface pressure:** Weighted by pressure-splitting functions:  

$$  p_{\frac{1}{2}} = \frac{1}{2}(p_L + p_R) + \frac{1}{2}(\beta_L^+ - \beta_L^-)(p_L - p_R)  $$

with 
$\beta^{+}(M) = \tfrac{1}{4}(2 - M)(M + 1)^2, \quad \beta^{-}(M) = \tfrac{1}{4}(2 + M)(M - 1)^2 \quad (|M| < 1)$

**Characteristics:** Clean formulation, excellent subsonic accuracy, simple implementation.

### SLAU2 (Kitamura & Shima, 2013)

Enhanced SLAU with velocity-based shock detection for better shock resolution.

**Interface pressure:** Improves dissipation scaling near normal shocks:  

$$ p_{\frac{1}{2}} = \frac{1}{2}(p_L + p_R) + \frac{1}{2}(\beta_L^+ - \beta_L^-)(p_L - p_R) + (\beta_L^+ + \beta_L^- - 1) \cdot f_v \cdot a_F \cdot \frac{1}{2}(\rho_L + \rho_R) $$

where the shock-detection term $f_v = \sqrt{\frac{1}{2}(|\mathbf{v}_L|^2 + |\mathbf{v}_R|^2)}$ represents a velocity-based shock indicator normalizing the dissipation strength.

**Key improvement:** The velocity weighting $f_v$ provides automatic strong dissipation near normal shocks (high velocity magnitude = shock present) while maintaining low-dissipation accuracy in smooth regions.

**Advantages:**  
- Best low-Mach performance in SLAU family (subsonic flows)
- Converges to exact solver near strong normal shocks
- Seamless all-speed handling from incompressible to hypersonic
- Conservative and naturally handles multi-species flows

---

## Shock Detection for Hybrid Solvers

Several solvers (HLLC+, HLLE++) employ adaptive dissipation controlled by
a shock-detection parameter $\beta \in [0,1]$. This parameter is computed from local pressure gradients using a **Jameson-type sensor**:

$$
s = \max\!\left( \left| \frac{p_E - 2p_C + p_W}{p_E + 2p_C + p_W} \right|,
       \left| \frac{p_N - 2p_C + p_S}{p_N + 2p_C + p_S} \right|,
       \left| \frac{p_T - 2p_C + p_B}{p_T + 2p_C + p_B} \right| \right)
$$

where subscripts denote center (C), east (E), west (W), north (N), south (S),
top (T), bottom (B) cells.

**Shock detection formula:**

$$
\beta = \begin{cases}
1 - \tanh(10\phi^3) & \text{if } s < 1/\Delta \quad (\text{smooth region}) \\
0 & \text{if } s \ge 1/\Delta \quad (\text{shock detected})
\end{cases}
$$

where $\phi = \max(s/\Delta, 0)$ and $\Delta = 20$ is a calibration parameter.

**Interpretation:**

| $\beta$ | Region Type | Solver Behavior |
|:-------:|:----------:|-----------------|
| $\beta \approx 1.0$ | Smooth flow, no shock | Use HLLC (contact-capturing) or HLLE++ (shear-capturing) |
| $0 < \beta < 1$ | Weak shock, transition zone | Blend toward dissipative scheme |
| $\beta \approx 0.0$ | Strong shock | Use HLLE (maximum stability) |

This adaptive approach **combines accuracy in smooth regions with robustness at shocks**,
making hybrid solvers like HLLC+ ideal for general-purpose simulations.

---

## Practical Solver Selection Guide

### Quick Reference

**For most compressible-flow simulations:** Start with **HLLC+** (or **SLAU/SLAU2** for all-speed flows).
HLLC+ balances accuracy, robustness, and cost.

**By application:**

| Flow Type | Recommended Solver | Reason |
|---|---|---|
| **Shock-dominated** (hypersonic, detonation) | HLLE, HLLC+ | Robust; HLLC+ adds contact resolution |
| **Smooth subsonic** (combustor, inlet) | SLAU2, AUSM+-up | Low-Mach accuracy; SLAU2 best |
| **Mixed subsonic/transonic** | HLLC+, AUSM+-up | HLLC+ adaptive; AUSM+-up is simpler |
| **Boundary layers + shocks** | HLLC+ | Shock detection aids shear resolution |
| **Very low Mach (incompressible limit)** | SLAU2, PLLF | Preconditioning essential |
| **Unsteady shock interaction** | HLLC+ | Shock detection tracks transients |
| **Emergency (solver divergence)** | HLLE, LLF | Maximum stability; accept extra diffusion |

### Convergence Tips

1. **Start robust, refine if needed:** Begin with HLLE or SLAU2, switch to HLLC+ after convergence behaves.
2. **Monitor divergence:** If solver crashes, switch to HLLE (fallback).
3. **Verify solution structure:** Check shock positions, contact discontinuities, and Mach numbers.
4. **Mesh refinement:** Coarse meshes tolerate more dissipative solvers (HLLE); fine meshes benefit from HLLC+.

---

## References

1. M.-S. Liou, "A sequel to AUSM, Part II: AUSM+-up for all speeds,"
   *J. Comput. Phys.*, 214(1), 2006.
2. E. F. Toro, M. Spruce, W. Speares, "Restoration of the contact surface
   in the HLL-Riemann solver," *Shock Waves*, 4, 1994.
3. P. Batten, N. Clarke, C. Lambert, D. M. Causon, "On the choice of
   wavespeeds for the HLLC Riemann solver," *SIAM J. Sci. Comput.*,
   18(6), 1997.
4. B. Einfeldt, C. D. Munz, P. L. Roe, B. Sjögreen, "On Godunov-type
   methods near low densities," *J. Comput. Phys.*, 92, 1991.
5. E. Shima, K. Kitamura, "Parameter-free simple low-dissipation
   AUSM-family scheme for all speeds," *AIAA J.*, 49(8), 2011.
6. K. Kitamura, E. Shima, "Towards shock-stable and accurate hypersonic
   heating computations: A new pressure flux for AUSM-family schemes,"
   *J. Comput. Phys.*, 245, 2013.
7. R. Tramel, R. Nichols, P. Buning, "Addition of improved shock-capturing
   schemes to OVERFLOW 2.1," *19th AIAA Computational Fluid Dynamics*,
   2009. AIAA Paper 2009-3988.
