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
| $M_{co}$ | Low-Mach cutoff Mach (shared dissipation floor, input `riemann-options-Mco`) |

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

### AUSM+ (Liou, 1996)

Mach-based splitting extending standard AUSM with improved behavior near $|M|=1$.

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
**Disadvantages:** Can exhibit odd–even decoupling at very low Mach; not all-speed.
**Role in MOSE:** kept as the simple AUSM **baseline** (teaching/reference). For
production all-speed work use **AUSM+M** below.

### AUSM+M (Chen et al., 2020)

An improved AUSM-family all-speed scheme (Chen, Cai, Xue, Wang & Yan, *Applied
Mathematical Modelling* 77, 2020). It fixes the two weaknesses of AUSM+-up — the
`Kp/fa` stagnation time-step restriction and carbuncle at high Mach — with three
ingredients:

- **Pressure-diffusion mass flux** whose denominator has **no Mach number**
  (Eq. 14), so no stagnation singularity and a larger stable time step than
  AUSM+-up at low speed: $M_p = -\tfrac12(1-f)\,\frac{\Delta p}{\rho_{1/2}c_{1/2}^2}(1-g)$,
  with $f = \tfrac12(1-\cos\pi M)$.
- **Pressure flux** with a scaling Mach function $f_o$ for low-Mach accuracy,
  plus a multidimensional **velocity-diffusion** term gated by a pressure-ratio
  shock sensor $g$ for carbuncle control.
- **AUSMPW+ numerical sound speed** (Kim et al.) for correct oblique shocks and
  no unphysical expansion shocks.

The scaling function is floored, $f_o = \min(1,\max(f, M_{co}))$, by the shared
low-Mach cutoff Mach $M_{co}$ — see [Low-Mach Cutoff Mach](#low-mach-cutoff-mach-m_co).
The paper uses $f_o=\max(f,M_\infty^2)$ with $M_\infty$ the *freestream* Mach; MOSE
floors it **linearly** in $M_{co}$ so a single input value is consistent across all
low-Mach solvers.

---

## HLL Family

**Physical principle:**
HLL solvers approximate the exact Riemann solution by **bracketing the entire wave fan between two bounding waves**.
Instead of resolving each acoustic/contact wave individually (as in Roe), HLL computes a single
intermediate state across the entire fan, leading to unconditional stability.

**Key variants:**  
- **HLLE:** Two-wave solver; very robust but highly dissipative  
- **HLLC:** Three waves (contact explicitly resolved); sharp contacts, may oscillate at shocks  
- **HLLC-PC:** Low-Mach preconditioned HLLC (steady)  
- **HLLC+Chen:** All-speed low-Mach shock-stable HLLC (Chen et al. 2020)  
- **HLLE++, HLLC+Tramel:** Tramel variants with improved eigenvalue handling and shock detection  

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

**Recommendation:** Use HLLC+Tramel instead for adaptive shock-aware fallback.

### HLLC+Chen (Chen et al., 2020)

A genuine **all-speed** HLLC (Chen, Lin, Li & Yan, *SIAM J. Sci. Comput.* 42(4),
2020). It is the standard HLLC star flux plus a single correction term
(their Eq. 3.25):

$$
\mathbf{F}^{\ast+}_K = \mathbf{F}^{\ast}_K +
\frac{\varphi_L \varphi_R}{\varphi_R - \varphi_L}
\bigl[\,0,\ (f^\ast\!-\!1)\,\Delta U\,\hat{\mathbf n} + \zeta\,\Delta\mathbf{u}_t,\ (f^\ast\!-\!1)\,\Delta U\,S^\ast\,\bigr]^T
$$

with $\varphi_K = \rho_K(S_K - U_K)$. Two effects, on the normal and transverse
velocity jumps respectively:

- **Low-Mach fix:** the $f(M)$ factor (Thornber interface Mach) scales away the
  $\mathcal O(a)$ pressure dissipation that makes plain HLLC over-diffuse at low
  Mach — added to **both** momentum and energy to recover the correct $M^2$
  pressure *and* density scaling. The interface Mach is floored,
  $M \to \max(M_{co}, M)$, by the shared low-Mach cutoff Mach
  (see [Low-Mach Cutoff Mach](#low-mach-cutoff-mach-m_co)).
- **Carbuncle cure:** a transverse shear dissipation, active only near strong
  shocks via the **multidimensional** pressure-ratio sensor $g = 1 - h^{M}$
  (where $h$ scans the neighbouring interfaces, shared with AUSM+M's sensor) and
  the factor $\zeta = g\,S_K/(S_K - S^\ast)$. The same $g$ also restrains the
  low-Mach $f^\ast$ near shocks (in place of the paper's separate sonic sensor).

It reduces to plain HLLC at $M\ge 1$ in smooth flow; its only tuning is the shared
low-Mach cutoff $M_{co}$ (the floor on the interface Mach).

**Advantages:** all-speed accuracy + carbuncle robustness in one solver; cheap
(a small add-on to HLLC). In MOSE testing it held the Gresho vortex at
$M=10^{-3}$ (peak velocity 0.995 vs exact 1.0, vs 0.91 for plain HLLC) and gave a
Sod profile essentially identical to HLLC. **Recommendation:** excellent default
for **reacting / low-Mach + shock** problems.

### HLLC+Tramel (Tramel, 2009)

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
While HLLC+Tramel uses shock detection in the original frame, the rotated variant
senses shock proximity through the velocity-difference orientation, providing
a complementary robustness mechanism for angled discontinuities.

---

## Roe Family

**Physical principle:**
Roe solvers linearise the Riemann problem about a special **Roe-averaged state** and decompose the jump between the left and right states onto the **eigenvectors** of the resulting flux Jacobian. Each wave (two acoustic, one entropy, two shear) is upwinded by the sign of its own eigenvalue, giving sharp resolution of every wave family.

$$
\mathbf{F}_{\frac{1}{2}} =
\tfrac{1}{2}(\mathbf{F}_L + \mathbf{F}_R)
- \tfrac{1}{2}\sum_i |\lambda_i|\,\alpha_i\,\mathbf{r}_i ,
\qquad
\mathbf{U}_R - \mathbf{U}_L = \sum_i \alpha_i\,\mathbf{r}_i
$$

with eigenvalues $\lambda_i = \{\tilde v_n - \tilde a,\ \tilde v_n,\ \tilde v_n,\ \tilde v_n,\ \tilde v_n + \tilde a\}$
and Roe-averaged eigenvectors $\mathbf r_i$.

The standard Roe flux is accurate and low-cost but suffers two well-known issues:
it admits expansion shocks without an **entropy fix**, and — crucially for this
code's target applications — it is **not accurate in the low-Mach limit**: its
momentum dissipation contains a term $\sim \tfrac{1}{2}\rho a\,\Delta v_n$ that
is $\mathcal{O}(1/M)$ too large. MOSE therefore exposes the low-Mach-corrected
variant directly.

### LMRoe (Rieper, 2011)

**Low-Mach fix for Roe.** The asymptotic analysis shows the offending term is the **jump in the normal velocity** $\Delta v_n$ carried by the acoustic waves. Rieper's remedy is a one-line change: multiply that jump by the local Mach number before forming the acoustic wave strengths,

$$
\Delta v_n \;\longrightarrow\; \phi\,\Delta v_n ,
\qquad
\phi = \min\!\bigl(1,\; \tilde M\bigr),
\quad
\tilde M = \frac{|\tilde v_n| + |\tilde v_t|}{\tilde a},
$$

applied **only** to the two acoustic characteristics. Eigenvalues, eigenvectors, and the entropy/shear wave strengths are **unchanged**. This removes the $\mathcal{O}(1/M)$ momentum dissipation while leaving density, shear, and energy dissipation at their correct level.

A floor $\phi = \max(M_{co}, \min(1,\tilde M))$ is retained, set by the shared
low-Mach cutoff Mach $M_{co}$ (input `riemann-options-Mco`): without it
$\phi \to 0$ at stagnation points removes the acoustic pressure–velocity coupling
and re-introduces odd–even (checkerboard) modes. See
[Low-Mach Cutoff Mach](#low-mach-cutoff-mach-m_co) for how to choose it.

**Advantages:**
- Accurate in the incompressible limit on a fixed mesh — preserves vortices and
  acoustic-scale structures that standard Roe/HLLC diffuse away.
- Trivial, parameter-light, and **time-accurate** (a pure flux modification — it
  does not touch the time integration).
- Robust: the eigen-structure is untouched, so there is no preconditioning
  singularity at stagnation.

**Disadvantages:**
- Mildly under-dissipative if the floor is too small (amplitude overshoot /
  checkerboard noise).
- For *steady-state convergence acceleration* it offers nothing — it is an
  accuracy fix, not a preconditioner.

### MiczekRoe (Miczek, Roepke & Edelmann, 2015)

**Preconditioned Roe** for the low-Mach regime. Where LMRoe rescales a single
jump, MiczekRoe applies a full **preconditioner** $\mathbf P$ to the Roe
dissipation, $\tfrac12\,\mathbf P^{-1}|\mathbf P\mathbf A|\,(\mathbf U_R-\mathbf U_L)$,
which rebalances *all* wave speeds so the acoustic and convective scales become
comparable as $M\to 0$. The preconditioner is built from a limited local Mach

$$
\mu = \min\!\bigl(1,\ \max(M_{\text{loc}},\,M_{co})\bigr),
\qquad \delta = \tfrac{1}{\mu}-1,
$$

so $\delta\to 0$ ($\mathbf P\to\mathbf I$) recovers the standard Roe flux for
$M\ge 1$. The floor $M_{co}$ (input `riemann-options-Mco`) is the **same**
cutoff Mach used by the other low-Mach solvers; here it prevents the
$\delta=1/\mu-1\to\infty$ singularity of the preconditioner at stagnation.

**Advantages:**
- Designed to also **accelerate steady convergence** at low Mach (a genuine
  preconditioner, not only an accuracy fix).

**Disadvantages:**
- The preconditioner couples the eigen-structure, so it is more intrusive and
  less time-accurate than LMRoe; validated in MOSE at $M=0.1$, but
  checkerboard-prone at $M=10^{-3}$ under explicit time stepping (the reference
  uses implicit integration).

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

Several solvers (HLLC+Tramel, HLLE++) employ adaptive dissipation controlled by
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
making hybrid solvers like HLLC+Tramel ideal for general-purpose simulations.

---

## Low-Mach Cutoff Mach $M_{co}$

All of MOSE's low-Mach-capable solvers (LMRoe, MiczekRoe, HLLC+Chen, AUSM+M)
share a single tuning parameter, the **cutoff Mach** $M_{co}$, set through the
input key `riemann-options-Mco`. This section explains what it does and how to
choose it.

### Why a floor is needed

A compressible upwind flux adds numerical dissipation scaled by the acoustic wave
speed. Written per unit of velocity jump, the offending acoustic (pressure↔velocity
coupling) dissipation is

$$
D_{\text{acoustic}} \;\sim\; \rho\,a\,\Delta v_n ,
$$

which is $\mathcal O(a) = \mathcal O(u/M)$ — a factor $1/M$ **too large** at low
Mach. It swamps the physical convective fluxes and destroys vortices, shear
layers and acoustic-scale structures. Every all-speed scheme cures this by
multiplying that term by a Mach-like factor (LMRoe's $\phi$, AUSM+M's $f_o$,
HLLC+Chen's $f(M)$, MiczekRoe's $\mu$) so the dissipation collapses back to the
convective scale $\mathcal O(u)$.

The catch: that factor $\to 0$ as $M\to 0$. On a collocated grid the
pressure–velocity system with **vanishing** acoustic dissipation has a non-trivial
null space — the classic **odd–even / checkerboard** decoupling. So the factor
cannot be allowed to reach zero; it must be **floored**. $M_{co}$ *is* that floor.

### One knob, consistent across schemes

Each scheme has its own low-Mach factor with its own parametrization, but $M_{co}$
maps to all of them so that the **effective acoustic-dissipation coefficient**
$D$ (the coefficient of $\rho a\,\Delta v_n$) is essentially the same number:

| Scheme | Native floor | Effective $D$ |
|--------|-------------|---------------|
| LMRoe        | $\phi_{\min} = M_{co}$            | $D \approx M_{co}$ |
| HLLC+Chen    | $M_{\text{loc}} \ge M_{co}$, $f(M)\!\approx\!\sqrt5\,M$, prefactor $\approx\tfrac12\rho a$ | $D \approx M_{co}$ |
| AUSM+M       | $f_o \ge M_{co}$ (pressure-flux coeff $\approx 0.67$) | $D \approx 0.67\,M_{co}$ |
| MiczekRoe    | $\mu \ge M_{co}$ (preconditioner)  | $D \approx \mathcal O(M_{co})$ |

The only subtlety is **AUSM+M**: the original paper floors its scaling function
*quadratically*, $f_o=\max(f, M_\infty^2)$ with $M_\infty$ the freestream Mach.
That makes the input number look very different (e.g. $M_\infty=0.1$ gives the
same dissipation as $\phi_{\min}=0.01$). MOSE instead floors $f_o$ **linearly** in
$M_{co}$, so the *same* $M_{co}$ value gives matched dissipation for every solver —
you set one number and they all behave consistently.

### Choosing $M_{co}$: the two-sided window

$M_{co}$ lives inside a window bounded below by stability and above by accuracy.

**Too low → odd–even / checkerboard instability.**
Below a (grid- and time-step-dependent) threshold the acoustic coupling is too
weak to suppress the spurious checkerboard pressure mode, and the solution
develops a growing high-frequency $\pm$ pattern that eventually crashes. The
floor must stay above the representative local Mach **everywhere**, including
regions where the physical Mach itself $\to 0$ (vortex cores, stagnation
points) — which is why a *constant* floor is needed rather than $\phi=M$.

**Too high → excessive dissipation / loss of low-Mach accuracy.**
This is the question often overlooked. The artificial velocity dissipation
introduced by the floor scales as

$$
\frac{D}{M_{\text{local}}} \;\approx\; \frac{M_{co}}{M_{\text{local}}},
$$

so raising $M_{co}$ *linearly* increases the over-dissipation. Concretely:

- vortices and shear layers are smeared, and **kinetic energy decays
  unphysically** — a steady vortex spins down, a free shear layer thickens;
- small pressure/acoustic fluctuations are over-damped;
- in the limit $M_{co}\to 1$ the low-Mach correction is switched off entirely and
  you recover the **baseline compressible scheme**: perfectly stable, but with
  the original $\mathcal O(1/M)$ dissipation — i.e. no low-Mach accuracy at all.

Note the asymmetry: too high is *stable but inaccurate* (more dissipation always
helps stability), whereas too low is *outright unstable*. When in doubt, err
slightly high.

**Practical rule.** Set $M_{co}$ to **a few times the flow's peak Mach**
($M_{co}\approx 5\,M_{\text{peak}}$ is a good default). For the Gresho vortex at
$M_{\text{peak}}=10^{-3}$, $M_{co}=0.005$ sits comfortably in the window: above
the checkerboard threshold, and only $\sim 5\times$ over-dissipative at the vortex
peak (vs $\sim 50\times$ for a generic $0.05$ floor calibrated for $M\sim 10^{-2}$
flows). Because the optimum tracks the problem's Mach number, $M_{co}$ should be
revisited per case rather than left at a single hard-coded value.

---

## Practical Solver Selection Guide

### Quick Reference

**For most compressible-flow simulations:** Start with **HLLC+Tramel** (or **SLAU/SLAU2** for all-speed flows).
HLLC+Tramel balances accuracy, robustness, and cost.

**By application:**

| Flow Type | Recommended Solver | Reason |
|---|---|---|
| **Shock-dominated** (hypersonic, detonation) | HLLE, HLLC+Tramel, AUSM+M | Robust; AUSM+M/HLLC+Chen are carbuncle-stable |
| **All-speed (low Mach + shocks)** | **HLLC+Chen**, **AUSM+M**, SLAU2 | single solver across the whole Mach range |
| **Smooth subsonic** (combustor, inlet) | SLAU2, AUSM+M | Low-Mach accuracy; SLAU2 best |
| **Mixed subsonic/transonic** | HLLC+Tramel, AUSM+M | HLLC+Tramel adaptive; AUSM+M all-speed |
| **Boundary layers + shocks** | HLLC+Tramel, HLLC+Chen | shock-aware; HLLC+Chen preserves shear |
| **Very low Mach, time-accurate** | **LMRoe**, SLAU2, HLLC+Chen | flux-only low-Mach fix; time-accurate |
| **Very low Mach, steady** | HLLC-PC, MiczekRoe, LMRoe | preconditioned / low-Mach; clean convergence |
| **Vortex / acoustic-near-field preservation** | LMRoe, HLLC+Chen | removes the $\mathcal{O}(1/M)$ momentum dissipation |
| **Unsteady shock interaction** | HLLC+Tramel | Shock detection tracks transients |
| **Emergency (solver divergence)** | HLLE, LLF | Maximum stability; accept extra diffusion |

### Convergence Tips

1. **Start robust, refine if needed:** Begin with HLLE or SLAU2, switch to HLLC+Tramel after convergence behaves.
2. **Monitor divergence:** If solver crashes, switch to HLLE (fallback).
3. **Verify solution structure:** Check shock positions, contact discontinuities, and Mach numbers.
4. **Mesh refinement:** Coarse meshes tolerate more dissipative solvers (HLLE); fine meshes benefit from HLLC+Tramel.

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
8. S.-s. Chen, B. Lin, Y. Li, C. Yan, "HLLC+: Low-Mach shock-stable HLLC-type
   Riemann solver for all-speed flows," *SIAM J. Sci. Comput.*, 42(4), 2020,
   B921–B950.
9. S.-s. Chen, F.-j. Cai, H.-c. Xue, N. Wang, C. Yan, "An improved AUSM-family
   scheme with robustness and accuracy for all Mach number flows,"
   *Appl. Math. Modelling*, 77, 2020, 1065–1081.
10. F. Rieper, "A low-Mach number fix for Roe's approximate Riemann solver,"
    *J. Comput. Phys.*, 230(13), 2011, 5263–5287.
11. F. Miczek, F. K. Röpke, P. V. F. Edelmann, "New numerical solver for flows
    at various Mach numbers," *Astron. Astrophys.*, 576, 2015, A50.
