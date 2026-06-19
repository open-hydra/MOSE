# Boundary Conditions for 1D Low-Mach Premixed Flame — Diagnosis Notes

## Case setup

- 1D laminar premixed flame, low Mach number.
- Solver: finite volume, structured grid, explicit time integration, Godunov-like flux with LMRoe.
- LMRoe rescales the Roe dissipation matrix on the entropy/shear fields to cure low-Mach inaccuracy, but **does not** change the acoustic eigenvalues (u±c) or relax the acoustic CFL restriction. The characteristic decomposition used by NSCBC-style boundary treatments is therefore unaffected by the LMRoe modification and can be applied as-is.
- Workflow: pre-compute S_L and unburnt-side thermochemical state with Cantera, then build inflow/outflow BCs around that target.
- Solver options available: time-accurate or local time-stepping. **Use time-accurate mode for all BC debugging** — local time-stepping decouples the physical phase relationships needed to diagnose acoustic effects.

## Symptom

Strong pressure and velocity oscillations; velocity oscillates between positive and negative values (i.e. flips sign). With a fixed back-pressure outflow, oscillations ring persistently. With pure zero-gradient (full extrapolation) outflow, the oscillation goes away but the mean pressure level drifts.

**Diagnostic check (not yet confirmed in this conversation, do this next):** compare the oscillation period to the acoustic round-trip time 2L/c of the domain. A match confirms a trapped acoustic mode rather than a chemistry/transport issue.

## Diagnosis trail

### Attempt 1 — fixed ρ (from p_ref) + u, hard inflow; hard back-pressure outflow
ρ at the inflow was computed once from a fixed reference pressure and target T, then held constant; u was imposed directly. Because ρ, u are both frozen, and T is also fixed, this pins p = ρRT = p_ref every step as well — i.e. **three quantities (ρ, u, p) are simultaneously over-constrained**, even though only two should be externally imposed for a subsonic inflow (one characteristic is genuinely outgoing and must be left free / extrapolated).

Combined with a hard back-pressure outflow (also fully reflecting on its one incoming characteristic), this builds a **closed resonant acoustic cavity**. Since LMRoe leaves acoustic-eigenvalue dissipation essentially at Roe level, the trapped mode rings with little damping. Any mismatch between the imposed S_L and the actual discrete equilibrium flame speed acts as a persistent low-level acoustic forcing — this is the likely source of the velocity sign flips (acoustic perturbation amplitude comparable to or larger than S_L itself, which is only cm/s).

### Attempt 2 — zero-gradient outflow
Removes the wall on the outflow side (waves can leave), which kills the oscillation, but removes any anchor on the mean pressure level. Any small imbalance between imposed inflow mass flux and what the flame actually consumes integrates in time as pressure drift, with nothing to counteract it.

### Attempt 3 — ρ computed from local/extrapolated (domain) pressure instead of fixed p_ref
This fixes the *literal* triple-constraint problem — p is no longer pinned to a fixed external value. However, a **second, more subtle issue remains** in how ρ and u are then assembled:

```
ρ' = p'/(R T_target) = (γ/c0²) p'        [since c0² = γ R T_target]
u' = −ṅ_target · ρ'/ρ0² = −u0 (ρ'/ρ0)
```

Combining these gives the implicit relation the boundary forces between velocity and pressure perturbations:

```
u' = −(γ Ma0 / (ρ0 c0)) · p'
```

Compare to the relation a genuinely non-reflecting boundary needs (the acoustic impedance relation):

```
u' = −p'/(ρ0 c0)
```

The recipe's slope is off from the correct one by a factor of **γ·Ma0** instead of 1. Since Ma0 ≪ 1 for a laminar flame, this boundary **nearly freezes u′ regardless of p′ — i.e. it behaves like a rigid wall on velocity, and gets *more* wall-like as Ma0 → 0.** This explains why fixing the over-constraint on paper did not fix the oscillations: deriving u from ṅ_target/ρ instantaneously, every step, still hard-clamps the acoustic channel — just through a less obvious algebraic path.

### Attempt 4 — switch to (u_target, T_target) instead of (ṅ_target, T_target)
This matters and is the *correct* variable choice, for a structural reason rooted in the LODI relations themselves. With wave amplitudes L1 (assoc. with u−c), L2 (entropy wave, assoc. with u), L5 (assoc. with u+c):

```
∂u/∂t = −(L5 − L1)/(2ρc)
∂p/∂t = −(L5 + L1)/2
∂ρ/∂t = −[L2 + (L5 + L1)/2]/c²
```

u and p depend only on the acoustic pair (L1, L5); L2 never enters their evolution. Density (and hence mass flux ṅ = ρu) is entangled with **both** the entropy wave and the acoustic pair. Targeting ṅ directly (with T also targeted separately) silently re-couples the entropy and acoustic channels — the same kind of cross-contamination this whole exercise is trying to remove, just hidden one layer deeper. Targeting u directly keeps the channels cleanly separated, which is why Yoo & Im (2007) and Poinsot & Lele (1992) formulate the inflow this way.

Practical implication: use `u_target = S_L` directly (already known from Cantera) as the **relaxed** acoustic-channel target, with T_target handling the entropy channel. Mass flow rate stops being imposed and becomes a diagnostic that fluctuates with genuine acoustic activity, settling on time-average to ρ_u·S_L once the domain reaches equilibrium — note ρ_u here is whatever the domain self-consistently produces given the floated pressure and T_target, not necessarily the exact Cantera value at p_ref (mean pressure is set by the outflow target, not the inflow).

### Attempt 5 — current implementation: no ghost cells, direct flux at interface
State at the inflow face: p_face from interior (extrapolated), ρ_face from EOS using bc_T, u_face = bc_u (hard), ṅ and other fluxes computed from that state directly into the flux formula (no Riemann solve at that face).

**This is still fully reflecting**, regardless of how well-behaved p_face and ρ_face are, for a reason independent of all the above: from u′ = (w₊ − w₋)/(2ρc), forcing u_face = bc_u exactly every step forces u′ = 0 identically, which algebraically requires w₊ = w₋ at every instant regardless of what w₋ (the wave actually arriving from the interior) is. **This is the textbook "fully reflecting inflow" case from Poinsot & Lele (1992)** — hard-imposed velocity is their canonical anti-example, behaving like a rigid acoustic wall (p reflects with +1 coefficient, u pinned to zero perturbation — same signature as a closed pipe end).

Using direct flux assembly instead of ghost cells + Riemann solver is **not** itself the problem — it's a legitimate, common implementation choice. The problem is purely the hard-vs-relaxed axis: whether u_face is reset from bc_u every step (reflecting) or evolved/relaxed toward bc_u over time (non-reflecting).

## Recommended fix

Do not reset the face state from algebra every step. Instead, give the face state memory and relax it toward targets on a timescale slow compared to the acoustic frequency — the same logic already applied (correctly, in spirit) to the outflow back-pressure idea.

**Wave-amplitude (rigorous) version**, at the inflow:
1. Compute the outgoing wave amplitude (associated with u−c) from a one-sided spatial derivative in the interior — already available from the existing stencil.
2. Replace the incoming wave amplitude (associated with u+c) not with the value that exactly cancels the outgoing one (the reflecting choice), but with a relaxed target:
   ```
   w₊ = w₋ − 2ρc · σ_in · (u_face − u_target)
   ```
3. Integrate p_face, u_face forward in time via the LODI relations rather than recomputing them from scratch each step.
4. Keep T_target (entropy channel) and Y_target (species channels) as hard targets — these are correctly externally-imposed incoming characteristics and don't suffer this problem.
5. Feed the evolved face state into the existing direct-flux assembly unchanged.

**Cheap stand-in**, before committing to full LODI bookkeeping:
```
du_face/dt = σ_in (u_target − u_face)
```
integrated alongside the domain; ρ_face still from T_target via EOS, p_face still extrapolated from interior. Not equivalent to the rigorous wave-amplitude version (σ isn't scaled correctly against the acoustic timescale), but removes the u′≡0 pathology and should noticeably reduce the oscillation as a first check.

**At the outflow**, apply the analogous relaxed pressure control (Rudy & Strikwerda, 1980) rather than a hard back-pressure clamp or pure zero-gradient:
```
L1 = σ_out (p − p_target)
```
with σ_out ≈ K(1 − Ma²)c / L_domain, K ~ 0.2 as a starting point (tune from there). This is the middle ground between the two outflow experiments already run: it controls the long-term mean pressure (fixing the drift seen with zero-gradient) while letting fast transients pass through largely unreflected (fixing the resonance seen with hard back-pressure).

Suggested starting point for σ_in: similar order of magnitude to σ_out, tune empirically; must stay well below the acoustic frequency of the domain or it degenerates back toward the hard/reflecting limit.

## Caveat on the direct-flux (no-Riemann-solve) approach

Evaluating the flux on a single, jump-free boundary state means **zero numerical dissipation is applied at that face** — none of LMRoe's tailored dissipation is present there, since there's no Riemann problem being solved. Usually harmless once the boundary state itself is well-behaved (i.e. after the relaxation fix above), but if high-frequency, grid-scale noise appears at that face specifically (as opposed to a smooth low-frequency oscillation) after implementing relaxation, this is the likely cause — a touch of artificial dissipation local to that face (not the whole domain) is the standard patch.

## Key references

- Poinsot, T., Lele, S., "Boundary conditions for direct simulations of compressible viscous flows," *J. Comput. Phys.*, 101(1), 1992, pp. 104–129. — Original NSCBC formulation; canonical reference for why hard-imposed velocity/pressure inflow/outflow is fully reflecting, and for the LODI relations.
- Rudy, D., Strikwerda, J., *J. Comput. Phys.*, 36, 1980, pp. 55–70. — Relaxed/partially-reflecting outflow pressure control, σ formulation.
- Yoo, C.S., Wang, Y., Trouvé, A., Im, H.G., "Characteristic boundary conditions for direct simulations of turbulent counterflow flames," *Combustion and Flame*, 2005. — Reacting-flow NSCBC inflow with imposed velocity + temperature.
- Yoo, C.S., Im, H.G., "Characteristic boundary conditions for simulations of compressible reacting flows with multi-dimensional, viscous and reaction effects," *Combustion Theory and Modelling*, 11(2), 2007, pp. 259–286. — Improved/generalized NSCBC for reacting flows; the (u, T) inflow pairing referenced in this conversation comes from here.
- Motheau, E., Almgren, A., Bell, J.B., "Navier–Stokes Characteristic Boundary Conditions Using Ghost Cells," *AIAA Journal*, 2017 (GC-NSCBC). — NSCBC reformulated specifically for finite-volume, ghost-cell, Godunov-type solvers; structurally the closest match to this solver's architecture. Built on the Yoo & Im reacting-flow wave-amplitude treatment.
- PeleC documentation (AMReX-Combustion, open source) — implements GC-NSCBC with relaxation parameters; useful as a working code reference for the ghost-cell update and relaxation coefficients in practice.

## RESOLUTION (2026-06-17) — diagnosis confirmed, oscillation eliminated

### Step 1: 2L/c diagnostic — CONFIRMED trapped acoustic mode
- Observed oscillation (probe pressure, cell 80, FFT): a single sharp tone at **3974 Hz (period 0.252 ms)**, with violent amplitude (p std ≈ 1686 Pa; extremes ~31→170 kPa about a 101.2 kPa mean; u flipping sign).
- Acoustic travel time integrated over the *real* field (T 289→2235 K, c 349→958 m/s, flame at x≈0.107 m): one-way τ = ∫dx/c ≈ 0.50 ms → fundamental c/(2L)-equivalent = 1/(2τ) ≈ **991 Hz**.
- **f_obs / f_fundamental = 4.01** → the oscillation is the n=4 acoustic harmonic of the domain (4·991 = 3966 Hz, match to 0.2%). It sits exactly on the integer harmonic ladder, ~3–4 orders of magnitude faster than any convective (94–880 ms) or chemical timescale. Pure acoustics, not chemistry/transport.

### Step 2–3: root cause in code (two distinct issues) + fix
- **Outflow bug** — `Lib_BC_Fluxes_Outflow.f90`: the relaxed `p_exit = rel·p_target + (1−rel)·p_bound` was computed but then **discarded**; the call passed the hard `BC_pexit` instead. So the `rel_fac` knob in `bc.txt` was inert and the outflow was always a hard back-pressure (reflecting). Fixed by passing `p_exit` into `Outflow`.
- **Inflow** — `Lib_BC_Fluxes_Inflow.f90` BC 404: `Un3 = BC_g` hard-imposed the velocity every step (`u′≡0`, the textbook fully-reflecting inflow); `BC_rel_fac` was passed in but never used. Replaced with `Un3 = rel·BC_g + (1−rel)·Un` (= explicit-Euler discretization of the relaxation ODE `du/dt = σ(u_target−u)` with α=rel·… via the boundary-cell normal velocity Un). T and species kept as hard targets (entropy/species channels untouched); mass flux ρ·Un3 is now a free diagnostic.
- **bc.txt**: inflow `rel_fac` 1.0→0.1, outflow `rel_fac` 0.0→0.1. Semantics now symmetric: rel=1 hard/reflecting, rel→0 transparent (but drifts). 0.1 = soft anchor, mostly transparent.

### Step 4: validation (time-accurate, resumed from the oscillating field, 30k iters)
Pressure decays monotonically to a flat steady state over ~49 ms:
| quarter | p std (Pa) |
|---|---|
| 1 | 111.9 |
| 2 | 14.1 |
| 3 | 1.6 |
| 4 | **0.67** |
- p std 1686 → 0.67 Pa (**~2500× reduction**); 3974 Hz tone gone; mean pressure pinned at **101325 Pa = outflow target (no drift)**; flame T steady at 2224 K (std 8 K).
- Both original failure modes resolved at once: resonance (was the hard-u + hard-p cavity) and mean-pressure drift (was zero-gradient). No further tuning needed.

## Remaining (optional) follow-ups
4. Switch back to local time-stepping now that the BC treatment is validated, for faster convergence.
5. Flame-position tracking: since S_L isn't known exactly at the discrete level, consider feedback control on u_target via an iso-T / iso-Y_fuel contour, so the flame stays anchored while u_target converges toward the true discrete S_L. Only needed if the flame drifts off its target location over long runs.
6. If high-frequency grid-scale boundary noise ever appears at the inflow face specifically, add local artificial dissipation at that face (see caveat above) — not observed in this run.
6. If switching to direct-flux boundary evaluation introduces high-frequency boundary noise after the relaxation fix, add local dissipation at that face specifically (see caveat above) rather than globally.