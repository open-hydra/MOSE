# MOSE Performance Optimization Notes

> Generated: 7 April 2026  
> Target file: `src/lib/numerics/fluxes/Mod_Fluxes.f90` and related flux routines

---

## Priority 0 — Bug: stray `stop` statement

There is a `stop` on line ~125 of `Mod_Fluxes.f90`, after the first-direction diffusive flux loop and before the direction-2 convective loop. This means **directions 2 and 3 (j- and k-fluxes) never execute**. If this is a debug leftover, remove it.

---

## Priority 1 — Replace `rlimiter` procedure pointer with static dispatch

`rlimiter` (in `MOSE_Lib_Limiters`) is a procedure pointer called inside the `elemental` subroutine `Reconstruction`, which is called **once per face per RK stage**. Since the limiter choice does not change during a run, the pointer prevents:

- **Inlining** — the limiter is a one-liner (e.g., Van Leer: 3 arithmetic ops), but the indirect call overhead dominates.
- **SIMD vectorization** — the compiler cannot vectorize through an indirect call.

### Suggested fix

Option A — Use a module-level integer `limiter_id` and a `select case` inside `Reconstruction`:

```fortran
select case (limiter_id)
  case (LIM_MINMOD)  ; slopel = rlimiter_MINMOD(slope1, slope0)
  case (LIM_VANLEER) ; slopel = rlimiter_VANLEER(slope1, slope0)
  ...
end select
```

Option B — Duplicate `Reconstruction` for each limiter using preprocessor macros or separate routines, and dispatch at the loop level in `Fluxes_blk`.

**Expected gain:** 10–20% on the full solver (limiter is in the innermost loop).

---

## Priority 2 — Eliminate array temporaries in `Convective_Flux` calls

### Already fixed: `dl` argument

The `dl(i-1:i+2,j,k) % c(1)` expression was replaced with an explicit `dl_loc` array.

### Still present: `Prim` and `Res` arguments

Each call passes:

| Argument | Example | Issue |
|---|---|---|
| `Prim(:,i-1:i+2,j,k)` | 4-column slice, dim-2 stride | Non-contiguous → runtime copy |
| `Res(:,i:i+1,j,k)` | 2-column slice, dim-2 stride | Non-contiguous → runtime copy back |
| `Prim(:,i,j-1:j+2,k)` | Dir-2: stride across dim-3 | Worse: stride across two dims |

### Suggested fix

Refactor `Convective_Flux` to receive the full `Prim` array plus index coordinates:

```fortran
subroutine Convective_Flux ( dl, normal, area, Prim, Res, beta, SD_limiter, &
                             i, j, k, dir )
  real(R8), intent(in)    :: dl(-1:2), normal(3), area, beta
  real(R8), intent(in)    :: Prim(nprim, 1-gc:*, 1-gc:*, 1-gc:*)
  real(R8), intent(inout) :: Res(nprim, 1-gc:*, 1-gc:*, 1-gc:*)
  integer, intent(in)     :: i, j, k, dir
  logical, intent(in)     :: SD_limiter
```

Then index into `Prim` and `Res` using offsets computed from `dir`. This eliminates all temporary copies. Alternatively, use local variables:

```fortran
real(R8) :: P_stencil(nprim, -1:2), R_face(nprim, 0:1)
! Manually copy the stencil before the call
P_stencil(:,-1) = Prim(:,i-1,j,k)
P_stencil(:, 0) = Prim(:,i  ,j,k)
P_stencil(:, 1) = Prim(:,i+1,j,k)
P_stencil(:, 2) = Prim(:,i+2,j,k)
```

**Expected gain:** 10–15%.

---

## Priority 3 — Precompute thermodynamic quantities per cell

`co_rotot_Rtot` (density + gas constant from species mass fractions) is called redundantly:

- 2× per face in `Convective_Flux` (left + right reconstructed states)
- 6× per face in `Diffusive_Flux` (2 in the main body + 4 in `Tangential_Gradient`)
- 1× more in `Compute_Diffusive_Flux`

Many of these calls operate on the **same cell-center state** but from different faces. For a cell with 6 neighbors, that's up to ~18 redundant evaluations.

### Suggested fix

Add per-cell precomputation before the face loops:

```fortran
real(R8) :: rho_cell(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
real(R8) :: T_cell  (1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
real(R8) :: Rgas_cell(1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)

!$omp do collapse(3)
do k = 1-gc, n(3)+gc
do j = 1-gc, n(2)+gc
do i = 1-gc, n(1)+gc
  call co_rotot_Rtot(Prim(1:nsc,i,j,k), rho_cell(i,j,k), Rgas_cell(i,j,k))
  T_cell(i,j,k) = Prim(np,i,j,k) / (rho_cell(i,j,k) * Rgas_cell(i,j,k))
end do; end do; end do
```

Then pass `rho_cell`, `T_cell` etc. to the flux routines instead of recomputing.

**Expected gain:** 5–15% depending on the cost of `co_rotot_Rtot`.

---

## Priority 4 — Fuse convective + diffusive loops

Currently, convective and diffusive fluxes for each direction are computed in separate triple loops. Each loop pass streams through the entire `Prim` and `Res` arrays, causing redundant cache misses.

### Suggested fix

Merge into a single loop per direction:

```fortran
!$omp do collapse (2)
do k = 1, n(3)
do j = 1, n(2)
do i = 1, n(1) - 1
  call Convective_Flux(...)   ! data already in L1/L2 cache
  if (model > 0) call Diffusive_Flux(...)  ! reuse the same data
enddo; enddo; enddo
```

**Expected gain:** 5–10% (cache locality).

---

## Priority 5 — Replace `Riemann` procedure pointer with static dispatch

Same issue as `rlimiter`. The `Riemann` solver is a procedure pointer called per face. While each call is heavier than `rlimiter` (so the relative overhead is smaller), the indirect call still blocks inlining and compiler optimization.

### Suggested fix

Dispatch at the loop level:

```fortran
select case (riemann_id)
  case (RIEMANN_HLLC)
    do k = 1, n(3); do j = 1, n(2); do i = 1, n(1)-1
      call riemann_HLLC(...)  ! compiler can inline
    enddo; enddo; enddo
  case (RIEMANN_AUSM)
    ...
end select
```

Or use a generic `Riemann` interface with a discriminating integer argument.

**Expected gain:** 5–10%.

---

## Priority 6 — AoS → SoA for metric types

The derived types used for mesh metrics store small arrays inside structs:

```fortran
type :: MOSE_vector_3D_type
  real(R8) :: c(3)         ! dl(i,j,k)%c(d) — Array of Structures
end type
```

This layout causes stride-3 access when sweeping over cells for a single component. The optimal layout for SIMD is Structure of Arrays:

```fortran
real(R8) :: dl(3, 1-gc:n1+gc, 1-gc:n2+gc, 1-gc:n3+gc)  ! dl(d,i,j,k)
```

Same applies to `MOSE_tensor_3D_type` (`c(3,3)`) and `MOSE_f_metrics_type` (`N(3)`, `A`).

### Impact

This is a large refactor affecting the mesh setup, I/O, boundary conditions, etc. Best done incrementally, starting with `dl` in the flux routines (which already required the `dl_loc` workaround).

**Expected gain:** 5–10% (better vectorization).

---

## Priority 7 — Compiler flags

### Currently used (Release)

`-O3`, `-funroll-loops`, `-ffast-math`, `-march=native`, `-finline-functions`

### Additions to consider

| Flag (gfortran) | Effect |
|---|---|
| `-flto` | Link-time optimization — enables cross-module inlining. Your cmake says it's "fragile" but modern GCC (≥10) handles it reliably. |
| `-fopt-info-vec-missed` | Diagnostic: shows what the compiler fails to vectorize — useful for identifying the remaining bottlenecks. |
| `-fopenmp-simd` | Enables `!$omp simd` directives for explicit vectorization hints without full OpenMP threading overhead. |

**Expected gain:** 5–10% from `-flto` alone (cross-module inlining of small routines like limiters).

---

## Summary

| # | Action | Effort | Expected gain |
|---|---|---|---|
| 0 | Remove stray `stop` | trivial | correctness |
| 1 | Static dispatch for `rlimiter` | small | 10–20% |
| 2 | Eliminate `Prim`/`Res` array temporaries | medium | 10–15% |
| 3 | Precompute `rho`, `T`, `a` per cell | medium | 5–15% |
| 4 | Fuse convective + diffusive loops | small | 5–10% |
| 5 | Static dispatch for `Riemann` | small | 5–10% |
| 6 | AoS → SoA for metric types | large | 5–10% |
| 7 | Enable `-flto` | trivial | 5–10% |
