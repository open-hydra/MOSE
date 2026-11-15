# Multirate Analysis

This page proves the order, stability and monotonicity statements of the [multirate integrator](multirate.md). Appendix A treats the macro step on a fixed mesh and the margin on condition (4); Appendix B the total variation at a face between levels.

## Appendix A: Order and stability of the multirate macro step

This appendix proves the two properties of the macro step that the section
[Recursive Heun scheme](multirate-method.md#recursive-heun-scheme) states: the second order on a
fixed mesh (Theorem 1) and the monotonicity for linear two-cell fluxes
(Theorem 2), with their hypotheses and their range, and what they give for a
map that changes from one macro step to the next (Corollary).

**Setting.**  The semi-discretisation is $y' = f(y) = \sum_k f_k(y)$, with $y$
the cell averages of the conservative variables (and of the turbulence
variables) and $k$ running over the faces (convective and diffusive fluxes),
the boundary records and the turbulence sources; $f_k$ depends only on the
cells of its stencil $S_k$.  The level map is nested ($|\ell_c - \ell_v| \le 1$
across every face, equal levels across connections), $h_\ell = H/2^\ell$, and
every stage is $U := U^0 + c_k\,(U - U^0 + h_\ell L(U))$ with $c = (1, \tfrac12)$:
Heun in the form of Shu and Osher [1], $U^{(1)} = U^0 + h\,L(Y^{(1)})$ and
$U^1 = U^0 + \tfrac h2\,[L(Y^{(1)}) + L(Y^{(2)})]$, where $Y^{(1)}$ and $Y^{(2)}$
are the stage vectors, the states that the whole stencil reads.

On a face
between a coarse level $m$ and a fine level $m+1$ the coarse side uses and
saves $G_s$ at its stage 1, the four fine stages accumulate
$A = \tfrac14\sum_{j=1}^4 G_j$, and the coarse side uses $2A - G_s$ at its
stage 2 without evaluating the face
([Conservative correction](multirate-method.md#conservative-correction)).

Two facts carry both proofs.

**(F1) Weights and first moments.**  Every stage is linear in the evaluations
of $L$, and the correction replaces a flux by a linear combination of fluxes,
so a macro step updates a cell $c$ as
$U^{n+1}_c = U^n_c + \sum_k \sum_j w_{kj}\, f_{k,c}(Y_{kj})$, with weights
$w_{kj} > 0$ and stage times $t_{kj} \in [0, H]$ counted from $t_n$.  For every
pair of a cell and a face, $\sum_j w_{kj} = H$ and $\sum_j w_{kj}\, t_{kj} = H^2/2$:

- a face without correction, seen by a cell of level $\ell$: the $2^\ell$ Heun
  sub-steps put the weight $h_\ell/2$ at the times $m h_\ell$ and $(m+1) h_\ell$,
  $m = 0, \dots, 2^\ell - 1$, so $\sum w = 2^\ell h_\ell = H$ and
  $\sum w t = \sum_m \tfrac{h_\ell^2}{2}(2m+1) = \tfrac{h_\ell^2}{2}\,4^\ell = \tfrac{H^2}{2}$;
- the coarse side of a corrected face: on each of its sub-steps
  $[\tau, \tau + h]$ it receives $\tfrac h2\,(G_s + 2A - G_s) = hA$, the weight
  $h/4$ at the four fine stage times $\tau$, $\tau + h/2$, $\tau + h/2$ and
  $\tau + h$: $\sum w = h$ and $\sum w t = h\tau + h^2/2$, the moments of a face
  without correction on the same sub-step;
- the fine side of a corrected face: its own four stages, as a face without
  correction.

**(F2) Consistency of the stage vectors.**  Every evaluation reads, on its whole
stencil, $Y_{kj} = y^n + t_{kj}\, f(y^n) + O(H^2)$, with a constant that depends
on the mesh and the map but not on the stage time.  By induction over the
recursion: the cells of the level are read at the start of their sub-step or
at the predictor, $y(\tau) + O(H^2)$ and $y(\tau + h_\ell) + O(H^2)$ (each
sub-step is first-order consistent, and at most $2^L$ of them precede $\tau$);
the finer cells are at the stage time, because the finer levels complete their
sub-steps before the coarser corrector; a coarser cell of level $m$ is read
through the interpolant,
$\hat U_m(t) = U_m(\tau_m) + (t - \tau_m)\, f_m(Y(\tau_m)) = y_m(t) + O(h_m^2)$ for
$t \in [\tau_m, \tau_m + h_m]$, also when it is two or more levels coarser; the
ghost cells are set at the stage time, from boundary data of that time and
from cells for which (F2) holds.

### Theorem 1 (order on a fixed mesh)

*Hypotheses.*  (H1) The mesh and the level map are fixed.  (H2) $f \in C^2$ in a
neighbourhood of the trajectory (no hypothesis for a linear $f$).  (H3) The
recursion, the interpolant and the correction are those above: steps in the
ratio 2 between levels, nesting, equal levels across connections, bands of one
cell allowed, the coarse corrector reading the final values of the finer
cells.  (H4) Time-dependent boundary data are evaluated at the stage times.

*Statement.*  The macro step $\Psi_H$ satisfies $\Psi_H(y) - \varphi_H(y) = O(H^3)$,
with $\varphi_H$ the exact flow of $y' = f(y)$ over $H$ and a constant that
depends on $f$, the mesh and the map, not on where a cell lies with respect to
the faces between levels.  The error of the macro steps up to a time $T$ is
$O(H^2)$, uniformly over the cells, the cells next to the faces between levels
included.

*Proof.*  By (F2), $f_k(Y_{kj}) = f_k(y) + t_{kj}\, f_k'(y)\, f(y) + O(H^2)$, where
the Jacobian $f_k'$ involves only the components of $S_k$, on which (F2)
holds.  With (F1),

$$
\Psi_H(y) - y = \sum_k \Bigl(\sum_j w_{kj}\Bigr) f_k(y)
+ \sum_k \Bigl(\sum_j w_{kj}\, t_{kj}\Bigr) f_k'(y)\, f(y) + O(H^3)
= H f + \frac{H^2}{2}\, f' f + O(H^3),
$$

the expansion of $\varphi_H(y) - y$; for a linear $f = Ly$ the same count gives
$\Psi_H = I + HL + \tfrac{H^2}{2}L^2 + O(H^3)$.  The macro step is a composition
of $2^{L+2} - 2$ maps $Y \mapsto Y + h\,L(Y)$ and of linear combinations with
bounded coefficients, so it is Lipschitz with a constant $1 + CH$, $C$
proportional to $\|f'\|$.  Then $|e_{n+1}| \le (1 + CH)\,|e_n| + DH^3$ and
$|e_n| \le DH^2\,(e^{CT} - 1)/C$. ∎

*Remark 1 (the two differences from [10]).*  (i) The coarse corrector reads the
finer cells at their final values, $y_F(\tau + h) + O(h^3)$, where the scheme of
[10] reads their predictor, $y_F(\tau + h) + O(h^2)$: both satisfy (F2), and the
choice enters only the $O(H^3)$ term.  With two-cell fluxes the two choices
give the same scheme, since the face between levels is not evaluated again at
the coarse stage 2.  (ii) With several levels and bands of one cell, (F1) and
(F2) hold level by level, and a halo cell two or more levels coarser changes
only the constant, through $O(h_m^2)$.  Neither difference adds an order
condition.

*Remark 2 (range).*  The constant of Theorem 1 contains $\|f'\| \sim 1/\Delta x$:
the theorem holds for $H \to 0$ on a fixed mesh, which is what the order tests
of [Multirate Verification](../vv/multirate.md) measure.  Where the limiter switches, and at the
switches of HLLC and of the turbulence models, (H2) fails locally, and the
local error of Heun itself is $O(H^2)$ there.  At a fixed Courant number
($H \propto \Delta x$) the theorem gives no bound (the Corollary below gives
one for linear two-cell fluxes).

At the coarse stage 2 the
coarse cells are at their forward-Euler predictor, with the error
$-\tfrac{h^2}{2}L^2 y$, and the finer cells at their final values, with an
error $O(h^3)$.  A flux that reads both sides of a face between levels (the
MUSCL slope of the coarse cell next to the face, the HLLC flux, the diffusive
flux) sees a jump $O(h^2)$ in the errors of its stencil, and its error,
$O(h^2)$, differs by as much from that of the next face: the error of the
state is $h \cdot O(h^2)/\Delta x = O(h^2)$ per macro step in the two cells of
the face, with opposite signs.

The local error is then $O(H^2)$ next to a face
between levels and $O(H^3)$ elsewhere: first order locally, as the
conventional boundary data of [7].

### Theorem 2 (monotonicity for linear two-cell fluxes)

*Hypotheses.*  $dU_i/dt = \sum_j \beta_{ij}\,(U_j - U_i)$ with $\beta_{ij} \ge 0$, the
flux of every face being a linear function of its two cells alone: the
first-order upwind flux ($\beta_{ij} = |a_f| A_f/V_i$ for the upwind neighbour,
with $a_f$ the normal velocity through the face of area $A_f$, and 0 for the
downwind one) and the central diffusion ($\beta_{ij} = D A_f/(\delta_f V_i)$, with
$\delta_f$ the distance between the two cell centres), on any mesh and in any
dimension (the diffusive flux of MOSE is this two-cell flux for a scalar on orthogonal cells; on skewed cells it also
reads the neighbours along the face).

With the upwind flux the velocity field has no discrete divergence
($\sum_f a_f A_f = 0$ over the faces of every cell, $a_f$ outward); this form
is then the conservative scheme of the face fluxes $|a_f| A_f U_j$, also in the
macro step, where every cell sends its outflow at its own stage values (on a
corrected face too: $2A - G_s = |a_f| A_f\,U_i^{(1)}$ when the coarse cell $i$
is upwind), and the conservation in the statement uses it.

The load of cell
$i$ with its step $h_i = H/2^{\ell_i}$ is
$x_i = h_i \sum_j \beta_{ij}$; the forward-Euler step
$(1 - x_i)\,U_i + h_i\sum_j \beta_{ij}U_j$ is positive for $x_i \le 1$, and on a
Cartesian mesh $x_i = h_i\sum_d\bigl(|a_d|/\Delta x_d + 2D/\Delta x_d^2\bigr)$,
the ratio of $h_i$ to the local step with CFL = 1 and VNN = 1/2
([Local time step](time-integration.md#local-time-step)).  Any nested map, any number of levels,
and:

- (C1) $x_i \le 1$ for every cell;
- (C2) $x_i \le 2\,c(x_j)$, with $c(x) = 1 - x + x^2/2 - x^3/8$, on every face
  where a finer cell $j$ feeds a coarser cell $i$ ($\beta_{ij} > 0$,
  $\ell_j = \ell_i + 1$).

*Statement.*  The matrix $M$ of the macro step has non-negative entries,
$M\mathbf 1 = \mathbf 1$ and $\mathbf 1^T V M = \mathbf 1^T V$, with $V$ the diagonal
of the cell volumes.  Hence $\|M\|_\infty = 1$ (maximum principle),
$\|M\|_{1,V} = 1$ (conservation and positivity) and $\|M\|_{2,V} \le 1$, and the
same holds for every product of macro steps, with different maps.

*Proof.*  By induction over the recursion, every stage vector, interpolant and
end-of-sub-step state is a combination of the values at $t_n$ with
non-negative coefficients that sum to 1; such combinations are preserved by
composition, and the only coefficient that can become negative arises in the
coarse corrector.  Take a cell $i$ on a sub-step $[\tau, \tau + h]$, $h = h_i$.
Stage 1 gives $U_i^{(1)} = (1 - x_i)\,U_i^0 + h\sum_j \beta_{ij}\,\hat U_j(\tau)$, a
convex combination by (C1), $\hat U_j$ being the neighbour at the stage time
(its state, or its interpolant, itself a convex combination).

On a corrected
face with a finer neighbour $j$ the four fine stages read the coarse cell at
$\eta = 0, \tfrac12, \tfrac12, 1$, whose mean is $(U_i^0 + U_i^{(1)})/2$, so
$2A - G_s = \beta_{ij}\,[\,2\bar U_j - U_j(\tau) - U_i^{(1)}\,]$, with $\bar U_j$ the
mean of the four fine stage vectors.  Stage 2 then gives

$$
U_i^1 = \Bigl(1 - \frac{x_i}{2}\Bigr)U_i^0 - \frac{x_i}{2}\,U_i^{(1)}
+ \frac h2 \sum_{j\ \text{not finer}} \beta_{ij}\,\bigl[\hat U_j(\tau) + \hat U_j(\tau + h)\bigr]
+ h \sum_{j\ \text{finer}} \beta_{ij}\,\bar U_j .
$$

With $U_i^{(1)}$ substituted, $U_i^0$ has the coefficient $1 - x_i + x_i^2/2 > 0$;
a neighbour that is not finer contributes
$\tfrac h2\,\beta_{ij}(1 - x_i)\,\hat U_j(\tau) + \tfrac h2\,\beta_{ij}\,\hat U_j(\tau + h)$,
non-negative by (C1); a finer neighbour contributes $h\beta_{ij}\,(c_j - x_i/2)$
on $U_j(\tau)$, with $c_j$ the coefficient of $U_j(\tau)$ in $\bar U_j$, plus
non-negative terms.  The four fine stage vectors carry $U_j(\tau)$ with the
coefficients 1, $1 - x_j$ (the predictor), at least $1 - x_j + x_j^2/2$ (the
state after the first fine sub-step: the diagonal of Heun,
$I + hL + h^2L^2/2$, plus non-negative products) and at least
$(1 - x_j)(1 - x_j + x_j^2/2)$ (the predictor of the second fine sub-step), so
$c_j \ge c(x_j)$, and (C2) makes the coefficient non-negative.

A cell of a
band of one cell is coarse on one side and fine on the other, and falls in the
cases above.  The rows sum to 1 because a constant is an exact solution (every
flux of a constant vanishes, and with it $A$, $G_s$ and the correction); the
columns weighted by the volumes sum to 1 because every face gives its two
sides the same impulse with opposite signs (in the conservative form of the
hypotheses), also through the correction ($h_m A$ on either side).

$M \ge 0$ with $M\mathbf 1 = \mathbf 1$ gives
$\|M\|_\infty = 1$; with $\mathbf 1^T V M = \mathbf 1^T V$ it gives
$\|M\|_{1,V} = 1$; and since the rows are convex,
$(MU)_i^2 \le \sum_j M_{ij}U_j^2$ (Jensen), so
$\sum_i V_i\,(MU)_i^2 \le \sum_j V_j U_j^2$. ∎

*Remark 3 (threshold).*  $c$ decreases on $[0, 1]$, with $c(1) = 3/8$, and
$2c(x) \ge x$ exactly when $x^3/4 - x^2 + 3x - 2 \le 0$, that is for
$x \le x^* = 0.8603$: $x_i \le 0.86$ for every cell implies (C1) and (C2).  The
threshold is exact for the pure upwind flux with the finer cell upstream.
With the additive local-step rule and $h_c \le \Delta t_c/\sigma$, writing
$r_c = \sum_d |a_d|/\Delta x_d$ and $r_v = \sum_d D/\Delta x_d^2$,

$$
x = h\,(r_c + 2r_v) \le \frac{r_c + 2r_v}{\sigma\,(r_c/\text{cfl} + r_v/\text{vnn})}
\le \frac{\max(\text{cfl},\, 2\,\text{vnn})}{\sigma},
$$

0.4 with $\sigma = 2$, cfl 0.8 and vnn 0.4.  A single-rate Heun step is
positive up to $x = 1$ (Hundsdorfer [8], §6.1, pp. 49–50): the faces between
levels cost at most 14 % of the positivity step.  On a mesh of equal rectangular cells with constant coefficients, with
cfl $\le 1$ and, in Navier–Stokes
runs, vnn $\le 1/2$, the bounds that the redo imposes when $L_{\max} > 0$
([Non-local CFL condition](multirate-method.md#non-local-cfl-condition)), and the default $\sigma = 2$, every cell starts the macro step with
$x \le 1/2$ (in Euler runs $r_v = 0$ and $x \le \text{cfl}/\sigma$), below $x^*$: (C1) and (C2) hold for every map that
the level rule builds, whatever the limiter of the run (the loads are those of the first-order fluxes of Theorem 2), for
every $\sigma \ge 2$, with or without the interface margin; with a smaller
$\sigma$ they hold while $\max(\text{cfl}, 2\,\text{vnn})/\sigma \le x^*$.

Within the macro
step the certificate of (4) (with the redo) accepts, on the state at every start $\tau$ of a sub-step, a cell with $x \le 1$, the load of
(C1), and the finer
cell of a face between levels, at a start that it shares with the coarser one ($\hat h = 2h$), with $x \le 1/2$, where
$2c(1/2) = 1.22 > 1$, the load of (C2).  Theorem 2 takes the coefficients fixed within the macro step; when they change,
these bounds at the starts of the sub-steps are what the certificate checks, not a theorem.

The bound
$x \le \max(\text{cfl}, 2\,\text{vnn})/\sigma$ is written for a mesh of equal rectangular cells with constant
coefficients; on other rectangular cells the load of the diffusion can exceed $2r_v$
([Stability Conditions](time-integration.md#stability-conditions)), and on skewed cells the diffusion of MOSE is not a two-cell flux and
is outside Theorem 2.

*Remark 4 (what Theorem 2 does not cover).*  With a limited reconstruction the
coarse corrector is
$U_i^1 = \tfrac12 U_i^0 + \tfrac12\,[\,U_i^{(1)} + h\,\tilde L_i(Y^{(2)})\,]$, where
$\tilde L$ has, on a face between levels, the extrapolation $2A - G_s$ (twice
the mean of the fine fluxes minus the flux at $\tau$) in place of a flux: it
is not a convex combination of forward-Euler steps.  For a monotone flux
$G = aU$, $2A - G_s = a\,(2\bar U_j - U_j(\tau))$ leaves the range of the values
when the fine value decreases fast.

In the linear two-cell case this is the
term $-\tfrac{x_i}{2}\,h\beta_{ij}\,U_j(\tau)$ that $c_j$ absorbs above; with a
limiter the corresponding coefficient depends on the limiter at different
times, and the argument by convex combinations, by which Heun inherits the
properties of its forward-Euler stages, does not go through at the faces
between levels.  Outside Theorem 2 are also HLLC on the Euler and
Navier–Stokes equations, the turbulence sources, the terms of the
axisymmetric wedge (faces 200, which are not two-cell fluxes) and shocks that
cross faces between levels.

### Corollary (every sequence of maps)

*Hypotheses.*  The macro steps $n = 0, \dots, N-1$ take the steps $H_n \le H$,
with $\sum_n H_n = T$, and nested maps $M_n$ with levels in
$\{0, \dots, L\}$, $L$ fixed.  $M_n$ and $H_n$ may depend on the numerical
solution $y_n$: the map rebuilt at every macro step, the levels raised by a
redo, a step shortened to land on an output instant.  The error
$e_n = y_n - y(t_n)$ is taken against the exact solution of $y' = f(y)$ on the
same mesh, as in Theorem 1.

*Statement.*  (a) On a fixed mesh, under (H2)–(H4) and for $H$ small enough,
$|e_N| \le DH^2\,(e^{CT} - 1)/C$, with $D$ and $C$ the largest constants of
Theorem 1 over the nested maps.  (b) In the setting of Theorem 2, with every
$M_n$ within (C1) and (C2) for its step $H_n$ (as the maps of the level rule
are, on a mesh of equal rectangular cells with constant coefficients, with cfl $\le 1$, vnn $\le 1/2$ and
$\sigma \ge 2$, Remark 3), let the
second and third time derivatives of the exact
solution be bounded on $[0, T]$ by $K_2$ and $K_3$ independently of the mesh
(since $y^{(k)}(t) = e^{tL} L^k y(0)$ and $e^{tL}$ is non-negative with unit
row sums, bounds on $L^2 y(0)$ and $L^3 y(0)$ suffice), and let $|B|$ bound the
volume of the cells within $2^{L+2}$ cells of a face between levels of every
$M_n$. Then

$$
\|e_N\|_\infty \le c_L K_2\, T H, \qquad
\|e_N\|_{1,V} \le \Bigl(c_L K_2\,\frac{|B|}{H} + \frac{K_3\,|\Omega|}{6}\Bigr)\, T H^2,
$$

with $c_L$ that depends only on $L$ and $|\Omega|$ the volume of the domain.
At a fixed Courant number ($H \propto \Delta x$), with faces between levels of
bounded total area on a quasi-uniform family of meshes, $|B| = O(\Delta x)$:
the error is $O(H)$ in the maximum norm and $O(H^2)$ in the volume-weighted
$L^1$ norm, for every sequence of maps, a face between levels that travels
with a wave included.

*Proof.*  In both parts
$e_{n+1} = [\Psi_{H_n}^{M_n}(y_n) - \Psi_{H_n}^{M_n}(y(t_n))] + \tau_n$, with
$\tau_n = \Psi_{H_n}^{M_n}(y(t_n)) - y(t_{n+1})$ the local error of the map
$M_n$ on the exact solution: every bound below holds for every admissible map,
so the choice of $M_n$ from the numerical solution does not enter.  (a) On a
fixed mesh the nested maps are finitely many, so the largest constants exist,
$|e_{n+1}| \le (1 + CH_n)\,|e_n| + DH_n^3$, and
$|e_N| \le DH^2 \sum_n H_n\, e^{C(T - t_{n+1})} \le DH^2\,(e^{CT} - 1)/C$; for
$H$ small enough, by induction, the numerical solution stays in the
neighbourhood of (H2).

(b) The macro step is linear, $y_{n+1} = M_n y_n$, so
$e_N = \sum_n M_{N-1} \cdots M_{n+1}\, \tau_n$, and every product has norm 1 in
$\|\cdot\|_\infty$ and in $\|\cdot\|_{1,V}$ (Theorem 2).  Take every value that
the macro step computes (start state, predictor, interpolant, end of a
sub-step) with its deviation from the exact solution at its time.  A predictor
is a convex combination, by (C1), of values read at $\tau$, and adds the
remainder $y_i(\tau + h) - y_i(\tau) - h\,y_i'(\tau)$, at most $h^2 K_2/2$; an
interpolant adds at most $h_m^2 K_2/8$; a corrector adds at most twice the
largest deviation of the values it reads (the weights of each face sum to $h$,
and $h \sum_j \beta_{ij} \le 1$) and a quadrature error at most $h^2 K_2$ (the
weights of each face have the moments of (F1), and
$|\beta_{ij}\,(y_j'' - y_i'')| \le 2\beta_{ij} K_2$).

Over the at most
$2^{L+2} - 2$ stages of a macro step every deviation, $\tau_n$ included, stays
within $c_L K_2 H_n^2$.  The value of a cell after the macro step depends only
on the cells within $2^{L+2}$ cells of it (each stage reads face neighbours);
if they all have its level $\ell$, it is the value after $2^\ell$ single-rate
Heun steps of $h_\ell$, whose rows there are non-negative with unit sums
(by (C1), since $I + hL + \tfrac{h^2}{2}L^2 = \tfrac12 I + \tfrac12 (I + hL)^2$),
each with the local error $y(\tau + h) - (I + hL + \tfrac{h^2}{2}L^2)\,y(\tau)$,
at most $h^3 K_3/6$: $|\tau_{n,i}| \le 2^\ell h_\ell^3 K_3/6 \le K_3 H_n^3/6$. With $\sum_n H_n^2 \le HT$ and $\sum_n H_n^3 \le H^2 T$ the bounds follow. ∎

Part (b) holds in the setting of Theorem 2 only: outside it (Remark 4) the
bound of the products of macro steps is missing, and at a fixed Courant number
no bound is proven.  Part (a) holds on a fixed mesh, with constants that grow
with $1/\Delta x$ (Remark 2).  In the maximum norm the order 1 is what the
errors of order $H^2$ of the $T/H$ macro steps give when they add up on the
same cells, as next to a face between levels that travels with a wave at its
speed.

### Proposition (the margin on condition (4))

*Setting.*  $\Delta t_c$ is the admissible step of cell $c$ for the attempt: the local step at $t_n$, a positive finite
number ([Local steps](multirate-method.md#levels-and-macro-step)), lowered by a redo to the smallest value that the cell showed;
$\Delta t_c(\tau)$ its value at the start $\tau$ of a sub-step, as the certificate evaluates it.  $v \sim c$ runs over the
face neighbours of $c$ in its block, every one of them on the process of $c$.  $\varepsilon = 10^{-12}$ and

$$
H_{\text{cap}} = \min_c\ 2^{\ell_c} \min\Bigl(\frac{\Delta t_c}{\sigma},\ \min_{v \sim c} \frac{\Delta t_v}{\kappa}\Bigr)
$$

over all the processes, for a nested map with equal levels across connections.

*Hypotheses.*  (M1) $\sigma \ge 2$, $L_{\max} > 0$, the interface margin on (with the redo or alone).
(M2) $\kappa = \max(\text{cfl}\,(1+\theta), \kappa_0)$ with $1 < \kappa_0 \le 2$ and $\text{cfl}\,(1+\theta) \le \sigma$ (from
cfl $\le 1$, which the redo imposes, with $\theta \le 1$ and $\sigma \ge 2$), so $1 < \kappa \le \sigma$.  (M3) $H \le (1 + \varepsilon)\, H_{\text{cap}}$: true of every
attempt with an adaptive macro step, since the redo halves $H$ while it is above that bound of the map it has rebuilt, a target and
a halving only shorten $H$, and a landing lengthens it by at most $\varepsilon H$ plus half a unit in the last place of the time
(a relative $1.1\cdot10^{-16}\,t/H$, above $\varepsilon$ when $H < 10^{-4}\,t$); (a) then holds with that larger $\varepsilon$, without
consequence next to $1 - 1/\kappa_0$. No hypothesis on the mesh, the state, the limiter, the fluxes, the number of dimensions, of blocks or of
processes.

*Statement.*

(a) If $\Delta t_v(\tau) \ge \Delta t_v/\kappa$ for every cell $v$ and every start $\tau$ of its sub-steps in the attempt,
the certificate finds no violation.  So every rejection of a violation of (4) comes from a cell whose admissible step fell by more
than $1 - 1/\kappa \ge 1 - 1/\kappa_0$ within the macro step.

(b) The margin raises a cell $c$ only where a face neighbour has $\Delta t_c/\Delta t_v > \sigma/\kappa$.  Where no face
of the mesh has that jump, the levels and $H$ are those of the rule without the margin, to the bit.

(c) Before nesting and connections the level of a cell is at most the largest level that its own term gives to the cell or to
one of its face neighbours.  So the margin does not propagate through the rule, and the deepest level does not change.

(d) In the adaptive step without a target, on a map built for $H_{\text{try}} = 2^{L}\Delta t_{\min}/\sigma$,
$2^{-\varepsilon} H_{\text{try}} \le H_{\text{cap}} \le H_{\text{try}}$: the first $H$ does not depend on $\kappa$, to a
relative $7\cdot10^{-13}$.

*Proof.*  (a) Take $v$ and a start $\tau$ of a sub-step.  If a face neighbour $g$ of $v$ in the block is one level coarser and
starts a sub-step at $\tau$, then $\hat h_v = h_{\ell_g} = H/2^{\ell_g}$, and $v \sim g$ gives in $H_{\text{cap}}$ the
term $2^{\ell_g} \Delta t_v/\kappa$, so $\hat h_v \le (1 + \varepsilon) \Delta t_v/\kappa$ by (M3).  Otherwise $\hat h_v =
h_{\ell_v} \le (1 + \varepsilon) \Delta t_v/\sigma \le (1 + \varepsilon) \Delta t_v/\kappa$ by (M2).

A neighbour across a connection has
the level of $v$, so it is never the coarser one, and a finer neighbour does not enter $\hat h_v$.  By hypothesis the
last member is at most $(1 + \varepsilon) \Delta t_v(\tau)$, the test of the certificate.  (b) The rule takes for $c$ the larger
of $\sigma H/\Delta t_c$ and $\kappa H/\Delta t_v$, $v \sim c$; the second exceeds the first only if
$\Delta t_c/\Delta t_v > \sigma/\kappa$.

Without that jump every cell has the first in the rule, and in
$H_{\text{cap}}$ the term of $v$ is not below that of $c$.  (c) $\kappa H/\Delta t_v \le \sigma H/\Delta t_v$ by (M2), and the
formula of the levels is monotone: the margin term of $c$ asks at most the level of the own term of $v$; it reads the local
steps of the neighbours, not their levels.  (d) By (c) the cell with $\Delta t_{\min}$ keeps level $L$, the rules cannot
take it higher (their fixed point does not exceed the largest level of the rule), and its term is $H_{\text{try}}$, so
$H_{\text{cap}} \le H_{\text{try}}$.

Every level satisfies $\log_2 x \le \ell_c + \varepsilon$ for both terms $x$ of the rule, so every
term of $H_{\text{cap}}$ is at least $2^{-\varepsilon} H_{\text{try}}$; no level is cut at $L_{\max}$, since $x \le 2^{L}$ by
(c).  ∎

*The extremes.*  cfl $\to 0$ gives $\kappa = \kappa_0$, the same relative fall.  At $\text{cfl}\,(1+\theta) = \kappa_0$ the two
terms of $\kappa$ meet, and at $\text{cfl}\,(1+\theta) = \sigma$, which (M1) and (M2) allow only with $\sigma = 2$, $\theta = 1$ and
cfl 1, the margin is the term of the cell itself.  With $\sigma > 2$ a coarse neighbour of a cell at the limit of its level is raised only if $\kappa > \sigma/2$;
(a) holds as it stands.

With $\sigma < \kappa$, outside (M1) and (M2), a cell at the limit of its own term has, at a start
that no coarser neighbour shares, $\hat h_v = \Delta t_v/\sigma$, and a fall of its admissible step by more than
$1 - 1/\sigma$ rejects the attempt: with $\sigma = 1$, every fall.  A cell with no neighbour in its block has no margin term.  At $L_{\max}$ a level can be cut, after a redo,
and (a) holds in that case too, because it rests on the term of $H_{\text{cap}}$ and not on the level.

With $L_{\max} = 0$ the margin
is off and there is no face between levels.  At $\kappa = 1$, $\sigma = 2$ (a) leaves no room: a cell at the limit of its level with a
coarser neighbour has $\hat h_v = \Delta t_v$.

*What it does not say.*  It does not say that (4) is sufficient (the source gives no proof either [13]); that a front is
not rejected (a fall by more than $1 - 1/\kappa$ is necessary for a rejection, not sufficient; where it happens, (4) asks for it); how much the raised cells cost, which depends on how many
faces have the jump of (b); or why $\kappa_0$ is 1.1, a value taken from another code (section
[Non-local CFL condition](multirate-method.md#non-local-cfl-condition)).

With a fixed macro step (M3) fails, and there the certificate only counts.
A map with the margin is a nested map, so Theorem 1, Theorem 2, the Corollary and Appendix B hold for it as for any other.

## Appendix B: Total variation at a face between levels

This appendix proves the bound of the total variation that the section
[Non-local CFL condition](multirate-method.md#non-local-cfl-condition) states for the limited reconstruction, in a restricted setting: the
advection equation in one dimension, two levels and one face between them (Theorem 3), and three levels with a band of one
cell (Theorem 3b), with the reconstruction of MOSE on cells of different widths.  It is a bound, not the TVD
property: one macro step can raise the total variation, by terms local to the faces between levels (Remark 6).

**Setting.**

- (B1) $u_t + a\,u_x = 0$, $a > 0$ constant, on cells $I_j$ of widths $\Delta x_j$, with cell averages $u_j$ and a constant
  state on each side far from the face between levels; $d_j = u_j - u_{j-1}$ is the difference across the face $j - \tfrac12$
  and $TV(u) = \sum_j |d_j|$.  The weights $w^+_j = \Delta x_j/(\Delta x_j + \Delta x_{j+1})$ and
  $w^-_j = \Delta x_j/(\Delta x_{j-1} + \Delta x_j)$ satisfy $w^+_{j-1} + w^-_j = 1$ and are $\tfrac12$ on cells of equal
  width.  The face between levels lies between the cells 0 and 1, the cell 1 downstream; $\rho = \Delta x_0/\Delta x_1$ is the
  ratio of the upstream to the downstream width there, and $r$ the ratio of the coarse to the fine width ($\rho = r$ with the
  coarse cell upstream, $\rho = 1/r$ with the fine cell upstream).
- (B2) The reconstruction of MOSE, with the distances between the cell centres
  $l^0_j = (\Delta x_{j-1} + \Delta x_j)/2$ and $l^1_j = (\Delta x_j + \Delta x_{j+1})/2$ and the half width $\Delta x_j/2$, and the
  upwind flux $a\,q_{j+1/2}$: $q_{j+1/2} = u_j + s_j$ with $s_j = \hat\phi_j\,w^+_j\,\Delta_+u_j$,
  $r'_j = (\Delta_-u_j/\Delta_+u_j)\,(l^1_j/l^0_j)$, $\hat\phi_j = \min(\phi(r'_j),\ 1/w^+_j,\ r'_j/w^-_j)$ for $r'_j > 0$ and
  $\hat\phi_j = 0$ for $r'_j \le 0$, and $s_j = 0$ when $\Delta_+u_j = 0$.  The limiter receives the bounds
  $l^1_j/(\Delta x_j/2) = 1/w^+_j$ and $l^0_j/(\Delta x_j/2) = 1/w^-_j$, which keep the face value between $u_j$ and $u_{j+1}$;
  with minmod and the first-order reconstruction they never act.  The small term ($10^{-30}$) added to the
  denominator of the ratio of the slopes is left out.  On cells of equal width the bounds are 2 and $2r'_j$, which no limiter of
  (B3) exceeds, and this is the form of the interface margin above.
- (B3) The limiter: $\phi(r) = 0$ for $r \le 0$, and $0 \le \phi(r) \le 2\theta$, $0 \le \phi(r)/r \le 2\theta$ for $r > 0$,
  with the $\theta$ of the interface margin: 0 with the first-order reconstruction, $\tfrac12$ for minmod, 1 for van Leer, MC
  and superbee.  Van Albada, whose $\phi$ is negative for $-1 < r < 0$, is outside.
- (B4) Two levels with the macro step of MOSE: the coarse cells advance with $H$ and the fine cells with $h = H/2$, with the
  Courant numbers $\lambda_j = aH/\Delta x_j$ and $\mu_j = ah/\Delta x_j$.  The coarse stage 1 at $t_n$; the two Heun
  sub-steps of the fine cells, whose four stages read the coarse cells through the interpolant at the fractions
  $0, \tfrac12, \tfrac12, 1$ of $H$; the coarse stage 2, which reads the fine cells at their final values.  On the face between
  levels the coarse cell uses and saves its flux $G_s$ at its stage 1, each fine stage adds a quarter of its flux to $A$, and
  the coarse cell uses $2A - G_s$ at its stage 2 ([Conservative correction](multirate-method.md#conservative-correction)).
- (B5) Every cell satisfies $c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$, with $c_j$ its own Courant number ($\lambda_j$ or
  $\mu_j$); on cells of equal width this is $c_j\,(1+\theta) \le 1$ (Remark 5).

A forward-Euler stage on a state $w$ (the values that the stage reads on its whole stencil) updates the cells of its level as
$w'_j = w_j - c_j\,(q_{j+1/2} - q_{j-1/2})$, with the face values of $w$.

**Lemma 1 (incremental form of a stage).**  For every state $w$ and cell $j$ there are
$\alpha_j \in [0, \min(2\theta w^+_j, 1)]$ and $\beta_j \in [0, \min(2\theta w^-_j, 1)]$ with
$s_j = \alpha_j\,(w_{j+1} - w_j) = \beta_j\,(w_j - w_{j-1})$, so that

$$
q_{j+1/2} - q_{j-1/2} = B_j\,(w_j - w_{j-1}), \qquad
B_j = 1 + \beta_j - \alpha_{j-1} \in [\,1 - \min(2\theta w^+_{j-1}, 1),\ 1 + \min(2\theta w^-_j, 1)\,] \subset [0, 2].
$$

*Proof.*  If $\Delta_+w_j = 0$ or $r'_j \le 0$, then $s_j = 0$: take $\alpha_j = \beta_j = 0$.  Otherwise take
$\alpha_j = w^+_j\,\hat\phi_j$ and, since $(l^1_j/l^0_j)\,w^+_j = w^-_j$, $\beta_j = w^-_j\,\hat\phi_j/r'_j$.  By (B2) and (B3),
$0 \le \hat\phi_j \le \phi(r'_j) \le 2\theta$ and $0 \le \hat\phi_j/r'_j \le \phi(r'_j)/r'_j \le 2\theta$, and by (B2)
$\hat\phi_j \le 1/w^+_j$ and $\hat\phi_j/r'_j \le 1/w^-_j$: so $\alpha_j \le \min(2\theta w^+_j, 1)$ and
$\beta_j \le \min(2\theta w^-_j, 1)$.  Then $q_{j+1/2} - q_{j-1/2} = (w_j - w_{j-1}) + s_j - s_{j-1}$, with
$s_j = \beta_j\,(w_j - w_{j-1})$ and $s_{j-1} = \alpha_{j-1}\,(w_j - w_{j-1})$, and $\alpha_{j-1} \le 1$, $\beta_j \le 1$ give
$0 \le B_j \le 2$. ∎

So $B_j \ge 0$ at every face, the face between levels included, on cells of any widths.  The bound of $\beta_j$ does not depend
on the width of the downstream cell: the slope of a cell that leaves through a face is written with its upstream difference,
whatever lies downstream, and the two sides of the face between levels can be treated apart.

**Lemma 2 (a region with a ghost upstream).**  Let $J = \{j_a, j_a + 1, \dots\}$ be consecutive cells that end in a constant
state or in an outflow cell $j_b$ (one whose slope reads a cell outside $J$), $g = w_{j_a - 1}$ the value upstream of $J$ (the
ghost), and $TV_J(w; g) = |w_{j_a} - g| + \sum_{j \in J,\ j > j_a} |w_j - w_{j-1}|$.  If (B5) holds in $J$, a forward-Euler stage on
$J$ with the ghost fixed gives

$$
TV_J(w'; g) \le TV_J(w; g) - c_{j_b}B_{j_b}\,|w_{j_b} - w_{j_b - 1}|,
$$

without the last term when $J$ ends in a constant state.

*Proof.*  By Lemma 1, $w'_{j_a} - g = (1 - c_{j_a}B_{j_a})(w_{j_a} - g)$ and
$d'_j = (1 - c_jB_j)\,d_j + c_{j-1}B_{j-1}\,d_{j-1}$ for $j > j_a$.  By Lemma 1 and (B5),
$0 \le c_jB_j \le c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$ in $J$.  Taking absolute values and summing, each
$c_{j-1}B_{j-1}|d_{j-1}|$ cancels the $-c_{j-1}B_{j-1}|d_{j-1}|$ of the row before; the outflow term remains. ∎

**Lemma 3 (a Heun sub-step).**  (a) Let $\tilde w = \tfrac12 w + \tfrac12 E_2(E_1(w))$, with $E_k$ a forward-Euler stage on $J$
with the ghost $g_k$, under the hypotheses of Lemma 2.  Then
$TV_J(\tilde w; g_2) \le TV_J(w; g_1) + |g_2 - g_1|$.
(b) If $J$ starts at a constant state and ends in an outflow cell $j_b$, then
$TV_J(\tilde w) \le TV_J(w) - D$, with
$D = \tfrac12\,c_{j_b}\,(B^{(1)}_{j_b}|d^{(1)}_{j_b}| + B^{(2)}_{j_b}|d^{(2)}_{j_b}|) \ge |\tilde w_{j_b} - w_{j_b}|$, the
superscripts naming the two stages.

*Proof.*  $TV_J(\cdot\,; g)$ is convex.  (a) $TV_J(\tilde w; g_2) \le \tfrac12\,TV_J(w; g_2) + \tfrac12\,TV_J(E_2E_1w; g_2)$, with
$TV_J(w; g_2) \le TV_J(w; g_1) + |g_2 - g_1|$ and, by Lemma 2 twice,
$TV_J(E_2E_1w; g_2) \le TV_J(E_1w; g_1) + |g_2 - g_1| \le TV_J(w; g_1) + |g_2 - g_1|$.  (b) Lemma 2 twice, and
the convexity; $\tilde w_{j_b} - w_{j_b} = -\tfrac12\,c_{j_b}(B^{(1)}_{j_b}d^{(1)}_{j_b} + B^{(2)}_{j_b}d^{(2)}_{j_b})$. ∎

### Theorem 3 (total variation at a face between levels)

*Statement.*  Under (B1)–(B5), one macro step with one face between levels gives

$$
TV(u^{n+1}) \le TV(u^n) + R,
$$

with $R$ as follows.

- *Coarse cell upstream* (cells $j \le 0$ coarse, $j \ge 1$ fine).  Let $p_0$ be the predictor of the cell 0,
  $d^p_0 = p_0 - p_{-1}$, $X = B^0_0\,d_0$ and $Y = B^1_0\,d^p_0$ the increments of Lemma 1 of the cell 0 at its two stages
  without the correction, and $\delta_0$ the change that the correction makes to the cell 0.  Then

  $$
  R = \tfrac12\lambda_0\,(|X| - |Y| + |X - Y|) + 2|\delta_0|
  \le \lambda_0\,\bigl[(1 + \tfrac32\theta)\,|d_0| + \theta\,|d^p_0|\bigr]
  \le \lambda_0\,\bigl[(1 + \tfrac52\theta)\,|d_0| + \theta\,|d_{-1}|\bigr],
  $$

- *Fine cell upstream* (cells $j \le 0$ fine, $j \ge 1$ coarse).  Let $\Phi = 2A - G_s$ and $\Phi_H$ be the face values that
  the cell 1 uses at its stage 2 with and without the correction (without it, from the fine cells at their final values),
  $d^{(s)}_0$ the difference of the last fine cell at the fine stage $s$ ($d^{(1)}_0 = d_0$) and $d^{n+1}_0$ at the end, and
  $m = \mu_0\,(1 + \theta) \le 1$.  Then

  $$
  R = \lambda_1\,|\Phi - \Phi_H| \le \lambda_1\,\bigl[(m + 3\theta)\,|d_0| + (m + 4\theta)\,|d_{-1}|
  + (\tfrac34 m + 3\theta)\,|d_{-2}| + (\tfrac14 m + \tfrac54\theta)\,|d_{-3}| + \tfrac14\theta\,|d_{-4}|\bigr],
  $$

The bounds of $R$ take the two cells upstream of the face of equal width (otherwise $\theta$ becomes $\min(2\theta w^-_0, 1)$ in
them).

*Proof, coarse cell upstream.*

1. *Regions.*  The coarse region $C = \{\dots, -1, 0\}$ reads the fine cells only in the slope $s_0$, which leaves through the
   face $\tfrac12$: written as $s_0 = \beta_0\,(w_0 - w_{-1})$ with $\beta_0 \in [0, \min(2\theta w^-_0, 1)]$, whatever the
   fine cell, it makes $C$ a region with the outflow cell 0, and $B_j \ge 0$ in $C$.  The fine region $F = \{1, 2, \dots\}$ reads
   the cells $-1$ and 0 through the interpolant: its ghost is $g(\tau) = u_0 + \tau\,(p_0 - u_0)$ at the fraction $\tau$ of $H$,
   and $w_{-1}$ enters only $\alpha_0 \in [0, \min(2\theta w^+_0, 1)]$, $w^+_0 = \rho/(\rho + 1)$; $B_j \ge 0$ for $j \ge 1$
   (Lemma 1).  $TV(u) = TV_C(u) + TV_F(u_F; u_0)$.
2. *Fine region.*  The four fine stages have the ghosts $g(0)$, $g(\tfrac12)$ and $g(\tfrac12)$, $g(1)$.  Lemma 3(a) on the
   two fine sub-steps (entry cell 1, $c = \mu_1$) gives $TV_F(u^{n+1}_F; p_0) \le TV_F(u_F; u_0) + |p_0 - u_0|$, the path
   of the ghost being monotone, of length $|p_0 - u_0|$; then
   $TV_F(u^{n+1}_F; u^{n+1}_0) \le TV_F(u^{n+1}_F; p_0) + |u^{n+1}_0 - p_0|$.
3. *Coarse region without the correction.*  The cells $j \le -1$ take the coarse Heun step, whose faces read no fine cell.
   With $v_0$ the Heun value of the cell 0 with its own flux at both stages, Lemma 3(b) gives
   $TV_C(\dots, u^{n+1}_{-1}, v_0) \le TV_C(u) - \tfrac12\lambda_0\,(|X| + |Y|)$.
4. *Correction.*  The cell 0 receives on the face $\tfrac12$ the mean of the four fine face values in place of the mean of its
   own two.  With $q_{1/2} = w_0 + \beta_0\,(w_0 - w_{-1})$ on each state, the parts $w_0$ have the same mean,
   $\tfrac14\,(u_0 + 2g(\tfrac12) + p_0) = \tfrac12\,(u_0 + p_0)$, and the interpolated differences are
   $(1 - \tau)\,d_0 + \tau\,d^p_0$, so that

   $$
   \delta_0 = -\lambda_0\,\bigl[(-\tfrac14\beta^{(1)} + \tfrac18\beta^{(2)} + \tfrac18\beta^{(3)})\,d_0
   + (\tfrac18\beta^{(2)} + \tfrac18\beta^{(3)} + \tfrac14\beta^{(4)} - \tfrac12\beta^{c})\,d^p_0\bigr],
   $$

   with $\beta^{(s)}$ that of the cell 0 at the fine stage $s$ and $\beta^c$ at its own stage 2: $|\delta_0| \le
   \lambda_0\,\theta\,(\tfrac14|d_0| + \tfrac12|d^p_0|)$ with the cells $-1$ and 0 of equal width.  Only the last difference of
   $C$ changes: $TV_C(u^{n+1}) \le TV_C(\dots, u^{n+1}_{-1}, v_0) + |\delta_0|$.
5. *Sum.*  With $|p_0 - u_0| = \lambda_0|X|$, $v_0 - p_0 = \tfrac12\lambda_0\,(X - Y)$ and
   $|u^{n+1}_0 - p_0| \le |v_0 - p_0| + |\delta_0|$, steps 2–4 give
   $TV(u^{n+1}) - TV(u^n) \le \lambda_0|X| + \tfrac12\lambda_0|X - Y| - \tfrac12\lambda_0(|X| + |Y|) + 2|\delta_0| = R$.
   Then $\tfrac12(|X| - |Y| + |X - Y|) \le |X| \le (1 + \theta)\,|d_0|$, and $|d^p_0| \le |d_0| + |d_{-1}|$ by Lemma 2 on the
   coarse cells, whose coefficients lie in $[0, 1]$. ∎

*Proof, fine cell upstream.*

1. *Regions.*  The fine region $F = \{\dots, -1, 0\}$ reads the cell 1 only in the slope $s_0 = \beta_0\,(w_0 - w_{-1})$,
   $\beta_0 \in [0, \min(2\theta w^-_0, 1)]$: it is a region with the outflow cell 0, with $B_j \ge 0$.  The coarse region
   $G = \{1, 2, \dots\}$ has the fine cell 0 as its ghost, and $B_1 = 1 + \beta_1 - \alpha_0$ lies in
   $[1 - \min(2\theta\rho/(\rho + 1), 1),\ 1 + \min(2\theta/(\rho + 1), 1)] \subset [0, 2]$; $B_j \ge 0$ for $j \ge 2$.  $TV(u) = TV_F(u_F) + TV_G(u_G; u_0)$.
2. *Fine region.*  Lemma 3(b) on the two fine sub-steps gives $TV_F(u^{n+1}_F) \le TV_F(u_F) - D_1 - D_2$, with
   $|u^{n+1}_0 - u_0| \le D_1 + D_2$.
3. *Coarse region without the correction.*  The coarse Heun step reads the ghost at $u_0$ at its stage 1 (face value $G_s$) and
   at $u^{n+1}_0$ at its stage 2 (face value $\Phi_H$); Lemma 3(a), with the entry cell 1 and $c = \lambda_1$, gives
   $TV_G(u^H; u^{n+1}_0) \le TV_G(u_G; u_0) + |u^{n+1}_0 - u_0|$.  Here (B5) of the cell 1 is used, since $B_1$ reaches
   $1 + \min(2\theta/(\rho + 1), 1)$.
4. *Correction.*  It changes the cell 1 only, by $\tfrac12\lambda_1\,(\Phi - \Phi_H)$, and the cell 1 enters two terms of
   $TV_G$: $TV_G(u^{n+1}_G; u^{n+1}_0) \le TV_G(u^H; u^{n+1}_0) + \lambda_1|\Phi - \Phi_H|$.
5. *Sum.*  Steps 2–4 give $TV(u^{n+1}) \le TV(u^n) + R$.
6. *The form of $R$.*  With $A = \tfrac14\sum_s q_{1/2}(S_s)$ over the states $S_s$ of the four fine stages, $G_s = q_{1/2}(S_1)$
   and $q_{1/2}(S) = w_0 + \beta_0\,(w_0 - w_{-1})$ on each state, $\Phi - \Phi_H = \tfrac12\,[q_{1/2}(S_2) + q_{1/2}(S_3) +
   q_{1/2}(S_4)] - \tfrac12\,q_{1/2}(S_1) - q_{1/2}(S^{n+1})$, with $S^{n+1}$ the fine cells at their final values.  The parts
   $w_0$, through the Heun stages of the fine cell 0, give $\tfrac12\mu_0\,(B^{(4)}_0d^{(4)}_0 - B^{(1)}_0d^{(1)}_0)$, and the
   slopes give $\tfrac12\,(\beta^{(2)}d^{(2)}_0 + \beta^{(3)}d^{(3)}_0 + \beta^{(4)}d^{(4)}_0) - \tfrac12\beta^{(1)}d^{(1)}_0 -
   \beta^{n+1}d^{n+1}_0$; with $B \le 1 + \theta$ and $\beta \le \theta$ on cells of equal width,
   $R \le \lambda_1\,[\tfrac12 m\,(|d^{(1)}_0| + |d^{(4)}_0|) + \theta\,(\tfrac12\sum_s |d^{(s)}_0| + |d^{n+1}_0|)]$.
7. *Locality.*  In $F$ a stage gives $|d'_j| \le |d_j| + |d_{j-1}|$ (Lemma 2, coefficients in $[0, 1]$): on the coefficients of
   $|d_0|, |d_{-1}|, \dots$ a stage adds the vector shifted by one cell upstream.  From $[1]$ this gives $[1, 1]$ after the
   first stage, $[1, 1, \tfrac12]$ at the end of the first sub-step, $[1, 2, \tfrac32, \tfrac12]$ after the third stage and
   $[1, 2, 2, 1, \tfrac14]$ at the end, which in step 6 give the bound of the statement. ∎

### Corollary (a solution Lipschitz continuous next to the face)

Let $\Delta x_*$ be the largest distance between the centres of two neighbouring cells among those whose differences $R$
reads (with cells of equal width within each level, $\Delta x_* = \Delta x_0$).  If $u$ is Lipschitz continuous with the constant
$K$, the cell averages have $|u_j - u_{j-1}| \le K\,(\Delta x_{j-1} + \Delta x_j)/2$.  If up to $T = NH$ every difference that
$R$ reads is at most $K\Delta x_*$, Theorem 3 summed over the $N$ macro steps gives, with the coarse cell upstream (the
first bound of $R$),

$$
TV(u^N) \le TV(u^0) + (1 + \tfrac52\theta)\,a K T\,\frac{\Delta x_*}{\Delta x_0},
$$

and with the fine cell upstream

$$
TV(u^N) \le TV(u^0) + (1 + 3\theta)\,a K T\,\frac{\Delta x_*}{\Delta x_1} .
$$

This is the form $TV(\bar u^n) \le TV(\bar u^0) + CT$ of Theorem 4.2 of Hoang et al. [11], here for the finite-volume scheme of
MOSE, and like theirs it asks the solution to be regular next to the face.  Each estimate is for one macro step and carries no
memory from one macro step to the next: it holds for every position of the face at every macro step, that is for every
sequence of two-level maps, with the regularity next to the face at every macro step.

Of the step, the proof uses only (B5), the own Courant number of every cell with the weights of the widths of (B1): not the coarse step measured on the first fine cell, which
the hypothesis (4.15) of [11] bounds and which in MOSE the interface margin bounds.

*Remark 5 (the step condition).*  With one time step for all the cells of a level, $c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$
follows from $c\,(1 + \theta) \le 1$, and at a face from a coarser to a finer cell (B5) asks nothing more.  At a face from a
finer to a coarser cell with $r > 1$ it asks, of the coarse cell downstream, $\lambda_1\,(1 + \min(2\theta r/(r + 1), 1)) \le 1$:
with $\theta = 1$ (van Leer, MC, superbee) this is $\lambda_1\,(1 + \theta) \le 1$; with $\theta = \tfrac12$ (minmod) it is
$\lambda_1\,(1 + r/(r + 1)) \le 1$, which $\lambda_1\,(1 + \theta) \le 1$ does not give.

Where the convective limit rules,
$\lambda_1 \le \text{cfl}/\sigma$, so with minmod $\text{cfl} \le \sigma/(1 + r/(r + 1))$ suffices, and $\text{cfl} \le \sigma/2$
for every $r$: with $\sigma = 2$, every $r$ with cfl 1, the largest cfl that the redo allows when $L_{\max} > 0$, so
the bound on cfl gives (B5) for every $r$ and every limiter.  It is a hypothesis of the proof, not a condition shown to be
necessary.

*Remark 6 (what the bound does not give).*  The total variation can grow at a face between levels, at first order too: with
the coarse cell upstream, $\lambda_0 = \mu_1 = 1$, cells of widths 1 and $\tfrac12$, $u_j = 0$ for $j \le -1$ and $u_j = 1$
for $j \ge 0$, the macro step gives 0.5 to the cell 0 and 0.375, 0.75, 0.875, 1 to the fine cells: the values stay in
$[0, 1]$, as Theorem 2 says, and the total variation goes from 1 to 1.25.

Away from the face between levels, on cells of any
widths, the stages of one level do not raise the total variation (Lemma 2).  More than three levels, three levels in another arrangement, several dimensions, systems, nonlinear fluxes, and the
positivity and the maximum principle at the faces between levels are outside this appendix.

### Theorem 3b (three levels, a band of one cell)

*Statement.*  Under (B1)–(B3) and (B5), let the map have three levels, with steps $H$, $H/2$ and $H/4$,
and one cell of the middle level (the band) between the coarse and the fine cells, in either order along the flow; the macro step is
the recursion above, every coarser cell read by a stage is taken at its interpolant at the stage time (the coarse cells also by the
stages of the fine cells, two levels finer), and the correction acts on both faces between levels, each with its own accumulator, the
band being the fine side of one face and, once per sub-step, the coarse side of the other.

Then one macro step gives
$TV(u^{n+1}) \le TV(u^n) + R_A + R_B$, with $R_A$ the remainder $R$ of Theorem 3 at the face between the coarse cells and the band
(the band in place of the fine cell or of the fine region) and $R_B$ the sum, over the two sub-steps of the band, of the remainder $R$
of Theorem 3 at the face between the band and the fine cells, with the band as the coarse cell and its own Courant number
$\mu_b = aH/(2\Delta x_b)$; with the fine cells upstream, $R_A$, the remainder of the downstream face, also carries one half of the remainder of the
upstream face in the second sub-step of the band, the correction that the band receives there.

With the coarse cells upstream every
difference read by $R_A + R_B$ lies at $t_n$ in the cells $-1$, $0$ and the band, and with $\theta = 1$ and a width ratio of 2,
$R_A + R_B \le \lambda_0\,(3.5\,|d_0| + |d_{-1}|) + \mu_b\,(7.34\,|d_b| + 5.34\,\lambda_0\,|d_0|)$, so that
$TV(u^N) \le TV(u^0) + 17.2\,aKT\,\Delta x_*/\Delta x_0$ when those differences stay below $K\Delta x_*$ up to $T = NH$.  With the fine
cells upstream the remainders read the differences of the last fine cell and of the band at the stages of the band; with $\theta = 1$
and a width ratio of 2, if every difference they read stays below $K\Delta x_*$, $TV(u^N) \le TV(u^0) + 13\,aKT\,\Delta x_*/\Delta x_2$,
$\Delta x_2$ the width of the coarse cells.

(B5) is required of the band with its own Courant number and, with the fine cells
upstream, of the first coarse cell after the band with $\min(2\theta \cdot 2/3, 1)$ (Remark 5).  As Theorem 3, this is a bound and not the
TVD property: the total variation can grow at a face between levels (Remark 6).

*Outline of the proof.*  Lemmas 1–3 hold as they stand.  The one new step is Lemma 3 written for a region of one cell: a band $b$ with
Courant number $c$ that takes a Heun sub-step with the ghosts $g_1$ and $g_2$, with $X$ and $Y$ the increments of Lemma 1 at its two
stages and $\tilde v_b$ the value of the sub-step, has $|\tilde v_b - g_2| \le |w_b - g_1| + |g_2 - g_1| - \tfrac12 c\,(|X| + |Y|)$.

The recursion is
nested: within each sub-step of the band, the band and the fine cells form the configuration of Theorem 3, with the band as the
coarse cell; on the scale of the macro step, the band and the fine cells together take the place of the fine region of Theorem 3
towards the coarse cells, with the same path of the ghost.  With the fine cells upstream the band reads the fine cell at its current
value and the coarse cells read the band at its final value; the correction that the band receives in its second sub-step reappears
in the remainder of the downstream face. Each face between levels then gives its own remainder.

---

## References

[1] C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comput. Phys.*, 77(2), 1988.

[7] M. H. Carpenter, D. Gottlieb, S. Abarbanel, W.-S. Don, "The theoretical accuracy of Runge–Kutta time discretizations for the initial boundary value problem: a study of the boundary error," *SIAM J. Sci. Comput.*, 16(6), 1241–1252, 1995, doi:10.1137/0916072 (ICASE Report 93-83, 1993).

[8] W. Hundsdorfer, *Numerical Solution of Advection–Diffusion–Reaction Equations*, lecture notes, Thomas Stieltjes Institute, CWI, Amsterdam, 2000.

[10] G. Jeanmasson, I. Mary, L. Mieussens, "On some explicit local time stepping finite volume schemes for CFD," *J. Comput. Phys.*, 397, 2019, 108818, doi:10.1016/j.jcp.2019.07.017 (Prop. RK2: section on the second-order scheme and its conservative correction).

[11] T.-T.-P. Hoang, L. Ju, W. Leng, Z. Wang, "High order explicit local time stepping methods for hyperbolic conservation laws," *Math. Comp.*, 89(324), 2020, 1807–1842, doi:10.1090/mcom/3507 (arXiv:1905.09705).

[13] N. Y. Gnedin, V. A. Semenov, A. V. Kravtsov, "Enforcing the Courant–Friedrichs–Lewy condition in explicitly conservative local time stepping schemes," *J. Comput. Phys.*, 359, 2018, 93–105, doi:10.1016/j.jcp.2018.01.008 (arXiv:1801.03108; section 2, conditions (3) and (4)).
