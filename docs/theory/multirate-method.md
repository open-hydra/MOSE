# Multirate: Method in Detail

This page details the [multirate integrator](multirate.md): how the levels are built, the checks on the local steps, the non-local
CFL condition, the recursion and the treatment of boundaries, turbulence and MPI. Read [Multirate Time Stepping](multirate.md) first.

In a time-accurate run with the global step every cell advances with the
smallest admissible step of the mesh, $\Delta t = \min_i \Delta t_i$.  When a
few cells (near walls, in a throat, in a refined band) set that minimum, most
of the work is spent advancing the other cells with steps far below their own
limit.  The **multirate** integrator lets each cell advance with a step close
to its own limit, keeping, on a fixed mesh, the second order in time, and the
exact conservation of the single-rate scheme.

The integrator has variants, which the sections below name.  The level map is
either rebuilt from the local steps every $k$ macro steps (a *recomputed map*,
for production runs) or built at the first macro step and kept (a *fixed
map*, for studies of the order).  The macro step is either adaptive or a
*fixed macro step*, a given multiple of the smallest local step of the first
macro step, for verification.  Condition (4) of the
[non-local CFL condition](#non-local-cfl-condition) is either imposed, by
rejecting and redoing the macro steps that violate it (the *redo*, the
default), only counted, or not checked; the *interface margin* comes with the
redo and can also be applied alone. The input keys of these choices are in
the [User Guide](../user/input.md#multirate-time-stepping).

## Levels and macro step

Each cell $c$ gets a **level** $\ell_c \in \{0,\dots,L\}$ and advances with
sub-steps $h_\ell = H/2^{\ell}$, where $H$ is the **macro step** and $L$ the
depth of the map, its deepest level.  The level
comes from the admissible local step $\Delta t^{\text{cap}}_c$ of the cell (the
[local time step](time-integration.md#local-time-step) above) and a safety factor $\sigma \ge 1$
(2 by default, see [Non-local CFL condition](#non-local-cfl-condition)), with $L_{\max}$ the deepest level
allowed, for the step $H$ for which the map is built:

$$
\ell_c = \min\Bigl(L_{\max},\ \max\bigl(0,\ \bigl\lceil \log_2(\sigma H/\Delta t^{\text{cap}}_c) - 10^{-12} \bigr\rceil\bigr)\Bigr),
$$

then two rules are imposed, iterated until nothing changes:

- **nesting**: the levels of two cells sharing a face differ by at most one;
- **equal levels across block connections** (direct and periodic): the two
  sides of a connection always advance together, so a connection is never a
  face between levels: its one flux per stage, which the record that owns the
  face computes for both cells, needs the two sides at the same stage times
  ([Blocks and boundary conditions](#blocks-and-boundary-conditions)), the
  correction stays within a process, and such a map is nested, so Appendix A
  holds for it.  SeisSol keeps the elements of a dynamic rupture face in one time
  cluster in the same way [18].

The macro step is then

$$
H = \min\Bigl(H_{\text{target}},\ \min_c \Delta t^{\text{cap}}_c\,2^{\ell_c}/\sigma\Bigr),
$$

so every cell advances with $h_{\ell_c} \le \Delta t^{\text{cap}}_c/\sigma$.  The
$10^{-12}$ in the levels, about 4500 units in the last place of a double, keeps a cell at the
limit of a level from being pushed one level finer by the rounding of
$\sigma H/\Delta t^{\text{cap}}_c$ and of its logarithm, which stays below $10^{-13}$; the other comparisons
of steps of the integrator use the same tolerance, as Castro does on the landing of its sub-steps [16].

**Local steps.**  The levels, the depth of the map and $H$ come from the
admissible local steps, so each of them must be a positive finite number.  At
the start of every macro step, in every variant, a cell
whose $\Delta t^{\text{cap}}_c$ is zero, negative, infinite or not a number stops
the run before the map is built: a negative step would make $H$ negative, and
the halving of $H$ while it exceeds a negative $H_{\text{cap}}$ in the redo of
[Non-local CFL condition](#non-local-cfl-condition) would never end; a step
that is not a number can make $H$ not a number, which the floor of the redo
never stops.

Such a step comes from a state that is not physical: a
pressure and a density of opposite signs make the temperature $p/(\rho R)$
negative and the speed of sound not a number; a density and a pressure both
negative make the diffusivity $\mu/\rho$ of a viscous run negative, and the
step negative where the diffusion rate outweighs the convective one; a negative
eddy viscosity can do the same where it makes the diffusivity of every equation
negative and the diffusion rate outweighs the convective rate and the bound of
the stiffness of the turbulence sources; a zero pressure at rest makes the step
of an inviscid run infinite, and a density so small that the diffusion rate
overflows makes the step zero.

A state that is not physical can also give a
positive finite step (a density and a pressure both negative in an inviscid
run): the check does not see it, and the run meets it at the first stage, in
the checks of the faces and of the cells.

**Depth of the map.**  The map of the adaptive step is built for
$H_{\hat L} = 2^{\hat L}\,\Delta t_{\min}/\sigma$ (at most $H_{\text{target}}$), with
$\Delta t_{\min} = \min_c \Delta t^{\text{cap}}_c$ and $\hat L$ the largest depth up to
$L_{\max}$ at which the cell with the largest admissible step,
$\Delta t^{\text{cap}}_{\max}$, keeps level 0 in the rule above (its own term,
without the interface margin below).  The test is monotone in $\hat L$ and holds for
$\hat L = 0$, so $\hat L$ is found from $L_{\max}$ down; without a target
$\hat L = \min(L_{\max}, \lfloor K + 10^{-12} \rfloor)$, with
$K = \log_2(\Delta t^{\text{cap}}_{\max}/\Delta t_{\min})$ the octaves spanned by
the local steps.

The cell with $\Delta t^{\text{cap}}_{\max}$ has level 0, so
$H_{\hat L} \le \Delta t^{\text{cap}}_{\max}/\sigma$ (to $10^{-12}$), and $H \le H_{\hat L}$: the
target bounds $H$ when it sets $H_{\hat L}$, otherwise the cell with $\Delta t_{\min}$
has level $\hat L$, which nesting and the connections do not raise, nor the margin while
$\kappa \le \sigma$ (the default $\sigma = 2$ with cfl $\le 1$), and
its term in $H$ is $H_{\hat L}$, so the depth of the map is $L = \hat L$.

The macro step
therefore does not exceed the admissible step of the slowest cells divided by
$\sigma$, and the rule leaves no coarse level empty; with $\kappa > \sigma$ the margin can raise that cell, and a coarse
level can stay empty.  Nesting, the margin and
the connections only raise levels, so in a strongly stretched mesh they can
raise every cell of the slowest octave as well; the run shows it as a first count
0 in its count of cells per level.

With a target at or below
$\Delta t^{\text{cap}}_{\max}/\sigma$ every depth passes the test and
$\hat L = L_{\max}$: the target already keeps $H$ within the step of the slowest
cells, and the map is the one built for $\min(2^{L_{\max}}\Delta t_{\min}/\sigma,
H_{\text{target}})$, whose depth $L$ is the one the target gives, at most $L_{\max}$.
On a uniform mesh ($K < 1$) the map has one level and $H =
\Delta t_{\min}/\sigma$.  The depth is that of the map of the adaptive step, also
the one that a fixed map keeps; where $H$ is given, with a fixed macro step
and a recomputed map and in the redo, which keeps
the $H$ of its attempt, the levels follow the cells up to $L_{\max}$: for a given
$H > 2^{\lfloor K\rfloor}\Delta t_{\min}/\sigma$ a cut at $\lfloor K \rfloor$
would take the cells with $\Delta t_{\min}$ beyond their limit. The depth
depends only on the state, as the map, and the largest admissible step is that
of all the MPI processes.

The default $L_{\max} = 10$ is a choice, not a bound of the scheme.  The ideal gain, the cell updates of the single-rate
step over those of the multirate step in the same simulated time, does not fall as the depth grows and stops growing at
$\lfloor K \rfloor$, while every level adds its stages, with their exchanges and synchronisations between the processes,
also when it holds few cells; and with $L_{\max}$ above 10 the floor of the redo is not within the step of a
single-rate run ([Non-local CFL condition](#non-local-cfl-condition)). Which depth costs least depends on how many
cells lie in each octave of the local steps and on the cost of a stage.

The map is computed at the first macro step and kept (fixed map) or
recomputed every $k$ macro steps (recomputed map), at the
iterations $1, 1 + k, \dots$ of the run counter.  With the redo a
rejected macro step rebuilds the map with levels that only rise
([Non-local CFL condition](#non-local-cfl-condition)); with $k > 1$ the raised
map stays in use until the next rebuild.  One iteration of
MOSE is one macro step.

**Which map for which run.**  For a fixed map the macro step is a one-step
method of order 2 in time: halving $H$ halves every $h_{\ell_c}$ and the error
drops by four.  If instead the map is rebuilt for each $H$, with every cell at
the coarsest level within its limit, halving $H$ lowers by one the level of
every cell finer than level 0, since every ratio that sets a level
($\sigma H/\Delta t^{\text{cap}}_c$, and $\kappa H/\Delta t^{\text{cap}}_v$ with
the interface margin) is proportional to $H$, so its step $h_{\ell_c}$ stays
the same: its error does not drop, and the measured order would be that of the
change of map, not of the scheme.

A study of the temporal order therefore uses
a fixed map: the map of the adaptive step of the first macro step, kept, also with
a fixed macro step.  A production run needs the opposite: as the state changes,
the map for the $H$ actually taken keeps every cell within its limit at the
coarsest admissible level, while a map kept from another $H$ either wastes
steps or exceeds the limits; so a recomputed map is built for the
$H$ taken, also with a fixed macro step.

A fixed map keeps
every cell within its limit only up to its own step, and it never rises: with an $H$
above that step the finest cells start the first macro step beyond $\Delta t^{\text{cap}}/\sigma$, and every later one
until their admissible step grows by the ratio of $H$ to the step of the map; the run reports it, as
with a recomputed map where $L_{\max}$ levels are not enough.

The macro step lands exactly on the output instants $t_k = k\,\Delta t_{out}$
($\Delta t_{out}$ the output interval, $k$ integer) and on the final time: when
$t_n + H$ reaches the next of them, or ends within $10^{-12} H$ of it, $H$
becomes $t_k - t_n$ and the time after the step is $t_k$ itself.  The map of
the step is kept: every $h_{\ell_c}$ only shrinks (up to the $10^{-12}$ of a
stretched step), so every cell stays within its limit, and the step is the
same scheme with a shorter $H$, conservative as any other.

A macro step whose
$H$ a redo halves does not land: a later macro step reaches the instant.  The
instants depend only on the time, and the map only on the state and on the run
counter, so a restart from a solution written with its counters goes on as the uninterrupted run (with $k > 1$, when
it falls on an iteration where the map is rebuilt).

## Non-local CFL condition

The CFL condition is not local: within a step a wave crosses the faces of a
cell and reaches its neighbours, which advance with their own steps.
Gnedin, Semenov and Kravtsov [13, section 2] state two requirements for a
local time stepping scheme:

- **(3)** the step of a cell satisfies its own admissible step and the
  admissible steps of all its neighbours;
- **(4)** if the admissible step of a cell at a time $t_0$ is $\Delta t$, all
  its immediate neighbours end their step at or before $t_0 + \Delta t$.

(3) is necessary but not sufficient [13, fig. 5]: when the admissible step of
a cell drops during the steps of its neighbours (a shock reaching it), a
neighbour that keeps its step delays the wave.  The source asks for (4) at any
time $t_0$ and imposes it in its code, where the tests show no delay; it gives
no proof that (4) is sufficient.

**(3) on the state at the start of the macro step.**  Let $c$ and $v$ share a
face in a block.  Nesting gives $\ell_c \ge \ell_v - 1$ and the levels give
$h_{\ell_v} \le \Delta t^{\text{cap}}_v/\sigma$, so

$$
h_{\ell_c} = \frac{H}{2^{\ell_c}} \le \frac{2H}{2^{\ell_v}} = 2\,h_{\ell_v}
\le \frac{2}{\sigma}\,\Delta t^{\text{cap}}_v ,
$$

and with $\sigma \ge 2$ every cell advances within the admissible step of each
of its face neighbours.  The two cells of a block connection have the same
level, so there $h_{\ell_c} = h_{\ell_p} \le \Delta t^{\text{cap}}_p/\sigma$.
Where the convective limit rules, the step of $c$ measured with the limit of
$v$ has a Courant number of at most $2\,\text{cfl}/\sigma$, which exceeds the
cfl of the limit of $v$ when $\sigma < 2$.  With the interface margin (below, always on with the redo) the levels and
$H$ also keep $h_{\ell_c} \le \Delta t^{\text{cap}}_v/\kappa$, $\kappa \ge 1.1$,
for every face neighbour $v$ in the block, so (3) holds there for every $\sigma \ge 1$.

So $\sigma \ge 2$ is what gives
(3) without the margin; $\sigma$ can take every value of at least 1.  Without the margin, a Sod tube with $\sigma = 1$ breaks down ([Multirate Verification](../vv/multirate.md#safety-factor)).  With
$L_{\max} = 0$ every cell advances with $H \le \min_c \Delta t^{\text{cap}}_c/\sigma$, (3) holds with $\sigma = 1$
and the macro step is the single-rate step.  The default is $\sigma = 2$, the smallest value that gives (3) without
the margin: the term $\sigma H/\Delta t^{\text{cap}}_c = 2^{L}\Delta t_{\min}/\Delta t^{\text{cap}}_c$ of the levels of
the adaptive step does not depend on $\sigma$, while $H$ is inversely proportional to it, so a larger $\sigma$ costs that
much more work per unit of simulated time.

With the margin, a $\sigma$ below $\kappa$ lets the margin set the levels of
every cell next to a neighbour with a step as small as its own, and a $\sigma$ near 1 makes the cell that sets $H$ start at
its limit, so that every fall of its admissible step within the macro step rejects the attempt (see the floor of
$\kappa$ below): more attempts, not a failure.

**Interface margin.**  On a face between levels the coarse cell $c$ next to
the finer cell $v$ advances with a Courant number, measured with the limit of
$v$, up to $2\,\text{cfl}/\sigma$.  In one dimension, for the advection
equation, the forward-Euler MUSCL step is positive and TVD for a Courant
number up to $1/(1+\theta)$: with
$u^-_{j+1/2} = \bar u_j + \tfrac12\phi(r_j)\,\Delta_+\bar u_j$ and
$r_j = \Delta_-\bar u_j/\Delta_+\bar u_j$,

$$
u^-_{j+1/2} - u^-_{j-1/2} = \Delta_-\bar u_j\,\Bigl[1 + \frac{\phi(r_j)}{2r_j}
- \frac{\phi(r_{j-1})}{2}\Bigr],
$$

so that, with $B_j$ the bracket, in $[0, 1+\theta]$, and $\nu\,(1+\theta) \le 1$,
the step $u^{n+1}_j = (1-\nu B_j)\,\bar u_j + \nu B_j\,\bar u_{j-1}$ is a
convex combination and
$u^{n+1}_j - u^{n+1}_{j-1} = (1-\nu B_j)\,\Delta_-\bar u_j + \nu B_{j-1}\,\Delta_-\bar u_{j-1}$
has non-negative coefficients, so the total variation does not grow (the
argument of Hundsdorfer [8], §6.4, p. 54, where the limiter condition (5.4) of his
§5.2, p. 39, bounds the bracket by $1+\mu$):
$\theta = 1$ for van Leer, MC and superbee ($\phi$ and $\phi/r$ in $[0, 2]$),
$\theta = 1/2$ for minmod ($[0, 1]$), $\theta = \sqrt2/2$ for van Albada
($[-0.207, 1.207]$), $\theta = 0$ with the first-order reconstruction.

These bounds hold on cells of equal width. The
reconstruction of MOSE limits the slopes divided by the
distances between the cell centres and bounds the limited slope so that the value of a cell at each of its faces lies between
the value of the cell and that of the neighbour across the face
([Spatial Discretization](numerics.md#muscl-reconstruction)), so on cells of different widths the face value is
$u^-_{j+1/2} = \bar u_j + \hat\phi_j\,w^+_j\,\Delta_+\bar u_j$, with $w^+_j = \Delta x_j/(\Delta x_j + \Delta x_{j+1})$,
$w^-_j = \Delta x_j/(\Delta x_{j-1} + \Delta x_j)$, $r'_j = r_j\,(\Delta x_j + \Delta x_{j+1})/(\Delta x_{j-1} + \Delta x_j)$ and
$\hat\phi_j = \min(\phi(r'_j),\ 1/w^+_j,\ r'_j/w^-_j)$ for $r'_j > 0$ and $\hat\phi_j = \phi(r'_j)$ for $r'_j \le 0$, which is 0
for every limiter but van Albada.

For a limiter whose $\phi$ and $\phi/r$ lie in $[0, 2\theta]$ (all but van Albada) the
bracket lies in $[1 - \min(2\theta w^+_{j-1}, 1),\ 1 + \min(2\theta w^-_j, 1)] \subset [0, 2]$ on every mesh (with van Albada
in $[w^-_j - 0.207,\ w^-_j + 1.207]$), and the step is a convex combination when $\nu_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$
([Appendix B](multirate-analysis.md#appendix-b-total-variation-at-a-face-between-levels), Lemma 1).  With one time step for all the cells this
follows from $\nu_j\,(1+\theta) \le 1$ on every cell.

Hoang
et al. [11, Theorem 4.2] prove a total-variation bound for a one-dimensional
scalar law, discretised by piecewise-linear discontinuous Galerkin elements
with the second-order strong-stability-preserving Runge–Kutta scheme in
two regions joined by one interface cell, with their limiter condition (4.14)
on the reconstructed values, for a solution regular at the interface and under
the hypothesis (4.15) that the coarse step meets this bound on the first fine
cell.  With the interface margin, always on with the redo, the levels and $H$ also satisfy
$h_{\ell_c} \le \Delta t^{\text{cap}}_v/\kappa$ for every face neighbour $v$ in
the block, $\kappa = \max(\text{cfl}\,(1+\theta),\ \kappa_0)$, $\kappa_0 = 1.1$, with the cfl
in force (when cfl changes during the run, $\kappa$ follows it from the next macro step): the level of a cell comes from
$\max(\sigma H/\Delta t^{\text{cap}}_c,\ \kappa H/\min_v \Delta t^{\text{cap}}_v)$
and $H_{\text{cap}}$ takes the term $\min_v \Delta t^{\text{cap}}_v\,2^{\ell_c}/\kappa$.

Here $\theta$ is that of the limiter: 1 for van Leer, MC and superbee, 1/2 for minmod, $\sqrt2/2$ for van Albada, and 0
for first order.
With $\sigma \ge \kappa$ this adds nothing on the faces between equal levels or
towards coarser cells, and nesting and equal levels across connections only
raise levels.  It raises a ring of cells on the coarse side of a face between
levels where the fine cell runs above $1/\kappa$ of its limit.

The coarse cells run with $\text{cfl}/\sigma$, within $1/(1+\theta)$ when
$\text{cfl} \le \sigma/(1+\theta)$: with the default $\sigma = 2$ and $\theta \le 1$, $\sigma/(1+\theta) \ge 1$, so every
cfl up to 1 meets it for every limiter (with $L_{\max} > 0$ and the redo, cfl is bounded by 1: the bounds of cfl and
vnn, below); a $\sigma$ below 2 meets it only for $\text{cfl} \le \sigma/(1+\theta)$, a property and not a constraint
that MOSE imposes.

This is the hypothesis of the theorem, not the theorem: the scheme is
a finite-volume MUSCL scheme, in several dimensions, for systems, with several
levels.  At the faces between levels the monotonicity is proven for linear
two-cell fluxes, the first-order upwind flux of a given velocity field without
discrete divergence and the central diffusion
([Appendix A](multirate-analysis.md#appendix-a-order-and-stability-of-the-multirate-macro-step),
Theorem 2).  For the advection equation in one dimension, with the limited reconstruction of MOSE, a limiter
whose $\phi$ and $\phi/r$ lie in $[0, 2\theta]$ (van Albada, whose $\phi$ is negative for $-1 < r < 0$, is outside) and two
levels, on cells of any widths, one macro step raises the total variation at most by a term local to the face between levels,
a multiple of the Courant number of the coarse cell of the face times the differences of the solution next to the face.

Summed
over a run, this bounds the total variation of a solution that is Lipschitz
continuous next to the face, as in Hoang et al. [11, Theorem 4.2], without their hypothesis (4.15). The derivation asks every
cell to satisfy $c_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$, with $c_j$ its own Courant number.  With
$\text{cfl}/\sigma \le 1/2$ (cfl $\le 1$ and the default $\sigma = 2$) every cell has $c_j \le 1/2$, and with it this condition, on every mesh
and for every such limiter: also, with minmod, at the coarse cell downstream of a face from a finer to a coarser cell, where
it reads $\lambda\,(1 + r/(r+1)) \le 1$ with $r$ the ratio of the coarse to the fine width, which $c\,(1+\theta) \le 1$ alone
does not give ([Appendix B](multirate-analysis.md#appendix-b-total-variation-at-a-face-between-levels)).

Several levels, several dimensions, systems and the
positivity at the faces between levels are outside this result: for them, and for the limited reconstruction outside it, the
stability at the faces between levels rests on the tests of [Multirate Verification](../vv/multirate.md).

**The floor of $\kappa$.**  $\text{cfl}\,(1+\theta)$ is the hypothesis above; $\kappa_0$ is the margin of condition (4)
below.  With $\kappa \le 1$ and $\sigma = 2$ the margin would raise no cell above the level that nesting gives it: the cell that sets $H$, next to a coarser cell,
would start with $\hat h_v = \Delta t^{\text{cap}}_v$ exactly, any fall of its admissible step at a later start of a
sub-step shared with that cell would reject the attempt, and the next macro step, whose map comes from the state alone, would
meet the same limit.

With $1 < \kappa \le \sigma$ (the default $\sigma = 2$ with cfl $\le 1$) a rejection for a violation of
(4) needs a fall of the admissible step of a cell of more
than $1 - 1/\kappa$ within the macro step, on every mesh and state, for every limiter and number of dimensions and of
processes ([Appendix A, Proposition](multirate-analysis.md#proposition-the-margin-on-condition-4)); a front is rejected only when
it makes the step fall by more than that, and then only where the pair of cells runs that close to its limit: a fall by more
than $1 - 1/\kappa$ is necessary for a rejection, not sufficient.

That proposition does not fix the value.  $\kappa_0 = 1.1$ is a practice taken from Castro, whose
check of the step tolerates the same quantity, the fall of the admissible step within the step, up to $1/1.1$
[16]; it is not a quantity that the scheme gives.  PETSc and DOPRI5 take a step of 0.9 of the one they estimate on the same
principle, that a rejection be rare [14, 17], with another quantity, the error of the step.

Castro accepts a step
above its limit by up to that factor; here the margin is built into the map, and condition (4) stays imposed as stated.  The
margin raises a cell only where the admissible step jumps by more than $\sigma/\kappa$ across one of its faces, and
the rule of a cell reads the local steps of its neighbours, not their levels, so it does not propagate beyond those
neighbours.  How many cells it raises depends on the mesh and on the state.

Below
$\text{cfl}\,(1+\theta) = 1.1$ the margin does not depend on cfl; a $\sigma$ above 2 gives more
margin at the price of a shorter $H$.

**The certificate of (4).**  A state exists only at the start of the
sub-steps, so (4) is checked there: for every cell $v$ and every start $\tau$
of one of its sub-steps,

$$
\hat h_v(\tau) \le \Delta t^{\text{cap}}_v(\tau), \qquad
\hat h_v = 2\,h_{\ell_v} \text{ if a face neighbour one level coarser starts a
sub-step at } \tau, \quad \hat h_v = h_{\ell_v} \text{ otherwise}
$$

(finer neighbours start and end within the sub-step of $v$, nesting limits the
jumps to one, connections have none).  At $t_n$ this is (3).  The certificate
evaluates $\Delta t^{\text{cap}}_v(\tau)$ with the local-step formula of the
single-rate step, cell by cell, in the first stage of every sub-step after the
turbulence sources, on the start state of the sub-step; a value that is not a
number is a violation, and $\hat h_v$ above $\Delta t^{\text{cap}}_v$ by less than
$10^{-12}$ of it is not, the tolerance of the formula of the levels: a cell that
starts at its limit gives a ratio of 1 within some units in the last place, which
no map and no $H$ would change.

Two cells that share only an edge and are two levels
apart satisfy neither (3) nor (4) ($h_c = 4h_d \le 2\Delta t^{\text{cap}}_d$
with $\sigma = 2$): a wave reaches a cell only through its faces, along which
(4) holds pair by pair.

- Not checked: (4) is not imposed; the run measures $r_{post}$ and $r_{nb}$ (below).
- Counted: the certificate counts the violations per level and the macro steps
  that would be rejected; the results do not change.
- Redo (the default): a macro step whose certificate fails is rejected.  The state goes back
  to $t_n$, the map is rebuilt with the same $H$ from the smallest admissible
  step that every cell showed in the attempt, so the levels only rise, up to
  $L_{\max}$ (the map can become deeper than that of the first attempt), and $H$
  is halved while a cell would need a level above $L_{\max}$ or a value was not a
  number.  An accepted macro step satisfies (4) at every start of a sub-step.
  With $L > 0$ an attempt is rejected also when a cell or a face takes a state
  that is not physical: a cell whose update gives a conservative variable that
  is not finite or a total density that is not positive (tested before the
  conversion to primitive variables, which divides by the density), or a
  primitive variable that is not finite, a negative pressure or a temperature
  outside the tables of the gas (after the floor of the partial densities, the
  test of the faces); an internal face of a block or a face between levels
  whose reconstructed left or right state has a negative pressure or a
  temperature outside the tables (the faces of block connections and of
  boundaries have no such test, as in the single-rate step).  The cell keeps the
  start state of its sub-step and the face gives no flux; from the next stage
  on, the process that took the state computes no internal face and updates no
  cell (its exchanges go on, so the processes stay together), and the attempt is
  rejected on every process.  The state goes back to $t_n$ and the attempt is
  redone on the same map and local steps with $H$ halved once: the admissible
  steps that the certificate sees after such a state come from states that are
  not the solution, so they do not lower the local steps.  A face that is not
  physical at the first stage of a level at $t_n$ sees only the state at $t_n$,
  which a shorter $H$ does not change: there the run stops with the lines of the
  face, as the single-rate step does.  In a macro step accepted after attempts
  that the redo can redo (a recomputed map, an adaptive macro step, $L > 0$)
  every cell ends with finite values, partial densities not below the floor, a
  pressure that is not negative and a temperature within the tables.
  The attempts end: the levels are bounded by $L_{\max}$, and as $H$ halves the states
  of the step tend to those at $t_n$, where (3) holds with the margins (the
  local step is continuous in the state where density and pressure are
  positive); every rejection for a state that is not physical halves $H$, so a
  macro step has at most 11 of them; below $2^{-10}$ of the first $H$ the run
  stops.  With a map of depth $L \le 10$, the default cap, the
  first $H$ is at most $2^{L}\Delta t_{\min}/\sigma$, so this floor is at most
  $\Delta t_{\min}/\sigma$: every cell then advances with at most the step of a single-rate
  run with the same margin.  Their number is bounded too, by a bound
  that grows with the mesh: $H$ takes at most 11 values, and with the same $H$
  the local steps only decrease, so the levels of the rebuilt maps only rise,
  and every rejection that keeps $H$, but the forced one of the test option
  below, raises a cell by at least one level; over the $N$ cells of the mesh a
  macro step has at most $1 + 11\,(N L_{\max} + 2)$ attempts (the term 1 is the
  forced rejection).  $H$ is halved, the factor of the retry of Castro when it has no suggested
  step [16], and not set to $H_{\text{cap}}$, so that it takes these 11 values, on which the
  bound rests; the macro step can then come out up to twice shorter than $H_{\text{cap}}$.
  A limit on the attempts (21 by default) bounds them: when that many
  attempts of a macro step are rejected (the forced rejection is not counted),
  the run stops.  The floor and the limit are rules that stop a run, not
  conditions of the scheme: the bound above ends the attempts
  without them.  The floor is 10 halvings of $H$, as many rejections as PETSc
  allows a step by default [14] and as many failures of the nonlinear solver or
  of the inequality constraints as end a step of the ARKODE integrators of
  SUNDIALS [15]; there every rejection shortens the step.  21 is one attempt for
  each of the 11 values of $H$ that the floor allows and 10 more for the
  rejections that keep $H$ and raise levels, which those integrators do not
  have.  The limit is checked after the floor, so with a limit of 11 or more a
  macro step whose rejections all halve $H$ ends as without it.
  With a fixed macro step or a fixed map the certificate only counts: there
  (4) is not imposed, and, as with one level, the stops of the single-rate step
  stay: a cell whose update gives a value that is not finite or a negative
  pressure, or a face that is not physical, stops the run.  With $k > 1$ the map raised by a redo
  stays in use until the next rebuild of the counter.  With $L = 0$ there is one
  level and the check is the own limit at $t_n$, true by construction.  A test
  option rejects the first attempt of every macro step that the certificate
  accepts and redoes it with the same $H$ and map: the state goes back to $t_n$
  to the bit, so the results do not change.

**The bounds of cfl and vnn.**  At $t_n$ the map runs every cell with $h \le \Delta t^{\text{cap}}/\sigma$; with
the redo the certificate accepts, at every later start $\tau$ of a sub-step, $\hat h_v$ up to
$\Delta t^{\text{cap}}_v(\tau)$ itself, without $\sigma$ (to the $10^{-12}$ of its test), so within the macro step a cell
can advance with the whole step that cfl and vnn admit on its state at $\tau$.  On a uniform mesh of equal cells, for the
frozen-coefficient advection–diffusion model with the first-order upwind flux and the central diffusion, the
checkerboard datum $u_j = (-1)^j$ gives a zero limited slope with every reconstruction of MOSE (first order, minmod,
superbee, van Leer, van Albada, MC), so on it every limited scheme is the first-order scheme, and the datum is an
eigenvector of Heun with the factor $R(z) = 1 + z + z^2/2$, $z = -2\,\text{cfl}\,a - 4\,\text{vnn}\,b$, where $a$ and $b$
($a + b = 1$) are the convective and diffusive shares of the local step of a cell at its limit; $|R(z)| > 1$ exactly when
$\text{cfl}\,a + 2\,\text{vnn}\,b > 1$.

Hence a sub-step with the whole admissible step is linearly stable on every cell
only if cfl $\le 1$ and vnn $\le 1/2$ (pure convection: $R(-2\,\text{cfl}) = 1 - 2\,\text{cfl} + 2\,\text{cfl}^2 > 1$ for
cfl $> 1$; pure diffusion: $R(-4\,\text{vnn}) = 1 - 4\,\text{vnn} + 8\,\text{vnn}^2 > 1$ for vnn $> 1/2$), for every
limiter; with a larger value it is not, on a cell where that number sets the step (with cfl $= 1$ and vnn $> 1/2$ on every
cell with $b > 0$; with cfl $< 1$ where $b > (1 - \text{cfl})/(2\,\text{vnn} - \text{cfl})$).

These bounds are also
sufficient for the first-order reconstruction and for the limiters with a linearisation (van Leer, MC, van Albada)
([Local time step](time-integration.md#local-time-step)), and are the bound of the other slopes for minmod and superbee, for which no linear
limit is proven: their one-sided linearisations are the fully upwind slope, whose Heun limit on the checkerboard mode is
1/2, and the central slope, which Heun does not integrate stably at any Courant number; what is proven for them is the
bound of the total variation of the forward-Euler step, 2/3 for minmod and 1/2 for superbee, on equal cells in one
dimension (interface margin, above).

In two and three dimensions, with the additive rule of the local step, the bounds
are those of the checkerboard in every direction and are the same numbers (on cells skewed in all three pairs of
directions the viscous rate takes a factor that raises the diffusive threshold on those cells).  For the Euler and
Navier–Stokes equations linearised about a uniform state the convective bound is attained by the fastest wave with the
Roe, Rusanov and HLL-type fluxes in one dimension, in several dimensions with Rusanov on every state and with Roe and HLLC
on a state at rest; the diffusive bound is attained by the scalar that owns the largest diffusivity (the energy with
gases of Prandtl number near 0.7, the species, the turbulence variables) in the limit of a cell where the diffusion sets
the step.

The bound vnn = 1/2 is proven on a uniform mesh of equal cells
of any skew and aspect ratio (in 3-D with the factor that the local step takes on cells skewed in all three pairs of
directions), with frozen coefficients and away from boundaries, for the scalars exactly and for the momentum
conservatively ([Local time step](time-integration.md#local-time-step)); on cells of varying size or shape, next to boundaries and at the
axis of the wedge the local-step rule applies the same bound cell by cell, and its support there is the evidence of
that section.

So with $L_{\max} > 0$ and the redo, where the certificate redoes the attempts
(a recomputed map and an adaptive macro step), MOSE refuses a cfl above 1 and, in Navier–Stokes runs, a vnn above
1/2: these bounds are necessary for every sub-step that the certificate accepts to be linearly stable there.  A larger value would not show at
$t_n$, where the cells run within $1/\sigma$ of their step.  That a run with a larger value fails is not proven: an
unstable sub-step occurs only where the admissible step of a cell falls within the macro step by more than
$\sigma/\text{cfl}$ or $\sigma/(2\,\text{vnn})$, for the sub-steps of that macro step that remain.

When (4) is not
checked or only counted, and with the redo and a fixed map or a fixed macro step, where the certificate only counts,
nothing bounds the sub-steps after the start of the macro step, where every
cell runs within cfl/$\sigma$ and vnn/$\sigma$, and neither bound is imposed.  The two conditions are distinct: at
the start of the macro step every cell runs within $\Delta t^{\text{cap}}/\sigma$, so every cell, the coarse cells of a
face between levels included, has on its own width a Courant number of at most $\text{cfl}/\sigma$, 1/2 with cfl $\le 1$
and the default $\sigma = 2$, within the bound $1/(1+\theta)$ of the total variation above for every limiter (measured
with the limit of the finer neighbour, the coarse step reaches $2\,\text{cfl}/\sigma$; the interface margin keeps it
within $1/(1+\theta)$); within the macro step, with the redo, the certificate keeps
every sub-step within the linear stability limit of Heun, on the state where the sub-step starts.

For a cell next to a wall where the diffusion sets the step, with 1 µs its step
at the limit of Heun (a diffusion number $h\,D_{\max}\sum_d|\nabla\xi_d|^2$ of 1/2 on its state):

| | step at the start | diffusion number at the start | the certificate rejects above a diffusion number of |
|---|---|---|---|
| single rate, vnn 0.5 | 1 µs | 0.5, the limit | — (the step is computed again at every step) |
| multirate, vnn 0.5 | 0.5 µs | 0.25 | 0.5, the limit |
| multirate, vnn 0.98, which the redo does not allow | 0.98 µs | 0.49 | 0.98: between 0.5 and 0.98 an accepted sub-step is not linearly stable |

For the convection, cfl $\le 1$ puts the threshold of the certificate at the linear limit of Heun in the same way.  With
$L_{\max} = 0$ the macro step is the single-rate step, whose limits are not imposed
([Local time step](time-integration.md#local-time-step)).

**Measures without the certificate.**  When (4) is not checked the levels and
$H$ come from the state at the start of the macro step and stay fixed during
it: if the admissible step of a cell drops within the macro step, its
neighbours keep their steps.  The run measures (4) on the state reached at
$t_{n+1} = t_n + H$, whose admissible steps the next macro step computes
anyway:

$$
r_{post} = \max_c \frac{h_{\ell_c}}{\Delta t^{\text{cap}}_c(t_{n+1})},
\qquad
r_{nb} = \max_c\, \max_{v \in N(c)} \frac{h_{\ell_c}}{\Delta t^{\text{cap}}_v(t_{n+1})},
$$

with $N(c)$ the face neighbours of $c$ in its block; across a block
connection $h_{\ell_c} = h_{\ell_p}$, so the ratio with the partner $p$ is the
$r_{post}$ of $p$.  The last sub-step of $c$ ends at $t_{n+1}$: if the
admissible step of $v$ does not fall below $\Delta t^{\text{cap}}_v(t_{n+1})$
during that sub-step, $r_{nb} \le 1$ gives (4) for it, and $r_{post} \le 1$
gives the own limit of $c$; a value above 1 flags a possible violation.

The
measure is partial: it sees only the state at the end of the macro step (not
the earlier sub-steps) and only the face neighbours in the block (not the
cells that share an edge or a vertex).

## Recursive Heun scheme

With $\mathbf{U}^0$ the state of a cell at the start of its sub-step, the
macro step is the recursion below.  It is the *Prop. RK2* scheme of
Jeanmasson, Mary and Mieussens [10] with two differences: the coarse corrector
reads the final values of the finer cells (their Table 16 reads the predictor
of the fine cells), and the recursion over several levels, with bands of one
cell, is ours (they write the scheme for two levels).  Both keep the
conservation, which needs only that the two sides of a face receive the same
flux with the same weights.

Both differences keep the second order on a fixed mesh
([Appendix A](multirate-analysis.md#appendix-a-order-and-stability-of-the-multirate-macro-step),
Theorem 1).  Every stage vector is first-order consistent on the stencil it
is read from (the coarse cells through the interpolant, the finer cells at
their final values, the cells of the level at the start state or at the
predictor), and every cell receives each face flux with the weights of Heun,
sum $H$ and first moment $H^2/2$, also on the corrected side of a face between
levels, where the coarse cell receives $h_m A = \tfrac{h_m}{4}\sum_j G_j$ at
the four fine stage times: the two differences enter only the $O(H^3)$ term,
so the macro step is a one-step method of order 2, uniformly at the faces
between levels, for every linear semi-discretisation and for $f \in C^2$ (the
limiter and the switches of the turbulence models break $C^2$ where they
act).

This is what the order tests of [Multirate Verification](../vv/multirate.md) measure: fixed mesh,
fixed map, $H \to 0$.  At a fixed Courant number ($H \propto \Delta x$) the
local error of a macro step is $O(H^2)$ in the cells next to a face between
levels and $O(H^3)$ elsewhere (first order locally, as the conventional
boundary data of [7]; Appendix A, Remark 2).  When the map is rebuilt
(a recomputed map) and a face between levels travels with a wave at the
speed of the wave, the errors of order $H^2$ that each macro step leaves next
to the face stay on the same part of the wave and can add up.
[Appendix A, Corollary](multirate-analysis.md#corollary-every-sequence-of-maps) bounds the error
for every sequence of nested maps, this one included: on a fixed mesh by
$O(H^2)$, under (H2)–(H4) of Theorem 1; at a fixed Courant number, for the
linear two-cell fluxes of Theorem 2, by $O(H)$ in the maximum norm and by
$O(H^2)$ in the volume-weighted $L^1$ norm (with faces between levels of
bounded total area).

Outside that linear setting no
bound is proven at a fixed Courant number: the moving-levels test of
[Multirate Verification](../vv/multirate.md) measures the $L^1$ order of this
case on a simple wave of the Euler equations (second order); its global
$L^\infty$ order, near 1.3 as for the single-rate step with
the same limiter, has its largest error at the maximum of the wave and does
not isolate the cells next to the moving faces, where the test reports,
without deciding on it, a lower $L^\infty$ order (1.67 on its last pair,
against 2.02 to 2.03 for the single-rate step on the same cells).

The
stability of the two differences is proven for the first-order upwind flux
of a velocity field without discrete divergence and the central diffusion:
the macro step is monotone (non-negative matrix,
maximum principle, conservation, norms 1, 2 and $\infty$ at most 1) for every
nested map and every number of levels when every cell advances within 0.86 of
its forward-Euler positivity step in the setting of Theorem 2 (linear two-cell fluxes, Appendix A): on a mesh of
equal rectangular cells with constant coefficients, where the diffusive flux of a scalar is the two-cell flux of
Theorem 2, the additive local-step rule
starts every cell within $\max(\text{cfl}, 2\,\text{vnn})/\sigma$ of that step, 0.4 with $\sigma = 2$, cfl 0.8 and
vnn 0.4, and at most 1/2 with cfl $\le 1$, vnn $\le 1/2$ and the default $\sigma = 2$.

The bounds that the redo
imposes are necessary conditions of linear stability
([Non-local CFL condition](#non-local-cfl-condition)), not of positivity: on skewed cells the face gradient also takes
the central differences of the two cells along the face, whose cross metric terms can give corner
neighbours a negative weight (on equal parallelograms they do, so that no step makes that diffusion positive), and
Theorem 2 does not cover it.  With a limited
reconstruction the coarse corrector is not a convex combination of
forward-Euler steps (Appendix A, Remark 4): for the advection equation in one
dimension, with two levels and with three levels and a band of one cell, Appendix B bounds the growth of the total variation
by terms local to the faces between levels, and elsewhere the stability at the
faces between levels rests on the tests of
[Multirate Verification](../vv/multirate.md).

```text
ADVANCE(l, tau):  U0 := U on the cells of level l
                  STAGE(l, 1, tau)                       predictor
                  if l < L:  ADVANCE(l+1, tau);  ADVANCE(l+1, tau + h_l/2)
                  STAGE(l, 2, tau + h_l)                 corrector
STAGE(l, k, t):   coarser cells read by the stage <- interpolant at time t
                  ghost cells read by the stage at time t
                  viscous runs: viscous variables of the cells and ghost cells read by the stage
                  fluxes of the faces with a cell of level l (direction by direction, the convective
                  then the diffusive faces, as the single-rate step), boundary fluxes
                  turbulence models: turbulence sources of the cells of level l
                  U := U0 + c_k (U - U0 + h_l L(U)),  c = (1, 1/2),  on the cells of level l only
```

A macro step with all levels $0\dots L$ populated has $2^{L+2}-2$ stages.
When a stage of level $\ell$ runs, every finer cell is already at the stage
time (finer levels complete their two sub-steps before the coarser
corrector), and every coarser cell $m<\ell$ that the stage reads receives the
**linear interpolant in conservative variables** between its start state and
its predictor (Hoang, Ju, Leng and Wang [11], eq. (3.2), with the predictor of
their App.

B.1, eq. (B.4): the full Heun stage of the coarse cell with its true
neighbours, as stage 1 above, not the step 1 of their section 3, which uses the
coarse cell alone; the same values as [10], Tables 13 and 16):

$$
\mathbf{U}_c(t) = \mathbf{U}^0_c + \eta\,(\mathbf{U}^{(1)}_c - \mathbf{U}^0_c),
\qquad \eta = \frac{t - \tau_m}{h_m},
$$

with exact copies at $\eta = 0$ and $\eta = 1$.  Stage times are dyadic
fractions of $H$ and are handled as integers, so $\eta$ is exact.  Boundary
conditions are evaluated at the stage time.

The halo can hold cells two or more levels coarser (bands of one cell, below).
The interpolant has an error $O(h_m^2)$ for every $\eta$ in $[0, 1]$, and the
equations (B.3)–(B.6) of [11] hold for any integer ratio of the steps, so these
cells are read with the same local order.  They exchange no flux with the
cells that read them: only face fluxes enter the balance of a cell, and with
the nesting and the equal levels across connections a face joins levels at
most one apart, so the correction below keeps the conservation.

For linear
two-cell fluxes (first-order upwind, central diffusion) Theorem 2 of
Appendix A covers them (every nested map, bands of one cell included); with the limited reconstruction
Appendix B bounds the total variation for the advection equation in one dimension with three levels and a band of one
cell (Theorem 3b), and their effect on the stability is tested in
[Multirate Verification](../vv/multirate.md) (bands of one cell, 2000 macro steps).

The coarser cells read by a stage (the **halo**) are the union of the true
stencils of the routines the stage calls: the four cells along the axis of
every convective face and the cells read through the ghost cells of the
boundary records.  With the laminar Navier–Stokes equations the halo also
holds the ten cells of every diffusive face (the two cells of the face and
their neighbours $\pm 1$ in the two other directions, diagonals included)
and the four tangential neighbours of the partner cell of every block
connection of the level, which the viscous connection flux reads through the
tangential ghost values refreshed at every stage after the ghost cells.

The
viscous fluxes take the viscous variables of these cells (mass fractions,
velocity, temperature and turbulence variables per unit mass) from an array
that the single-rate stage fills once for every cell after the ghost cells; a
stage of level $\ell$ fills it, after the interpolation and the ghost cells, on
the cells of the level, their face neighbours and their diagonal neighbours in
every pair of directions, ghost cells across one face included: the coarser
cells at the stage time, the finer ones at their current values, so that no
flux reads a value computed from an older state.

The
turbulence sources of a cell read its neighbours $\pm 1$ along the axes and its
own ghost cells, already in that union.  Bands of one cell are therefore
allowed.

## Conservative correction

On a face between a coarse cell (level $m$) and a fine cell (level $m+1$) the
two sides would integrate different fluxes.  The coarse side is corrected so
that it receives exactly the impulse of the fine side (the refluxing principle
of Berger and Colella [12], eq. (2), stated there for a one-step scheme; the
form for the two Heun stages is derived here, with the algebra of [11],
eq. (3.5)): at its stage 1 the coarse cell uses and saves its own flux $G_s$; each
of the four fine stages adds $\tfrac14$ of its flux to an accumulator $A$
(weight $b\,h_{m+1}/h_m = \tfrac12\cdot\tfrac12$); at its stage 2 the coarse cell
uses $2A - G_s$.

Over its step the coarse cell then receives
$h_m A = \sum_{\text{fine stages}} \tfrac12 h_{m+1} G$, the same as the fine cell,
and mass, momentum and energy are conserved to round-off.  The convective and
the diffusive flux of a face are corrected each with its own accumulator, so
every cell receives its contributions in the order of the single-rate step.

## Blocks and boundary conditions

The level map spans all the blocks.  Because the two sides of a block
connection (direct and periodic) always have the same level, a connection
is never a face between levels: at the stages of that level the record that
owns the face computes its one flux for both cells, as the single-rate step
does, from ghost cells that every stage of the level refreshes from the
partner block at the stage time.  The cells of
the partner block read through those ghost cells are part of the halo, so a
coarser one takes its interpolant.

Equal levels cost a few cells: a
connection on the edge between two bands raises the cells of the coarser side
next to it, and the nesting rule follows.

Every boundary record is evaluated at the stages of the level of its cell,
at the stage time $t_n + \tau$ of that level.  The time-dependent data of the
tabulated inlets and of the assigned boundary states are therefore read at the stage times of each level,
the conventional imposition of [7] applied level by level; with one level the
times coincide with those of the single-rate Heun step (section
[Stage times](time-integration.md#stage-times-and-time-dependent-boundary-data)).

## Turbulence models

With the Spalart–Allmaras or SST model the stage of a level ends, as the
single-rate stage, with the turbulence sources of its cells, evaluated on the
state of the stage; they read the ghost cells of the fluxes of the stage (at
the walls those of the wall values, one set for the fluxes and the sources).  In
time-accurate runs the sources are explicit and the admissible local step
includes a bound $\hat d_c$ of the rate of their destruction, $D_q/q$, on the state of the cell at the start of the
macro step alone (for the Spalart–Allmaras model twice the largest value of that rate, for SST 1.84 times the rate of
$k$ and at least twice that of $\omega$) ([turbulence sources](turbulence-numerics.md#time-accurate-runs-explicit-sources)):
$1/\Delta t^{\text{cap}}_c$ gains $\hat d_c$.

The levels and the macro step
inherit it, so every cell advances with
$h_{\ell_c}\,D_q/q \le h_{\ell_c}\,\hat d_c \le 1/\sigma$, and each stage of a level is a forward-Euler
step on its cells that keeps $\tilde\nu$, $k$ and $\omega$ positive against
their destruction (the coarser cells read by the stage take a convex
combination of two states, the interpolant with $0 \le \eta \le 1$).  The
correction on the faces between levels is not a convex combination, so the
positivity there rests on the step margin and is checked by the tests. With one level the macro step is the
single-rate step to the bit.

## Several processes (MPI)

With MPI every block belongs to one process.  The level map is the one of a
single process: at every sweep of the two rules the levels of the partner
cells of the connections between processes are exchanged, and the sweeps stop
when no process changes its map; $H$ and $L$ are global.  Every process runs
the same recursion, so the stages of a level are simultaneous everywhere.

A stage of level $\ell$ exchanges only what it reads from other processes:
the ghost cells of the connection records whose partner block is on another
process and whose ghost positions the stage reads, in viscous runs the
tangential ghost values of the connection records of its cells, and the share
of the flux of the connection faces whose owner record is on the other
process.  The process
of the partner block sends them at the stage time, as the single-rate
exchange does: cells of level $\ell$ or finer as they are, coarser cells with
their interpolant, which that process evaluates as part of its own halo of
level $\ell$.

The two processes agree on these lists, record by record and
level by level, every time the lists are built (first the tangential values,
whose sources the partner's process may read from a third process, then the
ghost cells).  Because both sides of a connection have the same level, no
face between levels crosses a connection and the conservative correction
stays within a process.

## Scope and limits of the analysis

The integrator advances time-accurate runs with Heun's scheme on the conservative variables, for the Euler and the
Navier–Stokes equations, laminar or with the Spalart–Allmaras or SST models. The configurations it runs, and those that stop at set-up,
are listed in the [User Guide](../user/input.md#scope).

The production of the rotation-corrected Spalart–Allmaras model and the dilatation part of the compressible SST production can
be a sink outside the bound of the stiffness that the step takes
([turbulence sources](turbulence-numerics.md#time-accurate-runs-explicit-sources)), as in the single-rate step. No multirate test
measures their accuracy.

What is proven, and where:

- **Stability** is proven only in the linear setting of Theorem 2 of
  [Appendix A](multirate-analysis.md#appendix-a-order-and-stability-of-the-multirate-macro-step) (first-order upwind flux,
  central diffusion) and, for the limited reconstruction, in the setting of
  [Appendix B](multirate-analysis.md#appendix-b-total-variation-at-a-face-between-levels): the advection equation in one
  dimension, two levels, and three levels with a band of one cell, a bound of the total variation.
- **The second order in time** is proven on a fixed mesh (Appendix A, Theorem 1, and for every sequence of maps its
  Corollary (a)). With the mesh and the macro step refined together at a fixed Courant number, the local error next to a
  face between levels is of first order (Appendix A, Remark 2). It is bounded only for linear two-cell fluxes (Corollary
  (b): $O(H)$ in the maximum norm and $O(H^2)$ in the volume-weighted $L^1$ norm); outside them it is measured
  ([Recursive Heun scheme](#recursive-heun-scheme)).

Outside the proofs are the limited reconstruction in several dimensions, more than three levels or three levels in another
arrangement, HLLC on the Euler and Navier–Stokes equations, the turbulence sources, the terms of the axisymmetric wedge (faces
`200`) and shocks that cross faces between levels. For them the stability rests on the tests of
[Multirate Verification](../vv/multirate.md).

!!! note "Stability of Heun with a limiter"
    Every cell advances within its own local step divided by $\sigma$.  In
    one dimension, for the advection equation with a limited reconstruction,
    the forward-Euler step is a convex combination of the values of the cell
    and of its upwind neighbour, so it keeps a local maximum principle and the
    positivity, for a Courant number up to $1/(1+\theta)$, where $\theta$
    bounds the face increment of the reconstruction relative to the upwind
    difference: $\theta = 1$ for the limiters whose slope can reach twice the
    smaller difference (van Leer, MC, superbee) and $\theta = 1/2$ for minmod.
    On cells of different widths the bounds of the slopes of the reconstruction keep the property on every mesh, for
    $\nu_j\,(1 + \min(2\theta w^-_j, 1)) \le 1$, which $\nu_j\,(1+\theta) \le 1$ gives with one time step for all the cells
    (see the interface margin in [Non-local CFL condition](#non-local-cfl-condition)).
    This is the bound $(1+\mu)\,\nu \le 1$ of Hundsdorfer [8] (§6.4, p. 54, with
    $\theta$ the $\mu$ of the limiter condition (5.4) of his §5.2, p. 39),
    stated there for the explicit trapezoidal rule, which is Heun: Heun keeps
    the property of its forward-Euler stages (Hundsdorfer [8], §6.2, Theorem 6.4
    and the example after it, pp. 51–52; Shu and Osher [1], Proposition 2.1 and Remark 2.1).  In several
    dimensions no second-order scheme is TVD
    ([Local time step](time-integration.md#local-time-step)): the forward-Euler stage is a convex
    combination of one-dimensional steps with the Courant number
    $\sum_d \nu_d$ (sum of the directional Courant numbers), so
    $\sum_d \nu_d \le 1/(1+\theta)$ keeps a local maximum principle.  The
    local step adds up the directions, so a cell at its limit has
    $\sum_d \nu_d \le \text{CFL}/\sigma$, 0.4 with $\sigma = 2$ and CFL 0.8, on
    every cell, isotropic ones included.  On a face between levels the coarse
    step measured with the limit of the fine cell can reach
    $2\,\text{CFL}/\sigma$, above $1/(1+\theta)$, without the interface
    margin; with it (always on with the redo, the default) it stays
    within $1/(1+\theta)$ on the state at the start of the macro step, where the
    levels are built ([Non-local CFL condition](#non-local-cfl-condition)),
    the bound that Hoang et al. [11] (Theorem 4.2, condition (4.15)) assume on
    the first fine cell.  Their theorem, for a one-dimensional scalar law
    discretised by piecewise-linear discontinuous Galerkin elements with the
    second-order strong-stability-preserving Runge–Kutta scheme in two
    regions, with a solution regular at the interface, is not a foundation of
    the stability at the faces between levels of a finite-volume scheme for
    systems, in several dimensions and with several levels: there the
    monotonicity is proven for linear two-cell fluxes (first-order upwind,
    central diffusion), for every cell within 0.86 of its forward-Euler
    positivity step (Appendix A, Theorem 2).  With a limited reconstruction the
    coarse corrector is not a convex combination of forward-Euler steps
    (Appendix A, Remark 4); for the advection equation in one dimension with two
    levels Appendix B bounds the total variation, without the condition (4.15),
    and elsewhere the stability there rests on the tests.

## References

[1] C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comput. Phys.*, 77(2), 1988.

[7] M. H. Carpenter, D. Gottlieb, S. Abarbanel, W.-S. Don, "The theoretical accuracy of Runge–Kutta time discretizations for the initial boundary value problem: a study of the boundary error," *SIAM J. Sci. Comput.*, 16(6), 1241–1252, 1995, doi:10.1137/0916072 (ICASE Report 93-83, 1993).

[8] W. Hundsdorfer, *Numerical Solution of Advection–Diffusion–Reaction Equations*, lecture notes, Thomas Stieltjes Institute, CWI, Amsterdam, 2000.

[10] G. Jeanmasson, I. Mary, L. Mieussens, "On some explicit local time stepping finite volume schemes for CFD," *J. Comput. Phys.*, 397, 2019, 108818, doi:10.1016/j.jcp.2019.07.017 (Prop. RK2: section on the second-order scheme and its conservative correction).

[11] T.-T.-P. Hoang, L. Ju, W. Leng, Z. Wang, "High order explicit local time stepping methods for hyperbolic conservation laws," *Math. Comp.*, 89(324), 2020, 1807–1842, doi:10.1090/mcom/3507 (arXiv:1905.09705).

[12] M. J. Berger, P. Colella, "Local adaptive mesh refinement for shock hydrodynamics," *J. Comput. Phys.*, 82(1), 1989.

[13] N. Y. Gnedin, V. A. Semenov, A. V. Kravtsov, "Enforcing the Courant–Friedrichs–Lewy condition in explicitly conservative local time stepping schemes," *J. Comput. Phys.*, 359, 2018, 93–105, doi:10.1016/j.jcp.2018.01.008 (arXiv:1801.03108; section 2, conditions (3) and (4)).

[14] PETSc, `TSSetMaxStepRejections` (option `-ts_max_step_rejections`) in `src/ts/interface/ts.c` and the default `max_reject = 10` of `TSCreate` in `src/ts/interface/tscreate.c`; a step rejected more often fails with `TS_DIVERGED_STEP_REJECTED` (`src/ts/impls/explicit/rk/rk.c`), https://github.com/petsc/petsc, commit febf3a0.

[15] SUNDIALS, ARKODE: the defaults `MAXNCF` and `MAXCONSTRFAILS` (10) and `MAXNEF` (7) of `src/arkode/arkode_impl.h`, the failures of the nonlinear solver, of the inequality constraints and of the error test at which a step ends (`src/arkode/arkode.c`), set by `ARKodeSetMaxConvFails`, `ARKodeSetMaxNumConstrFails` and `ARKodeSetMaxErrTestFails` (`src/arkode/arkode_io.c`), https://github.com/LLNL/sundials, commit 312fc0f.

[16] Castro: `change_max` (1.1 by default, `Source/driver/_cpp_parameters`) in the check before and after the advance of `Source/driver/Castro_advance.cpp`, the fall of the time step tolerated within a step; `retry_subcycle_factor` (0.5) in `Source/driver/Castro_advance_ctu.cpp`, with `dt_cutoff` (1e-12) on the landing of the sub-steps, https://github.com/AMReX-Astro/Castro, commit 371ba40.

[17] E. Hairer, G. Wanner, the code DOPRI5, `WORK(2)`, the safety factor in step size prediction, 0.9 by default; and PETSc [14], the default `safety = 0.9` of `TSAdaptCreate` in `src/ts/adapt/interface/tsadapt.c`.

[18] SeisSol, `docs/local-timestepping.rst` (the maximum difference property), https://github.com/SeisSol/SeisSol, commit a526776.
