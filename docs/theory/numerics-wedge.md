# Axisymmetric Wedge and 3-D Sector

This page details rotational periodicity on the $k$ faces: the axisymmetric wedge, with the pressure and viscous terms that its two
faces carry, and the 3-D sector. The summary is in [Spatial Discretization](numerics.md#rotational-periodicity-and-axisymmetric-wedge).

## Rotational Periodicity and Axisymmetric Wedge

Rotational periodicity on the $k$ faces closes a mesh that turns about the $x$ axis: an
axisymmetric wedge (one cell in $k$) or a 3-D sector (several). The flow
beyond face 5 is the flow next to face 6 turned back by the angle $\delta$
of the wedge or sector, and the flow beyond face 6 the flow next to face 5
turned forward. With $g = +1$ on face 5 and $g = -1$ on face 6, the ghost
cell takes the state of its source cell (on a wedge the cell itself, on a
sector the cell next to the opposite $k$ face of the same block, in the same
column) with the velocity turned about $x$ by $\varphi = -g\,\delta$:

$$
v_g = v\cos\varphi - w\sin\varphi = v\cos\delta + g\,w\sin\delta, \qquad
w_g = v\sin\varphi + w\cos\varphi = -g\,v\sin\delta + w\cos\delta .
$$

$u$ and the scalars do not change under a turn about $x$, and are copied.
The ghost cells and the fluxes of these faces take $\cos\delta$ and
$\sin\delta$ computed once, with $\delta$, when the mesh is read. A rotationally
periodic face joins a block to itself: a block closed so holds the whole wedge or
sector in $k$, and a column of cells with the condition on one of its $k$ faces
has it on the other ([Boundary Conditions](../user/boundary-conditions.md#euler-boundaries)).

MOSE finds the wedge and the sector from the nodes (the rule is in the note
*Mesh type* of [Using MOSE](../user/using.md#case-directory-structure)): the
last $k$ node plane of the block must be the first one turned or translated to
within 1 % of the move of the nodes per cell in $k$, measured along the $k$
lines: the largest, over the nodes of the first $k$ node plane, of the sum of the
moves of the node across the cells in $k$ (each its largest component), over the
number of those cells (the thickness of a 2-D block).

Along straight $k$ lines it
is the move from the first plane to the last per cell, and on a sector the
length of the arc per cell in the norm of the largest component (between 0.71
and 1 of the arc), which the chord from plane to plane understates past half a
turn (4.2 times at 270 degrees).

### Axisymmetric wedge

A 2-D axisymmetric case is a wedge of the 3-D mesh: one cell in $k$, whose two
$k$ faces lie in the meridional planes $\theta = \mp\delta/2$ and are rotationally
periodic. MOSE integrates the 3-D equations on the wedge; there are
no axisymmetric source terms. The ghost cell across each wedge face is the cell
itself rotated by $\mp\delta$ about the $x$ axis, so every face flux follows
from the state of the cell:

- inviscid: the pressure and the momentum carried by the cell state rotated to
  the face (no Riemann problem). In the radial momentum the two faces give the
  $+p/r$ and the centrifugal terms of the axisymmetric equations;
- viscous: the flux of a face between the cell and its rotated image, with the
  stencil of an interior face (the $i$ and $j$ neighbours of the image rotated
  in the same way) and, with the stencil metric, the face metric of
  the same differences of the positions, the row of the image the row of the
  cell turned (exact for the linear fields that the rotation allows, see the 3-D
  sector below); with the mean metric, the inverse-mean face metric of
  the cell and its image. Both are computed once at setup (they depend on the
  mesh only).

On a planar 2-D mesh (faces 5 and 6 parallel, no wedge angle) the condition makes
the image the cell itself, and the fluxes of the two faces cancel.

For an axisymmetric field the mass, energy and scalar fluxes of the two faces
are equal and opposite, and only the momentum fluxes are applied. Without swirl
the viscous traction of each face is the normal stress

$$
\tau_{\theta\theta} = 2\mu\,\frac{u_r}{r} + \lambda\,\nabla\!\cdot\mathbf{u},
\qquad \lambda = -\tfrac23\,\mu ,
$$

with $\mu$ including the eddy viscosity, and the two faces together give the
$-\tau_{\theta\theta}/r$ of the radial momentum equation [5, 6]; with swirl
they also carry $\tau_{r\theta}$ and $\tau_{x\theta}$. Near the axis
$\tau_{\theta\theta}/r$ is as large as the radial stress terms it balances: in
Homann's stagnation-point flow [7] $\tau_{rr} = \tau_{\theta\theta}$, and
$\tau_{rr}/r$ and $\tau_{\theta\theta}/r$ cancel exactly while each grows like
$1/r$ towards the axis.

The azimuthal derivatives of the $i$ and $j$ faces are differences across the
rotated ghost cells divided by the $k$ edge of the cell, which is a chord: they
carry a relative error of order $\delta^2$ that the mesh does not reduce
($\delta^2/12$ of the hoop stress $2\mu u_r/r^2$ in the radial balance,
$2.5\times10^{-5}$ for $\delta = 1^\circ$). Keep the wedge thin.

Put the axis nodes at a small radius, as MOSE's meshes do ($10^{-6}\,\Delta r$):
a viscous run on a wedge whose axis is exactly at $r = 0$ stops at setup, with
every face metric. On such an axis the faces have zero area and no
normal, so their ghost cells copy the cells. With the mean metric the
wedge faces and the faces $i$ of the first row of cells take the derivative in
$r$ across the axis with the ghost cell one cell away and would see half of
$\partial u_r/\partial r$: an error of the hoop stress of those cells that
refinement does not reduce ($1/6$ of $2\mu u_r/r^2$ in Homann's flow, on 32 and
64 cells).

With the stencil metric the copy sits at the cell and the faces take the
difference on one side of it: in Homann's flow, linear in $r$, the first row
gets the hoop stress of a small axis radius, so the error of $1/6$ is gone; a
field that is not linear in $r$ stays farther from the result of the small
radius (manufactured-solution tests). The stop stays with every metric, a declared
limit: the axis at $r = 0$ is not treated completely (ghost cells beyond the
axis, not copies).

The value of a cell is the mean of the solution over the cell. To compute a
flux it is convenient to read this mean as the value of the solution at some
average point of the cell [8, p. 18]: MOSE takes the cell centre, the mean of
the eight nodes, as cell-centred schemes typically do [9, Sec. III]. The mean
of a linear field over a cell is its value at the centroid of the cell volume
[13, Eq. (3.1) and (3.12)]. Vinokur, for whom a finite-volume method needs the
location of this point only in the output of the solution, notes that the
centroid "is strictly valid only if all the components of U vary linearly
throughout the cell" and that, "since the distribution of U is not known, the
centroid has no particular advantage over the cell center, defined as the
vectorial average of the cell vertices" [8, p. 18]. The two points coincide
on a parallelepiped, but not next to the axis of a wedge: in the first row of
cells, whose inner nodes are on the axis or at a small radius, the centroid
lies at two thirds of the height of the cell from the axis and the cell centre
at one half, at every refinement.

There the operators that take the cell
centres are not exact for a linear field given by its cell means: with uniform
radial spacing and the stencil metric, the radial gradient at the faces between the
first and the second row is $1/9$ too small ($1/45$ between the second and the
third), the cell gradient of the turbulence sources of the first row is off by
$1/9$, and the term $u_r/r$ of the hoop stress is $1/3$ too large in the first
row ($1/27$ in the second).

These errors do not decrease as the mesh is refined: with the exact solution
given by its cell means, the truncation error of the first rows does not
converge (in the one-dimensional model below it grows as $1/\Delta r$). The
error of the solution converges nonetheless. In a one-dimensional model of the
viscous terms of the radial momentum next to the axis, with these positions,
the error against the cell means decreases at first order in the first row and
at second order in its mean weighted by the volume (the $L_1$ norm); against
the exact values at the cell centres, where MOSE reads the value of a cell, its
largest error decreases at second order.

In MOSE, with a manufactured solution
on wedges of $1^\circ$ and $0.5^\circ$ with SST, uniform radial spacing and 32,
64 and 128 cells across the radius (the manufactured-solution driver with the
SST sources of the wedge added, which the repository does not contain), the
error of $k$ and $\omega$ in the first row decreases at orders 0.95 to 1.02
against either value, and their $L_1$ error at 1.8 to 2.1: the first row loses
one order of accuracy, a local loss.

The $L_1$ errors of the density, the
pressure and the axial and radial velocities of the same runs decrease at
orders 0.02 to 1.67 between 64 and 128 cells; at 128 cells those of the
density, the pressure and the axial velocity are 3.6 to 3.9 times smaller on
the $0.5^\circ$ wedge than on the $1^\circ$ one. On irregular grids the error of
the solution can keep the order of the scheme where the truncation error
converges at a lower order or does not converge at all [9, Sec. V.A]; a global
error better behaved than the local truncation error would indicate is the
property called supraconvergence [17, Sec. 1], here that of the $L_1$ error of
the model and of $k$ and $\omega$.

Next to the axis the hoop stress also stiffens the explicit step. In the radial
momentum of the first row of cells, at $r \approx \Delta r/2$, it is a damping
$(2\mu+\lambda)\,u_r/r^2 \approx \tfrac{16}{3}\,\mu\,u_r/\Delta r^2$ ($\mu$ with
the eddy viscosity), with couplings to the neighbours through $\lambda\,\nabla\cdot\mathbf{u}$.
The local time step includes them ([Time Integration](time-integration.md#local-time-step)),
so the VNN limits of the rule hold next to the axis too; without them, at a
Prandtl number of 10, the first rows broke down at the VNN the rule declares
stable.

### 3-D sector

On a 3-D sector (rotational periodicity with more than one cell in $k$) faces 5 and 6
make one interior face of the whole annulus, computed once, by the record of
face 6, for the inviscid and for the viscous flux. The image beyond face 6 is
the cell next to face 5 turned forward by $R$, the turn about $x$ by the sector
angle: the $i$ and $j$ differences of the image row are those of that row turned, the
metric of the image is that of its cell turned, $M R^{\mathsf T}$ (each row of
$M = \partial\xi/\partial\mathbf{x}$ is a vector), and the wall distance and
roughness of the face are the mean of the two cells, as on an interior face.

The cell next to face 6 takes the flux and the cell next to face 5 the same
flux with its momentum turned back, $R^{\mathsf T}F$: mass, energy, species,
turbulence variables, axial momentum and angular momentum about $x$ are
conserved. A translated sector (faces 5 and 6 parallel) is the same with
$R = I$. The viscous face of the wedge and of the sector takes, with
the stencil metric, the metric of the same differences of the
positions as the gradient takes of the values: across the face the image, at
the position of its ghost cell, minus the cell; along the face the differences
of the row of the cell and of the row of the image, which are those of the row
of the source turned by $R$, positions and values (a ghost cell of the row of
the source, behind a wall too, turned with it).

A linear field that the
periodicity allows, $\mathbf{u}(R\mathbf{x}) = R\,\mathbf{u}(\mathbf{x})$, is so
exact at these faces on any mesh; where a wall ghost cell enters the row of the
cell or of the source, the face takes one metric per group of quantities, as an
interior face. With the mean metric the face takes the inverse of the
mean of the inverse metrics of the cell and of its image, whose columns along
the face are the mean edges of the two cells, not the differences of the
centres: where the cell size grows along the face by a ratio $q$ the
derivatives along the face carry the factor $(1+q)^2/(4q)$ of the note
*Tangential derivatives* above, a relative error of $(q-1)^2/(4q)$
($8.3\times10^{-3}$ at $q = 1.2$), and in the first row next to an axis with
symmetry, whose ghost cell is the mirror image of the cell, $(3+q)/4$, a
relative error of $(q-1)/4$ ($5\times10^{-2}$ at $q = 1.2$); refinement at a fixed
ratio does not reduce them.

For $u_r = U r$, whose viscous force is zero, the
radial viscous force of a cell of the wedge is, with the inverse-mean metric on
every face, that relative error times the hoop term $\tau_{\theta\theta}/r$ on top of the error
of order $\delta^2$ above, and with the stencil metric that error alone,
$\sin^2(\delta/2)$ of the hoop term in every cell (manufactured-solution tests).
A viscous run on a sector whose axis nodes are
exactly at $r = 0$ stops at setup, as on the wedge, with every
face metric: behind faces of zero area the ghost cells copy the cells,
and with the mean metric the periodic faces of the first row of cells
would see half of $\partial u_r/\partial r$ ($1/6$ of the hoop term in Homann's
flow; with the stencil metric the error of a small axis radius, next to the periodic faces
too). Put the axis nodes at a small radius.

## References

[5] R. Vandenhoeck, A. Lani, "Implicit high-order flux reconstruction solver for high-speed compressible flows," arXiv:1811.08666v1, 2018, Appendix A.2, eqs. (A.6)–(A.7) (axisymmetric source $(0, 0, p - \tau_{\theta\theta}, 0)^T$).

[6] Eilmer 4, Gas Dynamics Toolkit, `src/eilmer/fvcell.d`, lines 1213–1228 (viscous source of axisymmetric flow), https://github.com/gdtk-uq/gdtk, commit 826496a; the code is described in N. N. Gibbons, K. A. Damm, P. A. Jacobs, R. J. Gollan, "Eilmer: an open-source multi-physics hypersonic flow solver," *Comput. Phys. Commun.*, 282, 108551, 2023, doi:10.1016/j.cpc.2022.108551.

[7] F. Homann, "Der Einfluss grosser Zähigkeit bei der Strömung um den Zylinder und um die Kugel," *Z. Angew. Math. Mech.*, 16, 153–164, 1936 (not read; as cited by J. Magnaudet and M. Abbas, arXiv:2011.10755v1, 2020). The equation of the flow, f''' + 2ff'' − f'² + 1 = 0 with f(0) = f'(0) = 0, f'(∞) = 1, u_r = a r f'(η), u_x = −2(aν)^{1/2} f(η), as written by P. Rajamanickam, arXiv:2508.17982v4, Eq. (12), and by H. Guo, X. Zhuang, X. Meng, T. Rabczuk, arXiv:2010.05682v1, Eq. (7) with β₀ = 2, β = 1, whose Table 8 gives f''(0) = 1.31194.

[8] M. Vinokur, "An analysis of finite-difference and finite-volume formulations of conservation laws," NASA CR-177416, 1986: pp. 15–16, Eq. (60a) (the area vector of a face of a hexahedron as half the cross product of its diagonals; the shape of each face "precisely defined and consistently used by the two neighboring cells" so that the volumes sum to the total volume; the doubly ruled face); p. 18 (the value of a cell taken as the value at some average point of the cell to compute a flux; the centroid and the vectorial average of the vertices); p. 26, Eq. (101) ("the two difference approximations must be the same"); p. 23, Eq. (90)–(91) (the transport terms on an auxiliary cell, with means of four cells along the face); p. 40, Eq. (153) (the difference across a wall). Journal version: *J. Comput. Phys.*, 81, 1989.

[9] B. Diskin, J. L. Thomas, E. J. Nielsen, H. Nishikawa, J. A. White, "Comparison of Node-Centered and Cell-Centered Unstructured Finite-Volume Discretizations: Viscous Fluxes," *AIAA J.*, 48(7), 1326–1338, 2010, doi:10.2514/1.44940, Sec. III, III.A and V.A.

[13] H. Jasak, *Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows*, PhD thesis, Imperial College, University of London, 1996, Eq. (3.1), (3.12) and (3.32).

[17] B. Wendroff, A. B. White, Jr., "Some supraconvergent schemes for hyperbolic equations on irregular grids," Los Alamos National Laboratory, LA-UR-88-1384, 1988 (2nd International Conference on Hyperbolic Problems, Aachen, 1988).
