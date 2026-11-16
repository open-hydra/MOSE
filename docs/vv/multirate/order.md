# Multirate: Temporal Order

The multirate macro step is second order in time when the level map is fixed ([theory](../../theory/multirate-method.md#recursive-heun-scheme)). When the map changes with the flow, so that the faces between levels move, the analysis guarantees second order only for simple linear fluxes. This page measures both.

## Fixed level map

Each case runs with its levels fixed from the initial state and with four macro steps, $H_0$, $H_0/2$, $H_0/4$ and $H_0/8$, up to $T = 64 H_0$. The error is taken against a single-rate run on the same mesh with a step 64 times smaller than the smallest cell step, so it is the error of the time integration alone.

The eight cases cover inviscid, laminar and turbulent flow:

- acoustic and entropy waves in 1D, on bands of cells of different size (levels 0 to 3);
- an isentropic vortex crossing two block connections (levels 0 and 1);
- the laminar flat plate from a uniform start (levels 0 to 3, by rows);
- a hot square carried by a viscous flow on skewed cells (levels 0 and 1);
- a no-slip wall on skewed cells whose level changes along the wall (levels 0 to 2, figure below);
- wall channels with Spalart–Allmaras and SST (levels 0 to 3).

<figure>
  --8<-- "docs/vv/images/multirate-order.svg"
  <figcaption>Density error of the eight cases against the macro step <span class="arithmatex">\(H\)</span>, relative to the error with <span class="arithmatex">\(H_0\)</span>. Left: <span class="arithmatex">\(L^1\)</span> error over all cells. Right: <span class="arithmatex">\(L^\infty\)</span> error over the cells within 4 cells of a face between levels. Dashed: second order. The legend gives the number of cells and the observed order of the last pair.</figcaption>
</figure>

**Result.** The observed order of the last pair is 1.93–2.15 in $L^1$ and 1.91–2.19 next to the faces between levels. The test requires at least 1.9 in $L^1$ (1.8 for the turbulent channels) and 1.8 next to the faces.

The smallest margin is on the entropy wave (1.93 against 1.9). The light form of the test, which a pull request runs, uses only $H_0$ and $H_0/2$ with the same thresholds; a change that lowers that ratio by about 2% fails it, and the full form then decides.

<figure>
  --8<-- "docs/vv/images/multirate-wall-s-levels.svg"
  <figcaption>Level map of the skewed-wall case (48 × 8 cells): bands of 8, 2.8, 1.28, 2.8 and 8 mm along the no-slip wall (black line) take levels 0, 1, 2, 1, 0 (grey, blue, green). Four faces between levels cut the wall, on cells skewed by 14 to 37 degrees.</figcaption>
</figure>

## Moving faces between levels

A simple wave of the Euler equations is refined in space and time together ($\Delta x$ and $H$ halved at each step, 256 to 2048 cells) up to $T = 0.5$, before the wave breaks, and compared with its exact solution. With two levels, the faces between them sit where the wave speed $u + c$ is halfway between its extremes, so they travel with the wave.

This is the case the analysis does not fully cover: a cell next to a face between levels has a local error of first order ([theory](../../theory/multirate-analysis.md#appendix-a-order-and-stability-of-the-multirate-macro-step), Remark 2), and when the face moves with the flow these errors accumulate on the same part of the wave. For linear fluxes the analysis still bounds the error as $O(H^2)$ in $L^1$ and $O(H)$ in $L^\infty$ ([theory](../../theory/multirate-analysis.md#corollary-every-sequence-of-maps)); for the limited reconstruction and HLLC of this test there is no proof.

<figure>
  --8<-- "docs/vv/images/multirate-moving-levels-order.svg"
  <figcaption>Density error of the simple wave against the number of cells <span class="arithmatex">\(n\)</span>. Left: <span class="arithmatex">\(L^1\)</span> error. Right: <span class="arithmatex">\(L^\infty\)</span> error. Thick: multirate with the level map rebuilt at every macro step. Also shown: the map of the first macro step kept, and the single rate. Dashed: second order; dotted: first order.</figcaption>
</figure>

**Result.**

- In $L^1$ the multirate with moving faces is second order, as the analysis predicts for linear fluxes. The test requires at least 1.9.
- In $L^\infty$ every run, the single rate included, stays near 1.3, because the largest error is at the crest of the wave, away from the faces.
- On the cells next to the moving faces alone, the $L^\infty$ order of $\rho$, $u$ and $p$ is 1.80, 1.75 and 1.67, against 2.0 for the single rate and 2.1 for a fixed map on the same cells. This is the slower local convergence that the analysis expects; it is reported, not tested.

The test is one-dimensional, inviscid, periodic, in one block and with two levels. Moving faces in several dimensions, with more levels, boundaries or viscous terms are not measured.
