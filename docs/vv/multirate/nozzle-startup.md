# Multirate: Nozzle Start-Up

A quasi-1D convergent-divergent nozzle starts as its inlet total pressure rises from $10^5$ to $2.6\cdot10^5$ Pa in 8 ms against an outlet at $10^5$ Pa. The flow chokes at the throat, and a shock forms in the divergent and travels down it, crossing a face between levels. The cells shrink from 10 mm at the ends to 0.64 mm at the throat, so the throat cells take the finest levels.

Six runs on the same mesh (371 cells):

| Run | Step | Role |
|---|---|---|
| single rate | $dt_0/4$ | reference |
| single rate | $dt_0$, $dt_0/2$ | baseline error |
| multirate | $H = dt_0$ | throat at level 1 once the flow is fast |
| multirate | production $H$ and $H/2$ | levels set by the flow, map moving with it |

The step $dt_0$ is small enough that every cell of the fixed-step runs stays within its stability limit.

**Criteria.** In the throat band, the multirate run with $H = dt_0$ must have an error within 1.05 times the single rate's; the production runs must converge with the step (order at least 0.8, since their maps move with the flow, where second order is not proven: see [Temporal order](order.md#moving-faces-between-levels)); and the moving shock must stay within a face or so of the reference, judged at the outputs where it is a real shock (density jump of at least 1.5).

The test is sensitive to the coupling between levels: a run in which the finer levels read the coarser cells at the start of their step, instead of at their own stage time, fails the 1.05 criterion.

<figure>
  --8<-- "docs/vv/images/multirate-nozzle-short-shock.svg"
  <figcaption>Position of the moving shock (the face of largest compression in the divergent) at the outputs where it is a shock. Grey lines: the faces between the bands of 1.6 and 4 mm and of 4 and 10 mm, which are faces between levels.</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-nozzle-short-errors.svg"
  <figcaption>Left: density error in the throat band of every run against the reference. Right: error of the multirate run with <span class="arithmatex">\(H = dt_0\)</span> over that of the single rate at the outputs where it has levels, bounded by 1.05.</figcaption>
</figure>
