# Multirate: Viscous Shock Tube of Daru and Tenaud

A shock moves along a closed tube at Re = 200, reflects from the end wall and interacts with the boundary layer it has created, producing a complex pattern of shocks, vortices and a separated wall flow (Daru and Tenaud, *Computers & Fluids* 30, 2001). The lower half of the tube is computed. Near the walls the viscous stability limit is much smaller than in the core, so the multirate integrator uses several levels.

The same mesh is run four times:

| Run | Step | Role |
|---|---|---|
| single rate | $dt_0/4$ | reference |
| single rate | $dt_0$ | baseline error |
| multirate | $H = dt_0$ | core cells at the single-rate step, near-wall cells finer |
| multirate | $H = 4\,dt_0$ | up to 3 levels, the production setting |

The multirate run with $H = dt_0$ must have an error within 1.05 times that of the single rate with $dt_0$ at every output time. The test runs on three meshes: 50 × 25 cells in every pull request, and 100 × 50 and 500 × 250 cells by hand. Only the finest mesh is compared with the published solution.

## Medium mesh (100 × 50)

<figure>
  --8<-- "docs/vv/images/multirate-daru-medium-errors.svg"
  <figcaption>Left: density error against the reference of the single rate with <span class="arithmatex">\(dt_0\)</span> and of the two multirate runs, at the ten output times. Right: error of the multirate run with <span class="arithmatex">\(H = dt_0\)</span> over that of the single rate; the test requires at most 1.05.</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-medium-levels.svg"
  <figcaption>Share of the cells on levels 1 to <span class="arithmatex">\(L\)</span> at every macro step of the two multirate runs (level 0 holds the rest), and their work: cell updates relative to the single rate with <span class="arithmatex">\(dt_0\)</span>.</figcaption>
</figure>

## Fine mesh (500 × 250)

- The multirate error with $H = dt_0$ stays within 1.0016 times the single-rate error at every output time.
- With $H = 4\,dt_0$ the multirate run costs 0.27 times the cell updates of the single rate.
- At the final time the wall densities of the four runs differ by less than $10^{-5}$, and they agree with table 1 of Zhou, Xu and Liu (*Physics of Fluids* 30, 2018; 1500 × 750 cells) within 0.79% on average and 4.2% at most, at $x = 0.6123$.

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-errors.svg"
  <figcaption>Left: density error of the single rate with <span class="arithmatex">\(dt_0\)</span> and of the two multirate runs at the ten output times. Right: error of the multirate run with <span class="arithmatex">\(H = dt_0\)</span> over that of the single rate, within the limit of 1.05.</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-levels.svg"
  <figcaption>Share of the cells on levels 1 to <span class="arithmatex">\(L\)</span> at every macro step, and the work relative to the single rate: 0.27 for <span class="arithmatex">\(H = 4\,dt_0\)</span> and 1.01 for <span class="arithmatex">\(H = dt_0\)</span>.</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-wall.svg"
  <figcaption>Density along the bottom wall at the final time: the four runs lie on one another, with the 20 points of table 1 of Zhou, Xu and Liu.</figcaption>
</figure>

<figure>
  --8<-- "docs/vv/images/multirate-daru-full-density.svg"
  <figcaption>Density field of the multirate run with <span class="arithmatex">\(H = 4\,dt_0\)</span> at the final time, <span class="arithmatex">\(0.3 \le x \le 1\)</span>, <span class="arithmatex">\(y \le 0.3\)</span>.</figcaption>
</figure>
