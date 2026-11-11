#!/usr/bin/env bash
#===============================================================================
#  Multirate order with moving level boundaries: the 1D periodic
#  simple wave mr1d-simple-<n> of gen_mr_cases.py (c = 1 + 0.05 phi, phi = (1 - cos 2 pi x) / 2, u = 2 (c - 1) / (gamma - 1),
#  non-dimensional gas), n = 256, 512, 1024, 2048, refined jointly in space and time up to T = 0.5, before the wave breaks
#  (t = 1.06); Heun, HLLC, MUSCL with van Leer, cfl 0.8. A probe run of one step of the one-level path gives dt0, the
#  step of the single rate at t = 0; then for every n:
#    mr    multirate, max-level 1, safety 2, levels recomputed at every macro step, macro step H = T / N with
#          N = ceil(T / (0.565 dt0)): the level boundaries lie where u + c is halfway between its extremes (phi = 1/2) and
#          travel with the wave, on the characteristics of x0 = 1/4 and 3/4;
#    mrf   the same with the map of the first macro step kept (levels = fixed);
#    srH   single rate with the same step H;
#    srdt  single rate with the step T / M, M = ceil(T / dt0);
#  fixed-macro-step with macro-step-target = H in mr and mrf, so that the map is built for the macro step taken; shell-diter
#  = the number of macro steps in mr and mrf, so that the log has the map of the first and of the last macro step.
#  Decided on mr (verify_moving_levels.py): on the last pair (1024 -> 2048) the L1 order of rho, u and p against the exact solution
#  by characteristics at least 1.9; every run of mr and mrf valid (final time T within 1e-12, cells at both levels in
#  the first and in the last map); mass drift of mr at most 1e-12. Reported: the Linf orders, global and on the cells
#  within 4 cells of the two faces between levels at the final time (the faces placed from the exact solution), mrf, srH
#  and srdt, the differences of mr from srdt and from srH.
#  Figure for the V&V page (plot_moving_levels.py, skipped without matplotlib; the test does not depend on it).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
NS=(256 512 1024 2048)
mr_cases cases $(printf 'mr1d-simple-%s ' "${NS[@]}") || { echo "mr-moving-levels: FAIL - cannot build the cases"; exit 1; }
status=0
BASE=("MOSE-Parameters|time-threshold|1.0e9" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|recompute-every|1")
for n in "${NS[@]}"; do
  mr_run runs/probe-n$n cases/mr1d-simple-$n "${BASE[@]}" "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
    "MOSE-Multirate|fixed-macro-step|1.0" "MOSE-Parameters|iter-threshold|1" || { echo "mr-moving-levels: FAIL - probe run n $n"; exit 1; }
  read -r N FMR M FSR HMR < <(python3 -B - runs/probe-n$n <<'EOF'
import math
import sys
sys.path.insert(0, '../mr-common')
import mrcheck as mc
dt0 = mc.read_run(sys.argv[1])[0]['time']
T = 0.5
N, M = math.ceil(T / (0.565 * dt0)), math.ceil(T / dt0)
print(N, repr(T / (N * dt0)), M, repr(T / (M * dt0)), repr(T / N))
EOF
)
  [[ -n ${HMR:-} ]] || { echo "mr-moving-levels: FAIL - dt0 of the probe run n $n not read"; exit 1; }
  echo "  n $n: mr, mrf, srH: $N macro steps of H = $HMR; srdt: $M steps"
  for spec in "mr|1|2.0|$FMR|$N|recompute|$HMR" "mrf|1|2.0|$FMR|$N|fixed|$HMR" "srH|0|1.0|$FMR|$N|fixed|0.0" \
              "srdt|0|1.0|$FSR|$M|fixed|0.0"; do
    IFS='|' read -r lab L sig f it lev tgt <<< "$spec"
    # mr and mrf print the cells per level at the first and at the last macro step (shell-diter = the macro steps)
    sd=1000000000; [[ $L == 1 ]] && sd=$it
    mr_run runs/$lab-n$n cases/mr1d-simple-$n "${BASE[@]}" "MOSE-Multirate|max-level|$L" "MOSE-Multirate|safety|$sig" \
      "MOSE-Multirate|fixed-macro-step|$f" "MOSE-Parameters|iter-threshold|$it" "MOSE-Multirate|levels|$lev" \
      "MOSE-Multirate|macro-step-target|$tgt" "MOSE-IO|shell-diter|$sd" || { echo "  run $lab n $n failed"; status=1; }
  done
done
[[ $status -eq 0 ]] || { echo "mr-moving-levels: FAIL"; exit 1; }
V=(python3 -B verify_moving_levels.py cases runs)
"${V[@]}" mr "${NS[@]}" --multirate --decide --ref srdt --json runs/mr.json || status=1
"${V[@]}" mr "${NS[@]}" --ref srH | command grep -a "difference\|order of the difference"
"${V[@]}" mrf "${NS[@]}" --multirate --json runs/mrf.json || status=1
for lab in srH srdt; do "${V[@]}" $lab "${NS[@]}" --json runs/$lab.json || status=1; done
python3 -B plot_moving_levels.py runs || echo "  (figure not drawn)"
if [[ $status -eq 0 ]]; then echo "mr-moving-levels: ok"; else echo "mr-moving-levels: FAIL"; fi
exit $status
