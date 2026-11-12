#!/usr/bin/env bash
#===============================================================================
#  Start-up of a nozzle with the multirate integrator: a quasi-1D convergent-divergent nozzle (case mr1d-nozzle-su of
#  gen_mr_cases.py, 371 cells) whose inlet total pressure rises from 1e5 to 2.6e5 Pa in 8 ms against an outlet at 1e5 Pa.
#  The flow starts and chokes; the shock that forms in the divergent moves down it across a face between levels. The
#  test has one form (CTest Multirate.Nozzle, a few minutes).
#  dt_min(0) from one step on the one-level path; dt0 = (t_fin/10)/M with M the smallest multiple of 32 such that
#  dt0 <= 0.45 dt_min(0) (dt_min falls to about 0.5 dt_min(0) as the flow speeds up through the throat). Six runs on the
#  multirate path, the level map rebuilt at every macro step, ten solutions at t_k = k t_fin/10:
#    ref  max-level 0, safety 1, step dt0/4 (the reference);   sr, sr2  max-level 0, safety 1, steps dt0 and dt0/2;
#    mrb  max-level 3, safety 2, H = dt0 (the cells of the throat band on level 1 once the flow is fast);
#    mra  max-level 3, safety 2, the production macro step H = 2^L dt_min / safety (no fixed step, output instants hit);
#    mra2 the same with safety 4 (H/2).
#  Checks: verify_nozzle.py; figures: plot_nozzle.py (skipped without matplotlib).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
CASE=mr1d-nozzle-su
rm -rf cases runs
mr_cases cases "$CASE" || { echo "nozzle: FAIL - cannot build the case"; exit 1; }
C=cases/$CASE
TFIN=$(python3 -B -c "import re, sys; print(re.search(r'time-threshold\s*=\s*(\S+)', open(sys.argv[1]).read()).group(1))" "$C/input.ini")

# dt_min(0): one step with H = dt_min (one level, fixed-macro-step 1); its time is dt_min(0) exactly
mr_run runs/probe "$C" "MOSE-Parameters|iter-threshold|1" "MOSE-Parameters|time-threshold|1e9" \
  "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
  "MOSE-Multirate|fixed-macro-step|1.0" || { echo "nozzle: FAIL - the run of one step"; exit 1; }
read -r M DT0 FREF FSR FSR2 < <(python3 -B - "$TFIN" runs/probe <<'EOF'
import math, os, sys
sys.path.insert(0, os.path.join('..', 'mr-common'))
import mrcheck as mc
tfin = float(sys.argv[1])
dtf = mc.read_run(sys.argv[2])[0]['time']
M = 32 * math.ceil((tfin / 10.0) / (32 * 0.45 * dtf))
dt0 = (tfin / 10.0) / M
print(M, repr(dt0), repr(dt0 / (4.0 * dtf)), repr(dt0 / dtf), repr(dt0 / (2.0 * dtf)))
EOF
)
[[ -n ${FSR2:-} ]] || { echo "nozzle: FAIL - dt0"; exit 1; }
echo "  M = $M steps of sr per t_fin/10, dt0 = $DT0 s"
COM=("MOSE-IO|sol-overwrite|false" "MOSE-IO|res-diter|$M" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|levels|recompute"
     "MOSE-Multirate|recompute-every|1")
FIX=("MOSE-Parameters|time-threshold|1e9")
ONE=("MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0")
PROD=("MOSE-Parameters|time-threshold|$TFIN" "MOSE-Parameters|iter-threshold|1000000000" "MOSE-IO|sol-diter|1000000000"
      "MOSE-IO|sol-dtime|$(python3 -c "print(repr($TFIN / 10.0))")" "MOSE-IO|shell-diter|1" "MOSE-Multirate|max-level|3")
status=0
mr_run runs/ref "$C" "${COM[@]}" "${FIX[@]}" "${ONE[@]}" "MOSE-Parameters|iter-threshold|$((40 * M))" \
  "MOSE-IO|sol-diter|$((4 * M))" "MOSE-Multirate|fixed-macro-step|$FREF" || status=1
mr_run runs/sr "$C" "${COM[@]}" "${FIX[@]}" "${ONE[@]}" "MOSE-Parameters|iter-threshold|$((10 * M))" \
  "MOSE-IO|sol-diter|$M" "MOSE-Multirate|fixed-macro-step|$FSR" || status=1
mr_run runs/sr2 "$C" "${COM[@]}" "${FIX[@]}" "${ONE[@]}" "MOSE-Parameters|iter-threshold|$((20 * M))" \
  "MOSE-IO|sol-diter|$((2 * M))" "MOSE-Multirate|fixed-macro-step|$FSR2" || status=1
mr_run runs/mrb "$C" "${COM[@]}" "${FIX[@]}" "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" \
  "MOSE-Parameters|iter-threshold|$((10 * M))" "MOSE-IO|sol-diter|$M" "MOSE-IO|shell-diter|1" \
  "MOSE-Multirate|fixed-macro-step|$FSR" "MOSE-Multirate|macro-step-target|$DT0" || status=1
mr_run runs/mra "$C" "${COM[@]}" "${PROD[@]}" "MOSE-Multirate|safety|2.0" || status=1
mr_run runs/mra2 "$C" "${COM[@]}" "${PROD[@]}" "MOSE-Multirate|safety|4.0" || status=1
[[ $status -eq 0 ]] && { python3 -B verify_nozzle.py runs "$C" "$M" --json runs/nozzle.json || status=1; }
[[ $status -eq 0 ]] && { python3 -B plot_nozzle.py runs/nozzle.json runs || echo "  plot_nozzle: figures not written"; }
if [[ $status -eq 0 ]]; then echo "mr-nozzle: ok"; else echo "mr-nozzle: FAIL"; fi
exit $status
