#!/usr/bin/env bash
#===============================================================================
#  Viscous shock tube of Daru and Tenaud with the multirate integrator: a shock and a contact running into a closed
#  box, reflected by the end wall onto the boundary layer of the bottom wall, Re = 200, in three forms on cases of
#  gen_mr_cases.py: short (CTest Multirate.Daru), ./run.sh, mr2d-daru-coarse, 50 x 25 cells; medium (Multirate.DaruMedium),
#  ./run.sh medium, mr2d-daru, 100 x 50 cells; full (Multirate.DaruFull), ./run.sh full, mr2d-daru-full, 500 x 250
#  cells, the resolution of the reference solutions (hours; OMP_NUM_THREADS as set by the caller).
#  dt_min(0) from one step on the one-level path; dt0 = (t_fin/10)/M with M the smallest multiple of 32 such that
#  dt0 <= 0.9 dt_min(0). Four runs on the multirate path, the level map rebuilt at every macro step, fixed-macro-step
#  together with macro-step-target, ten numbered solutions at t_k = k t_fin/10:
#    ref  max-level 0, safety 1, step dt0/4 (the reference);
#    sr   max-level 0, safety 1, step dt0 (the single-rate step);
#    mra  max-level 3, safety 2, H = 4 dt0: the largest power of two of dt0 within 2^L dt_min(0)/safety = 4.44 dt0;
#    mrb  max-level 3, safety 2, H = dt0: the core cells take the step of sr, the cells next to the walls, where the
#         viscous limit rules, smaller ones.
#  Checks: verify_daru.py. Figures (docs/vv/multirate.md): plot_daru.py writes runs/multirate-daru-<form>-*.svg from the
#  runs and from the errors that verify_daru.py writes with --json (mra and mrb log their cells per level at every macro
#  step, shell-diter 1); every figure names its mesh and form, and only the full form is drawn next to the published
#  solutions. Without matplotlib the figures are skipped and the test is unchanged.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
case ${1:-short} in
  short)  CASE=mr2d-daru-coarse ;;
  medium) CASE=mr2d-daru ;;
  full)   CASE=mr2d-daru-full ;;
  *) echo "usage: $0 [short|medium|full]"; exit 2 ;;
esac
rm -rf cases runs
mr_cases cases "$CASE" || { echo "daru-tenaud: FAIL - cannot build the case"; exit 1; }
C=cases/$CASE

# dt_min(0): one step with H = dt_min (one level, fixed-macro-step 1); its time is dt_min(0) exactly
mr_run runs/probe "$C" "MOSE-Parameters|iter-threshold|1" "MOSE-Parameters|time-threshold|1e9" \
  "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
  "MOSE-Multirate|fixed-macro-step|1.0" || { echo "daru-tenaud: FAIL - the run of one step"; exit 1; }
read -r M DT0 FREF FSR FA < <(python3 -B - "$C/input.ini" runs/probe <<'EOF'
import math, os, re, sys
sys.path.insert(0, os.path.join('..', 'mr-common'))
import mrcheck as mc
tfin = float(re.search(r'time-threshold\s*=\s*(\S+)', open(sys.argv[1]).read()).group(1))
dtf = mc.read_run(sys.argv[2])[0]['time']
M = 32 * math.ceil((tfin / 10.0) / (32 * 0.9 * dtf))
dt0 = (tfin / 10.0) / M
print(M, repr(dt0), repr(dt0 / (4.0 * dtf)), repr(dt0 / dtf), repr(4.0 * dt0 / dtf))
EOF
)
[[ -n ${FA:-} ]] || { echo "daru-tenaud: FAIL - dt0"; exit 1; }
echo "  M = $M steps of sr per t_fin/10, dt0 = $DT0 s"
COM=("MOSE-Parameters|time-threshold|1e9" "MOSE-IO|sol-overwrite|false" "MOSE-IO|res-diter|$M"
     "MOSE-Multirate|enabled|.true." "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1")
ONE=("MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0")
LEV=("MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0")
status=0
mr_run runs/ref "$C" "${COM[@]}" "${ONE[@]}" "MOSE-Parameters|iter-threshold|$((40 * M))" "MOSE-IO|sol-diter|$((4 * M))" \
  "MOSE-Multirate|fixed-macro-step|$FREF" || status=1
mr_run runs/sr "$C" "${COM[@]}" "${ONE[@]}" "MOSE-Parameters|iter-threshold|$((10 * M))" "MOSE-IO|sol-diter|$M" \
  "MOSE-Multirate|fixed-macro-step|$FSR" || status=1
mr_run runs/mra "$C" "${COM[@]}" "${LEV[@]}" "MOSE-Parameters|iter-threshold|$((10 * M / 4))" "MOSE-IO|sol-diter|$((M / 4))" \
  "MOSE-Multirate|fixed-macro-step|$FA" "MOSE-IO|shell-diter|1" \
  "MOSE-Multirate|macro-step-target|$(python3 -c "print(repr(4 * $DT0))")" || status=1
mr_run runs/mrb "$C" "${COM[@]}" "${LEV[@]}" "MOSE-Parameters|iter-threshold|$((10 * M))" "MOSE-IO|sol-diter|$M" \
  "MOSE-Multirate|fixed-macro-step|$FSR" "MOSE-Multirate|macro-step-target|$DT0" "MOSE-IO|shell-diter|1" || status=1
[[ $status -eq 0 ]] && { python3 -B verify_daru.py runs "$C" "$M" --json runs/daru.json || status=1; }
[[ $status -eq 0 ]] && { python3 -B plot_daru.py runs "$C" runs/daru.json || echo "  plot_daru: figures not written"; }
if [[ $status -eq 0 ]]; then echo "mr-daru: ok"; else echo "mr-daru: FAIL"; fi
exit $status
