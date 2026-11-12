#!/usr/bin/env bash
#===============================================================================
#  Multirate smoke test (a few seconds): six checks of the multirate integrator, the first five on the 1-D
#  cases of gen_mr_cases.py (284 cells in bands of 10 / 4 / 1.6 / 0.64 / 1.6 / 4 / 10 mm):
#   1. one level (max-level 0, safety 1) is the single-rate Heun step: the acoustic wave, 20 steps, fields and residual
#      history identical byte by byte;
#   2. four levels (max-level 3, safety 2, map rebuilt at every macro step) conserve: the Sod tube closed by slip
#      walls, 20 macro steps, mass and energy drift <= 1e-12;
#   3. the restart is exact: the same run against 10 macro steps + restart (vtk raw) + 10, the two output trees
#      identical byte by byte (numbered and final solutions, residual history);
#   4. a configuration outside the scope stops at set-up: cfl 1.2 with max-level 3 and nonlocal-cfl = redo (the
#      default), a non-zero exit code, the message "[MOSE-Multirate] cfl" and no solution;
#   5. a run that diverges stops at the state check of the stage, as the single-rate step does: the Sod tube with one
#      level and cfl 20, exit code 1, the message of the check of the multirate (a value that is not finite, or a
#      negative pressure) and no solution;
#   6. the lines of the gas entering through 400 faces, which the boundary fluxes print at the stages of the
#      single-rate step, come at the stages of the multirate: a 1-D tube of 8 cells with 400 ends and the gas entering
#      through face 1 at 0.5 a (the case count of test/checks/boundary-conditions/extrapolation-400/entering_gas.py),
#      20 Heun steps at first order; one level gives the lines and the end line of the single-rate step word for
#      word: one line (iteration 1, faces: 1) and 40 flux evaluations, one a stage.
#  The complete checks are the tests Multirate.OneLevel, Multirate.Conservation and Multirate.Reject.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
CHECK=../../checks/restart/check_restart_vtk.py
rm -rf cases runs
mr_cases cases mr1d-acu mr1d-sod || { echo "mr-smoke: FAIL - cannot build the cases"; exit 1; }
status=0
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0"
    "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1")
OUT=("MOSE-IO|res-diter|1" "MOSE-IO|sol-overwrite|false")

# part DIR -> one run of MOSE in DIR (input.ini set), appended to solve.log; 0 if exit 0 on the threads asked
part() {
  ( cd "$1" && "$MR_MASTER" >> solve.log 2>&1 ) || { echo "  run $1: exit non-zero"; tail -3 "$1/solve.log"; return 1; }
  mr_threads_ok "$1" "$OMP_NUM_THREADS"
}

# 1. one level = single-rate step
K=("MOSE-IO|res-diter|1" "MOSE-Parameters|iter-threshold|20")
if mr_run runs/one-sr cases/mr1d-acu "${K[@]}" &&
   mr_run runs/one-mr cases/mr1d-acu "${K[@]}" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" \
     "MOSE-Multirate|safety|1.0" &&
   mr_same runs/one-sr runs/one-mr; then
  echo "  one level = single-rate step: identical (fields and residual history)"
else
  echo "  one level = single-rate step: FAIL"; status=1
fi

# 2. conservation with four levels: 20 macro steps in one go
d=runs/restart
if mr_prepare "$d/full" cases/mr1d-sod "${MR[@]}" "${OUT[@]}" "MOSE-Parameters|iter-threshold|20" "MOSE-IO|sol-diter|8" &&
   part "$d/full" &&
   python3 -B ../mr-common/verify_conservation.py cases/mr1d-sod "$d/full" --quantities mass,energy --max 1e-12 &&
   [[ $(command grep -a -m1 'Multirate: cells per level' "$d/full/solve.log" | sed 's/.*level//' | wc -w) -eq 4 ]]; then
  echo "  conservation with four levels: ok ($(command grep -a -m1 'Multirate: cells per level' "$d/full/solve.log" | sed 's/^ *//'))"
else
  echo "  conservation with four levels: FAIL"; status=1
fi

# 3. exact restart: the run of 2. against 10 macro steps + restart + 10
if [[ -f $d/full/OUTPUT/field.vtm ]] &&
   mr_prepare "$d/split" cases/mr1d-sod "${MR[@]}" "${OUT[@]}" "MOSE-Parameters|iter-threshold|10" "MOSE-IO|sol-diter|8" &&
   part "$d/split" &&
   python3 -B "$MR_COMMON/set_ini.py" "$d/split/input.ini" "MOSE-Parameters|newrun|false" "MOSE-Parameters|iter-threshold|10" &&
   part "$d/split" &&
   [[ -s $d/full/OUTPUT/residual-history.dat && -f $d/full/OUTPUT/field.vtm ]] &&
   python3 -B "$CHECK" "$d/full/OUTPUT" "$d/split/OUTPUT" > "$d/compare.log" 2>&1; then
  echo "  restart (20 = 10 + restart + 10): identical ($(cd "$d/full/OUTPUT" && find . -type f ! -name '*.log' | wc -l) files)"
else
  echo "  restart (20 = 10 + restart + 10): FAIL"; [[ -f $d/compare.log ]] && sed 's/^/    /' "$d/compare.log" | head -6
  status=1
fi

# 4. a configuration outside the scope stops at set-up
r=runs/reject
if mr_prepare "$r" cases/mr1d-acu "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Numerics|cfl|1.2" \
     "MOSE-Parameters|iter-threshold|2"; then
  ( cd "$r" && "$MR_MASTER" > solve.log 2>&1 ); rc=$?
  if [[ $rc -ne 0 && ! -f $r/OUTPUT/field.vtm ]] && command grep -a -q -F "[MOSE-Multirate] cfl" "$r/solve.log"; then
    echo "  cfl 1.2 with max-level 3: stopped at set-up (exit $rc)"
  else
    echo "  cfl 1.2 with max-level 3: NOT stopped as expected (exit $rc)"; tail -3 "$r/solve.log"; status=1
  fi
else
  echo "  cfl 1.2 with max-level 3: cannot prepare"; status=1
fi

# 5. a run that diverges stops at the state check of the stage
r=runs/diverge
if mr_prepare "$r" cases/mr1d-sod "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
     "MOSE-Numerics|cfl|20" "MOSE-Parameters|iter-threshold|20"; then
  ( cd "$r" && "$MR_MASTER" > solve.log 2>&1 ); rc=$?
  if [[ $rc -eq 1 && ! -f $r/OUTPUT/field.vtm ]] &&
     command grep -a -q -F "NaN, infinity or p<0 detected (multirate)" "$r/solve.log"; then
    echo "  cfl 20 with one level: stopped by the state check (exit 1)"
  else
    echo "  cfl 20 with one level: NOT stopped as expected (exit $rc; exit 1, the message of the state check of the"
    echo "  multirate and no solution required)"; tail -3 "$r/solve.log"; status=1
  fi
else
  echo "  cfl 20 with one level: cannot prepare"; status=1
fi

# 6. the lines of the gas entering through 400 faces at the stages: one level against the single-rate step
EG=$HERE/../../checks/boundary-conditions/extrapolation-400/entering_gas.py
d=runs/enter400
ok=1
for m in sr mr; do
  r=$d/$m
  rm -rf "$r"; mkdir -p "$r/INPUT" "$r/OUTPUT"
  ln -s "$(cd "$MR_AIR" && pwd)/phase.txt" "$r/INPUT/phase.txt"
  ln -s "$(cd "$MR_AIR" && pwd)/thermo.dat" "$r/INPUT/thermo.dat"
  ( cd "$r" && python3 -B "$EG" case count ) || { ok=0; continue; }
  printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 20' 'newrun = true' 'res-threshold = 1e-300' '' \
    '[MOSE-Physics]' 'equations = euler' '' \
    '[MOSE-Numerics]' 'cfl = 0.5' 'time-scheme = RK2' 'space-reconstruction = first-order' 'time-accurate = true' \
    'integration-variables = cons' 'riemann-solver = HLLC' '' \
    '[MOSE-IO]' 'sol-format = tecplot ascii' > "$r/input.ini"
  [[ $m == mr ]] && python3 -B "$MR_COMMON/set_ini.py" "$r/input.ini" "MOSE-Multirate|enabled|.true." \
    "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0"
  ( cd "$r" && "$MR_MASTER" > solve.log 2>&1 ) && [[ -f $r/OUTPUT/field.tec ]] || { echo "  run $r: failed"; tail -3 "$r/solve.log"; ok=0; }
  command grep -a -E '^ 400 boundary of block|400 faces with gas entering' "$r/solve.log" > "$r/lines400.txt"
done
if [[ $ok -eq 1 ]] && cmp -s "$d/sr/lines400.txt" "$d/mr/lines400.txt" &&
   [[ $(command grep -a -c -E '^ 400 boundary of block 1, face 1, gas entering at iteration 1, faces: 1, ' "$d/mr/lines400.txt") -eq 1 &&
      $(command grep -a -c -E '^ 400 boundary of block' "$d/mr/lines400.txt") -eq 1 &&
      $(command grep -a -c -E '400 faces with gas entering, flux evaluations: 40, ' "$d/mr/lines400.txt") -eq 1 ]]; then
  echo "  400 faces with gas entering, one level: the lines of the single-rate step (one line, 40 evaluations)"
else
  echo "  400 faces with gas entering, one level: FAIL (single rate, then multirate)"
  sed 's/^/    /' "$d/sr/lines400.txt" "$d/mr/lines400.txt" 2>/dev/null | head -8; status=1
fi

if [[ $status -eq 0 ]]; then echo "mr-smoke: ok"; else echo "mr-smoke: FAIL"; fi
exit $status
