#!/usr/bin/env bash
#===============================================================================
#  Multirate with a fixed macro step: with fixed-macro-step = f and levels = recompute (the default) the level map is
#  built for the macro step taken, H = f dt_min of the first macro step, and not for 2^L dt_min / safety. With
#  levels = fixed the map is that of the adaptive step, kept for the studies of the temporal order (mr-order, mr-bc,
#  mr-deep), and with a larger H the finest cells run beyond their limit, which the run reports with a warning
#  (mr-reject, runs 34 and 34b). 1D acoustic wave and 2D box with a pressure pulse,
#  max-level 3, safety 2, one macro step with trace = .true.:
#   A  fixed-macro-step f = 1 and f = 2, no macro-step-target;
#   C  the adaptive macro step with macro-step-target = f dt_first, where dt_first, the dt_min of the first step, is
#      the TIME of a one-step run on the one-level path (max-level 0, safety 1, fixed-macro-step 1), written with
#      17 significant digits (exact round trip).
#  The level map of the first macro step (MR-LEVEL lines) and the lines "Multirate: H" and "cells per level" of A and
#  C must be identical, and A must report no cell beyond its limit.
#  D  the depth of the map of the adaptive step: L = min(max-level, floor(K)), K the octaves between the largest and the
#     smallest local step of the mesh, so max-level above floor(K) changes nothing (20 macro steps, outputs and residual
#     history byte by byte, every line "Multirate: H" with that L and H = 2^L/safety dt_min, and every line "cells per
#     level" with cells at level 0): mr1d-acu (K = 3.97) with max-level 10 and with max-level left out (the default, 10)
#     against max-level 3; mr1d-deep12 (K above 12) with max-level left out against max-level 10; mr1d-simple-64 (K below
#     1) with max-level 10 against max-level 0 (one level), safety 2 in both.
#  E  the same depth with levels = fixed and with macro-step-target, mr1d-acu, 20 macro steps: levels = fixed with
#     fixed-macro-step 1, max-level 10 against max-level 3 (the map kept is the one with depth 3, L = 3); a target of
#     7 dt_first, within the step of the slowest cells over safety (7.8 dt_min), leaves the map of the target (the cap of
#     the depth does not act): max-level 10 against max-level 4 with the same target (L = 4, the depth the target gives);
#     a target of 20 dt_first, above it: L = 3 and H = 4 dt_min, as max-level 3 without a target. Outputs and residual
#     history byte by byte. With levels = fixed a fixed macro step above the step of that map runs on with the finest
#     cells beyond their limit, which it reports (mr-reject, runs 34 and 34b).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-acu mr2d-box || { echo "mr-fixedstep: FAIL - cannot build the cases"; exit 1; }
status=0
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|trace|.true."
    "MOSE-Parameters|iter-threshold|1")
# lines of the first macro step that define the map and the step
map_lines() { command grep -a -E '^MR-LEVEL|Multirate: H =|Multirate: cells per level' "$1/solve.log"; }
for c in mr1d-acu mr2d-box; do
  mr_run "runs/$c-dt" "cases/$c" "MOSE-Parameters|iter-threshold|1" "MOSE-Multirate|enabled|.true." \
    "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" "MOSE-Multirate|fixed-macro-step|1.0" || { status=1; continue; }
  dt=$(python3 -B -c "import sys; sys.path.insert(0, '../mr-common'); import mrcheck; print('%.17e' % mrcheck.read_run('runs/$c-dt')[0]['time'])") \
    || { echo "  $c: cannot read dt_first"; status=1; continue; }
  for f in 1 2; do
    tgt=$(python3 -c "print('%.17e' % ($f * float('$dt')))")
    mr_run "runs/$c-A$f" "cases/$c" "${MR[@]}" "MOSE-Multirate|fixed-macro-step|$f.0" || { status=1; continue; }
    mr_run "runs/$c-C$f" "cases/$c" "${MR[@]}" "MOSE-Multirate|macro-step-target|$tgt" || { status=1; continue; }
    nA=$(map_lines "runs/$c-A$f" | command grep -a -c '^MR-LEVEL [0-9]')
    if [[ $nA -gt 0 ]] && cmp -s <(map_lines "runs/$c-A$f") <(map_lines "runs/$c-C$f"); then
      echo "  $c f = $f: map of the fixed step = map of the adaptive step with H = $tgt s ($nA cells): identical"
    else
      echo "  $c f = $f: map of the fixed step and of the adaptive step DIFFERENT"
      diff <(map_lines "runs/$c-A$f" | command grep -a 'Multirate') <(map_lines "runs/$c-C$f" | command grep -a 'Multirate') | sed 's/^/    /'
      status=1
    fi
    if command grep -a -q 'exceeds the limit of some cell' "runs/$c-A$f/solve.log"; then
      echo "  $c f = $f: a cell beyond its limit with the fixed step"; status=1
    fi
  done
done
echo "  D: the depth of the adaptive map, min(max-level, floor(K))"
mr_cases cases mr1d-deep12 mr1d-simple-64 || { echo "mr-fixedstep: FAIL - cannot build the cases of D"; exit 1; }
DK=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|safety|2.0" "MOSE-Parameters|iter-threshold|20" "MOSE-IO|res-diter|1"
    "MOSE-IO|shell-diter|1")
# depth LABEL RUN L HDT -> 0 if every line "Multirate: H" of the run (one per macro step) has L and H = HDT dt_min (HDT as
# the line prints it, ES10.3)
depth() {
  local n bad hdt
  hdt=$(printf '%s' "$4" | sed 's/[.+]/\\&/g')
  n=$(command grep -a -c 'Multirate: H =' "$2/solve.log")
  bad=$(command grep -a 'Multirate: H =' "$2/solve.log" | command grep -a -v -E "= +$hdt dt_min, L = +$3," | wc -l)
  if [[ $n -gt 0 && $bad -eq 0 ]]; then
    echo "  $1: $n macro steps, L = $3, H = $4 dt_min"
  else
    echo "  $1: $bad of $n macro steps without L = $3 and H = $4 dt_min: $(command grep -a -m1 'Multirate: H =' "$2/solve.log" | sed 's/^ *//')"
    return 1
  fi
}
# pair LABEL CASE KEYS_A... -- KEYS_B... -> the runs a and b of CASE, their outputs identical byte by byte
pair() {
  local lab=$1 case=$2 ka=() kb=()
  shift 2
  while [[ $# -gt 0 && $1 != -- ]]; do ka+=("$1"); shift; done
  shift
  kb=("$@")
  mr_run "runs/$lab-a" "cases/$case" "${DK[@]}" "${ka[@]}" || return 1
  mr_run "runs/$lab-b" "cases/$case" "${DK[@]}" "${kb[@]}" || return 1
  if mr_same "runs/$lab-a" "runs/$lab-b"; then echo "  $lab: identical (fields and residual history)"; else echo "  $lab: DIFFERENT"; return 1; fi
}
# level0 LABEL RUN -> 0 if every line "Multirate: cells per level" of the run has cells at level 0: no coarse level left
# empty, the step of the map within the step of its slowest cells
level0() {
  local n z
  n=$(command grep -a -c 'Multirate: cells per level' "$2/solve.log")
  z=$(command grep -a 'Multirate: cells per level' "$2/solve.log" | awk '$5 + 0 == 0' | wc -l)
  if [[ $n -gt 0 && $z -eq 0 ]]; then
    echo "  $1: level 0 occupied in the $n macro steps"
  else
    echo "  $1: $z of $n macro steps with level 0 empty: $(command grep -a -m1 'Multirate: cells per level' "$2/solve.log" | sed 's/^ *//')"
    return 1
  fi
}
# lev LABEL RUN L -> 0 if every line "Multirate: H" of the run has L
lev() {
  local n bad
  n=$(command grep -a -c 'Multirate: H =' "$2/solve.log")
  bad=$(command grep -a 'Multirate: H =' "$2/solve.log" | command grep -a -v -E "dt_min, L = +$3," | wc -l)
  if [[ $n -gt 0 && $bad -eq 0 ]]; then
    echo "  $1: $n macro steps, L = $3"
  else
    echo "  $1: $bad of $n macro steps without L = $3: $(command grep -a -m1 'Multirate: H =' "$2/solve.log" | sed 's/^ *//')"
    return 1
  fi
}
pair D-acu10 mr1d-acu "MOSE-Multirate|max-level|3" -- "MOSE-Multirate|max-level|10" || status=1
depth D-acu10 runs/D-acu10-b 3 4.000E+00 || status=1
level0 D-acu10 runs/D-acu10-b || status=1
pair D-acudef mr1d-acu "MOSE-Multirate|max-level|3" -- || status=1
depth D-acudef runs/D-acudef-b 3 4.000E+00 || status=1
level0 D-acudef runs/D-acudef-b || status=1
pair D-deep12 mr1d-deep12 "MOSE-Multirate|max-level|10" -- || status=1
depth D-deep12 runs/D-deep12-b 10 5.120E+02 || status=1
level0 D-deep12 runs/D-deep12-b || status=1
pair D-simple mr1d-simple-64 "MOSE-Multirate|max-level|0" -- "MOSE-Multirate|max-level|10" || status=1
depth D-simple runs/D-simple-b 0 5.000E-01 || status=1
level0 D-simple runs/D-simple-b || status=1
echo "  E: the depth of the map with levels = fixed and with macro-step-target"
dt=$(python3 -B -c "import sys; sys.path.insert(0, '../mr-common'); import mrcheck; print('%.17e' % mrcheck.read_run('runs/mr1d-acu-dt')[0]['time'])") \
  || { echo "  E: cannot read dt_first of mr1d-acu"; dt=; status=1; }
if [[ -n $dt ]]; then
  # E-fixed: levels = fixed keeps the map of the adaptive step of the first macro step, with its depth: with fixed-macro-step 1
  #     max-level 10 gives the map and the outputs of max-level 3
  pair E-fixed mr1d-acu "MOSE-Multirate|levels|fixed" "MOSE-Multirate|fixed-macro-step|1.0" "MOSE-Multirate|max-level|3" -- \
    "MOSE-Multirate|levels|fixed" "MOSE-Multirate|fixed-macro-step|1.0" "MOSE-Multirate|max-level|10" || status=1
  lev E-fixed runs/E-fixed-b 3 || status=1
  # E-target7: a target within the step of the slowest cells (7 dt_first, below dtcap_max/safety = 7.8 dt_min) already keeps H
  #     there, so the cap of the depth does not act: the map is that of max-level 4 with the same target (L = 4, H = the target)
  t7=$(python3 -c "print('%.17e' % (7 * float('$dt')))")
  pair E-target7 mr1d-acu "MOSE-Multirate|max-level|4" "MOSE-Multirate|macro-step-target|$t7" -- \
    "MOSE-Multirate|max-level|10" "MOSE-Multirate|macro-step-target|$t7" || status=1
  lev E-target7 runs/E-target7-b 4 || status=1
  # E-target20: a target above it (20 dt_first): the depth from the octaves of the local steps, L = 3 and H = 4 dt_min, as
  #     max-level 3 without a target
  t20=$(python3 -c "print('%.17e' % (20 * float('$dt')))")
  pair E-target20 mr1d-acu "MOSE-Multirate|max-level|3" -- "MOSE-Multirate|max-level|10" "MOSE-Multirate|macro-step-target|$t20" \
    || status=1
  depth E-target20 runs/E-target20-b 3 4.000E+00 || status=1
fi
if [[ $status -eq 0 ]]; then echo "mr-fixedstep: ok"; else echo "mr-fixedstep: FAIL"; fi
exit $status
