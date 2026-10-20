#!/usr/bin/env bash
#===============================================================================
#  Time-dependent boundary data with levels: 1D duct with the inlet on the
#  finest cells (level 3 of 0..3 with safety 2) and a back-pressure outlet; inlet 402 with p0(t) = 1.2e5 Pa
#  (1 + 0.1 sin(2 pi t / Tp)) (mr1d-402), or inlet 410 with the assigned state of a right-running acoustic wave of
#  relative amplitude 0.1 tabulated in a Q2D file (mr1d-410); Tp = 6.4e-5 s, about 16 H0. Levels fixed from the initial
#  state, H0 = 2^3 dt_min / 2, fixed macro steps H0 ... H0/8 up to T = 64 H0, reference single rate with dt_min/64 (its
#  boundary data at the stage times too). Decided on the last pair: Linf over the 4 inlet cells (floor: 5 % of the
#  global Linf) >= 1.8 for rho, u, p; L1 reported. The boundary data of a stage are read at the stage time of the level
#  of the inlet cells: with the time t^(n+1) at every stage the order would drop to 1.
#  The order compares runs that share the data, so verify_inlet.py checks that the data act: u > 0 in the 4 inlet cells,
#  and the first cell against the table at the same time, at the peak of the data (4 macro steps H0, about Tp/4) within
#  half of the distance the table has moved from t = 0, the midpoint between a cell that follows the data (a lag of half
#  a cell over u + a, second order at a peak) and one that stays at their first value, and at T within the amplitude of
#  the table, since an error above the signal would not measure the data.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-402 mr1d-410 || { echo "mr-bc: FAIL - cannot build the cases"; exit 1; }
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|fixed")
status=0
for c in mr1d-402 mr1d-410; do
  st=0
  runs=()
  for k in 0 1 2 3; do
    f=$(python3 -c "print(4.0 / 2**$k)"); n=$((64 * 2**k))
    tr=(); [[ $k -eq 0 ]] && tr=("MOSE-Multirate|trace|.true.")      # level map in the log (MR-LEVEL lines)
    mr_run "runs/$c-h$k" "cases/$c" "${MR[@]}" "MOSE-Multirate|fixed-macro-step|$f" "MOSE-Parameters|iter-threshold|$n" "${tr[@]}" || st=1
    runs+=("runs/$c-h$k")
  done
  mr_run "runs/$c-ref" "cases/$c" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
         "MOSE-Multirate|fixed-macro-step|0.015625" "MOSE-Parameters|iter-threshold|16384" || st=1
  mr_run "runs/$c-peak" "cases/$c" "${MR[@]}" "MOSE-Multirate|fixed-macro-step|4.0" "MOSE-Parameters|iter-threshold|4" || st=1
  echo "  $c:"
  if [[ $st -eq 0 ]]; then
    python3 -B ../mr-common/verify_order.py --runs "${runs[@]}" --ref "runs/$c-ref" --levels "runs/$c-h0" \
            --vars "rho(1),u,p" --inlet 4 --l1-informative --min-int 1.8 || st=1
    python3 -B ../mr-common/verify_slab.py "cases/$c" "runs/$c-h0" --h-max || st=1
    python3 -B verify_inlet.py "cases/$c" "runs/$c-peak" "runs/$c-h0" || st=1
    # information: with a fixed macro step the waves that enter can take a cell above its limit dtcap/sigma (warning)
    for r in "${runs[@]}"; do
      w=$(grep -a -m1 'exceeds the limit' "$r/solve.log" | sed 's/.*= *//')
      echo "  $(basename "$r"): largest h/(dtcap/sigma) above 1: ${w:-no}"
    done
  fi
  [[ $st -eq 0 ]] || status=1
done
if [[ $status -eq 0 ]]; then echo "mr-bc: ok"; else echo "mr-bc: FAIL"; fi
exit $status
