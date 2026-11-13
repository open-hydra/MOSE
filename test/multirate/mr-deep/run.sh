#!/usr/bin/env bash
#===============================================================================
#  Deep hierarchies: 1D periodic acoustic wave on meshes whose bands have
#  a step ratio of 2.1 (the level of every band equals its index with safety 2), safety 2:
#  - mr1d-deep: 9 levels (L = 8), bands of 4 cells: order, conservation and stability;
#  - mr1d-deep1: 9 levels, bands of one cell, so the halo holds levels m <= l - 2: order and
#    conservation;
#  - mr1d-deep12-3b: 13 levels (L = 12), bands of 4 cells, in three blocks, one connection on the edge of two bands and
#    one inside the deepest band: order (interface mask with the cells next to the connections), conservation and
#    stability.
#  Order: levels fixed from H0 = 2^L dt_min / 2, macro steps H0 ... H0/8 up to T = 16 H0, reference single rate with
#  dt_min/64: last pair L1 >= 1.9 and interface Linf (floor 5 % of the global Linf) >= 1.8 (rho, u, p). Conservation after 200 macro steps
#  (levels recomputed): mass, momentum, energy drift <= 1e-12. Stability over 2000 macro steps: the run completes and the
#  final amplitude of p - p0 is <= 1.01 x the initial one. Every step stays below slab thickness / a.
#  Two forms. Full (MR_FORM=full, CTest Multirate.DeepFull): the three cases as above. Short (the default, Multirate.Deep):
#  the order of mr1d-deep and mr1d-deep1 on the macro steps H0 and H0/2 against the reference step dt_min/16, as in the
#  short form of mr-order, the conservation of the three cases and the stability of mr1d-deep; the order and the
#  stability with 13 levels are in the full form only.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
case ${MR_FORM:-short} in
  short) FORM=short; KS=(0 1); REFF=0.0625; REFN=16 ;;        # H0, H0/2; reference step dt_min/16
  full)  FORM=full; KS=(0 1 2 3); REFF=0.015625; REFN=64 ;;  # H0 ... H0/8; reference step dt_min/64
  *) echo "mr-deep: MR_FORM must be short or full"; exit 2 ;;
esac
rm -rf cases runs
mr_cases cases mr1d-deep mr1d-deep1 mr1d-deep12-3b || { echo "mr-deep: FAIL - cannot build the cases"; exit 1; }
status=0

# deep_case CASE L ORDER(yes|no) STABILITY(yes|no) [extra arguments of verify_order.py]
deep_case() {
  local c=$1 L=$2 ord=$3 long=$4 k f n nref st=0; shift 4
  local MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" "MOSE-Multirate|safety|2.0")
  local runs=()
  local h0=$(( 2 ** (L - 1) ))                                            # H0 = 2^L dt_min / 2
  if [[ $ord == yes ]]; then
    echo "  $c (L = $L), order ($FORM form):"
    for k in "${KS[@]}"; do
      f=$(python3 -c "print($h0 / 2**$k)"); n=$((16 * 2**k))
      local tr=(); [[ $k -eq 0 ]] && tr=("MOSE-Multirate|trace|.true.")        # level map in the log (MR-LEVEL lines)
      mr_run "runs/$c-h$k" "cases/$c" "${MR[@]}" "MOSE-Multirate|levels|fixed" "MOSE-Multirate|fixed-macro-step|$f" \
             "MOSE-Parameters|iter-threshold|$n" "${tr[@]}" || st=1
      runs+=("runs/$c-h$k")
    done
    nref=$((16 * h0 * REFN))
    mr_run "runs/$c-ref" "cases/$c" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
           "MOSE-Multirate|fixed-macro-step|$REFF" "MOSE-Parameters|iter-threshold|$nref" || st=1
    if [[ $st -eq 0 ]]; then
      python3 -B ../mr-common/verify_order.py --runs "${runs[@]}" --ref "runs/$c-ref" --levels "runs/$c-h0" \
              --vars "rho(1),u,p" --min-l1 1.9 --min-int 1.8 "$@" || st=1
      python3 -B ../mr-common/verify_slab.py "cases/$c" "runs/$c-h0" --h-max || st=1
    fi
  fi
  echo "  $c, conservation, 200 macro steps:"
  if mr_run "runs/$c-cons" "cases/$c" "${MR[@]}" "MOSE-Parameters|iter-threshold|200"; then
    python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-cons" --quantities mass,momx,energy --max 1e-12 || st=1
  else st=1; fi
  if [[ $long == yes ]]; then
    echo "  $c, stability, 2000 macro steps:"
    if mr_run "runs/$c-long" "cases/$c" "${MR[@]}" "MOSE-Parameters|iter-threshold|2000"; then
      python3 -B - "cases/$c" "runs/$c-long" <<'PY' || st=1
import os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, '../mr-common')
import numpy as np, mrcheck as mc
p0 = np.concatenate([b['cell']['p'].ravel() for b in mc.read_ic_tec(os.path.join(sys.argv[1], 'INPUT', 'ic.tec'))])
p1 = np.concatenate([b['cell']['p'].ravel() for b in mc.read_run(sys.argv[2])])
a0, a1 = np.abs(p0 - 1.0e5).max(), np.abs(p1 - 1.0e5).max()
print('  amplitude of p - p0: initial %.6e Pa, final %.6e Pa, ratio %.5f (<= 1.01)' % (a0, a1, a1 / a0))
sys.exit(0 if a1 <= 1.01 * a0 else 1)
PY
    else st=1; fi
  fi
  return $st
}

deep_case mr1d-deep 8 yes yes || status=1
deep_case mr1d-deep1 8 yes no || status=1
if [[ $FORM == full ]]; then
  deep_case mr1d-deep12-3b 12 yes yes --case cases/mr1d-deep12-3b || status=1
else
  deep_case mr1d-deep12-3b 12 no no || status=1
fi
if [[ $status -eq 0 ]]; then echo "mr-deep ($FORM form): ok"; else echo "mr-deep ($FORM form): FAIL"; fi
exit $status
