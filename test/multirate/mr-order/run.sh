#!/usr/bin/env bash
#===============================================================================
#  Multirate temporal order with levels. Levels fixed from the initial state with H0 = 2^L dt_min / 2; fixed macro steps
#  up to T = 64 H0 against a single rate (max-level 0, safety 1) with a fixed step. Two forms:
#  - short (the default, validation tier): the macro steps H0 and H0/2, the reference step dt_min/16;
#  - full (MR_FORM=full, label verification): the macro steps H0, H0/2, H0/4 and H0/8, the reference step dt_min/64.
#  In both the reference step is a quarter of the finest sub-step of the run with the smallest H. Decided on the last
#  pair of the form: L1 >= 1.9; interface Linf >= 1.8 (the cells within 4 cells of a face between levels, floor 5 % of
#  the global Linf). Flux limiter: van Leer (MOSE has no unlimited MUSCL: without flux-limiter, MUSCL runs van Leer).
#  - 1D periodic acoustic and entropy waves on the banded mesh (levels 0..3, L = 3): rho, u, p; rho only for the entropy
#    wave, whose u and p must stay uniform to 1e-12;
#  - isentropic vortex carried through three blocks in x, the central one finer by 4 (levels 0 and 1, L = 1; the
#    connections take the level of their finer side): rho, u, v, p, interface Linf with the cells next to the
#    connections;
#  - laminar Navier-Stokes: flat plate started from the uniform state, levels from the wall rows (0..3, L = 3), and
#    the skewed periodic hot square with viscosity 0.3 Pa s, square of level 1 with sharp corners (L = 1): rho, u, v,
#    p; on the hot square the Linf over the cells next to the faces between levels (1-cell mask,
#    same 5 % floor) decides too, last pair >= 1.8: there the faces between levels read their coarser neighbours,
#    diagonals included (a halo without the diagonals drops it to about 1).
#  - a no-slip adiabatic wall on smoothly skewed cells whose levels change along the wall (mr2d-wall-s: x periodic,
#    bands of 8 / 2.8 / 1.28 / 2.8 / 8 mm on the levels 0 | 1 | 2 | 1 | 0 along the wall, L = 2, skew angle at the
#    wall 14 to 37 degrees, viscosity 0.3 Pa s): rho, u, v, p; the Linf on the cells next to the faces between levels
#    (1-cell mask, same 5 % floor) decides too, last pair >= 1.8: a wall face of a finer cell reads the ghost cells of
#    its coarser neighbour along the wall, which the finer stages refresh at their own times (without them
#    the order there drops).
#  - Spalart-Allmaras and SST: the wall channel of the test TimeIntegration/TurbulenceSources with each model, from the
#    turbulence profiles of that test (twice the near-wall equilibrium for SA), levels from the wall rows (0..3, L = 3):
#    rho, u, v, p and the turbulence variables (rhoNuTilde; rhoK, rhoOmega), L1 >= 1.8 and interface Linf >= 1.8.
#  Every step stays below slab thickness / a (wall mode across the slab), checked from the log of the H0 runs.
#  Every case is verified even when an earlier one failed.
#  Figures for the V&V page, from the full form (plot_order.py, skipped without matplotlib; the test does not depend on them):
#  runs/multirate-order.svg, the density error of every case against H with the slope of second order, and
#  runs/multirate-wall-s-levels.svg, the level map of mr2d-wall-s on its skewed mesh.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
case ${MR_FORM:-short} in
  short) FORM=short; KS=(0 1); REFF=0.0625; REFN=16 ;;        # H0, H0/2; reference step dt_min/16
  full)  FORM=full; KS=(0 1 2 3); REFF=0.015625; REFN=64 ;;  # H0 ... H0/8; reference step dt_min/64
  *) echo "mr-order: MR_FORM must be short or full"; exit 2 ;;
esac
rm -rf cases runs
CASES=(mr1d-acu mr1d-ent mr2d-vortex-3b mr2d-plate mr2d-hot mr2d-wall-s mr2d-chan-sa mr2d-chan-sst)
[[ $# -gt 0 ]] && CASES=("$@")           # a subset of the cases (development runs); ctest runs them all
mr_cases cases "${CASES[@]}" || { echo "mr-order: FAIL - cannot build the cases"; exit 1; }
status=0
# case L H0/dt_min vars extra-arguments-of-verify_order
for spec in "mr1d-acu 3 4.0 rho(1),u,p" "mr1d-ent 3 4.0 rho(1)" "mr2d-vortex-3b 1 1.0 rho(1),u,v,p --case cases/mr2d-vortex-3b" \
            "mr2d-plate 3 4.0 rho(1),u,v,p" "mr2d-hot 1 1.0 rho(1),u,v,p --case cases/mr2d-hot --min-int1 1.8" \
            "mr2d-wall-s 2 2.0 rho(1),u,v,p --min-int1 1.8" \
            "mr2d-chan-sa 3 4.0 rho(1),u,v,p,rhoNuTilde --min-l1 1.8" "mr2d-chan-sst 3 4.0 rho(1),u,v,p,rhoK,rhoOmega --min-l1 1.8"; do
  set -- $spec; c=$1; L=$2; h0=$3; vars=$4; shift 4
  [[ " ${CASES[*]} " == *" $c "* ]] || continue
  MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|fixed")
  st=0
  runs=()
  for k in "${KS[@]}"; do
    f=$(python3 -c "print($h0 / 2**$k)"); n=$((64 * 2**k))
    tr=(); [[ $k -eq 0 ]] && tr=("MOSE-Multirate|trace|.true.")      # level map in the log (MR-LEVEL lines)
    mr_run "runs/$c-h$k" "cases/$c" "${MR[@]}" "MOSE-Multirate|fixed-macro-step|$f" "MOSE-Parameters|iter-threshold|$n" "${tr[@]}" || st=1
    runs+=("runs/$c-h$k")
  done
  nref=$(python3 -c "print(int(round(64 * $h0 * $REFN)))")
  mr_run "runs/$c-ref" "cases/$c" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
         "MOSE-Multirate|fixed-macro-step|$REFF" "MOSE-Parameters|iter-threshold|$nref" || st=1
  echo "  $c ($FORM form):"
  if [[ $st -eq 0 ]]; then
    python3 -B ../mr-common/verify_order.py --runs "${runs[@]}" --ref "runs/$c-ref" --levels "runs/$c-h0" \
            --vars "$vars" --min-l1 1.9 --min-int 1.8 "$@" --json "runs/$c-order.json" || st=1
    # the entropy wave keeps u and p uniform: their error is round-off and has no order; they are checked to stay uniform
    # instead
    if [[ $c == mr1d-ent ]]; then
      python3 -B ../mr-common/verify_uniform.py "cases/$c" "${runs[@]}" "runs/$c-ref" --vars u,p --max 1e-12 || st=1
    fi
    python3 -B ../mr-common/verify_slab.py "cases/$c" "runs/$c-h0" --h-max || st=1
  fi
  [[ $st -eq 0 ]] || status=1
done
[[ $FORM == full ]] && { python3 -B plot_order.py runs || echo "  plot_order: figures not written"; }
if [[ $status -eq 0 ]]; then echo "mr-order ($FORM form): ok"; else echo "mr-order ($FORM form): FAIL"; fi
exit $status
