#!/usr/bin/env bash
#===============================================================================
#  Multirate conservation: safety 2, max-level 3, 400 macro
#  steps. Periodic 1D acoustic and entropy waves and the acoustic wave in three blocks (mass, momentum, energy); closed
#  1D Sod tube, 2D box and 2D box in 2 x 2 blocks (mass, energy): relative drift <= 1e-12 (the single-rate drift over
#  400 steps is reported). Positive control: Sod tube with the conservative correction off -> drift >= 1e-9 on mass or
#  energy, otherwise the criterion cannot see a defect. Two more runs (drift <= 1e-12): levels
#  recomputed every 5 macro steps (Sod tube) and a macro-step target below the adaptive macro step (acoustic wave).
#  The steps of the first macro step stay below slab thickness / a (wall mode across the slab), checked from the
#  logs. Laminar Navier-Stokes: the skewed periodic hot square in one block and in 2 x 2 blocks (max-level 1, the levels
#  recomputed every 5 macro steps; mass, x and y momentum, energy), with the conservative correction of the diffusive
#  fluxes; its 3D twin, the hot cube (mass, x, y and z momentum, energy). Spalart-Allmaras and
#  SST: the wall channel of the test TimeIntegration/TurbulenceSources with each model, and with Spalart-Allmaras on
#  a rough wall (max-level 3, levels recomputed every 5 macro steps; mass and energy: periodic along the flow, adiabatic
#  wall, symmetry plane). diffusive-metric = mean: the channel with Spalart-Allmaras (max-level 3) and the hot cube
#  (max-level 1), with the levels recomputed every 5 macro steps.
#  The 2-D axisymmetric wedge with levels (mr2d-wedge-band: laminar, max-level 2, levels 2 | 2 | 1 | 0 from the axis,
#  closed; mass and energy), with its own positive control (the conservative correction off -> drift >= 1e-9). A depth
#  of the map that changes from one macro step to the next: the blast tube mr1d-blast with max-level 10 (the depth
#  follows the octaves of the local steps, between 2 and 3 as the front moves), 200 macro steps, mass and energy, at
#  least two depths printed, and cells at level 0 in every macro step. A frozen mixture of four species, laminar with a
#  Schmidt number (mr1d-mix, max-level 3, the levels recomputed every 5 macro steps, 400 macro steps): the mass of each
#  species, the total mass and the momentum, with its own positive control; the energy of the checks holds for a perfect gas.
#  Turbulence options and a safety below 2, with levels (max-level 3, levels recomputed every 5 macro steps): SA-R and
#  sst-production = compressible on the wall channels (mass, energy), k-energy-coupling on the SST channel (mass; with k
#  in the total energy the sources of k enter the energy row, and its drift is reported), safety 1.5 on the acoustic
#  wave in three blocks (mass, momentum, energy); each option acts: the run of SA-R differs from that of SA on the same
#  channel with the same keys, and those of sst-production = compressible and k-energy-coupling from that of SST.
#  An SRM grain in its inflow mode (502, mode 2) at the end of a duct
#  closed by a slip wall, the gas at rest at the start (mr1d-502c, max-level 3, 400 macro steps): with the burn-rate
#  exponent n = 0 the grain lets in mdot_ref at every pressure, so the mass of the domain is M(0) + mdot_ref A t to
#  1e-12 (verify_grain.py), with its own positive control (the conservative correction off -> difference >= 1e-9).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-acu mr1d-ent mr1d-sod mr2d-box mr1d-acu-3b mr2d-box-4b mr2d-hot mr2d-hot-4b mr3d-hot \
               mr2d-chan-sa mr2d-chan-sst mr2d-chanrough-sa mr2d-wedge-band mr1d-blast mr1d-mix mr1d-502c \
  || { echo "mr-conservation: FAIL - cannot build the cases"; exit 1; }
MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-Parameters|iter-threshold|400")
status=0
for cq in "mr1d-acu mass,momx,energy" "mr1d-ent mass,momx,energy" "mr1d-sod mass,energy" "mr2d-box mass,energy" \
          "mr1d-acu-3b mass,momx,energy" "mr2d-box-4b mass,energy"; do
  set -- $cq; c=$1; q=$2
  mr_run "runs/$c-mr" "cases/$c" "${MR[@]}" || { status=1; continue; }
  mr_run "runs/$c-sr" "cases/$c" "MOSE-Parameters|iter-threshold|400" || { status=1; continue; }
  python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-mr" --quantities "$q" --max 1e-12 || status=1
  python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-sr" --quantities "$q" --report
  python3 -B ../mr-common/verify_slab.py "cases/$c" "runs/$c-mr" --h-max || status=1
done
for cq in "mr2d-hot 1 mass,momx,momy,energy" "mr2d-hot-4b 1 mass,momx,momy,energy" "mr3d-hot 1 mass,momx,momy,momz,energy" \
          "mr2d-chan-sa 3 mass,energy" "mr2d-chan-sst 3 mass,energy" "mr2d-chanrough-sa 3 mass,energy"; do
  set -- $cq; c=$1; L=$2; q=$3
  MRV=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" "MOSE-Multirate|safety|2.0" "MOSE-Parameters|iter-threshold|400"
       "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|5")
  mr_run "runs/$c-mr" "cases/$c" "${MRV[@]}" || { status=1; continue; }
  mr_run "runs/$c-sr" "cases/$c" "MOSE-Parameters|iter-threshold|400" || { status=1; continue; }
  python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-mr" --quantities "$q" --max 1e-12 || status=1
  python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-sr" --quantities "$q" --report
  python3 -B ../mr-common/verify_slab.py "cases/$c" "runs/$c-mr" --h-max || status=1
done
echo "  diffusive-metric = mean with levels:"
for cq in "mr2d-chan-sa 3 mass,energy" "mr3d-hot 1 mass,momx,momy,momz,energy"; do
  set -- $cq; c=$1; L=$2; q=$3
  if mr_run "runs/$c-mean" "cases/$c" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" "MOSE-Multirate|safety|2.0" \
       "MOSE-Parameters|iter-threshold|400" "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|5" \
       "MOSE-Numerics|diffusive-metric|mean"; then
    python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-mean" --quantities "$q" --max 1e-12 || status=1
  else status=1; fi
done
echo "  turbulence options and a safety below 2, with levels:"
# SA-R, sst-production = compressible and k-energy-coupling on the wall channels, and safety 1.5 on the acoustic wave in
# three blocks: max-level 3, levels recomputed every 5 macro steps, 400 macro steps. With k in the total energy the
# sources of k enter the energy row, so there the energy drift is reported, not checked
for spec in "sar mr2d-chan-sa mass,energy - MOSE-Physics|turbulence|SA-R" \
            "sstcomp mr2d-chan-sst mass,energy - MOSE-Turbulence|sst-production|compressible" \
            "sstkec mr2d-chan-sst mass energy MOSE-Turbulence|k-energy-coupling|.true." \
            "safety15 mr1d-acu-3b mass,momx,energy - MOSE-Multirate|safety|1.5"; do
  read -r lab c q qr kv <<< "$spec"
  if mr_run "runs/$c-$lab" "cases/$c" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" \
       "MOSE-Parameters|iter-threshold|400" "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|5" "$kv"; then
    python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-$lab" --quantities "$q" --max 1e-12 || status=1
    if [[ $qr != - ]]; then
      python3 -B ../mr-common/verify_conservation.py "cases/$c" "runs/$c-$lab" --quantities "$qr" --report
    fi
    # the option acts: the run without it, the channel of the loop above with the same keys, gives another result; with
    # k-energy-coupling the stress of the fluxes changes too (Stress_Vector_2eq), so there the check does not isolate the
    # copy of the sources, which mr-onelevel compares with its original (verify_copies.py)
    if [[ $c != mr1d-acu-3b ]]; then
      if mr_same "runs/$c-mr" "runs/$c-$lab" > /dev/null 2>&1; then
        echo "  $c-$lab: the same result as $c without the option: the option does not act"; status=1
      else
        echo "  $c-$lab: a result other than that of $c without the option (the option acts)"
      fi
    fi
  else status=1; fi
done
echo "  the SRM grain in its inflow mode with levels (mr1d-502c, max-level 3, 400 macro steps):"
if mr_run runs/mr1d-502c-mr cases/mr1d-502c "${MR[@]}"; then
  python3 -B verify_grain.py cases/mr1d-502c runs/mr1d-502c-mr --max 1e-12 || status=1
  grep -a -m1 'Multirate: cells per level' runs/mr1d-502c-mr/solve.log | sed 's/^ */  /'
else status=1; fi
if mr_run runs/mr1d-502c-nocorr cases/mr1d-502c "${MR[@]}" "MOSE-Multirate|conservative-correction|.false."; then
  echo "  positive control on the grain (correction off):"
  python3 -B verify_grain.py cases/mr1d-502c runs/mr1d-502c-nocorr --min 1e-9 \
    || { echo "  positive control on the grain: difference below 1e-9, the criterion cannot see a missing correction"; status=1; }
else status=1; fi
mr_run runs/mr1d-sod-nocorr cases/mr1d-sod "${MR[@]}" "MOSE-Multirate|conservative-correction|.false." || status=1
echo "  positive control (correction off):"
python3 -B ../mr-common/verify_conservation.py cases/mr1d-sod runs/mr1d-sod-nocorr --quantities mass,energy --min 1e-9 --any \
  || { echo "  positive control: drift below 1e-9, the criterion cannot see a missing correction"; status=1; }
echo "  levels recomputed every 5 macro steps (Sod tube) and macro-step target 3e-6 s (acoustic wave, adaptive H about 5.2e-6 s):"
if mr_run runs/mr1d-sod-every5 cases/mr1d-sod "${MR[@]}" "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|5"; then
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-sod runs/mr1d-sod-every5 --quantities mass,energy --max 1e-12 || status=1
else status=1; fi
if mr_run runs/mr1d-acu-target cases/mr1d-acu "${MR[@]}" "MOSE-Multirate|macro-step-target|3.0e-6"; then
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-acu runs/mr1d-acu-target --quantities mass,momx,energy --max 1e-12 || status=1
  grep -a -m1 'Multirate: H =' runs/mr1d-acu-target/solve.log | sed 's/^ */  /'
else status=1; fi
echo "  the 2-D axisymmetric wedge with levels (mr2d-wedge-band, max-level 2, 400 macro steps):"
MRW=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|2" "MOSE-Multirate|safety|2.0" "MOSE-Parameters|iter-threshold|400")
if mr_run runs/mr2d-wedge-band-mr cases/mr2d-wedge-band "${MRW[@]}"; then
  python3 -B ../mr-common/verify_conservation.py cases/mr2d-wedge-band runs/mr2d-wedge-band-mr --quantities mass,energy \
    --max 1e-12 || status=1
  grep -a -m1 'Multirate: cells per level' runs/mr2d-wedge-band-mr/solve.log | sed 's/^ */  /'
else status=1; fi
if mr_run runs/mr2d-wedge-band-sr cases/mr2d-wedge-band "MOSE-Parameters|iter-threshold|400"; then
  python3 -B ../mr-common/verify_conservation.py cases/mr2d-wedge-band runs/mr2d-wedge-band-sr --quantities mass,energy --report
else status=1; fi
if mr_run runs/mr2d-wedge-band-nocorr cases/mr2d-wedge-band "${MRW[@]}" "MOSE-Multirate|conservative-correction|.false."; then
  echo "  positive control on the wedge (correction off):"
  python3 -B ../mr-common/verify_conservation.py cases/mr2d-wedge-band runs/mr2d-wedge-band-nocorr --quantities mass,energy \
    --min 1e-9 --any || { echo "  positive control on the wedge: drift below 1e-9, the criterion cannot see a missing correction"; status=1; }
else status=1; fi
echo "  a depth that changes from one macro step to the next (mr1d-blast, max-level 10, 200 macro steps):"
if mr_run runs/mr1d-blast-depth cases/mr1d-blast "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|10" \
     "MOSE-Multirate|safety|2.0" "MOSE-Parameters|iter-threshold|200" "MOSE-IO|shell-diter|1"; then
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-blast runs/mr1d-blast-depth --quantities mass,energy --max 1e-12 \
    || status=1
  depths=$(command grep -a 'Multirate: H =' runs/mr1d-blast-depth/solve.log | sed -E 's/.*L = +([0-9]+),.*/\1/' | sort -n | uniq | tr '\n' ' ')
  if [[ $(wc -w <<< "$depths") -ge 2 ]]; then echo "  mr1d-blast-depth: depths printed: $depths"
  else echo "  mr1d-blast-depth: one depth only ($depths): the depth does not change in the run"; status=1; fi
  # no coarse level left empty: cells at level 0 in every macro step
  n0=$(command grep -a -c 'Multirate: cells per level' runs/mr1d-blast-depth/solve.log)
  z0=$(command grep -a 'Multirate: cells per level' runs/mr1d-blast-depth/solve.log | awk '$5 + 0 == 0' | wc -l)
  if [[ $n0 -gt 0 && $z0 -eq 0 ]]; then echo "  mr1d-blast-depth: level 0 occupied in the $n0 macro steps"
  else echo "  mr1d-blast-depth: $z0 of $n0 macro steps with level 0 empty"; status=1; fi
else status=1; fi
echo "  a frozen mixture of four species with levels (mr1d-mix, max-level 3, 400 macro steps):"
MRM=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-Parameters|iter-threshold|400"
     "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|5")
if mr_run runs/mr1d-mix-mr cases/mr1d-mix "${MRM[@]}"; then
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-mix runs/mr1d-mix-mr --quantities mass,momx,mass1,mass2,mass3,mass4 \
    --max 1e-12 || status=1
  grep -a -m1 'Multirate: cells per level' runs/mr1d-mix-mr/solve.log | sed 's/^ */  /'
else status=1; fi
if mr_run runs/mr1d-mix-sr cases/mr1d-mix "MOSE-Parameters|iter-threshold|400"; then
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-mix runs/mr1d-mix-sr --quantities mass,momx,mass1,mass2,mass3,mass4 --report
else status=1; fi
if mr_run runs/mr1d-mix-nocorr cases/mr1d-mix "${MRM[@]}" "MOSE-Multirate|conservative-correction|.false."; then
  echo "  positive control on the mixture (correction off):"
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-mix runs/mr1d-mix-nocorr --quantities mass1,mass2,mass3,mass4 \
    --min 1e-9 || { echo "  positive control on the mixture: a species with a drift below 1e-9, the criterion cannot see a missing correction for it"; status=1; }
else status=1; fi
if [[ $status -eq 0 ]]; then echo "mr-conservation: ok"; else echo "mr-conservation: FAIL"; fi
exit $status
