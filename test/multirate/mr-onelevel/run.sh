#!/usr/bin/env bash
#===============================================================================
#  Multirate with one level = single-rate step: with max-level 0
#  and safety 1 one macro step is one single-rate Heun step; fields, wall output and residual history identical byte
#  by byte to the single-rate run (100 steps): 1D acoustic wave, 1D Sod tube, 2D box with a pressure pulse; the
#  acoustic wave in three blocks and the box in 2 x 2 blocks; the 1D duct with the time-dependent inlets 402 and
#  410, whose data are read at the stage times (the same stage times as the single-rate step), with the inlets 401,
#  403, 405, 408 and 420, and with an SRM grain in its inflow mode (mr1d-502: a 502 of mode 2, the mass-flux inlet 404
#  at the mass flux of Saint-Robert's law and the static temperature of total T_af at the pressure of the boundary
#  cell), so that every boundary type that the multirate integrator accepts runs here; the wall channel with SA-RC and
#  with SST-RC, whose rotation and curvature terms the block routine of the single-rate stage computes, which the stage
#  calls with one level. Laminar
#  Navier-Stokes: the flat plate in one block and in two blocks (connection 101 across the boundary layer),
#  each with diffusive-metric = mean (and the plate with stencil, the default), and the skewed hot
#  square in 2 x 2 blocks.
#  Its 3D twin, the hot cube on an oblique periodic mesh. Spalart-Allmaras and SST: the
#  wall channel of the test TimeIntegration/TurbulenceSources with each model and with Spalart-Allmaras on a rough wall
#  (SA-rough), the
#  turbulent flat plate with each model, and the plate with Spalart-Allmaras in 2 x 2 blocks (connections 101 across
#  the turbulent boundary layer). With one level the stage calls the block routines of the turbulence sources; the
#  per-cell copies used with more levels are compared with them statement by statement (verify_copies.py), as is the
#  copy of the loop of Connection_Fluxes (one flux per connection face), and the MPI tags of the multirate are checked
#  against those of Mod_GhostExchange (no value shared). The 2-D axisymmetric wedge with the field of part 4 of the
#  test MMSResiduals (viscous stress on the wedge faces). The laminar plate in 2 x 2 blocks (a vertex of four
#  blocks in the boundary layer, the case of mr-hybrid) with the two diffusive metrics. A 1-D slab
#  between symmetry j and k faces, planes normal to y and z, whose i faces are inclined
#  (mr1d-skew: the y and z momentum frozen, with the residual the pressure on the inclined faces gives them), a 1-D duct
#  whose top j and k faces are inclined (mr1d-nozzle: not planes normal to y and z, so the y and z momentum are
#  integrated), each with the set-up lines that say so in both runs, and the plate with the
#  records of bc.txt in reverse order (mr2d-plate-rev: the boundary fluxes of a corner cell summed in the order of the
#  faces 1..6, as BC_Fluxes), and the plate with an isothermal wall (mr2d-plate-302). The copy of the select of
#  BC_Fluxes (BC_Flux_Record) is compared with it branch by branch by verify_copies.py. Every run writes its residual at every step
#  (res-diter = 1), so the residual histories compared are not empty.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
CASES=(mr1d-acu mr1d-sod mr2d-box mr1d-acu-3b mr2d-box-4b mr1d-402 mr1d-410 mr2d-plate mr2d-plate-2b mr2d-hot-4b
       mr3d-hot mr2d-chan-sa mr2d-chan-sst mr2d-chanrough-sa mr2d-plate-sa mr2d-plate-sst mr2d-plate-sa-4b mr2d-wedge-visc
       mr2d-plate-4b mr1d-nozzle mr1d-skew mr2d-plate-rev mr1d-401 mr1d-403 mr1d-405 mr1d-408 mr1d-420 mr1d-502
       mr2d-plate-302)
mr_cases cases "${CASES[@]}" || { echo "mr-onelevel: FAIL - cannot build the cases"; exit 1; }
status=0
# run label, case, extra keys
for spec in "mr1d-acu" "mr1d-sod" "mr2d-box" "mr1d-acu-3b" "mr2d-box-4b" "mr1d-402" "mr1d-410" \
            "mr2d-plate mean" "mr2d-plate-2b mean" "mr2d-hot-4b" \
            "mr3d-hot" "mr2d-chan-sa" "mr2d-chan-sst" "mr2d-chanrough-sa" "mr2d-plate-sa" "mr2d-plate-sst" \
            "mr2d-plate-sa-4b" "mr2d-wedge-visc" "mr2d-plate-4b mean" \
            "mr2d-plate stencil" "mr2d-plate-4b stencil" "mr1d-nozzle" "mr1d-skew" "mr2d-plate-rev" \
            "mr1d-401" "mr1d-403" "mr1d-405" "mr1d-408" "mr1d-420" "mr1d-502" "mr2d-plate-302"; do
  set -- $spec; c=$1; met=${2:-}
  lab=$c; keys=("MOSE-IO|res-diter|1")
  if [[ -n $met ]]; then lab=$c-$met; keys+=("MOSE-Numerics|diffusive-metric|$met"); fi
  mr_run "runs/$lab-sr" "cases/$c" "${keys[@]}" || { status=1; continue; }
  mr_run "runs/$lab-mr" "cases/$c" "${keys[@]}" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" \
    "MOSE-Multirate|safety|1.0" || { status=1; continue; }
  if mr_same "runs/$lab-sr" "runs/$lab-mr"; then echo "  $lab: identical (fields and residual history)"; else echo "  $lab: DIFFERENT"; status=1; fi
done
# the 1-D slabs between symmetry j and k faces: the set-up freezes the y and the z velocity of mr1d-skew (j and k faces
# planes normal to y and z) and keeps those of mr1d-nozzle (inclined top faces), and says so in the log of both runs
for spec in "mr1d-skew|j faces: the y velocity is set to zero|k faces: the z velocity is set to zero|frozen" \
            "mr1d-nozzle|j faces that are not planes normal to y|k faces that are not planes normal to z|integrated"; do
  IFS='|' read -r c wj wk what <<< "$spec"; ok=1
  for r in "runs/$c-sr" "runs/$c-mr"; do
    for w in "$wj" "$wk"; do
      command grep -a -q -F "1-D case with symmetry (300) $w" "$r/solve.log" 2>/dev/null \
        || { echo "  $r: no set-up line \"1-D case with symmetry (300) $w\""; ok=0; }
    done
  done
  if [[ $ok -eq 1 ]]; then echo "  $c: the y and z velocity $what (set-up lines in both runs)"; else status=1; fi
done
# the rotation and curvature correction with one level: the stage calls the block routine of the single-rate stage,
# which computes its terms, so the macro step is the single-rate step to the bit
for spec in "mr2d-chan-sa SA-RC" "mr2d-chan-sst SST-RC"; do
  set -- $spec; c=$1; m=$2; lab=$c-$m
  keys=("MOSE-IO|res-diter|1" "MOSE-Physics|turbulence|$m")
  mr_run "runs/$lab-sr" "cases/$c" "${keys[@]}" || { status=1; continue; }
  mr_run "runs/$lab-mr" "cases/$c" "${keys[@]}" "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" \
    "MOSE-Multirate|safety|1.0" || { status=1; continue; }
  if mr_same "runs/$lab-sr" "runs/$lab-mr"; then echo "  $lab: identical (fields and residual history)"; else echo "  $lab: DIFFERENT"; status=1; fi
done
python3 -B ../mr-common/verify_copies.py ../../../src || status=1
if [[ $status -eq 0 ]]; then echo "mr-onelevel: ok"; else echo "mr-onelevel: FAIL"; fi
exit $status
