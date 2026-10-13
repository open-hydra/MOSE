#!/usr/bin/env bash
#===============================================================================
#  Patch test of the viscous fluxes (fast tier)
#
#  Linear fields have a uniform viscous stress and a uniform heat flux, so
#  the viscous force and the heat conduction of every cell are zero: a face
#  gradient exact for linear fields gives zero to round-off on any mesh
#  (make_case.py). Meshes with nodes moved at random, two blocks joined along
#  a fold, jumps of cell size with shear (1.5, and 2 next to the faces with
#  dyadic coordinates, where the ghost cells must not degenerate), a fan of
#  index lines at a face (the edges along it growing 7/4 times across the
#  first cell: the second layer of extrapolated ghost cells, which no flux
#  reads, turned over), points shared by 3 and by 5 blocks,
#  a 3-D pair of blocks, a symmetry plane under skewed cells. One explicit
#  step with the Navier-Stokes equations and one reference step (Euler, or
#  Navier-Stokes with mu = k = 1e-30 next to walls) from the same state: the
#  difference of r/dt in OUTPUT/diagnostic.tec is the viscous right-hand side.
#  It must be below 1e-11 of mu |grad u| (face areas)/volume, or of
#  k |grad T| (face areas)/volume, in the cells inside the blocks, next to the
#  faces between blocks, at the points shared by 3 and 5 blocks, next to
#  symmetry planes and next to walls with a prescribed heat flux (301); also
#  at every stage of a three-stage step at the points shared by 3 and 5 blocks.
#  The cells next to walls of prescribed temperature (302) are the test
#  ViscousWallPatch, also on a wall that does not turn with a rotating frame,
#  whose velocity in the frame varies along it (with every face metric), with
#  a T linear at uniform pressure (the density not uniform) next to walls 301
#  and 302 and 400 faces, whose ghost cells must extrapolate T itself, and a
#  run where the limits of those ghost cells act (counted at the end of the
#  run). Controls that every face metric passes: walls 301 with q = 0 (here)
#  and 302 with T_w uniform along the wall (ViscousWallPatch) on parallelogram
#  cells. Plus the face metrics inverse-mean and mean: they must give the
#  viscous force of an independent replica of their formulas.
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PATCH_TEST_DIR:-$HERE}"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../common/Air-viscous
NAME=${PATCH_TEST_NAME:-ViscousLinearPatch}
MC=${PATCH_TEST_MAKE_CASE:-./make_case.py}

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$VISC/phase.txt"  INPUT/phase.txt
ln -sf "../$VISC/thermo.dat" INPUT/thermo.dat

# step <tag> <mesh> <run> [metric]: one step with input-<tag>.ini and INPUT/bc-<tag>.txt; keep the diagnostic output
step() {
  local tag=$1 mesh=$2 run=$3 metric=${4:-}
  cp "input-$tag.ini" input.ini
  cp "INPUT/bc-$tag.txt" INPUT/bc.txt
  rm -f INPUT/transport.dat
  if [[ $tag == ref && ( $run == ene || $run == enep || $run == momw || $run == momr ) ]]; then
    ln -sf transport-ref.dat INPUT/transport.dat
  else
    ln -sf "../$VISC/transport.dat" INPUT/transport.dat
  fi
  rm -f OUTPUT/diagnostic.tec
  local out="$mesh-$run-$tag${metric:+-$metric}"
  "$LOCAL" > "OUTPUT/$out.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$out.log"; then
    echo "  [FAIL] the step $out did not complete (see OUTPUT/$out.log)"
    echo "$NAME: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diag-$out.tec"
}

status=0
if [[ $NAME == ViscousWallPatch ]]; then
  echo "-- viscous patch test: cells next to walls of prescribed temperature (302), T varying along the wall --"
  for mesh in fold2 jump jump2 star3 star5 box3d sym paraT; do
    python3 -B "$MC" case "$mesh" ene > /dev/null || exit 1
    step ns "$mesh" ene
    step ref "$mesh" ene
    python3 -B "$MC" check "$mesh" ene "OUTPUT/diag-$mesh-ene-ns.tec" "OUTPUT/diag-$mesh-ene-ref.tec" wall302 || status=1
  done
  echo "-- velocity linear and zero on a no-slip wall (302), wall cells with nodes moved at random: their momentum --"
  for mesh in rand1 box3d; do
    python3 -B "$MC" case "$mesh" momw > /dev/null || exit 1
    step ns "$mesh" momw
    step ref "$mesh" momw
    python3 -B "$MC" check "$mesh" momw "OUTPUT/diag-$mesh-momw-ns.tec" "OUTPUT/diag-$mesh-momw-ref.tec" wall302 || status=1
  done
  echo "-- the same on a wall that does not turn with a rotating frame: its velocity there, -Omega x r, varies along it --"
  for mesh in rand1 box3d; do
    python3 -B "$MC" case "$mesh" momr > /dev/null || exit 1
    step ns "$mesh" momr
    step ref "$mesh" momr
    python3 -B "$MC" check "$mesh" momr "OUTPUT/diag-$mesh-momr-ns.tec" "OUTPUT/diag-$mesh-momr-ref.tec" wall302 || status=1
  done
  echo "-- the same on parallelograms, inverse-mean (wall 302) and mean (adiabatic 301); also the work of the stress --"
  for case in paraT:inverse-mean:wall302 paraQ:mean:wall301; do
    IFS=: read -r mesh metric group <<< "$case"
    echo "    $mesh, diffusive-metric = $metric:"
    python3 -B "$MC" case "$mesh" momr "$metric" > /dev/null || exit 1
    step ns "$mesh" momr "$metric"
    step ref "$mesh" momr "$metric"
    python3 -B "$MC" check "$mesh" momr "OUTPUT/diag-$mesh-momr-ns-$metric.tec" "OUTPUT/diag-$mesh-momr-ref-$metric.tec" \
      "$group" || status=1
    python3 -B "$MC" work "$mesh" momr "OUTPUT/diag-$mesh-momr-ns-$metric.tec" "OUTPUT/diag-$mesh-momr-ref-$metric.tec" \
      "interior,$group" || status=1
  done
  echo "-- T linear at uniform pressure (density not uniform): walls 301, 302 and 400 faces extrapolate T itself --"
  for mesh in fold2 jump jump2 star3 star5 box3d sym paraT; do
    python3 -B "$MC" case "$mesh" enep > /dev/null || exit 1
    step ns "$mesh" enep
    step ref "$mesh" enep
    python3 -B "$MC" check "$mesh" enep "OUTPUT/diag-$mesh-enep-ns.tec" "OUTPUT/diag-$mesh-enep-ref.tec" \
      interior,connection,corner,symmetry,wall301,wall302 || status=1
  done
  echo "-- limits: a 302 wall colder than 2/3 of its first cell, the pressure doubling at the 400 faces: counted --"
  python3 -B "$MC" case paraT hot > /dev/null || exit 1
  step ns paraT hot
  python3 -B "$MC" limits OUTPUT/paraT-hot-ns.log || status=1
else
  echo "-- viscous patch test: linear velocity (mom) and linear temperature (ene) --"
  for mesh in fold2 jump jump2 fan star3 star5 box3d sym; do
    for run in mom ene; do
      python3 -B "$MC" case "$mesh" "$run" > /dev/null || exit 1
      step ns "$mesh" "$run"
      step ref "$mesh" "$run"
      groups=interior,connection,corner,symmetry
      [[ $run == ene ]] && groups=$groups,wall301
      python3 -B "$MC" check "$mesh" "$run" "OUTPUT/diag-$mesh-$run-ns.tec" "OUTPUT/diag-$mesh-$run-ref.tec" "$groups" \
        || status=1
    done
  done
  echo "-- control: walls 301 with q = 0 on parallelogram cells, exact with every face metric --"
  python3 -B "$MC" case paraQ ene > /dev/null || exit 1
  step ns paraQ ene
  step ref paraQ ene
  python3 -B "$MC" check paraQ ene OUTPUT/diag-paraQ-ene-ns.tec OUTPUT/diag-paraQ-ene-ref.tec interior,wall301 || status=1
  echo "-- at every stage of a three-stage step, around the points shared by 3 and 5 blocks --"
  for mesh in star3 star5; do
    python3 -B "$MC" case "$mesh" rk3 > /dev/null || exit 1
    step ns "$mesh" rk3
    python3 -B "$MC" check "$mesh" rk3 "OUTPUT/diag-$mesh-rk3-ns.tec" none centre || status=1
  done
  echo "-- the face metrics inverse-mean and mean against an independent replica of their formulas --"
  for metric in inverse-mean mean; do
    python3 -B "$MC" case rand1 mom "$metric" > /dev/null || exit 1
    step ns rand1 mom "$metric"
    step ref rand1 mom "$metric"
    python3 -B "$MC" replica rand1 "OUTPUT/diag-rand1-mom-ns-$metric.tec" "OUTPUT/diag-rand1-mom-ref-$metric.tec" "$metric" \
      || status=1
  done
fi

if [[ $status -eq 0 ]]; then
  echo "$NAME: all checks passed"
  exit 0
fi
echo "$NAME: FAILED"
exit 1
