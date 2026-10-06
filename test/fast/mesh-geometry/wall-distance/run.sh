#!/usr/bin/env bash
#===============================================================================
#  Wall distance: exact minimum distance to the wall (fast tier)
#
#  The distance to the nearest wall of the Spalart-Allmaras model must be the
#  minimum distance from the cell centre to the wall surface, which near a
#  convex corner is the distance to the corner (NASA Turbulence Modeling
#  Resource, Spalart-Allmaras page). With the fluid at rest and nu~ uniform,
#  one explicit Euler step gives, in OUTPUT/diagnostic.tec, the right-hand
#  side of rho nu~, which away from the walls is rho (nu~/d)^2 times a
#  function of nu~/nu alone, and so the distance d of every cell
#  (make_case.py). Three cases: a convex corner of 90 degrees between two
#  walls, a curved wall under a leaning mesh, and the lip of a nozzle on an
#  axisymmetric wedge (distance in the meridional plane); and five cases
#  with periodic boundaries, across which the wall goes on: a 2-D channel and a
#  pipe on a wedge periodic in x, a 3-D sector closed by the rotational
#  periodicity, the channel with the two sides of the period matched to 1e-8
#  only, and a 3-D slab periodic in k by type 200, each with the wall on part
#  of the period only. On the cells that do not touch the wall, d must equal
#  the exact distance to the wall of the mesh and its periodic images within
#  1e-9 (1e-6 for the channel matched to 1e-8).
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../common/Air-viscous

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$VISC/phase.txt"     INPUT/phase.txt
ln -sf "../$VISC/thermo.dat"    INPUT/thermo.dat
ln -sf "../$VISC/transport.dat" INPUT/transport.dat

cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = SA
schmidt = 0.7

[MOSE-Turbulence]
point-implicit = false

[MOSE-Numerics]
cfl = 0.5
vnn = 0.3
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF

echo "-- distance to the wall read back from the Spalart-Allmaras source --"
errors=()
for name in corner arc lip periodic pipe sector noisy span; do
  python3 -B make_case.py case "$name" || exit 1
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/$name.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$name.log"; then
    echo "  [FAIL] $name: the run did not complete (see OUTPUT/$name.log)"
    echo "WallDistance: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$name.tec"
  read -r err cells corner_cells < <(python3 -B make_case.py error "$name" "OUTPUT/diagnostic-$name.tec")
  echo "  $name: largest relative error of d $err over $cells cells ($corner_cells nearest to a wall node)"
  errors+=("$name=$err")
done

if python3 -B make_case.py check "${errors[@]}"; then
  echo "  [ok  ] exact distance to the wall"
  echo "WallDistance: all checks passed"
  exit 0
fi
echo "  [FAIL] the distance is not the minimum distance to the wall"
echo "WallDistance: FAILED"
exit 1
