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
#  axisymmetric wedge (distance in the meridional plane); and eleven cases
#  with periodic boundaries, across which the wall goes on: a 2-D channel and a
#  pipe on a wedge periodic in x, a 3-D sector closed by the rotational
#  periodicity, the channel with the two sides of the period matched to 1e-8
#  only, and a 3-D slab periodic in k by type 200, each with the wall on part
#  of the period only; a 3-D channel periodic in x and z, its periods a
#  hundredth of its height, with the wall on the whole period: its images are
#  one plane, many times over; the 2-D channel ten periods away from the
#  origin; a staggered row of plane blades, each 1.7 periods long along
#  the periodic direction; and a 3-D block periodic along three nearly
#  dependent translations, whose sum is fifty times shorter than each of
#  them and than any pairwise combination; and the 2-D channel turned by 30
#  degrees about z, its node coordinates rounded to 5 and to 4 significant
#  digits, whose periodic faces give translations that differ by up to
#  1.1e-4 and 8.8e-4 of the cells across them: translations that differ by
#  at most 1e-2 of the thinnest cell across the periodic boundaries are the
#  same. On the cells that do not touch the wall, d must
#  equal the exact distance to the wall of the mesh and its periodic images
#  within 1e-9 (1e-6 for the channel matched to 1e-8; for the turned
#  channel, it must lie within 1e-9 between the distances with the
#  translation of each periodic face, of which the run takes one). Rounded
#  to 2 digits, the turned channel gives translations that differ by up to
#  0.15 of its cells: the run must stop at set-up with the message on the
#  periodic translations, which asks for sides that match. On a 3-D wall of
#  twisted faces, where the four triangles of a face are not the face, the
#  run must report at start-up the largest bound of their gap over the
#  distance of the cell next to the face; on plane faces it must report
#  nothing.
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
for name in corner arc lip periodic pipe sector noisy span thin far staggered oblique turned5 turned4; do
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

fail=0
if python3 -B make_case.py check "${errors[@]}"; then
  echo "  [ok  ] exact distance to the wall"
else
  echo "  [FAIL] the distance is not the minimum distance to the wall"
  fail=1
fi

echo "-- periodic sides that do not match to 1e-2 of the cells: stop at set-up --"
python3 -B make_case.py case turned2 || exit 1
rm -f OUTPUT/diagnostic.tec
"$LOCAL" > OUTPUT/turned2.log 2>&1
if [[ ! -f OUTPUT/diagnostic.tec ]] && python3 -B make_case.py turnstop OUTPUT/turned2.log; then
  echo "  [ok  ] stopped at set-up with the message on the periodic translations"
else
  echo "  [FAIL] turned2: no stop at set-up with the message on the periodic translations (see OUTPUT/turned2.log)"
  fail=1
fi

echo "-- twisted wall faces in 3-D: the start-up report of the error bound --"
python3 -B make_case.py case twist || exit 1
"$LOCAL" > OUTPUT/twist.log 2>&1
if python3 -B make_case.py twist OUTPUT/twist.log; then
  echo "  [ok  ] reported"
else
  echo "  [FAIL] the twist of the wall faces is not reported, or not right (see OUTPUT/twist.log)"
  fail=1
fi
for name in corner arc lip periodic pipe sector noisy span thin far staggered oblique turned5 turned4; do
  if grep -a -q 'twisted wall faces' "OUTPUT/$name.log"; then
    echo "  [FAIL] $name: a twist reported on plane wall faces"
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "WallDistance: FAILED"
  exit 1
fi
echo "WallDistance: all checks passed"
exit 0
