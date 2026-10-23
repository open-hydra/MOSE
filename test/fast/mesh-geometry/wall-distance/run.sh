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
#  periodic translations, which asks for sides that match. Two 3-D walls of
#  twisted faces, each the bilinear surface through its four nodes, whose four
#  triangles between the centre and the sides are not the face: a block whose
#  wall is y = 2 (x - 0.5)(z - 0.5) (twist), where in some cells the nearest
#  face by the triangles is not the nearest face, and a block periodic in x
#  whose wall is y = 0.08 sin(2 pi x) sin(2 pi z) (twistper): d must equal the
#  exact distance to the bilinear faces and their periodic images within 1e-9
#  (references twist_exact.txt and twistper_exact.txt, written by make_case.py
#  exact NAME). No run may print a warning on twisted wall faces.
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
fan_differs=0
for name in corner arc lip periodic pipe sector noisy span thin far staggered oblique turned5 turned4 twist twistper; do
  python3 -B make_case.py case "$name" || exit 1
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/$name.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$name.log"; then
    echo "  [FAIL] $name: the run did not complete (see OUTPUT/$name.log)"
    echo "WallDistance: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$name.tec"
  read -r err cells third < <(python3 -B make_case.py error "$name" "OUTPUT/diagnostic-$name.tec")
  if [[ $name == twist* ]]; then
    echo "  $name: largest relative error of d $err over $cells cells ($third whose nearest face by the four triangles is not the nearest face)"
    [[ $name == twist ]] && fan_differs=$third
  else
    echo "  $name: largest relative error of d $err over $cells cells ($third nearest to a wall node)"
  fi
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

echo "-- twisted wall faces: the exact distance, no warning --"
if [[ $fan_differs -ge 1 ]]; then
  echo "  [ok  ] twist: $fan_differs cells whose nearest face by the four triangles is not the nearest face"
else
  echo "  [FAIL] twist: no cell whose nearest face by the four triangles differs: the case does not test the search"
  fail=1
fi
for name in corner arc lip periodic pipe sector noisy span thin far staggered oblique turned5 turned4 twist twistper; do
  if grep -a -q 'twisted wall faces' "OUTPUT/$name.log"; then
    echo "  [FAIL] $name: a warning on twisted wall faces"
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "WallDistance: FAILED"
  exit 1
fi
echo "WallDistance: all checks passed"
exit 0
