#!/usr/bin/env bash
#===============================================================================
#  Sign of the wall heat flux (fast tier)
#
#  q_wall of a wall 301 is positive into the fluid (docs/user/boundary-
#  conditions.md), whatever face the wall is on; the heat flux qw written in
#  OUTPUT/wall.tec is k grad T . n along the face normal n toward growing
#  index: into the fluid on faces 2, 4, 6, into the wall on faces 1, 3, 5.
#  A planar 2-D channel of air at rest with the same wall on both j faces,
#  which mirror each other (make_case.py), runs one step:
#    1. type 301 with q_wall = 1e4 W/m2 on both walls: the cells next to the
#       two walls must warm up alike, and qw must be -1e4 W/m2 on face 3 and
#       1e4 W/m2 on face 4;
#    2. type 302 at 400 K on both walls: qw must be opposite on the two walls,
#       positive on face 4 (the walls heat the gas).
#  Gas: test/common/Air-viscous with the transport of test/common/Air-sutherland.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../common/Air-viscous
SUTH=../../../common/Air-sutherland

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$VISC/phase.txt"     INPUT/phase.txt
ln -sf "../$VISC/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SUTH/transport.dat" INPUT/transport.dat

cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes

[MOSE-Numerics]
cfl = 0.1
time-accurate = true
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
INI

fail=0
for kind in flux temp; do
  echo "-- the same wall (type $([[ $kind == flux ]] && echo 301 || echo 302)) on both j faces, one step --"
  python3 -B make_case.py case "$kind" || exit 1
  rm -f OUTPUT/field.tec OUTPUT/wall.tec
  "$LOCAL" > "OUTPUT/$kind.log" 2>&1
  if [[ $? -ne 0 || ! -f OUTPUT/field.tec || ! -f OUTPUT/wall.tec ]]; then
    echo "  [FAIL] the run did not complete (see OUTPUT/$kind.log)"
    fail=1
    continue
  fi
  cp OUTPUT/field.tec "OUTPUT/field-$kind.tec"
  cp OUTPUT/wall.tec "OUTPUT/wall-$kind.tec"
  if [[ $kind == flux ]]; then
    if python3 -B make_case.py check-flux "OUTPUT/field-$kind.tec"; then
      echo "  [ok  ] both walls heat the gas alike"
    else
      echo "  [FAIL] the walls do not heat the gas alike"
      fail=1
    fi
    if python3 -B make_case.py check-qw-flux "OUTPUT/wall-$kind.tec"; then
      echo "  [ok  ] qw is k grad T . n along the face normal: -q_wall on face 3, q_wall on face 4"
    else
      echo "  [FAIL] qw is not -q_wall on face 3 and q_wall on face 4"
      fail=1
    fi
  else
    if python3 -B make_case.py check-qw-temp "OUTPUT/wall-$kind.tec"; then
      echo "  [ok  ] qw is opposite on the two walls, positive on face 4: the walls heat the gas"
    else
      echo "  [FAIL] qw is not opposite on the two walls with the sign of k grad T . n"
      fail=1
    fi
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "WallHeatSign: FAILED"
  exit 1
fi
echo "WallHeatSign: all checks passed"
exit 0
