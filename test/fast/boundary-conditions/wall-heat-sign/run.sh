#!/usr/bin/env bash
#===============================================================================
#  Sign of the wall heat flux and of the wall stress (fast tier)
#
#  q_wall of a wall 301 is positive into the fluid (docs/user/boundary-
#  conditions.md), whatever face the wall is on, and so are the heat flux qw
#  of OUTPUT/wall.tec and its stress tauX, tauY, tauZ, the traction of the gas
#  on the wall (tau . n, n the normal into the gas): the same on every face.
#  A channel with the same wall on the two faces of one index direction, one
#  face of each parity, which mirror each other (make_case.py): a planar 2-D
#  channel with the walls on the j faces (3 and 4), and a 3-D channel with the
#  walls on the i faces (1 and 2) and on the k faces (5 and 6). Each runs one
#  step:
#    1. type 301 with q_wall = 1e4 W/m2 on both walls, air at rest: the cells
#       next to the two walls must warm up alike, and qw must be 1e4 W/m2 on
#       both walls;
#    2. type 302 at 400 K on both walls, air at rest: qw must be the same on
#       the two walls and positive (the walls heat the gas);
#    3. type 302 at the temperature of the air, which slides along x at
#       10 m/s: tauX must be the same on the two walls and positive (the gas
#       drags the walls along), tauY opposite and tauZ the same (the traction
#       of mirror walls).
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
for dir in j i k; do
  for kind in flux temp shear; do
    tag=$kind-$dir
    echo "-- the same wall (type $([[ $kind == flux ]] && echo 301 || echo 302)) on both $dir faces, one step ($kind) --"
    python3 -B make_case.py case "$kind" "$dir" || exit 1
    rm -f OUTPUT/field.tec OUTPUT/wall.tec
    "$LOCAL" > "OUTPUT/$tag.log" 2>&1
    if [[ $? -ne 0 || ! -f OUTPUT/field.tec || ! -f OUTPUT/wall.tec ]]; then
      echo "  [FAIL] the run did not complete (see OUTPUT/$tag.log)"
      fail=1
      continue
    fi
    cp OUTPUT/field.tec "OUTPUT/field-$tag.tec"
    cp OUTPUT/wall.tec "OUTPUT/wall-$tag.tec"
    if [[ $kind == flux ]]; then
      if python3 -B make_case.py check-flux "$dir" "OUTPUT/field-$tag.tec"; then
        echo "  [ok  ] both walls heat the gas alike"
      else
        echo "  [FAIL] the walls do not heat the gas alike"
        fail=1
      fi
      if python3 -B make_case.py check-qw-flux "$dir" "OUTPUT/wall-$tag.tec"; then
        echo "  [ok  ] qw is q_wall on both walls: the heat flux into the gas"
      else
        echo "  [FAIL] qw is not q_wall on both walls"
        fail=1
      fi
    elif [[ $kind == temp ]]; then
      if python3 -B make_case.py check-qw-temp "$dir" "OUTPUT/wall-$tag.tec"; then
        echo "  [ok  ] qw is the same on the two walls and positive: the walls heat the gas"
      else
        echo "  [FAIL] qw is not the same positive heat flux into the gas on the two walls"
        fail=1
      fi
    else
      if python3 -B make_case.py check-tau "$dir" "OUTPUT/wall-$tag.tec"; then
        echo "  [ok  ] tauX is the same on the two walls and positive (the gas drags them along), tauY opposite, tauZ the same"
      else
        echo "  [FAIL] the traction of the gas is not the same positive tauX on the two walls, or not mirrored in tauY, tauZ"
        fail=1
      fi
    fi
  done
done

if [[ $fail -ne 0 ]]; then
  echo "WallHeatSign: FAILED"
  exit 1
fi
echo "WallHeatSign: all checks passed"
exit 0
