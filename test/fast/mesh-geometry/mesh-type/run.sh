#!/usr/bin/env bash
#===============================================================================
#  Mesh type: plane mesh, axisymmetric wedge, rotational periodicity (fast)
#
#  MOSE tells a plane 2-D mesh from an axisymmetric wedge, and a 3-D sector
#  from a translated periodic slab, by the nodes of the mesh (make_case.py):
#  - slab: a plane slab with y falling along j and z from -0.5 to 0.5 must be
#    plane; the wall distance read back from the Spalart-Allmaras source must
#    be the plane distance within 1e-9;
#  - flipwedge: a right-handed wedge with the axis on its last j face must be a
#    wedge; the distance must be the meridional one within 1e-9;
#  - wedge5: two blocks of a wedge of 5 degrees, their coordinates written with
#    5 significant digits, must be a wedge; the distance must be the
#    meridional one to the wall of the nodes as written within 1e-9;
#  - farslab: a plane slab far from the x axis (y from 1 to 1.001, z from 0 to
#    1e-4), which a small turn nearly fits, must be plane; the distance must be
#    the plane one within 1e-9;
#  - the set-up summary of these four names the type of the mesh;
#  - coarsewedge: a wedge of 1 degree written with 2 significant digits, too
#    few to tell it from a plane block, must stop the run at set-up with a
#    non-zero exit code and a message that asks for more digits;
#  - tiltslab: a plane slab whose k planes are not parallel (2 % of its
#    thickness), nearer a translation than a turn but not to 1 %, must stop
#    the run at set-up with a non-zero exit code and a message on the mesh
#    type;
#  - sectorshift100, sectorshift20: a 3-D sector of 90 degrees closed by type
#    200 whose last k node plane is turned by half a cell (100 cells in k) or
#    a tenth of a cell (20 cells) more than the first, but at two nodes, must
#    stop the run at set-up with a non-zero exit code and a message on the
#    mesh type: the plane does not match the turn to 1 % of a cell;
#  - mixed: a wedge block and a plane block in one mesh must stop the run at
#    set-up with a non-zero exit code and a message on the mesh type;
#  - slab3d: a 3-D slab periodic in k (type 200, no rotation) must keep a
#    uniform flow uniform within 1e-12 for five steps.
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

sa_input() {
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
}

euler_input() {
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = false
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF
}

fail=0
echo "-- wall distance on plane slabs and on wedges, read back from Spalart-Allmaras --"
for name in slab flipwedge wedge5 farslab; do
  python3 -B make_case.py case "$name" || exit 1
  sa_input
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/$name.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$name.log"; then
    echo "  [FAIL] $name: the run did not complete (see OUTPUT/$name.log)"
    fail=1
    continue
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-$name.tec"
  read -r err cells < <(python3 -B make_case.py error "$name" "OUTPUT/diagnostic-$name.tec")
  if python3 -c "import sys; sys.exit(0 if float('$err') <= 1e-9 else 1)"; then
    echo "  [ok  ] $name: largest relative error of the distance $err over $cells cells"
  else
    echo "  [FAIL] $name: largest relative error of the distance $err over $cells cells"
    fail=1
  fi
  case $name in
    slab|farslab) want='2-D, plane' ;;
    flipwedge)    want='2-D, axisymmetric wedge of -?(0\.99|1\.00)[0-9]* degrees' ;;
    wedge5)       want='2-D, axisymmetric wedge of -?(4\.99|5\.00)[0-9]* degrees' ;;
  esac
  if grep -a -q -E "Mesh type +$want" "OUTPUT/$name.log"; then
    echo "  [ok  ] $name: the set-up summary names the mesh type ($(grep -a -o -E 'Mesh type .*' "OUTPUT/$name.log" | tr -s ' '))"
  else
    echo "  [FAIL] $name: the set-up summary does not name the mesh type as '$want' (see OUTPUT/$name.log)"
    fail=1
  fi
done

echo "-- a wedge written with 2 significant digits, too few to tell it --"
python3 -B make_case.py case coarsewedge || exit 1
euler_input 1
"$LOCAL" > OUTPUT/coarsewedge.log 2>&1
rc=$?
if [[ $rc -ne 0 ]] && grep -a -q 'nearer a wedge of angle' OUTPUT/coarsewedge.log; then
  echo "  [ok  ] coarsewedge: stopped at set-up, exit code $rc, more digits asked"
else
  echo "  [FAIL] coarsewedge: exit code $rc, the run was not stopped on the digits of the wedge (see OUTPUT/coarsewedge.log)"
  fail=1
fi

echo "-- a plane slab with tilted k planes, and 3-D sectors whose periodic plane does not match by a fraction of a cell --"
for name in tiltslab sectorshift100 sectorshift20; do
  python3 -B make_case.py case "$name" || exit 1
  euler_input 1
  "$LOCAL" > "OUTPUT/$name.log" 2>&1
  rc=$?
  want='neither the first one turned about the x axis nor translated'
  [[ $name == tiltslab ]] && want='is nearer a plane block than a wedge'
  if [[ $rc -ne 0 ]] && grep -a -q "$want" "OUTPUT/$name.log"; then
    echo "  [ok  ] $name: stopped at set-up, exit code $rc"
  else
    echo "  [FAIL] $name: exit code $rc, the run was not stopped on the mesh type (see OUTPUT/$name.log)"
    fail=1
  fi
done

echo "-- a wedge block and a plane block in one mesh --"
python3 -B make_case.py case mixed || exit 1
euler_input 1
"$LOCAL" > OUTPUT/mixed.log 2>&1
rc=$?
if [[ $rc -ne 0 ]] && grep -a -q 'Mesh type' OUTPUT/mixed.log; then
  echo "  [ok  ] mixed: stopped at set-up, exit code $rc"
else
  echo "  [FAIL] mixed: exit code $rc, the run was not stopped on the mesh type (see OUTPUT/mixed.log)"
  fail=1
fi

echo "-- 3-D slab periodic in k (type 200, no rotation), uniform flow --"
python3 -B make_case.py case slab3d || exit 1
euler_input 5
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/slab3d.log 2>&1
if [[ ! -f OUTPUT/field.tec ]] || grep -a -qE 'NaN|Integration failed' OUTPUT/slab3d.log; then
  echo "  [FAIL] slab3d: the run did not complete (see OUTPUT/slab3d.log)"
  fail=1
else
  dev=$(python3 -B make_case.py uniform OUTPUT/field.tec)
  if python3 -c "import sys; sys.exit(0 if float('$dev') <= 1e-12 else 1)"; then
    echo "  [ok  ] slab3d: largest change of the uniform flow $dev"
  else
    echo "  [FAIL] slab3d: largest change of the uniform flow $dev"
    fail=1
  fi
fi

if [[ $fail -ne 0 ]]; then
  echo "MeshType: FAILED"
  exit 1
fi
echo "MeshType: all checks passed"
exit 0
