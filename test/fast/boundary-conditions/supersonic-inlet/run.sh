#!/usr/bin/env bash
#===============================================================================
#  Supersonic inlet (405) with the flow direction normal (fast tier)
#
#  With alpha and beta both normal the 405 imposes its Mach number along the
#  unit normal of the face into the domain; with one of them normal, that
#  angle is the one of that normal, alpha its azimuth and beta its elevation
#  (docs/user/boundary-conditions.md). A 1-D tube of air, one forward Euler
#  step with first-order reconstruction and HLLC: the 405 at Mach 2 at one
#  end, the gas of the tube moving in through it at 0.1 a, an extrapolation
#  face (400) at the other end (supersonic_normal.py). The tube along x, its
#  mirror image (i growing along -x), the tube turned so that the normal of
#  its end faces is (1, 1, 1)/sqrt(3), and the tube along z:
#    normal, normal on the faces 1 and 2 of the tube and of its mirror image
#    and on face 1 of the turned tube; normal, 0 on face 2 of the tube; 0.3,
#    normal on the faces 1 and 2 of the tube along z: the mass flux into the
#    domain is that of 0, 0 on face 1 of the tube within 1e-9 of rho a, and
#    the velocity of the inlet cell along the normal into the domain is that
#    of 0, 0 within 1e-9 of a, with a component across it below 1e-9 of a;
#    0, normal and 0.2, normal on the faces 1 and 2 of the tube, and normal,
#    0.3 on face 1 of the tube along z, whose normal has no azimuth: the mass
#    flux and the velocity of the inlet cell of the angles the token stands
#    for on that face (0, 0; 0.2, 0; 0, 0.3), within 1e-9 of rho a and of a.
#  Gas: test/common/Air.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
COMMON=../../../common/Air

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"  INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat" INPUT/thermo.dat

cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.1
time-accurate = true
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC
INI

fail=0
for name in ref x1n x2n m1n m2n r1n x2a z1b z2b x1c x2c x2e x1d x2d x1f x2f z1e z1g; do
  python3 -B supersonic_normal.py case "$name" || exit 1
  rm -f OUTPUT/field.tec "OUTPUT/$name.tec"
  timeout 60 "$LOCAL" > "OUTPUT/$name.log" 2>&1
  if [[ $? -ne 0 || ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] $name: the run did not complete (see OUTPUT/$name.log)"
    fail=1
    continue
  fi
  cp OUTPUT/field.tec "OUTPUT/$name.tec"
done

python3 -B supersonic_normal.py check || fail=1

if [[ $fail -ne 0 ]]; then
  echo "SupersonicInlet: FAILED"
  exit 1
fi
echo "SupersonicInlet: all checks passed"
exit 0
