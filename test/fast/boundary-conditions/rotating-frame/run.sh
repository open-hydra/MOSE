#!/usr/bin/env bash
#===============================================================================
#  Rotating frame (fast tier)
#
#  The walls that do not rotate with the frame (casing, hub) are given in
#  [MOSE-Rotating-Frame] as stationary-face-1, stationary-face-2, ... = block
#  face, numbered from 1. The Sod tube of test/1D/Sod79 (100 cells) runs 5
#  time-accurate steps in a rigid-body frame (omega = 1 rad/s about x) with
#  two such entries: the run must end within 60 s with exit code 0, and its
#  set-up summary must list the two patches, block 1 faces 3 and 4.
#  Then a stationary wall (make_case.py): a short tube along x, the frame
#  turning at 100 rad/s about x, the gas at rest in the absolute frame and
#  face 1 an adiabatic wall that does not rotate (stationary-face-1 = 1 1).
#  In the frame the wall and the cell next to it move alike, with -Omega x r:
#  after one step, the wall shear in the wall output must stay below 1e-3 of
#  mu Omega y_c / (dx/2), the shear of the same wall taken as rotating.
#  Reuses the Sod inputs (test/1D/Sod79) and the gas of test/common/Air.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../../1D/Sod79
COMMON=../../../common/Air
VISC=../../../common/Air-viscous
SUTH=../../../common/Air-sutherland

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin INPUT OUTPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/transport.dat          # left by the viscous case of an earlier run
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"  INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"  INPUT/bc.txt

cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 5
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler
rotational-frame = rigid-body

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-Rotating-Frame]
omega = 1.0
axis = 1.0 0.0 0.0
origin = 0.0 0.0 0.0
stationary-face-1 = 1 3
stationary-face-2 = 1 4
INI

echo "-- rigid-body frame with two stationary faces, 5 steps --"
timeout 60 "$LOCAL" > OUTPUT/setup.log 2>&1
rc=$?
fail=0
if [[ $rc -eq 124 ]]; then
  echo "  [FAIL] the run did not end within 60 s (set-up loop, see OUTPUT/setup.log)"
  fail=1
elif [[ $rc -ne 0 ]]; then
  echo "  [FAIL] exit code $rc (see OUTPUT/setup.log)"
  fail=1
else
  echo "  [ok  ] the run ended with exit code 0"
  if command grep -q 'stationary patches = 2:' OUTPUT/setup.log && \
     command grep -q 'block 1  face-dir 3' OUTPUT/setup.log && \
     command grep -q 'block 1  face-dir 4' OUTPUT/setup.log; then
    echo "  [ok  ] the set-up summary lists the two stationary patches, block 1 faces 3 and 4"
  else
    echo "  [FAIL] the set-up summary does not list the two stationary patches (see OUTPUT/setup.log)"
    fail=1
  fi
fi

# stationary wall: the wall output must use the wall velocity of the boundary fluxes
rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
ln -sf "../$VISC/phase.txt"     INPUT/phase.txt
ln -sf "../$VISC/thermo.dat"    INPUT/thermo.dat
ln -sf "../$SUTH/transport.dat" INPUT/transport.dat
python3 -B make_case.py case || exit 1
cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
schmidt = 0.7
rotational-frame = rigid-body

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii

[MOSE-Rotating-Frame]
omega = 100.0
axis = 1.0 0.0 0.0
origin = 0.0 0.0 0.0
stationary-face-1 = 1 1
INI
echo "-- stationary wall in the rotating frame, gas at rest, one step --"
rm -f OUTPUT/wall.tec
timeout 60 "$LOCAL" > OUTPUT/wall.log 2>&1
rc=$?
if [[ $rc -ne 0 || ! -f OUTPUT/wall.tec ]]; then
  echo "  [FAIL] exit code $rc or no wall output (see OUTPUT/wall.log)"
  fail=1
else
  r=$(python3 -B make_case.py check OUTPUT/wall.tec)
  if python3 -c "import sys; sys.exit(0 if float('$r') < 1e-3 else 1)"; then
    echo "  [ok  ] wall shear on the stationary wall $r of mu Omega y_c / (dx/2)"
  else
    echo "  [FAIL] wall shear on the stationary wall $r of mu Omega y_c / (dx/2): the wall taken as rotating"
    fail=1
  fi
fi

if [[ $fail -ne 0 ]]; then
  echo "RotatingFrame: FAILED"
  exit 1
fi
echo "RotatingFrame: all checks passed"
exit 0
