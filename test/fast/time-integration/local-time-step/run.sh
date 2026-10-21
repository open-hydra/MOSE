#!/usr/bin/env bash
#===============================================================================
#  Local time step (fast tier)
#
#  The explicit step must be stable in the cases where one direction or one
#  process does not dominate:
#    box   sound waves crossing square 2-D cells in both directions, cfl 0.8;
#    skew  an entropy wave at uniform pressure and velocity on a stretched and
#          skewed periodic 2-D mesh, cfl 0.8;
#    visc  heat conduction in a 1-D Navier-Stokes case at cell Peclet number
#          0.034, with the default cfl 0.5 and vnn 0.3;
#    box3d sound waves crossing cubic cells in the three directions, cfl 0.8;
#    species  species diffusion (schmidt 0.2, diffusivity 5 nu) in a 1-D
#          four-species case at cell Peclet number below 0.1, with the default
#          cfl 0.5 and vnn 0.3;
#  each with Heun (RK2) and RK3, and
#    sod   the 1-D shock tube of test/1D/Sod79 at cfl 0.99 with Heun;
#    rhombo  heat conduction on a periodic 3-D block of 8 x 8 x 8
#          rhombohedra, every edge angle 115 degrees, with Heun at vnn 0.5 for
#          600 steps: the temperature wave of 135 degrees per cell along the
#          three indices must decay (with the step of the rate
#          4 D sum_d |grad xi_d|^2 alone it grows by 1.076 per step).
#  Cases and checks: make_case.py; gases of test/common/Air and Air-viscous
#  and the four-species gas of test/1D/diff1.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
AIR=../../../common/Air
VISC=../../../common/Air-viscous
SOD=../../../1D/Sod79
MIX=../../../1D/diff1/INPUT

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi

fail=0

# gas <directory>: link the phase, thermo (and transport, composition) files
gas() {
  rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat INPUT/composition.txt
  ln -sf "../$1/phase.txt"  INPUT/phase.txt
  ln -sf "../$1/thermo.dat" INPUT/thermo.dat
  [[ -f "$1/transport.dat" ]] && ln -sf "../$1/transport.dat" INPUT/transport.dat
  [[ -f "$1/composition.txt" ]] && ln -sf "../$1/composition.txt" INPUT/composition.txt
  return 0
}

# solve <case> <equations> <time scheme> <cfl> <vnn> <steps> [schmidt]: time-accurate run, then the check of
# make_case.py; the final field is kept as OUTPUT/<case>-<scheme>.tec and the solver output as OUTPUT/<case>-<scheme>.log
solve() {
  local tag="$1-$3"
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $6
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = $2
schmidt = ${7:-0.7}

[MOSE-Numerics]
cfl = $4
vnn = $5
time-accurate = true
time-scheme = $3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  echo "-- $1, $3, cfl $4, vnn $5, $6 steps --"
  if [[ ! -f OUTPUT/field.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$tag.log"; then
    echo "  [FAIL] the run blew up (see OUTPUT/$tag.log)"
    fail=1
    return 1
  fi
  cp OUTPUT/field.tec "OUTPUT/$tag.tec"
  local done_steps
  done_steps=$(grep -a -oE 'Iter = +[0-9]+' "OUTPUT/$tag.log" | tail -1 | grep -a -oE '[0-9]+$')
  if [[ "$done_steps" != "$6" ]]; then
    echo "  [FAIL] the run stopped after ${done_steps:-0} of $6 steps (see OUTPUT/$tag.log)"
    fail=1
    return 1
  fi
  if python3 -B make_case.py check "$1" "OUTPUT/$tag.tec"; then
    echo "  [ok  ] stable"
  else
    echo "  [FAIL] unstable"
    fail=1
  fi
}

gas "$AIR"
rm -f INPUT/ic.tec INPUT/bc.txt   # never write through the links of the sod case
python3 -B make_case.py box || exit 1
solve box euler RK2 0.8 0.3 100
solve box euler RK3 0.8 0.3 100

rm -f INPUT/ic.tec INPUT/bc.txt   # never write through the links of the sod case
python3 -B make_case.py skew || exit 1
solve skew euler RK2 0.8 0.3 100
solve skew euler RK3 0.8 0.3 100

gas "$VISC"
rm -f INPUT/ic.tec INPUT/bc.txt   # never write through the links of the sod case
python3 -B make_case.py visc || exit 1
solve visc navier-stokes RK2 0.5 0.3 300
solve visc navier-stokes RK3 0.5 0.3 300

rm -f INPUT/ic.tec INPUT/bc.txt   # never write through the links of the sod case
python3 -B make_case.py rhombo || exit 1
solve rhombo navier-stokes RK2 0.5 0.5 600

gas "$AIR"
rm -f INPUT/ic.tec INPUT/bc.txt   # never write through the links of the sod case
python3 -B make_case.py box3d || exit 1
solve box3d euler RK2 0.8 0.3 100
solve box3d euler RK3 0.8 0.3 100

gas "$MIX"
rm -f INPUT/ic.tec INPUT/bc.txt   # never write through the links of the sod case
python3 -B make_case.py species || exit 1
solve species navier-stokes RK2 0.5 0.3 300 0.2
solve species navier-stokes RK3 0.5 0.3 300 0.2

gas "$AIR"
rm -f INPUT/ic.tec INPUT/bc.txt
ln -sf "../$SOD/INPUT/ic_x1.tec" INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt" INPUT/bc.txt
solve sod euler RK2 0.99 0.3 50

if [[ $fail -ne 0 ]]; then
  echo "LocalTimeStep: FAILED"
  exit 1
fi
echo "LocalTimeStep: all checks passed"
exit 0
