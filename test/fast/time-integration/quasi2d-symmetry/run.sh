#!/usr/bin/env bash
#===============================================================================
#  Quasi-2D case with symmetry planes on the k faces (fast tier)
#
#  A planar case meshed with one cell in z, thickness LZ much smaller than the
#  x and y cells, and symmetry (300) on both k faces: the two mirror planes
#  give w the mode dw/dt = -2 a w / LZ. At the time step of the x and y cells
#  (cfl 0.5) 2 a dt / LZ is about 2.2: outside the real stability interval of
#  Heun, [-2, 0], inside that of RK3, [-2.51, 0]. A flow at uniform pressure
#  and velocity carrying an entropy wave, with a z velocity of 1e-8 m/s (the
#  size of round-off), must keep p, u and v uniform and w must not grow,
#    1. with Heun (RK2) and with RK3, k faces 300;
#    2. with Heun, null (0) k faces (no flux through them): the control case;
#    3. with Heun, k faces 300 and a z velocity of 1 m/s at the start (a field
#       restarted from a run hit by the mode, or from a 3-D field).
#    4. a channel whose section varies along x in y (j faces 300, the top one
#       inclined), in a layer between parallel planes inclined to z (k faces
#       300): here w is not set to zero, as symmetry on the inclined k faces
#       asks w = c v; steady subsonic flow from a reservoir, p must follow the
#       quasi-1-D isentropic flow of the same section within 0.03 p0.
#  The same in a 1-D tube, one cell in j and one in k, both 1 mm thick against
#  x cells of 5 mm, with y and z velocities of 1e-8 m/s: symmetry on all the j
#  faces gives v the mode dv/dt = -2 a v / dy, on all the k faces w the mode
#  dw/dt = -2 a w / dz. p and u must stay uniform and v and w must not grow,
#    4. with Heun and with RK3, j and k faces 300;
#    5. with Heun, j faces 300 and null k faces;
#    6. with Heun, extrapolation (400) j faces and null k faces: the control.
#    7. a channel whose section varies along x, one cell in y with the top j
#       face inclined (j faces 300, null k faces), steady subsonic flow from a
#       reservoir: here v is not set to zero; p must follow the quasi-1-D
#       isentropic flow within 0.03 p0.
#  Box, inputs and checks: make_case.py; gas of test/common/Air.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
COMMON=../../../common/Air
N=50          # time steps

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"  INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat" INPUT/thermo.dat
fail=0

# solve <tag> <k-face type> <time scheme> <initial w>: N time-accurate steps at cfl 0.5;
# the final field is kept as OUTPUT/<tag>.tec and the solver output as OUTPUT/<tag>.log
solve() {
  local tag="$1"
  python3 -B make_case.py ic "$4" || return 1
  python3 -B make_case.py bc "$2" || return 1
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $N
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = $3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  echo "-- k faces $2, $3, initial w $4 m/s, $N steps --"
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] no solution written (see OUTPUT/$tag.log)"
    fail=1
    return 1
  fi
  cp OUTPUT/field.tec "OUTPUT/$tag.tec"
  python3 -B make_case.py check "OUTPUT/$tag.tec" "$4"
  if [[ $? -eq 0 ]]; then
    echo "  [ok  ] p, u, v stay uniform and w does not grow"
  else
    echo "  [FAIL] p, u, v do not stay uniform or w grows"
    fail=1
  fi
}

solve sym-RK2    300 RK2 1e-8
solve sym-RK3    300 RK3 1e-8
solve null-RK2   0   RK2 1e-8
solve sym-RK2-w1 300 RK2 1

# solve_channel <tag> <kind> <warning>: steady flow in the channel of variable section (make_case.py channel), 6000
# iterations of RK3 at cfl 0.5; the faces along the walls are symmetry planes, one of them inclined: MOSE must say that
# it integrates the velocity across the cells (<warning>), and p must follow the quasi-1-D isentropic flow
solve_channel() {
  local tag="$1"
  python3 -B make_case.py channel "$2" || return 1
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 6000
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
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  echo "-- channel of variable section ($2), symmetry faces along its walls, one of them inclined, steady --"
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] no solution written (see OUTPUT/$tag.log)"
    fail=1
    return 1
  fi
  cp OUTPUT/field.tec "OUTPUT/$tag.tec"
  if grep -q "$3" "OUTPUT/$tag.log"; then
    echo "  [ok  ] the velocity across the cells is integrated, with the warning"
  else
    echo "  [FAIL] no warning that the velocity across the cells is integrated (see OUTPUT/$tag.log)"
    fail=1
  fi
  python3 -B make_case.py checkchannel "OUTPUT/$tag.tec" "$2"
  if [[ $? -eq 0 ]]; then
    echo "  [ok  ] p follows the quasi-1-D flow of the same section"
  else
    echo "  [FAIL] p does not follow the quasi-1-D flow of the same section"
    fail=1
  fi
}

solve_channel layer layer "k faces that are not planes normal to z"

# solve1d <tag> <j-face type> <k-face type> <time scheme>: the same N steps on the 1-D tube
solve1d() {
  local tag="$1"
  python3 -B make_case.py ic1d || return 1
  python3 -B make_case.py bc1d "$2" "$3" || return 1
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $N
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = $4
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOF
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  echo "-- 1-D tube, j faces $2, k faces $3, $4, initial v and w 1e-8 m/s, $N steps --"
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] no solution written (see OUTPUT/$tag.log)"
    fail=1
    return 1
  fi
  cp OUTPUT/field.tec "OUTPUT/$tag.tec"
  python3 -B make_case.py check1d "OUTPUT/$tag.tec"
  if [[ $? -eq 0 ]]; then
    echo "  [ok  ] p, u stay uniform and v, w do not grow"
  else
    echo "  [FAIL] p, u do not stay uniform or v, w grow"
    fail=1
  fi
}

solve1d sym1d-RK2     300 300 RK2
solve1d sym1d-RK3     300 300 RK3
solve1d jsym1d-RK2    300 0   RK2
solve1d null1d-RK2    400 0   RK2
solve_channel cd1d cd1d "j faces that are not planes normal to y"

if [[ $fail -ne 0 ]]; then
  echo "Quasi2DSymmetry: FAILED"
  exit 1
fi
echo "Quasi2DSymmetry: all checks passed"
exit 0
