#!/usr/bin/env bash
#===============================================================================
#  Extrapolation faces (400) and the lines of the gas entering through them
#  (fast tier)
#
#  A 400 face takes the flux of the state of its boundary cell, F(q1), for
#  either direction of the gas: where the gas of the boundary cell moves into
#  the domain (w1 < 0, w1 its velocity along the normal out of the domain) it
#  enters with the state of the cell. The run prints a line the first time the
#  gas enters through the faces of a 400 boundary (the 400 faces of one face of
#  a block, on each grid level), with the iteration, the faces and the largest
#  -w1/a1 (a1 the speed of sound of the cell), one more each time that largest
#  value at a stage is more than 10 times that of the last line of the
#  boundary, and at the end of a run in which the gas entered the flux
#  evaluations with gas entering and their largest -w1/a1 (entering_gas.py).
#  The line of the ghost cells of walls and 400 faces comes only when the cut
#  of their extrapolation acted, which a uniform state never makes it do.
#  1. One forward Euler step with first-order reconstruction, a 1-D tube of air
#     (j faces 300, k faces 0, i faces 400) at a uniform velocity (u0, 0, w0):
#     the boundary cells do not change (within 1e-12 of rho, rho a, rho a^2)
#     and the flux of each 400 face from the change of its cell is F(q0) . S
#     within 1e-9 of rho a, rho a^2, rho a^3:
#       rest       u0 = 0;
#       sub, subm  u0 = +-0.5 a;
#       sup, supm  u0 = +-1.327 a, w0 = 2 a;
#       vac, vacm  u0 = +-6 a;
#     one line for the face where the gas enters (iteration 1, faces 1,
#     |u0|/a), none at rest; the end line with 1 evaluation, none at rest; no
#     line of the ghost cells.
#  2. The tube with the gas entering through face 1 at 0.5 a, 20 RK3 steps at
#     first order: the state stays uniform, one line, and the end line counts
#     60 evaluations, one a stage (count).
#  3. The tube with the gas at 1e-3 a and the pressure 10 % lower in its right
#     half, 40 forward Euler steps at first order: the expansion speeds up the
#     gas entering through face 1, which has a second line, more than 10 times
#     the first (tenfold).
#  4. A run that stops: the gas at 6 a moving away from face 1 of type 0, the
#     j faces 400 with the gas entering through face 3 at 1e-3 a, one forward
#     Euler step at cfl 1.13: the boundary cell gets a negative pressure at the
#     first stage and the run stops; the line of face 3 comes before the line
#     of the stop, and no end line is printed (stop).
#  5. A tube of 130 cells with all j faces 400, the gas entering through
#     face 1 at 0.5 a and through face 3 at 1e-3 a, one step with 1 and with 4
#     threads (130 records on face 3, more than two pieces of the dynamic
#     schedule of 64 of the loops of the boundary fluxes): the same lines and
#     the same end line, faces: 130 on face 3 (threads). A MOSE built without
#     OpenMP writes "Serial execution" where it writes the number of threads:
#     the part then stops after its run with 1 thread and is skipped, with a
#     line that says so.
#  6. The tube of part 1 with slip walls (300) at its ends, no 400 face, the
#     gas at 0.5 a, one step: no line of the gas entering, no end line and no
#     line of the ghost cells (no400).
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

# write_ini <steps> <scheme>: <steps> steps at first order, forward Euler at cfl 0.1 (euler), at cfl 0.5 (euler05) or at cfl
# 1.13 (stop), or RK3 at cfl 0.5 (rk3)
write_ini() {
  local numerics
  case $2 in
    euler)   numerics=('cfl = 0.1' 'time-scheme = euler') ;;
    euler05) numerics=('cfl = 0.5' 'time-scheme = euler') ;;
    stop)    numerics=('cfl = 1.13' 'time-scheme = euler') ;;
    rk3)     numerics=('cfl = 0.5' 'time-scheme = RK3') ;;
  esac
  printf '%s\n' '[MOSE-Parameters]' "iter-threshold = $1" 'newrun = true' 'res-threshold = 1e-300' '' \
    '[MOSE-Physics]' 'equations = euler' '' \
    '[MOSE-Numerics]' "${numerics[@]}" 'space-reconstruction = first-order' 'time-accurate = true' \
    'integration-variables = cons' 'riemann-solver = HLLC' '' \
    '[MOSE-IO]' 'sol-format = tecplot ascii' > input.ini
}

fail=0

# solve <name> <log>: write case <name> and solve it, the output in <log>; 0 when the run completes with a field
solve() {
  python3 -B entering_gas.py case "$1" || exit 1
  rm -f OUTPUT/field.tec
  "$LOCAL" > "$2" 2>&1
  if [[ $? -ne 0 || ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] $1: the run did not complete (see $2)"
    fail=1
    return 1
  fi
}

# run_case <name>: solve case <name> and check it
run_case() {
  local name=$1
  solve "$name" "OUTPUT/$name.log" || return
  cp OUTPUT/field.tec "OUTPUT/$name.tec"
  if python3 -B entering_gas.py check "$name" "OUTPUT/$name.tec" "OUTPUT/$name.log"; then
    echo "  [ok  ] $name"
  else
    echo "  [FAIL] $name"
    fail=1
  fi
}

echo "-- 1. one step: the flux of a 400 face and its lines --"
write_ini 1 euler
for name in rest sub subm sup supm vac vacm; do
  echo "   $name"
  run_case "$name"
done

echo "-- 2. 20 RK3 steps with the gas entering: one line, every evaluation in the end line --"
write_ini 20 rk3
echo "   count"
run_case count

echo "-- 3. the gas entering faster: a line at each factor of 10 --"
write_ini 40 euler05
echo "   tenfold"
run_case tenfold

echo "-- 4. a run that stops at its first stage: the line before the stop --"
write_ini 1 stop
echo "   stop"
python3 -B entering_gas.py case stop || exit 1
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/stop.log 2>&1
rc=$?
if python3 -B entering_gas.py stop stop OUTPUT/stop.log "$rc"; then
  echo "  [ok  ] stop"
else
  echo "  [FAIL] stop"
  fail=1
fi

echo "-- 5. 130 records on a face, 1 and 4 threads: the same lines --"
write_ini 1 euler
echo "   threads"
ok=1
skip=0
for nt in 1 4; do
  OMP_NUM_THREADS=$nt solve threads "OUTPUT/threads_$nt.log" || { ok=0; continue; }
  # a MOSE built without OpenMP writes "Serial execution" in place of the number of threads: it has one thread only
  if [[ $nt -eq 1 ]] && grep -q "^ Serial execution\$" OUTPUT/threads_1.log; then
    skip=1
    break
  fi
  if ! grep -q "Number of threads --> *$nt\$" "OUTPUT/threads_$nt.log"; then
    echo "  [FAIL] threads: the run did not use $nt threads (see OUTPUT/threads_$nt.log)"
    ok=0
  fi
done
if [[ $skip -eq 1 ]]; then
  echo "  [skip] threads: MOSE built without OpenMP (OUTPUT/threads_1.log: Serial execution)"
elif [[ $ok -eq 1 ]] && python3 -B entering_gas.py same OUTPUT/threads_1.log OUTPUT/threads_4.log; then
  echo "  [ok  ] threads"
else
  echo "  [FAIL] threads"
  fail=1
fi

echo "-- 6. a case without 400 faces: no line and no end line --"
write_ini 1 euler
echo "   no400"
run_case no400

if [[ $fail -ne 0 ]]; then
  echo "Extrapolation400: FAILED"
  exit 1
fi
if [[ $skip -eq 1 ]]; then
  echo "Extrapolation400: all checks passed; threads not run: MOSE built without OpenMP"
else
  echo "Extrapolation400: all checks passed"
fi
exit 0
