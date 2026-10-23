#!/usr/bin/env bash
#===============================================================================
#  Time control (fast tier)
#
#  How a run uses its time and when it ends: a 1-D subsonic duct (100 cells,
#  stagnation inlet, back-pressure outlet) for 1 to 3, a 1-D supersonic tube
#  for 4; inputs written by make_case.py:
#    1. order in time with a time-dependent inlet (BC 402, p0(t) =
#       p0 (1 + 1e-4 sin(2 pi t / Tp))): every Runge-Kutta stage must see the
#       boundary data of its own time t^n + c_k dt, which gives second order
#       for RK2 and RK3 (Carpenter, Gottlieb, Abarbanel, Don, ICASE 93-83,
#       1993); data of the end of the step at every stage give first order.
#       Runs with dt0/2^k (k = 0..3) over one period, reference dt0/64;
#       observed order of the last pair >= 1.8 for p and u;
#    2. a time table is held at its first value before its first instant:
#       BC 402 with a constant table gives the same field as BC 401 with the
#       same p0, in a steady run (time -1) and in a time-accurate run that
#       starts before the table;
#    3. res-threshold does not end a time-accurate run: a run started from
#       the steady state, with res-threshold = 2, which would end a steady
#       run at its first iteration, reaches its time-threshold;
#    4. res-threshold ends a steady run when the density residual has fallen
#       by that factor from the first iteration, whatever the level of
#       density, and a restart ends where the run in one go ends. A tube of
#       air at Mach 2 (40 cells, extrapolation on the i faces) carries an
#       entropy wave out: steady, forward Euler, res-threshold = 1e-6. The
#       Euler equations do not change when density and pressure are scaled
#       by the same factor at fixed velocity and temperature, and with an
#       even power of 2 the scaling is exact in floating point, the square
#       roots of the density in the Roe averages included: the same run at
#       p = 1 bar and at p = 1 bar x 2^-24 takes the same iterations, every
#       residual scaled. The run at low pressure must not end before its
#       residual has fallen by 1e-6, and the two runs must end at the same
#       iteration (a residual compared with res-threshold as it is, in
#       kg/m3, ends the run at low pressure after its first iteration). The
#       run at 1 bar in two parts, 60 iterations and a restart from the vtk
#       raw solution of the first part, which holds the residual of the
#       first iteration (MOSE_RES_REFERENCE), must end at the same iteration
#       as the run in one go, with the same ratio of the last residual to
#       the first.
#  First-order reconstruction keeps the semi-discrete operator smooth.
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

N0=40         # steps per period of the inlet signal at dt0
CFL0=0.5      # Courant number of dt0
EPS=1.0e-4    # relative amplitude of p0(t)
DT0=$(python3 -B make_case.py dt0 "$CFL0") || exit 1
TP=$(python3 -c "print(repr($N0 * $DT0))")

fail=0

# solve <tag> <time-accurate> <scheme> <cfl> <MOSE-Parameters lines...>
solve() {
  local tag="$1" accurate="$2" scheme="$3" cfl="$4"
  shift 4
  {
    printf '[MOSE-Parameters]\n'
    local line
    for line in "$@"; do
      printf '%s\n' "$line"
    done
    printf '\n[MOSE-Physics]\nequations = euler\n\n[MOSE-Numerics]\n'
    printf 'cfl = %s\ntime-accurate = %s\ntime-scheme = %s\n' "$cfl" "$accurate" "$scheme"
    printf 'integration-variables = cons\nspace-reconstruction = first-order\nriemann-solver = HLLC\n'
  } > input.ini
  rm -f OUTPUT/field.tec "OUTPUT/$tag.tec"
  "$LOCAL" > "OUTPUT/$tag.log" 2>&1
  if [[ -f OUTPUT/field.tec ]]; then
    mv OUTPUT/field.tec "OUTPUT/$tag.tec"
  else
    echo "  ($tag: no OUTPUT/field.tec, see OUTPUT/$tag.log)"
  fi
}

python3 -B make_case.py uniform || exit 1

echo "-- 1. order in time with a time-dependent inlet --"
python3 -B make_case.py bc402 "$TP" "$EPS" || exit 1
for scheme in RK2 RK3; do
  for k in 0 1 2 3; do
    f=$((1 << k))
    solve "${scheme}_$k" true "$scheme" "$(python3 -c "print($CFL0 / $f)")" \
          "iter-threshold = $((N0 * f))" "res-threshold = 1e-300"
  done
  solve "${scheme}_ref" true "$scheme" "$(python3 -c "print($CFL0 / 64)")" \
        "iter-threshold = $((N0 * 64))" "res-threshold = 1e-300"
done
python3 -B check.py order 1.8 RK2 RK3 || fail=1

echo "-- 2. time table before its first instant --"
python3 -B make_case.py bc401 || exit 1
solve steady-401   false RK3 0.5 "iter-threshold = 50" "res-threshold = 1e-300"
solve accurate-401 true  RK3 0.5 "iter-threshold = 20" "res-threshold = 1e-300"
python3 -B make_case.py bc402 "$TP" 0.0 || exit 1
solve steady-402   false RK3 0.5 "iter-threshold = 50" "res-threshold = 1e-300"
python3 -B make_case.py bc402 "$TP" 0.0 "$TP" || exit 1
solve accurate-402 true  RK3 0.5 "iter-threshold = 20" "res-threshold = 1e-300"
python3 -B check.py same steady-401 steady-402 1e-12 || fail=1
python3 -B check.py same accurate-401 accurate-402 1e-12 || fail=1

echo "-- 3. res-threshold in a time-accurate run --"
python3 -B make_case.py bc401 || exit 1
TEND=$(python3 -c "print(repr(20.5 * $DT0))")
solve accurate-steady true RK3 0.5 "time-threshold = $TEND" "res-threshold = 2"
python3 -B check.py reached accurate-steady "$TEND" || fail=1

echo "-- 4. res-threshold in a steady run --"
SPLIT=60
# solve_steady <tag> <iter-threshold> <newrun>: the tube, steady, with vtk raw solutions; the residual
# history of the run is kept in OUTPUT/residual-history-<tag>.dat (removed first: no check reads an older run)
solve_steady() {
  rm -f "OUTPUT/$1.log" "OUTPUT/residual-history-$1.dat"
  printf '%s\n' '[MOSE-Parameters]' "iter-threshold = $2" "newrun = $3" 'res-threshold = 1e-6' '' \
                '[MOSE-Physics]' 'equations = euler' '' \
                '[MOSE-Numerics]' 'cfl = 0.5' 'time-accurate = false' 'time-scheme = euler' \
                'integration-variables = cons' 'space-reconstruction = first-order' 'riemann-solver = HLLC' '' \
                '[MOSE-IO]' 'sol-format = vtk raw' > input.ini
  "$LOCAL" > "OUTPUT/$1.log" 2>&1 || echo "  ($1: exit code $?, see OUTPUT/$1.log)"
  cp OUTPUT/residual-history.dat "OUTPUT/residual-history-$1.dat" 2>/dev/null || echo "  ($1: no residual-history.dat)"
}
for tag in bar tiny; do
  python3 -B make_case.py tube "$tag" || exit 1
  rm -rf OUTPUT/residual-history.dat OUTPUT/field.vtm OUTPUT/vtk
  solve_steady "$tag" 3000 true
done
# The run at 1 bar in two parts: the restart reads the final solution of the first part, and no other
python3 -B make_case.py tube bar || exit 1
rm -rf OUTPUT/residual-history.dat OUTPUT/field.vtm OUTPUT/vtk
solve_steady split-first "$SPLIT" true
solve_steady split-restart 3000 false
python3 -B check.py threshold 1e-6 "$SPLIT" || fail=1

if [[ $fail -ne 0 ]]; then
  echo "TimeControl: FAILED"
  exit 1
fi
echo "TimeControl: all checks passed"
exit 0
