#!/usr/bin/env bash
#===============================================================================
#  Inlet boundary conditions (fast tier)
#
#  1. Far-field inlet with the flow direction "normal". The flow direction of
#     an inlet record is an angle or the token normal, the face normal
#     (docs/user/boundary-conditions.md, common fields). A planar duct with an
#     ambient boundary (407) where the gas enters at 50 m/s runs 20 steps twice
#     (farfield_normal.py): with alpha = beta = normal, and with the angles the
#     token stands for on that face, alpha = beta = 0. The two runs must give
#     the same field to the bit, with no transverse velocity: |v| below the
#     bound of rounding, 64 eps (|u| + c) (the estimate is in
#     farfield_normal.py).
#  2. Inlets at total conditions: when the gas enters, and the face state. A
#     1-D tube of air, one forward Euler step with first-order reconstruction,
#     so that the density change of a boundary cell gives the mass flux of its
#     face (total_conditions.py). Reservoir at T0 = T unless said otherwise,
#     flow direction normal:
#       farfield36 ambient (407) on both ends, p0 = p_back = 36 p, gas at rest:
#                  the reservoir drives the gas in, through a choked face
#                  (normal Mach 1);
#       farfield2  the same with p0 = p_back = 2 p: subsonic inflow;
#       inlet2     inlet (401) on both ends, p0 = 2 p, gas at rest: subsonic
#                  inflow;
#       inlet36    inlet (401) on both ends, p0 = 36 p, gas at rest: choked
#                  inflow;
#       nofeed     inlet (401) with p0 = 0.2 p and the gas moving in at 0.1 a:
#                  p0 is below the pressure the gas of the cell reaches at
#                  rest, the reservoir cannot feed the face, which is a slip
#                  wall;
#       outlet05   ambient (407) with p0 = p_back = 0.4 p and the gas moving in
#                  at 0.5 a: p0 is below p_stop = 0.478 p, the face is an
#                  outlet at p_back and lets gas out;
#       outlet12   the same at 1.2 a with p0 = p_back = 0.1 p (p_stop = 0.147 p);
#       outletstop as outlet05 with p_back = 0.6 p, above p_stop: the face is at
#                  rest at p_stop, no mass through it;
#       inlet6     inlet (401) with p0 = p and the gas moving in at 6 a, at
#                  2a/(gamma-1) or faster: no pressure brings it to rest
#                  (p_stop = 0), the reservoir feeds the face, which is choked;
#       farfield6  the same with an ambient face (407), p0 = p_back = p;
#       inlethot   inlet (401) with T0 = 25 T, p0 = 0.1 p and the gas moving in
#                  at 5.5 a: a subsonic face (normal Mach 0.83);
#       inletvacuum inlet (401) with p0 = p and the gas moving in at 20 a,
#                  faster than the gas of the reservoir can follow: choked;
#       riemplus, riemminus  as inlethot with the invariant a + (gamma-1)/2 w
#                  of the cell at +1e-9 a and -1e-9 a.
#     Each mass flux must equal that of the face state within 1e-9 of rho a,
#     and the mass fluxes of riemplus and riemminus must agree within 1e-9 of
#     rho a: the face state is continuous where the invariant changes sign.
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

fail=0

echo "-- 1. far-field inlet with the flow direction normal --"
cat > input.ini <<'INI'
[MOSE-Parameters]
iter-threshold = 20
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
INI

ffn=0
for dir in normal angle; do
  python3 -B farfield_normal.py case "$dir" || exit 1
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$dir.log" 2>&1
  if [[ $? -ne 0 || ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] the run with the direction $dir did not complete (see OUTPUT/$dir.log)"
    ffn=1
    continue
  fi
  cp OUTPUT/field.tec "OUTPUT/field-$dir.tec"
  if vline=$(python3 -B farfield_normal.py check "OUTPUT/field-$dir.tec"); then
    echo "  [ok  ] inlet direction $dir: $vline"
  else
    echo "  [FAIL] inlet direction $dir, a transverse velocity above rounding: $vline"
    ffn=1
  fi
done

if [[ $ffn -eq 0 ]]; then
  if cmp -s OUTPUT/field-normal.tec OUTPUT/field-angle.tec; then
    echo "  [ok  ] normal and the equivalent angles give the same field to the bit"
  else
    echo "  [FAIL] normal and the equivalent angles give different fields"
    ffn=1
  fi
fi
[[ $ffn -eq 0 ]] || fail=1

echo "-- 2. inlets at total conditions --"
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

pairok=1
for name in farfield36 farfield2 inlet2 inlet36 nofeed outlet05 outlet12 outletstop \
            inlet6 farfield6 inlethot inletvacuum riemplus riemminus; do
  echo "   $name"
  python3 -B total_conditions.py case "$name" || exit 1
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$name.log" 2>&1
  if [[ $? -ne 0 || ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] the run did not complete (see OUTPUT/$name.log)"
    fail=1
    [[ $name == riem* ]] && pairok=0
    continue
  fi
  cp OUTPUT/field.tec "OUTPUT/$name.tec"
  if python3 -B total_conditions.py check "$name" "OUTPUT/$name.tec"; then
    echo "  [ok  ] $name: the mass fluxes are those of the face state"
  else
    echo "  [FAIL] $name: the mass fluxes are not those of the face state"
    fail=1
  fi
done

echo "   riemplus and riemminus"
if [[ $pairok -eq 1 ]] && python3 -B total_conditions.py pair riemplus OUTPUT/riemplus.tec riemminus OUTPUT/riemminus.tec; then
  echo "  [ok  ] the mass flux is continuous where the invariant of the cell changes sign"
else
  echo "  [FAIL] the mass flux jumps where the invariant of the cell changes sign"
  fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "Inlets: FAILED"
  exit 1
fi
echo "Inlets: all checks passed"
exit 0
