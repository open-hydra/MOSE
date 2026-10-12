#!/usr/bin/env bash
#===============================================================================
#  Spalart-Allmaras: right-hand side on constructed fields (fast tier)
#
#  1. Variable density: a constructed periodic field at rest with variable
#     temperature and nu~ and no wall (make_case.py): one explicit Euler step
#     with local time stepping gives, in OUTPUT/diagnostic.tec, the right-hand
#     side of the equation for rho nu~ that the solver evaluates. It must
#     converge to the exact right-hand side of the standard conservation form
#     (NASA Turbulence Modeling Resource, Spalart-Allmaras page) with second
#     order, on 32 x 32 and 64 x 64 cells.
#  2. Rotation tensor of SA-RC and SST-RC: W_ij = (du_i/dx_j - du_j/dx_i)/2
#     (same page) is a tensor, so the right-hand side of rho nu~, and of rho k
#     and rho omega with SST-RC, must not depend on how the flow is turned in
#     space. A cube of 8^3 cells, one block (-RC needs one), with a 3-D
#     velocity whose vorticity has three components and whose strain varies
#     in space, so that the term r~ of the rotation-curvature correction is
#     not zero, and the same field turned by the rotation that takes x to y,
#     y to z and z to x; one explicit step each, with SA-RC and with SST-RC
#     (k = 1 m2/s2, omega = 3000 1/s). The change of rho nu~ (rho k and rho
#     omega) must be the same in corresponding cells within 1e-9 of its
#     largest value. A rotation tensor whose W_23 has the sign opposite to
#     that of W_12 and W_13 (+omega_x/2 instead of -omega_x/2) makes the two
#     changes differ by their own size.
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../../bin/MOSE
LOCAL=./bin/MOSE
VISC=../../../../common/Air-viscous

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

fail=0
errors=()
for n in 32 64; do
  python3 -B make_case.py case "$n" || exit 1
  rm -f OUTPUT/diagnostic.tec
  "$LOCAL" > "OUTPUT/n$n.log" 2>&1
  if [[ ! -f OUTPUT/diagnostic.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/n$n.log"; then
    echo "  [FAIL] the run on $n x $n cells did not complete (see OUTPUT/n$n.log)"
    echo "MMSResiduals: FAILED"
    exit 1
  fi
  cp OUTPUT/diagnostic.tec "OUTPUT/diagnostic-n$n.tec"
  errors+=("$(python3 -B make_case.py error "OUTPUT/diagnostic-n$n.tec")")
done

echo "-- 1. right-hand side of rho nu~ against the standard conservation form --"
if python3 -B make_case.py check "${errors[0]}" "${errors[1]}"; then
  echo "  [ok  ] second-order agreement"
else
  echo "  [FAIL] the solver does not evaluate the standard form"
  fail=1
fi

for model in SA-RC SST-RC; do
  echo "-- 2. $model: a 3-D flow and the same flow turned, one step --"
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = $model

[MOSE-Numerics]
cfl = 0.1
time-accurate = true
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
EOF
  for name in a b; do
    python3 -B make_case.py rotation "$name" "$model" || exit 1
    rm -f OUTPUT/field.tec "OUTPUT/field-$model-$name.tec"
    "$LOCAL" > "OUTPUT/rotation-$model-$name.log" 2>&1
    if [[ -f OUTPUT/field.tec ]]; then
      cp OUTPUT/field.tec "OUTPUT/field-$model-$name.tec"
    else
      echo "  (rotation $model $name: no OUTPUT/field.tec, see OUTPUT/rotation-$model-$name.log)"
    fi
  done
  if python3 -B make_case.py rotation-check "OUTPUT/field-$model-a.tec" "OUTPUT/field-$model-b.tec"; then
    echo "  [ok  ] $model: the change of the turbulence variables does not depend on how the flow is turned"
  else
    echo "  [FAIL] $model: the change of the turbulence variables depends on how the flow is turned"
    fail=1
  fi
done

if [[ $fail -ne 0 ]]; then
  echo "MMSResiduals: FAILED"
  exit 1
fi
echo "MMSResiduals: all checks passed"
exit 0
