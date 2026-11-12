#!/usr/bin/env bash
#===============================================================================
#  Preconditioning: inviscid flow past a cylinder at low Mach number
#  (validation tier)
#
#  The case of make_case.py: an O-grid of 96 x 48 cells, the free stream of air
#  at p = 1e5 Pa, T = 300 K, first-order steady runs of 3000 iterations
#  (forward Euler, local time steps, cfl 0.8):
#    prec-m2  integration-variables = prec, HLLC-PC, Mach 0.01, eps-min 1e-3
#    prec-m3  the same at Mach 0.001, eps-min 1e-4 (Ur below 1e-3 a: the
#             update takes the Ur of the flux and of the time step)
#    cons-m2  integration-variables = cons, HLLC, Mach 0.01
#  Checks:
#    1. the two preconditioned runs converge: the density residual falls below
#       1e-2 of its first value (about 2e-3 here; it levels off there, the
#       small unsteadiness of the wake);
#    2. they give the same Cp within 2e-3 near the cylinder (r < 1.5 R): the
#       incompressible limit, Mach-independent to O(M^2);
#    3. on the front half of the cylinder the preconditioned Cp is within 0.3
#       of potential flow (the dissipation of the scheme separates the flow
#       at the back, and the wake changes the front by up to about 0.2);
#    4. the unpreconditioned run (cons, HLLC) differs from it by more than 1
#       in Cp: its dissipation scales with the speed of sound, and it loses the
#       incompressible limit (pressure errors of order 1/M in Cp).
#  plot_vv.py draws the figures of docs/vv/2D-low-mach-cylinder.md from the runs.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
AIR=../../../common/Air

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin INPUT runs
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
rm -f INPUT/phase.txt INPUT/thermo.dat
ln -sf "../$AIR/phase.txt"  INPUT/phase.txt
ln -sf "../$AIR/thermo.dat" INPUT/thermo.dat

fail=0

# check <status> <description>: status 0 = pass
check() {
  if [[ $1 -eq 0 ]]; then
    echo "  [ok  ] $2"
  else
    echo "  [FAIL] $2"
    fail=1
  fi
}

# run <tag> <Mach number> <integration variables> <riemann solver> <eps-min>: solve, keep OUTPUT as runs/<tag>
run() {
  local tag="$1" mach="$2"
  MACH="$mach" python3 -B make_case.py case || exit 1
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 3000
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.8
time-accurate = false
time-scheme = euler
integration-variables = $3
space-reconstruction = first-order
riemann-solver = $4
preconditioning-eps-min = $5

[MOSE-IO]
sol-format = tecplot ascii
sol-diter = 1000000000
shell-diter = 1000000000
EOF
  rm -rf OUTPUT "runs/$tag"
  mkdir -p OUTPUT
  "$LOCAL" > "runs/$tag.log" 2>&1
  if [[ ! -f OUTPUT/field.tec ]] || grep -a -qE 'NaN|Integration failed|ERROR' "runs/$tag.log"; then
    echo "  [FAIL] the run $tag did not complete (see runs/$tag.log)"
    fail=1
    return 1
  fi
  mv OUTPUT "runs/$tag"
}

echo "-- preconditioned runs at Mach 0.01 and 0.001 --"
run prec-m2 0.01  prec HLLC-PC 1e-3
run prec-m3 0.001 prec HLLC-PC 1e-4
for tag in prec-m2 prec-m3; do
  [[ -f runs/$tag/residual-history.dat ]] || continue
  python3 -B make_case.py converged runs/$tag/residual-history.dat 1e-2
  check $? "$tag converges"
done
if [[ -f runs/prec-m2/field.tec && -f runs/prec-m3/field.tec ]]; then
  python3 -B make_case.py same runs/prec-m2/field.tec 0.01 runs/prec-m3/field.tec 0.001 2e-3
  check $? "the Cp field does not depend on the Mach number (incompressible limit)"
  python3 -B make_case.py front runs/prec-m2/field.tec 0.01 0.3
  check $? "the front half is within 0.3 of potential flow"
fi

echo "-- unpreconditioned run at Mach 0.01 --"
run cons-m2 0.01 cons HLLC 0.05
if [[ -f runs/cons-m2/field.tec && -f runs/prec-m2/field.tec ]]; then
  python3 -B make_case.py differs runs/cons-m2/field.tec 0.01 runs/prec-m2/field.tec 0.01 1
  check $? "cons with HLLC loses the incompressible limit (Cp off by more than 1)"
fi

if [[ $fail -ne 0 ]]; then
  echo "Preconditioning: FAILED"
  exit 1
fi
echo "Preconditioning: all checks passed"
exit 0
