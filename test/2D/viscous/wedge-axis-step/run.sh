#!/usr/bin/env bash
#===============================================================================
#  Local time step next to the axis of an axisymmetric wedge (validation tier)
#
#  The viscous stress on the wedge faces adds to the cells next to the axis a
#  hoop term, -(2 mu + lambda) u_r / (rho r^2), and couplings to the
#  neighbours; the time step must include them. The diffusion-dominated wedge
#  of make_case.py (prandtl = 10: momentum diffusion limits the step), local
#  time stepping, 2000 iterations at the largest vnn the step rule declares
#  stable: RK3 at vnn 0.62 and Heun at vnn 0.5 on square cells, RK3 at vnn
#  0.62 on cells twice as long in x as in r, and RK3 at vnn 0.62 on square
#  cells of a mesh with the axis on the last j face (negative wedge angle).
#  Then the vnn the guide recommends (Heun 0.4, RK3 0.5) on cells four times
#  as long in x as in r, with the Prandtl number of air (0.72) and with 0.9,
#  as where the eddy viscosity dominates. The flow swirls a little, so the
#  row of u_theta takes part too. Every run must complete with no NaN.
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

fail=0
# solve <scheme> <vnn> <NX> <NR> <LX> <LR> <tag> [flip or -] [prandtl, 10 by default]
solve() {
  python3 -B make_case.py case "$3" "$4" "$5" "$6" $([[ ${8:-} == flip ]] && echo flip) || exit 1
  cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 2000
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
schmidt = 0.7
prandtl = ${9:-10}

[MOSE-Numerics]
cfl = 0.5
vnn = $2
cfl-rise-threshold = 0
time-accurate = false
time-scheme = $1
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-diter = 2000
shell-diter = 500
EOI
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$7.log" 2>&1
  local done_steps
  done_steps=$(grep -a -oE 'Iter = +[0-9]+' "OUTPUT/$7.log" | tail -1 | grep -a -oE '[0-9]+$')
  if grep -a -qE 'NaN|Integration failed' "OUTPUT/$7.log" || [[ "$done_steps" != "2000" ]]; then
    echo "  [FAIL] $7: unstable, stopped after ${done_steps:-0} of 2000 iterations (see OUTPUT/$7.log)"
    fail=1
  else
    echo "  [ok  ] $7: 2000 iterations, no NaN"
  fi
}

echo "-- wedge next to the axis, prandtl 10, local time step --"
solve RK3 0.62 32 32 1.0e-4 1.0e-4 square-RK3
solve RK2 0.5  32 32 1.0e-4 1.0e-4 square-Heun
solve RK3 0.62 16 32 1.0e-4 1.0e-4 long-RK3
solve RK3 0.62 32 32 1.0e-4 1.0e-4 flip-RK3 flip
echo "-- cells four times as long in x, the recommended vnn, prandtl 0.72 and 0.9 --"
solve RK2 0.4 32 32 4.0e-4 1.0e-4 long4-air-Heun - 0.72
solve RK3 0.5 32 32 4.0e-4 1.0e-4 long4-air-RK3 - 0.72
solve RK2 0.4 32 32 4.0e-4 1.0e-4 long4-pr0.9-Heun - 0.9
solve RK3 0.5 32 32 4.0e-4 1.0e-4 long4-pr0.9-RK3 - 0.9

if [[ $fail -ne 0 ]]; then
  echo "WedgeAxisStep: FAILED"
  exit 1
fi
echo "WedgeAxisStep: all checks passed"
exit 0
