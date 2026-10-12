#!/usr/bin/env bash
#===============================================================================
#  Turbulence source terms in time-accurate runs (fast tier)
#
#  A 2-D channel with a wall, periodic along the flow (make_case.py), with
#  Spalart-Allmaras and with SST, each with Heun (RK2) and RK3:
#    order  nu~, k and omega converge with the order of the time scheme:
#           100 steps at cfl 0.8, 200 at 0.4 and 400 at 0.2 end at the same
#           time (the step is set by cells whose state does not change), and
#           the observed order must be at least 1.8 with Heun and 2.7 with RK3;
#    warn   one step of SST-RC with RK3 on the case of SST, with no
#           cfl-rise-threshold: the start-up warnings of the k-omega models and
#           of the rotation-curvature correction must both be printed;
#    pos    near-wall turbulence far above its equilibrium decays faster than
#           the transport acts: no stage may make nu~, k or omega negative
#           (count printed by the solver at the end of the run), and the step
#           must have been limited by the turbulence sources.
#  and, on its own case,
#    wilcox the vortex-stretching function of Wilcox 2006: a cube of 6^3
#           cells, one block, a velocity whose every component varies
#           linearly along the other two coordinates (divergence-free: with
#           uniform density and pressure the flux of rho omega through the
#           faces of a cell adds up to zero), uniform k and omega, one
#           explicit step. The change of rho omega must be
#           dt [2 gamma rho S_ij S_ij - beta_0 f_beta rho omega^2], with
#           f_beta from chi_w = |W_ij W_jk S^_ki|/(beta* omega)^3 and
#           S^_ki = S_ki - (1/2) du_m/dx_m delta_ki (NASA Turbulence Modeling
#           Resource, Wilcox 2006 page), within 1e-6 of it. Here S^ is zero
#           on the diagonal: a S^ without its off-diagonal part gives
#           f_beta = 1 instead of 0.85 and a change 39 % larger.
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

# solve <case> <time scheme> <cfl> <vnn> <steps> <tag>: time-accurate run of the case written by make_case.py;
# the final field is kept as OUTPUT/<tag>.tec and the solver output as OUTPUT/<tag>.log
solve() {
  local model=${1#*-} wall=asymptotic
  [[ $1 == pos-SST ]] && wall=practical
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = $5
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = $model
schmidt = 0.7

[MOSE-Turbulence]
omega-wall-bc = $wall

[MOSE-Numerics]
cfl = $3
vnn = $4
time-accurate = true
time-scheme = $2
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-diter = 1000000000
shell-diter = 1
EOF
  rm -f OUTPUT/field.tec
  "$LOCAL" > "OUTPUT/$6.log" 2>&1
  if [[ ! -f OUTPUT/field.tec ]] || grep -a -qE 'NaN|Integration failed' "OUTPUT/$6.log"; then
    echo "  [FAIL] the run blew up (see OUTPUT/$6.log)"
    fail=1
    return 1
  fi
  local done_steps
  done_steps=$(grep -a -oE 'Iter = +[0-9]+' "OUTPUT/$6.log" | tail -1 | grep -a -oE '[0-9]+$')
  if [[ "$done_steps" != "$5" ]]; then
    echo "  [FAIL] the run stopped after ${done_steps:-0} of $5 steps (see OUTPUT/$6.log)"
    fail=1
    return 1
  fi
  cp OUTPUT/field.tec "OUTPUT/$6.tec"
}

echo "-- start-up warnings, SST-RC, RK3, no cfl-rise-threshold, one step --"
python3 -B make_case.py order-SST || exit 1
if solve order-SST-RC RK3 0.8 0.4 1 warn-SST-RC; then
  if grep -a -q '^\[WARNING\] k-omega, time-accurate, no cfl-rise-threshold' OUTPUT/warn-SST-RC.log &&
     grep -a -q '^\[WARNING\] The rotation-curvature correction (-RC)' OUTPUT/warn-SST-RC.log; then
    echo "  [ok  ] the k-omega and the -RC warnings are both printed"
  else
    echo "  [FAIL] a start-up warning is missing (see OUTPUT/warn-SST-RC.log)"
    fail=1
  fi
fi

for model in SA SST; do
  python3 -B make_case.py "order-$model" || exit 1
  for scheme in RK2 RK3; do
    echo "-- order, $model, $scheme: 100 steps at cfl 0.8, 200 at 0.4, 400 at 0.2 --"
    tag="order-$model-$scheme"
    solve "order-$model" "$scheme" 0.8 0.4 100 "$tag-1" || continue
    solve "order-$model" "$scheme" 0.4 0.2 200 "$tag-2" || continue
    solve "order-$model" "$scheme" 0.2 0.1 400 "$tag-4" || continue
    if python3 -B make_case.py order "$model" "$scheme" "OUTPUT/$tag-1.tec" "OUTPUT/$tag-2.tec" "OUTPUT/$tag-4.tec"; then
      echo "  [ok  ] order of the time scheme"
    else
      echo "  [FAIL] order below that of the time scheme"
      fail=1
    fi
  done
done

for model in SA SST; do
  python3 -B make_case.py "pos-$model" || exit 1
  for scheme in RK2 RK3; do
    echo "-- positivity, $model, $scheme, cfl 0.8, vnn 0.4, 200 steps --"
    tag="pos-$model-$scheme"
    solve "pos-$model" "$scheme" 0.8 0.4 200 "$tag" || continue
    if python3 -B make_case.py pos "pos-$model" 0.8 0.4 "OUTPUT/$tag.log"; then
      echo "  [ok  ] positive"
    else
      echo "  [FAIL] not shown positive"
      fail=1
    fi
  done
done

echo "-- Wilcox 2006, vortex stretching: a 3-D flow, one step --"
python3 -B make_case.py wilcox || exit 1
cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 1
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = navier-stokes
turbulence = Wilcox2006
schmidt = 0.7

[MOSE-Numerics]
cfl = 0.1
vnn = 0.1
time-accurate = true
time-scheme = euler
integration-variables = cons
space-reconstruction = first-order
riemann-solver = HLLC
EOF
rm -f OUTPUT/field.tec OUTPUT/wilcox.tec
"$LOCAL" > OUTPUT/wilcox.log 2>&1
[[ -f OUTPUT/field.tec ]] && cp OUTPUT/field.tec OUTPUT/wilcox.tec
if python3 -B make_case.py wilcox-check OUTPUT/wilcox.tec; then
  echo "  [ok  ] the change of rho omega is the one of f_beta"
else
  echo "  [FAIL] the change of rho omega is not the one of f_beta"
  fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "TurbulenceSources: FAILED"
  exit 1
fi
echo "TurbulenceSources: all checks passed"
exit 0
