#!/usr/bin/env bash
#===============================================================================
#  Riemann smoke-matrix (fast tier, F5)
#
#  Runs a short coarse Sod shock tube once per Riemann solver and asserts that
#  density and pressure stay positive and finite.  This is a crash/robustness
#  smoke test that guarantees every solver at least executes — not an accuracy
#  check.  It reuses the coarse 1-D Sod inputs (test/1D/Sod79/INPUT) so no new
#  mesh/IC/BC data is stored here.  The two outer branches of HLLC: the Sod
#  tube moved by u = +3, where the flux of every face is that of its left state
#  (S_L >= 0), and its mirror image moved by u = -3, where it is that of its
#  right state (S_R <= 0); the second must be the mirror image of the first
#  (mirror_sod.py).
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../1D/Sod79
COMMON=../../common/Air

# Shock-capturing / all-speed solvers for which a strong Sod shock tube is an
# appropriate robustness vehicle.  Each must run without producing NaN/p<0.
SOLVERS=(
  HLLC HLLC+Tramel
  HLLE HLLE++
  SLAU SLAU2
  AUSM+ AUSM+M
  LLF Rusanov
  exact
)

# NOTE — low-Mach / preconditioned solvers are intentionally NOT in the Sod
# smoke set: they are designed for smooth low-Mach flow and diverge on a strong
# shock tube by construction (and HLLC-PC additionally requires
# integration-variables = prec).  They are exercised by the low-Mach Gresho
# vortex case instead.  See test/TESTING-PLAN.md (F5) for the plan to add a
# dedicated low-Mach smoke covering:  HLLC+Chen  HLLC-PC  LMRoe  MiczekRoe

# --- environment -----------------------------------------------------------
ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# --- stage inputs (coarse Sod) ---------------------------------------------
mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"        INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"       INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec"     INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt"     INPUT/bc.txt

# --- sweep -----------------------------------------------------------------
fail=0
for s in "${SOLVERS[@]}"; do
  cat > input.ini <<EOF
[MOSE-Parameters]
time-threshold = 0.04

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = $s
riemann-options-Mco = 0.2
EOF

  rm -f OUTPUT/field.tec
  if ! "$LOCAL" > "OUTPUT/solve_${s//[^A-Za-z0-9]/_}.log" 2>&1; then
    echo "  [$s] FAIL: solver exited non-zero"
    fail=1
    continue
  fi
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [$s] FAIL: no OUTPUT/field.tec produced"
    fail=1
    continue
  fi
  if ! python3 -B check_positivity.py "$s"; then
    fail=1
  fi
done

# --- HLLC: the left state (S_L >= 0) and the right state (S_R <= 0) ---------
echo "  [HLLC] the Sod tube moved by u = +3 and its mirror image moved by u = -3"
for sign in plus minus; do
  python3 -B mirror_sod.py case "$sign" 3.0 || { fail=1; break; }
  sed 's/^riemann-solver = .*/riemann-solver = HLLC/' input.ini > input-mirror.ini && mv input-mirror.ini input.ini
  rm -f OUTPUT/field.tec
  if ! "$LOCAL" > "OUTPUT/solve_HLLC_$sign.log" 2>&1 || [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [HLLC $sign] FAIL: solver exited non-zero or wrote no field"
    fail=1
    break
  fi
  cp OUTPUT/field.tec "OUTPUT/field-$sign.tec"
done
rm -f INPUT/ic.tec
ln -sf "../$SOD/INPUT/ic_x1.tec" INPUT/ic.tec
if [[ -f OUTPUT/field-plus.tec && -f OUTPUT/field-minus.tec ]]; then
  python3 -B mirror_sod.py compare OUTPUT/field-plus.tec OUTPUT/field-minus.tec || fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "RiemannSmoke: FAILED"
  exit 1
fi
echo "RiemannSmoke: all $(( ${#SOLVERS[@]} )) solvers passed; HLLC gives the mirror image on its two outer branches"
exit 0
