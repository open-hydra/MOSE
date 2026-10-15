#!/usr/bin/env bash
#===============================================================================
#  Limiter defaults (fast tier)
#
#  The reconstruction that runs must be the one the input declares, and the
#  set-up warnings must describe it:
#    1. MUSCL without `flux-limiter` falls back to van Leer and warns; the
#       field is bit-identical to MUSCL with `flux-limiter = vanleer`;
#    2. first-order with a `flux-limiter` ignores the limiter and warns; the
#       field is bit-identical to first-order without a limiter;
#    3. consistent inputs (MUSCL with a limiter, first-order without one)
#       print no limiter warning;
#    4. the set-up summary names the limiter that runs (van Leer for MUSCL
#       without flux-limiter, none for first-order);
#    5. sanity: the first-order and MUSCL/van Leer fields differ.
#  Short coarse Sod shock tube; inputs shared with test/1D/Sod79.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../../1D/Sod79
COMMON=../../../common/Air

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"    INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"   INPUT/thermo.dat
ln -sf "../$SOD/INPUT/ic_x1.tec" INPUT/ic.tec
ln -sf "../$SOD/INPUT/bc_x1.txt" INPUT/bc.txt

# Set-up warnings printed by Setup_Space_Scheme (src/lib/numerics/space/Mod_Space.f90)
WARN_DEFAULT="MUSCL reconstruction but not a flux limiter"
WARN_IGNORED="Ignoring flux limiter"

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

# run_case <tag> <space-reconstruction> [flux-limiter]
#   keeps the final field as OUTPUT/<tag>.tec and the solver output as OUTPUT/<tag>.log
run_case() {
  local tag="$1" recon="$2" limiter="${3:-}"
  {
    printf '[MOSE-Parameters]\ntime-threshold = 0.04\n\n'
    printf '[MOSE-Physics]\nequations = euler\n\n'
    printf '[MOSE-Numerics]\ncfl = 0.5\ntime-accurate = true\ntime-scheme = RK3\n'
    printf 'integration-variables = cons\nriemann-solver = HLLC\n'
    printf 'space-reconstruction = %s\n' "$recon"
    if [[ -n "$limiter" ]]; then
      printf 'flux-limiter = %s\n' "$limiter"
    fi
  } > input.ini
  rm -f OUTPUT/field.tec "OUTPUT/$tag.tec"
  if ! "$LOCAL" > "OUTPUT/$tag.log" 2>&1; then
    echo "  [FAIL] $tag: solver exited non-zero (see OUTPUT/$tag.log)"
    fail=1
    return
  fi
  if [[ ! -f OUTPUT/field.tec ]]; then
    echo "  [FAIL] $tag: no OUTPUT/field.tec produced (see OUTPUT/$tag.log)"
    fail=1
    return
  fi
  mv OUTPUT/field.tec "OUTPUT/$tag.tec"
}

same()  { cmp -s "OUTPUT/$1.tec" "OUTPUT/$2.tec"; }
warns() { [[ "$(< "OUTPUT/$1.log")" == *"$2"* ]]; }
# summary_limiter <tag>: print the limiter named in the set-up summary line
# "   Flux limiter    <name>"; status 1 if the summary names none
summary_limiter() {
  local line
  while IFS= read -r line; do
    if [[ "$line" =~ ^\ +Flux\ limiter\ +(.*[^\ ])\ *$ ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}"
      return 0
    fi
  done < "OUTPUT/$1.log"
  return 1
}

run_case muscl-vanleer MUSCL       vanleer
run_case muscl-default MUSCL
run_case first         first-order
run_case first-vanleer first-order vanleer

echo "-- MUSCL without flux-limiter --"
same muscl-default muscl-vanleer
check $? "runs van Leer: field bit-identical to MUSCL + flux-limiter = vanleer"
warns muscl-default "$WARN_DEFAULT"
check $? "warns that van Leer is used by default"
[[ "$(summary_limiter muscl-default)" == "Van Leer" ]]
check $? "set-up summary names the van Leer limiter"

echo "-- first-order with flux-limiter = vanleer --"
same first-vanleer first
check $? "runs first order: field bit-identical to first-order without a limiter"
warns first-vanleer "$WARN_IGNORED"
check $? "warns that the flux limiter is ignored"
! summary_limiter first-vanleer > /dev/null
check $? "set-up summary names no flux limiter"

echo "-- consistent inputs --"
! warns first "$WARN_IGNORED" && ! warns first "$WARN_DEFAULT"
check $? "first-order without flux-limiter: no limiter warning"
! warns muscl-vanleer "$WARN_IGNORED" && ! warns muscl-vanleer "$WARN_DEFAULT"
check $? "MUSCL with flux-limiter: no limiter warning"

echo "-- sanity --"
[[ -f OUTPUT/first.tec && -f OUTPUT/muscl-vanleer.tec ]] && ! same first muscl-vanleer
check $? "first-order and MUSCL/van Leer fields differ"

if [[ $fail -ne 0 ]]; then
  echo "LimiterDefaults: FAILED"
  exit 1
fi
echo "LimiterDefaults: all checks passed"
exit 0
