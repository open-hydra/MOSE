#!/usr/bin/env bash
#===============================================================================
#  One flux on the faces between blocks (fast tier)
#
#  A periodic 1-D domain written as one block or as two blocks (make_case.py):
#  a shock tube cut at its discontinuities, so that a shock leaves the block cut
#  and the periodic connection, or a low Mach number wave crossing them. A face
#  between two blocks has one flux, identical on the two sides and equal to
#  that of an interior face, so:
#    - the two-block run gives the one-block run bit for bit;
#    - mass, momentum and energy are conserved to round-off (time-accurate runs).
#  Riemann solvers: HLLC (no shock sensor), HLLC+Tramel and HLLC+Chen (shock
#  sensor) on the shock tube; HLLC-PC (preconditioned, steady pseudo-time
#  with local steps) on a low Mach number entropy wave and pressure pulse:
#  bit-for-bit check only. Gas: test/common/Air. Then the shock tube with the
#  Navier-Stokes equations and HLLC+Tramel (test/common/Air-viscous), which
#  also carries the viscous flux of the faces: conservation only (a face
#  between blocks has its own viscous stencil, not that of an interior face).
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
AIR=../../../common/Air
VISC=../../../common/Air-viscous

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin INPUT runs
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi

# gas <directory>: link the phase, thermo (and transport) files
gas() {
  rm -f INPUT/phase.txt INPUT/thermo.dat INPUT/transport.dat
  ln -sf "../$1/phase.txt"  INPUT/phase.txt
  ln -sf "../$1/thermo.dat" INPUT/thermo.dat
  [[ -f "$1/transport.dat" ]] && ln -sf "../$1/transport.dat" INPUT/transport.dat
  return 0
}

# write_ini <riemann solver> [equations]: 40 RK3 steps, time-accurate (HLLC-PC: steady, local steps,
# preconditioned variables)
write_ini() {
  local accurate=true vars=cons extra="" eq=${2:-euler}
  if [[ $1 == HLLC-PC ]]; then accurate=false; vars=prec; fi
  if [[ $1 == HLLC+Chen ]]; then extra="riemann-options-Mco = 0.2"; fi
  cat > input.ini <<EOF
[MOSE-Parameters]
iter-threshold = 40
newrun = true
res-threshold = 1e-300

[MOSE-Physics]
equations = $eq
schmidt = 0.7

[MOSE-Numerics]
cfl = 0.5
time-accurate = $accurate
time-scheme = RK3
integration-variables = $vars
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = $1
$extra

[MOSE-IO]
sol-format = vtk raw
sol-diter = 1000000000
shell-diter = 1000000000
ini-diter = 1000000000
EOF
}

# run <tag>: solve the case in INPUT/, keep OUTPUT as runs/<tag>
run() {
  rm -rf OUTPUT "runs/$1"
  mkdir -p OUTPUT/vtk
  "$LOCAL" > "runs/$1.log" 2>&1
  if [[ ! -f OUTPUT/field.vtm ]] || grep -a -qE 'NaN|Integration failed|ERROR' "runs/$1.log"; then
    echo "  [FAIL] the run $1 did not complete (see runs/$1.log)"
    return 1
  fi
  mv OUTPUT "runs/$1"
}

fail=0
gas "$AIR"
for rs in HLLC HLLC+Tramel HLLC+Chen HLLC-PC; do
  echo "-- $rs --"
  write_ini "$rs"
  ok=1
  kind=sod; [[ $rs == HLLC-PC ]] && kind=wave
  for layout in one two; do
    python3 -B make_case.py case "$layout" "$kind" || exit 1
    run "$rs-$layout" || { ok=0; fail=1; }
  done
  [[ $ok -eq 1 ]] || continue
  echo "  two blocks against one block:"
  python3 -B make_case.py same "runs/$rs-one" "runs/$rs-two" || fail=1
  if [[ $rs != HLLC-PC ]]; then
    for layout in one two; do
      echo "  conservation, $layout block(s):"
      python3 -B make_case.py conserved "runs/$rs-$layout" || fail=1
    done
  fi
done

echo "-- HLLC+Tramel, Navier-Stokes --"
gas "$VISC"
write_ini HLLC+Tramel navier-stokes
ok=1
for layout in one two; do
  python3 -B make_case.py case "$layout" sod || exit 1
  run "ns-$layout" || { ok=0; fail=1; }
done
if [[ $ok -eq 1 ]]; then
  for layout in one two; do
    echo "  conservation, $layout block(s):"
    python3 -B make_case.py conserved "runs/ns-$layout" || fail=1
  done
fi

if [[ $fail -ne 0 ]]; then
  echo "ConnectionFlux: FAILED"
  exit 1
fi
echo "ConnectionFlux: all checks passed"
exit 0
