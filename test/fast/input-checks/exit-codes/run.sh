#!/usr/bin/env bash
#===============================================================================
#  Exit codes of the error stops (fast tier)
#
#  A run that MOSE stops on an error must end with a non-zero exit code and
#  its message, and write no solution. Reuses the inputs of the 100-cell Sod
#  tube (test/1D/Sod79) and the gas of test/1D/diff1; no new mesh or IC data
#  is stored here.
#
#  Default Schmidt number. schmidt defaults to 0.0, i.e. mixture-averaged
#  multicomponent diffusion, which needs the binary diffusion table
#  INPUT/diffusion.dat when the gas has more than one species (one species
#  does not diffuse: AblatingWall runs a single-species column without
#  schmidt). The Sod tube as a viscous mixture of four species (those of
#  test/1D/diff1, with the mass fractions of air), without schmidt and
#  without a table, must stop at set-up with the diffusion-table error and
#  write no solution.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../../bin/MOSE
LOCAL=./bin/MOSE
SOD=../../../1D/Sod79/INPUT
MIX=../../../1D/diff1/INPUT
ERR='Sc<=0 (multicomponent diffusion) requires a binary diffusion table'

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

mkdir -p bin OUTPUT INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi

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

echo "-- default schmidt: four species, no diffusion table --"
rm -f INPUT/*
for f in phase.txt thermo.dat transport.dat; do
  ln -sf "../$MIX/$f" "INPUT/$f"
done
ln -sf "../$SOD/bc_x1.txt" INPUT/bc.txt
python3 -B make_case.py mixture "$SOD/ic_x1.tec" INPUT/ic.tec || exit 1
cat > input.ini <<EOI
[MOSE-Parameters]
iter-threshold = 2
newrun = true

[MOSE-Physics]
equations = navier-stokes

[MOSE-Numerics]
cfl = 0.5
time-accurate = true
time-scheme = RK3
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC
EOI
rm -f OUTPUT/field.tec
"$LOCAL" > OUTPUT/schmidt-mixture.log 2>&1
rc=$?
[[ $rc -ne 0 ]]
check $? "the solver stops with an error (exit status $rc)"
grep -a -qF "$ERR" OUTPUT/schmidt-mixture.log
check $? "it reports the missing diffusion table"
[[ ! -f OUTPUT/field.tec ]]
check $? "no solution is written"

if [[ $fail -ne 0 ]]; then
  echo "ExitCodes: FAILED"
  exit 1
fi
echo "ExitCodes: all checks passed"
exit 0
