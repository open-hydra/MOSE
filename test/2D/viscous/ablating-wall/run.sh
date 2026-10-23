#!/usr/bin/env bash
#===============================================================================
#  Ablative wall boundary conditions -- BC 503 (melting) and 504 (pyrolysis)
#
#  Runs the blown-column case once per wall model and checks each result against
#  the analytic solution (see build_case.py and verify.py). Every configuration
#  is generated from scratch, so nothing but the scripts needs to be stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

CONFIGS=${CONFIGS:-"melting htpb hdpe pp melting-top hdpe-top melting-inert srm default-schmidt"}
status=0

for cfg in $CONFIGS; do
  [[ "$cfg" == default-schmidt || "$cfg" == melting-inert ]] && continue   # below, in their own directories
  python3 -B build_case.py "$cfg" > /dev/null || { echo "AblatingWall: FAIL - could not build '$cfg'"; exit 1; }

  rm -f OUTPUT/wall.tec OUTPUT/field.tec
  if ! ./MOSE.sh solve > OUTPUT/solve-"$cfg".log 2>&1; then
    echo "AblatingWall [$cfg]: FAIL - solver exited non-zero"
    tail -20 OUTPUT/solve-"$cfg".log
    status=1
    continue
  fi
  # The burning-grain BC is not a viscous wall, so it writes no wall file; every
  # other configuration must produce both.
  need_wall=1
  [[ "$cfg" == srm ]] && need_wall=0
  if [[ ! -f OUTPUT/field.tec ]] || [[ $need_wall -eq 1 && ! -f OUTPUT/wall.tec ]]; then
    echo "AblatingWall [$cfg]: FAIL - solver produced no wall/field output"
    tail -20 OUTPUT/solve-"$cfg".log
    status=1
    continue
  fi

  python3 -B verify.py || status=1
done

# melting-inert: the inert column runs 8000 iterations, in melting-inert/ with
# the input.ini of this case at 8000 iterations. Its last check (no blowing)
# looks at the start-up flow while it decays: in these cells the additive local
# time step is about half the smallest of the direction limits, so the column
# needs 8000 iterations to reach the pseudo-time of 4000 steps of that limit;
# the other columns pass their checks with the 4000 of input.ini.
if [[ " $CONFIGS " == *" melting-inert "* ]]; then
  D=melting-inert
  mkdir -p "$D/INPUT" "$D/OUTPUT"
  sed -E 's/^(iter-threshold[[:space:]]*=[[:space:]]*)[0-9]+/\18000/' input.ini > "$D/input.ini"
  for f in phase.txt thermo.dat transport.dat; do
    ln -sfn "../../../../../common/Air-viscous/$f" "$D/INPUT/$f"
  done
  rm -f "$D/OUTPUT/wall.tec" "$D/OUTPUT/field.tec"
  if ! (cd "$D" && python3 -B ../build_case.py melting-inert > /dev/null); then
    echo "AblatingWall [melting-inert]: FAIL - could not build 'melting-inert'"
    status=1
  elif ! (cd "$D" || exit 1; ulimit -s unlimited 2>/dev/null; KMP_STACKSIZE=100M OMP_NUM_THREADS=1 \
          ../../../../../bin/MOSE) > "$D/OUTPUT/solve.log" 2>&1; then
    echo "AblatingWall [melting-inert]: FAIL - solver exited non-zero"
    tail -20 "$D/OUTPUT/solve.log"
    status=1
  elif [[ ! -f "$D/OUTPUT/field.tec" || ! -f "$D/OUTPUT/wall.tec" ]]; then
    echo "AblatingWall [melting-inert]: FAIL - solver produced no wall/field output"
    tail -20 "$D/OUTPUT/solve.log"
    status=1
  else
    (cd "$D" && python3 -B ../verify.py) || status=1
  fi
fi

# default-schmidt: the melting column again, in default-schmidt/ with the
# input.ini of this case less its schmidt line. schmidt then takes its default,
# 0 (mixture-averaged multicomponent diffusion), which needs the binary
# diffusion table INPUT/diffusion.dat only with more than one species: this gas
# is one species, which does not diffuse, so the run must go without the table
# and give the melting column.
if [[ " $CONFIGS " == *" default-schmidt "* ]]; then
  D=default-schmidt
  mkdir -p "$D/INPUT" "$D/OUTPUT"
  grep -v -E '^[[:space:]]*schmidt[[:space:]]*=' input.ini > "$D/input.ini"
  for f in phase.txt thermo.dat transport.dat; do
    ln -sfn "../../../../../common/Air-viscous/$f" "$D/INPUT/$f"
  done
  rm -f "$D/OUTPUT/wall.tec" "$D/OUTPUT/field.tec"
  if ! (cd "$D" && python3 -B ../build_case.py melting > /dev/null); then
    echo "AblatingWall [default-schmidt]: FAIL - could not build 'melting'"
    status=1
  elif ! (cd "$D" || exit 1; ulimit -s unlimited 2>/dev/null; KMP_STACKSIZE=100M OMP_NUM_THREADS=1 \
          ../../../../../bin/MOSE) > "$D/OUTPUT/solve.log" 2>&1; then
    echo "AblatingWall [default-schmidt]: FAIL - solver exited non-zero"
    tail -20 "$D/OUTPUT/solve.log"
    status=1
  elif [[ ! -f "$D/OUTPUT/field.tec" || ! -f "$D/OUTPUT/wall.tec" ]]; then
    echo "AblatingWall [default-schmidt]: FAIL - solver produced no wall/field output"
    tail -20 "$D/OUTPUT/solve.log"
    status=1
  else
    echo
    echo "AblatingWall [default-schmidt]: the melting column without schmidt"
    (cd "$D" && python3 -B ../verify.py) || status=1
  fi
fi

echo
if [[ $status -ne 0 ]]; then
  echo "AblatingWall: FAILED"
else
  echo "AblatingWall: PASSED"
fi
exit $status
