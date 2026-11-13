#!/usr/bin/env bash
#===============================================================================
#  Ablative wall boundary conditions -- BC 503 (melting), 504 (pyrolysis) and the SRM grain 502
#
#  Runs the blown-column case once per wall model and checks each result against
#  the analytic solution (see build_case.py and verify.py). Every configuration
#  is generated from scratch, so nothing but the scripts needs to be stored here.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# The walls of surface reactions (505) and the pyrolysing ablator (506) run on the gas CHO-N2-viscous, each in gsi/<cfg>/
GSI_CONFIGS="carbon carbon-top carbon-o2 ablator ablator-top ablator-phi0 ablator-hv0"
GSI_REFUSED="ablator-nofile ablator-pyro2 ablator-cond"
CONFIGS=${CONFIGS:-"melting htpb hdpe pp melting-top hdpe-top melting-inert srm srm-top srm-inflow srm-inflow-top default-schmidt $GSI_CONFIGS $GSI_REFUSED"}
status=0

for cfg in $CONFIGS; do
  [[ "$cfg" == default-schmidt || "$cfg" == melting-inert ]] && continue   # below, in their own directories
  [[ " $GSI_CONFIGS $GSI_REFUSED " == *" $cfg "* ]] && continue             # below, in gsi/<cfg>/
  python3 -B build_case.py "$cfg" > /dev/null || { echo "AblatingWall: FAIL - could not build '$cfg'"; exit 1; }

  rm -f OUTPUT/wall.tec OUTPUT/field.tec
  if ! ./MOSE.sh solve > OUTPUT/solve-"$cfg".log 2>&1; then
    echo "AblatingWall [$cfg]: FAIL - solver exited non-zero"
    tail -20 OUTPUT/solve-"$cfg".log
    status=1
    continue
  fi
  # The burning-grain BC is not a viscous wall (nor an inlet, its inflow mode), so it
  # writes no wall file; every other configuration must produce both.
  need_wall=1
  [[ "$cfg" == srm* ]] && need_wall=0
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

# The walls of surface reactions (505) and of the pyrolysing ablator (506): the gas test/common/CHO-N2-viscous, each
# configuration in gsi/<cfg>/ with the input.ini of this case at 1000 iterations. The checks are relations of the wall
# model at the state reached (verify.py, verify_gsi_reactions): the column has no steady state, its oxidisers are
# consumed, and 1000 iterations are enough for a blown column with its oxidisers still there.
run_gsi() {   # run_gsi <cfg>: build gsi/<cfg>, run it at one thread, verify it; status=1 on a failure
  local cfg=$1 D=gsi/$1
  mkdir -p "$D/INPUT" "$D/OUTPUT"
  sed -E 's/^(iter-threshold[[:space:]]*=[[:space:]]*)[0-9]+/\11000/' input.ini > "$D/input.ini"
  for f in phase.txt thermo.dat transport.dat; do
    ln -sfn "../../../../../../common/CHO-N2-viscous/$f" "$D/INPUT/$f"
  done
  rm -f "$D/OUTPUT/wall.tec" "$D/OUTPUT/field.tec"
  if ! (cd "$D" && python3 -B ../../build_case.py "$cfg" > /dev/null); then
    echo "AblatingWall [$cfg]: FAIL - could not build '$cfg'"
    status=1
  elif ! (cd "$D" || exit 1; ulimit -s unlimited 2>/dev/null; KMP_STACKSIZE=100M OMP_NUM_THREADS=1 \
          ../../../../../../bin/MOSE) > "$D/OUTPUT/solve.log" 2>&1; then
    echo "AblatingWall [$cfg]: FAIL - solver exited non-zero"
    tail -20 "$D/OUTPUT/solve.log"
    status=1
  elif [[ ! -f "$D/OUTPUT/field.tec" || ! -f "$D/OUTPUT/wall.tec" ]]; then
    echo "AblatingWall [$cfg]: FAIL - solver produced no wall/field output"
    tail -20 "$D/OUTPUT/solve.log"
    status=1
  else
    (cd "$D" && python3 -B ../../verify.py) || status=1
  fi
}
for cfg in $GSI_CONFIGS; do
  [[ " $CONFIGS " == *" $cfg "* ]] && run_gsi "$cfg"
done

# 506 with phi = 0 and h_v0 = 0 is 505 bit for bit (one thread): ablator-phi0 against carbon, the same column; the control,
# ablator-hv0 (h_v0 = -8.44e5 J/kg), must differ from carbon, or the comparison could not see a change
if [[ " $CONFIGS " == *" carbon "* && " $CONFIGS " == *" ablator-phi0 "* && " $CONFIGS " == *" ablator-hv0 "* ]]; then
  if cmp -s gsi/carbon/OUTPUT/field.tec gsi/ablator-phi0/OUTPUT/field.tec && \
     cmp -s gsi/carbon/OUTPUT/wall.tec gsi/ablator-phi0/OUTPUT/wall.tec; then
    echo "   [ok  ] 506 with phi = 0 and h_v0 = 0: the field and wall files of 505 (carbon) byte for byte"
  else
    echo "   [FAIL] 506 with phi = 0 and h_v0 = 0 differs from 505 (gsi/carbon, gsi/ablator-phi0)"
    status=1
  fi
  if cmp -s gsi/carbon/OUTPUT/wall.tec gsi/ablator-hv0/OUTPUT/wall.tec; then
    echo "   [FAIL] 506 with h_v0 = -8.44e5 J/kg gives the wall file of 505: the comparison above could not see a change"
    status=1
  else
    echo "   [ok  ] 506 with h_v0 = -8.44e5 J/kg moves the wall file of 505 (control of the comparison above)"
  fi
fi

# A 506 that blows the pyrolysis gas (phi > 0) and whose table cannot serve it stops at set-up, exit code 1, with the
# message and no solution: no INPUT/pyrolysis.dat, pyro_id beyond its tables, a table with a condensed-carbon column
refuse_gsi() {   # refuse_gsi <cfg> <message>
  local cfg=$1 D=gsi/$1
  mkdir -p "$D/INPUT" "$D/OUTPUT"
  sed -E 's/^(iter-threshold[[:space:]]*=[[:space:]]*)[0-9]+/\11000/' input.ini > "$D/input.ini"
  for f in phase.txt thermo.dat transport.dat; do
    ln -sfn "../../../../../../common/CHO-N2-viscous/$f" "$D/INPUT/$f"
  done
  rm -f "$D/INPUT/pyrolysis.dat" "$D/OUTPUT/field.tec"
  if ! (cd "$D" && python3 -B ../../build_case.py "$cfg" > /dev/null); then
    echo "AblatingWall [$cfg]: FAIL - could not build '$cfg'"
    status=1
    return
  fi
  (cd "$D" || exit 1; ulimit -s unlimited 2>/dev/null; KMP_STACKSIZE=100M OMP_NUM_THREADS=1 timeout 60 \
   ../../../../../../bin/MOSE) > "$D/OUTPUT/solve.log" 2>&1
  local rc=$?
  if [[ $rc -eq 1 && ! -f "$D/OUTPUT/field.tec" ]] && grep -a -q -F -- "$2" "$D/OUTPUT/solve.log"; then
    echo "   [ok  ] $cfg: stopped at set-up, exit code 1"
  else
    echo "   [FAIL] $cfg: exit code $rc, not refused with the message (see $D/OUTPUT/solve.log)"
    status=1
  fi
}
[[ " $CONFIGS " == *" ablator-nofile "* ]] && refuse_gsi ablator-nofile \
  "[ERROR] INPUT/pyrolysis.dat: the pyrolysis-gas tables of the walls 506 with phi > 0 cannot be read"
[[ " $CONFIGS " == *" ablator-pyro2 "* ]] && refuse_gsi ablator-pyro2 \
  "(record 25, type 506): pyro_id is 2, but INPUT/pyrolysis.dat holds 1 table(s)"
[[ " $CONFIGS " == *" ablator-cond "* ]] && refuse_gsi ablator-cond \
  "has a condensed-carbon column (treatment condensed): the wall 506 blows the gas of a table of treatment gas-only or strip"

echo
if [[ $status -ne 0 ]]; then
  echo "AblatingWall: FAILED"
else
  echo "AblatingWall: PASSED"
fi
exit $status
