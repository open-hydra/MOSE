#!/bin/bash -
#===============================================================================
#          FILE: MOSE.sh
#   DESCRIPTION: run the 2D viscous Method-of-Manufactured-Solutions order study (N2)
#
#   Uses the dedicated MOSE_mms driver (src/app/mms.f90), which injects the
#   analytic manufactured source term making a smooth periodic field an EXACT
#   steady Navier-Stokes solution.  The meshes/ICs/BCs for each resolution are
#   pre-generated and committed (INPUT/ic-<N>.tec, INPUT/bc-<N>.txt).  To
#   regenerate them (needs the Hydra GRIB + ATLAS tools) set nx=ny=<N> in
#   input.ini, then:
#       GRIB meshgen            # -> mesh.tec
#       python3 build_ic.py     # -> ic.tec   (manufactured field)
#       ATLAS BCB               # -> fromATLAStoSolver/bc.txt (periodic)
#   and copy ic.tec/bc.txt to INPUT/ic-<N>.tec / INPUT/bc-<N>.txt.
#
#   The flat transport table (INPUT/transport.dat, constant mu=10) and the
#   Air phase/thermo (linked from ../../../../common/Air) fix the constants that
#   the baked source in mms.f90 must match.
#===============================================================================
function print_usage {
  echo
  echo "Tasks"
  echo " solve                 -->     run the middle grid (N=64) only"
  echo " test                  -->     run all grids (32/64/128) for the order study"
  echo " testfast              -->     run the fast subset (32/64) for CTest"
  echo " kill                  -->     kill the background process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

MASTERDIR=../../../../
MASTER=$MASTERDIR/bin/MOSE_mms
LOCAL=./bin/MOSE_mms

BG=0
NTHREADS=1
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then shift; break; fi
  case $1 in
    -b | --background) BG=1; shift ;;
    -p | --parallel)   NTHREADS=$2; shift 2 ;;
    -h | --help)       print_usage ;;
    * ) break ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

if [[ $1 == solve || $1 == test || $1 == testfast ]]; then
  cd INPUT
  ln -sf ../../../../common/Air/phase.txt phase.txt
  ln -sf ../../../../common/Air/thermo.dat thermo.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi

  if [[ $1 == solve ]]; then
    GRIDS="64"
  elif [[ $1 == testfast ]]; then
    GRIDS="32 64"
  else
    GRIDS="32 64 128"
  fi

  for N in $GRIDS; do
    cd INPUT
    ln -sf ic-$N.tec ic.tec
    ln -sf bc-$N.txt bc.txt
    cd ..
    if [[ $BG == 0 ]]; then
      $LOCAL
    else
      $LOCAL 2>errors_file >logfile
    fi
    mv OUTPUT/field.tec OUTPUT/field_$N.tec
  done
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
