#!/bin/bash -
#===============================================================================
#          FILE: MOSE.sh
#   DESCRIPTION: manufactured-solution order study on skewed and non-smooth 2D
#                grids
#
#   Uses the MOSE_mms driver (src/app/mms.f90), which adds the analytic source
#   of the field chosen with MOSE_MMS_FIELD, so that the field is an exact
#   steady solution of the Navier-Stokes equations.  make_case.py builds the
#   mesh, INPUT/ic.tec and INPUT/bc.txt of each case before its run (numpy
#   only): the grid family zigzag (every cell skewed, its neighbours the other
#   way) and iface (a cell-size jump of 1.5 and kinked grid lines at y = 1/2),
#   each made of two blocks that meet with permuted indices.  Cases:
#     zigzag, periodic2d   periodic on every side;
#     iface,  periodic2d   periodic on every side;
#     iface,  wall2d       periodic in x, walls at y = 0 (type 302, T = T_w)
#                          and y = 1 (type 301, the exact heat flux).
#   The flat transport table of ../INPUT/transport.dat (mu = 10) and the Air
#   phase and thermodynamics of test/common/Air fix the constants of the
#   sources in mms.f90.
#===============================================================================
function print_usage {
  echo
  echo "Tasks"
  echo " solve                 -->     run the three cases on N = 20"
  echo " test                  -->     run the three cases on N = 20, 40, 80 (full form, V&V figures)"
  echo " testfast              -->     run the three cases on N = 10, 20 (short form, for CTest)"
  echo " kill                  -->     kill the background process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

MASTERDIR=../../../../../
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

if [[ $1 == solve || $1 == test || $1 == testfast ]]; then
  mkdir -p INPUT OUTPUT bin
  cd INPUT
  ln -sf ../../../../../common/Air/phase.txt phase.txt
  ln -sf ../../../../../common/Air/thermo.dat thermo.dat
  ln -sf ../../INPUT/transport.dat transport.dat
  cd ..
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi

  if [[ $1 == solve ]]; then
    GRIDS="20"
  elif [[ $1 == testfast ]]; then
    GRIDS="10 20"
  else
    GRIDS="20 40 80"
  fi
  echo "$1" > OUTPUT/form.txt

  # verify.py reads every OUTPUT/field_<case>.tec present: drop those of earlier runs
  rm -f OUTPUT/field.tec OUTPUT/wall.tec OUTPUT/field_*.tec OUTPUT/wall_*.tec
  for CASE in zigzag:periodic2d iface:periodic2d iface:wall2d; do
    FAM=${CASE%%:*}
    FIELD=${CASE##*:}
    for N in $GRIDS; do
      python3 make_case.py --fam $FAM --dim 2 --n $N --field $FIELD || exit 1
      export MOSE_MMS_FIELD=$FIELD
      if [[ $BG == 0 ]]; then
        $LOCAL
      else
        $LOCAL 2>errors_file >logfile
      fi
      mv OUTPUT/field.tec OUTPUT/field_$FAM-$FIELD-$N.tec || exit 1
      if [[ $FIELD == wall2d ]]; then
        mv OUTPUT/wall.tec OUTPUT/wall_$FAM-$FIELD-$N.tec || exit 1
      fi
    done
  done
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
