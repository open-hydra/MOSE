#!/bin/bash -
#===============================================================================
#          FILE: MOSE.sh
#   DESCRIPTION: turbulent flow through a planar converging-diverging nozzle
#                on grid families that differ only in the direction of their
#                grid lines (O, H, ZZ, SB2), with the Spalart-Allmaras and SST
#                models
#
#   make_case.py builds the mesh, INPUT/ic.tec and INPUT/bc.txt of each case
#   before its run (numpy only); the number of iterations of the form is
#   written into a copy of input-<model>.ini.  The Air-viscous phase and
#   thermodynamics and the Sutherland viscosity of test/common fix the gas.
#
#   testfast (short form, for CTest): 300 iterations on the coarse level of
#     the O and H grids, SA and SST, from the quasi-1D initial field: the
#     solution must stay mirror-symmetric about the centre line, which the
#     four blocks with different index orientations do not impose; and the
#     ZZ grid with SA, which must give a finite flow;
#   test (full form): O, H, ZZ and SB2 on the three levels, SA and SST, up to
#     the residual threshold of input-<model>.ini (long: use -p).
#===============================================================================
function print_usage {
  echo
  echo "Tasks"
  echo " solve                 -->     run the O grid, level 0, SA, to convergence"
  echo " test                  -->     full form: O, H, ZZ, SB2 on levels 0-2, SA and SST, to convergence"
  echo " testfast              -->     short form: symmetry of O and H (SA, SST) and ZZ (SA), 300 iterations, level 0"
  echo " kill                  -->     kill the background process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

MASTERDIR=../../../../
MASTER=$MASTERDIR/bin/MOSE
LOCAL=./bin/MOSE

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
  ln -sf ../../../../common/Air-viscous/phase.txt phase.txt
  ln -sf ../../../../common/Air-viscous/thermo.dat thermo.dat
  ln -sf ../../../../common/Air-sutherland/transport.dat transport.dat
  cd ..
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi

  if [[ $1 == solve ]]; then
    RUNS="O:0:SA"; ITER=200000
  elif [[ $1 == testfast ]]; then
    RUNS="O:0:SA H:0:SA O:0:SST H:0:SST ZZ:0:SA"; ITER=300
  else
    RUNS=""
    for MODEL in SA SST; do for LEV in 0 1 2; do for FAM in O H ZZ SB2; do RUNS="$RUNS $FAM:$LEV:$MODEL"; done; done; done
    ITER=200000
  fi
  echo "$1" > OUTPUT/form.txt

  # verify.py reads every OUTPUT/field_<case>.tec present: drop those of earlier runs
  rm -f OUTPUT/field.tec OUTPUT/wall.tec OUTPUT/field_*.tec OUTPUT/wall_*.tec
  for RUN in $RUNS; do
    IFS=: read -r FAM LEV MODEL <<< "$RUN"
    python3 make_case.py --fam $FAM --level $LEV --model $MODEL || exit 1
    sed "s/^iter-threshold.*/iter-threshold = $ITER/" input-$MODEL.ini > input.ini
    if [[ $BG == 0 ]]; then
      $LOCAL
    else
      $LOCAL 2>errors_file >logfile
    fi
    mv OUTPUT/field.tec OUTPUT/field_$FAM-L$LEV-$MODEL.tec || exit 1
    mv OUTPUT/wall.tec OUTPUT/wall_$FAM-L$LEV-$MODEL.tec || exit 1
  done
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
