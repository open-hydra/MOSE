#!/bin/bash -
#===============================================================================
#          FILE: MOSE.sh
#   DESCRIPTION: compile / run MOSE for the Shu-Osher case
#===============================================================================
function print_usage {
  echo
  echo "Tasks"
  echo " compile               -->     recompile MOSE"
  echo " solve                 -->     run MOSE"
  echo " genic [N]             -->     regenerate INPUT/ic.tec & bc.txt (N cells)"
  echo " kill                  -->     kill the background process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

MASTERDIR=../../..
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

DIR=$(pwd)

if [[ $1 == genic ]]; then
  python3 gen_ic.py "${2:-200}"
fi

if [[ $1 == compile ]]; then
  mkdir -p bin
  rm -f $LOCAL
  cd $MASTERDIR
  ./install.sh compile
  cd $DIR
  cp $MASTER $LOCAL
fi

if [[ $1 == solve ]]; then
  cd INPUT
  ln -sf ../../../common/Air/phase.txt phase.txt
  ln -sf ../../../common/Air/thermo.dat thermo.dat
  ln -sf ic_x1.tec ic.tec        # coarse grid (N=200) for the fast single-grid check
  ln -sf bc_x1.txt bc.txt
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi
  if [[ $BG == 0 ]]; then
    $LOCAL
  else
    $LOCAL 2>errors_file >logfile &
    echo $! > .ID
  fi
fi

if [[ $1 == test ]]; then
  # Grid-convergence study: run the staged meshes (ic_x1/x2/x4 -> N=200/400/800)
  # and leave OUTPUT/field_x{r}.tec for convergence.py to compare vs the
  # resolved reference (reference/reference.dat, N=1600).
  cd INPUT
  ln -sf ../../../common/Air/phase.txt phase.txt
  ln -sf ../../../common/Air/thermo.dat thermo.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  if [[ "$MASTER" -nt "$LOCAL" ]]; then cp $MASTER $LOCAL; fi
  for r in 1 2 4; do
    cd INPUT
    ln -sf ic_x$r.tec ic.tec
    ln -sf bc_x$r.txt bc.txt
    cd ..
    $LOCAL
    mv OUTPUT/field.tec OUTPUT/field_x$r.tec
  done
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
