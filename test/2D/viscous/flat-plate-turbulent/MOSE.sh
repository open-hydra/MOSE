#!/bin/bash -
#===============================================================================
#
#          FILE: MOSE.sh
#
#         USAGE: run "./MOSE.sh [options]" in the current shell
#
#   DESCRIPTION: A script to compile and run MOSE
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " compile               -->     compile accordingly with the preset file"
  echo " solve                 -->     run MOSE"
  echo " test                  -->     run all turbulence models (coarse + fine SA)"
  echo " kill                  -->     kill the process" 
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -p | --parallel <n>   -->     launch solver with <n> threads"
  echo
  exit 1
}

# Directories and files definition
MASTERDIR=../../../../
MASTER=$MASTERDIR/bin/MOSE
LOCAL=./bin/MOSE

# Default Options
BG=0
NTHREADS=1

# Parse command-line options
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    # detect argument termination
    shift
    break
  fi
  case $1 in

    -b | --background)
        #echo " -> Background"
        BG=1
        shift
        ;;
    -p | --parallel)
        #echo " -> OpenMP($2)"
        NTHREADS=$2
        shift 2
        ;;
    -h | --help)
        print_usage
        shift
        ;;
    * )
      break
      ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

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
  ln -sfn ../../../../common/Air-viscous/phase.txt phase.txt
  ln -sfn ../../../../common/Air-viscous/thermo.dat thermo.dat
  ln -sfn ../../../../common/Air-viscous/transport.dat transport.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  # Check executable
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
  # Run the solver
  if [[ $BG == 0 ]]; then
    $LOCAL
  else
    $LOCAL 2>errors_file >logfile &
    echo $! > .ID
  fi
fi
   
if [[ $1 == test ]]; then
  cd INPUT
  ln -sfn ../../../../common/Air-viscous/phase.txt phase.txt
  ln -sfn ../../../../common/Air-viscous/thermo.dat thermo.dat
  ln -sfn ../../../../common/Air-viscous/transport.dat transport.dat
  cd ..
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  # Check executable
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
  # Input file
  input_file="input.ini"

  # IC and BC files per model (coarse)
  declare -A coarse_ic=([SA]="ic-coarse-SA.tec" [SST]="ic-coarse-kw.tec" [Wilcox2006]="ic-coarse-kw.tec")
  declare -A coarse_bc=([SA]="bc-coarse-SA.txt" [SST]="bc-coarse-kw.txt" [Wilcox2006]="bc-coarse-kw.txt")

  # --- Coarse runs ---
  models=(SA SST Wilcox2006)
  for scheme in "${models[@]}"; do
      echo "Running coarse $scheme"

      # Link IC and BC to coarse files
      ln -sf "${coarse_ic[$scheme]}" INPUT/ic.tec
      ln -sf "${coarse_bc[$scheme]}" INPUT/bc.txt

      # Set turbulence model
      perl -i -pe "s/^turbulence\s*=\s*.*/turbulence = $scheme/" "$input_file"

      # Run the solver
      if [[ $BG == 0 ]]; then
        $LOCAL
      else
        $LOCAL 2>errors_file >logfile &
        echo $! > .ID
      fi

      mv OUTPUT/wall.tec OUTPUT/wall-coarse-$scheme.tec

      echo "Completed coarse run for $scheme"
      echo "----------------------------------"
  done

  # --- Fine SA run ---
  echo "Running fine SA"
  ln -sf "ic-fine-SA.tec" INPUT/ic.tec
  ln -sf "bc-fine-SA.txt" INPUT/bc.txt
  perl -i -pe "s/^turbulence\s*=\s*.*/turbulence = SA/" "$input_file"

  if [[ $BG == 0 ]]; then
    $LOCAL
  else
    $LOCAL 2>errors_file >logfile &
    echo $! > .ID
  fi

  mv OUTPUT/wall.tec OUTPUT/wall-fine-SA.tec

  echo "Completed fine SA run"
  echo "----------------------------------"
fi

if [[ $1 == kill ]]; then
read PID < .ID && kill $PID
fi
