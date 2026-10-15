#!/usr/bin/env bash
#===============================================================================
#  The same output bits on any layout of MPI ranks and OpenMP threads
#  (needs-mpi tier)
#
#  A run must write the same bits with one thread, with N threads, with N
#  ranks and with ranks and threads together. Five small cases, each solved on
#      1x1 (reference), 1x4, 2x1, 4x1 and 2x2   (MPI ranks x OpenMP threads);
#  every file the run writes in OUTPUT/ (field, diagnostic and wall solutions
#  in vtk raw with their counters, residual-history.dat, the probe file) must
#  be identical byte for byte to that of 1x1, and the list of files the same:
#    sa-steady    a 3-D box of 2 x 2 blocks of 16 x 16 x 4 cells around a
#                 vertex the four blocks share (make_case.py): connections
#                 (101) between the blocks, inlets 404 and 407, outlets 406
#                 and 400, walls 301 and 302, symmetry planes (300) in z, a
#                 probe in block 4; Spalart-Allmaras, steady, local time step,
#                 implicit residual smoothing (irs-beta = 0.5), 30 iterations;
#    sst-steady   the same with SST;
#    sa-unsteady  the same box, Spalart-Allmaras, time-accurate, RK3, 20 steps;
#    star3        the mesh of ViscousLinearPatch with a point shared by three
#                 blocks, star5 the one with a point shared by five blocks
#                 (../../fast/viscous-patch/viscous-linear-patch/make_case.py: nodes moved
#                 at random, a linear velocity with shear and dilatation at
#                 10 Pa, outer faces 400, symmetry planes in z): Navier-Stokes,
#                 time-accurate, RK3, 40 steps, so the state around the shared
#                 point changes at every stage.
#  On 2 and 4 ranks the blocks of each case sit on different ranks, so the
#  faces between them, and the shared points, are between ranks.
#  Gas: test/common/Air-viscous.
#===============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MASTER=../../../bin/MOSE
LOCAL=./bin/MOSE
COMMON=../../common/Air-viscous
PATCH=../../fast/viscous-patch/viscous-linear-patch/make_case.py
MPIEXEC=${MOSE_MPIEXEC:-mpiexec}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}
LAYOUTS="1x4 2x1 4x1 2x2"

ulimit -s unlimited 2>/dev/null || true
export KMP_STACKSIZE=100M

mkdir -p bin INPUT
if [[ ! -x "$LOCAL" || "$MASTER" -nt "$LOCAL" ]]; then
  cp "$MASTER" "$LOCAL"
fi
ln -sf "../$COMMON/phase.txt"     INPUT/phase.txt
ln -sf "../$COMMON/thermo.dat"    INPUT/thermo.dat
ln -sf "../$COMMON/transport.dat" INPUT/transport.dat

fail() { echo "MPIBitIdentity: FAIL — $*"; exit 1; }

# write_ini <case>
write_ini() {
  case $1 in
    star3|star5)
      printf '%s\n' '[MOSE-Parameters]' 'iter-threshold = 40' 'res-threshold = 1e-300' 'newrun = true' '' \
        '[MOSE-Physics]' 'equations = navier-stokes' 'schmidt = 0.7' '' \
        '[MOSE-Numerics]' 'cfl = 0.5' 'vnn = 0.3' 'time-accurate = true' 'time-scheme = RK3' \
        'integration-variables = cons' 'space-reconstruction = MUSCL' 'flux-limiter = vanleer' 'riemann-solver = HLLC' '' \
        '[MOSE-IO]' 'sol-format = vtk raw' 'sol-diter = 1000000000' 'shell-diter = 1000000000' \
        'ini-diter = 1000000000' > input.ini
      return ;;
  esac
  local model=SA accurate=false iters=30 irs=0.5
  case $1 in
    sst-steady)  model=SST ;;
    sa-unsteady) accurate=true; iters=20; irs=0.0 ;;
  esac
  printf '%s\n' '[MOSE-Parameters]' "iter-threshold = $iters" 'res-threshold = 1e-300' 'newrun = true' '' \
    '[MOSE-Physics]' 'equations = navier-stokes' "turbulence = $model" '' \
    '[MOSE-Numerics]' 'cfl = 0.5' "time-accurate = $accurate" 'time-scheme = RK3' \
    'integration-variables = cons' 'space-reconstruction = MUSCL' 'flux-limiter = vanleer' \
    'riemann-solver = HLLC' "irs-beta = $irs" '' \
    '[MOSE-IO]' 'sol-format = vtk raw' 'sol-diter = 1000000000' 'shell-diter = 1000000000' \
    'ini-diter = 1000000000' '' \
    '[MOSE-Probes]' 'probe1 = p1' '' '[p1]' 'variables = rho(1) u v w p' 'diter = 1' \
    'index-position = 4 9 9 2' > input.ini
}

# make_inputs <case>: INPUT/ic.tec and INPUT/bc.txt of the case
make_inputs() {
  case $1 in
    sa-steady|sa-unsteady) python3 -B make_case.py SA > /dev/null ;;
    sst-steady)            python3 -B make_case.py SST > /dev/null ;;
    star3|star5)           python3 -B "$PATCH" case "$1" mom > /dev/null && cp INPUT/bc-ns.txt INPUT/bc.txt ;;
  esac
}

# run <case> <ranks> <threads>: solve and keep OUTPUT as runs/<case>/<ranks>x<threads>
run() {
  local c="$1" np="$2" nt="$3" label="$2x$3"
  rm -rf OUTPUT "runs/$c/$label"
  mkdir -p OUTPUT/vtk "runs/$c"
  if ! OMP_NUM_THREADS="$nt" "$MPIEXEC" "$NPFLAG" "$np" "$LOCAL" > "runs/$c/solve_$label.log" 2>&1; then
    tail -20 "runs/$c/solve_$label.log"; fail "solver exited non-zero ($c $label)"
  fi
  ls OUTPUT/vtk/fieldBlock*.vts > /dev/null 2>&1 || fail "no field output ($c $label)"
  grep -q "Number of threads --> *$nt\$" "runs/$c/solve_$label.log" || fail "the run did not use $nt threads ($c $label)"
  grep -q "Number of ranks *--> *$np\$" "runs/$c/solve_$label.log" || fail "the run did not use $np ranks ($c $label)"
  mv OUTPUT "runs/$c/$label"
}

# same <case> <label>: every file of runs/<case>/1x1 identical to that of runs/<case>/<label>, and the same list of files
same() {
  local c="$1" b="$2" n=0 f
  diff <(cd "runs/$c/1x1" && find . -type f | sort) <(cd "runs/$c/$b" && find . -type f | sort) > /dev/null \
    || { echo "  different file lists: $c 1x1 vs $b"; return 1; }
  while IFS= read -r f; do
    n=$((n + 1))
    cmp -s "runs/$c/1x1/$f" "runs/$c/$b/$f" || { echo "  differs: $c ${f#./} (1x1 vs $b)"; return 1; }
  done < <(cd "runs/$c/1x1" && find . -type f | sort)
  [[ $n -gt 0 ]] || { echo "  nothing compared: $c 1x1 vs $b"; return 1; }
  echo "  identical: $c 1x1 vs $b ($n files)"
}

status=0
for c in sa-steady sst-steady sa-unsteady star3 star5; do
  make_inputs "$c" || fail "case generation ($c)"
  write_ini "$c"
  run "$c" 1 1
  for cfg in $LAYOUTS; do
    run "$c" "${cfg%x*}" "${cfg#*x}"
    same "$c" "$cfg" || status=1
  done
done

if [[ $status -ne 0 ]]; then
  echo "MPIBitIdentity: FAIL"
  exit 1
fi
echo "MPIBitIdentity: ok — SA and SST steady with walls and IRS, SA time-accurate, points shared by 3 and 5 blocks: 1x1, 1x4, 2x1, 4x1 and 2x2 bit-identical"
exit 0
