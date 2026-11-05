#!/usr/bin/env bash
#===============================================================================
#  Multirate with several MPI processes and OpenMP threads: single rate
#  (300 steps) and multirate (safety 2, max-level 3, 100 macro steps), 1 process x 1 thread against P processes x T
#  threads: fields, wall output and residual history identical byte by byte.
#  - mr2d-plate-4b, the laminar flat plate in 2 x 2 blocks (a vertex of four
#    blocks in the boundary layer, so viscous fluxes and tangential ghost values Pg cross the processes):
#    2x1, 4x1, 1x2, 1x4, 2x2, 4x2;
#  - mr2d-plate-sa-4b, the same with Spalart-Allmaras: 2x1, 4x1, 2x2; mr2d-plate-sst-4b, with SST: 2x1, 4x1, 2x2;
#  - mr3d-hot-8b, the hot cube in 2 x 2 x 2 blocks, a vertex of eight blocks inside the cells of level 1 (viscous fluxes and
#    tangential ghost values Pg cross the processes in the three directions): 2x1, 4x1, 2x2;
#  - mr2d-plate-2br, the plate in two blocks with block 2 turned by 180 degrees about z, so that the connection
#    across the boundary layer has the orientation -1 0 0 1 (the tangential ghost values Pg of its records, read
#    by the viscous fluxes, come from the other process at 2x1): 2x1, 2x2; and its 1x1 runs against those of
#    mr2d-plate-2b, the same plate without the turn: every cell value within 1e-10 of its scale (verify_turned.py:
#    the two differ only in the order of the floating-point operations of block 2);
#  - mr2d-box-4b and mr1d-acu-3b (Euler, connections inside bands and between bands): 2 and 4, 2 and 3 processes;
#  - mr2d-plate-sa-4b with SA-RC and max-level 0 (safety 1, 100 macro steps): the stage calls the block routine of the
#    single-rate stage, whose rotation and curvature terms take the velocity gradients across the processes: 1x1
#    against 2x1;
#  - mr2d-box-4bx: the box cut one cell from the edge of a band, so a face between levels is one cell from a
#    process boundary and the coarser cells read through the connection are interpolated by the other process:
#    2 and 4 processes, and conservation of mass and energy within 1e-12 on 4 processes;
#  - mr1d-acu-3b with max-level 10, whose local steps span 3.97 octaves (depth 3), with cells of the slowest octave (the
#    10 mm band) in every block: 1 and 3 processes identical, with L = 3 in every macro step;
#  - mr1d-deep12-3b with max-level left out (10 by default): its local steps span 12.85 octaves (depth 10), those of its
#    middle block 8.56, so the depth must come from the largest local step of all the processes: 1 and 3 processes
#    identical, with L = 10 and H = 512 dt_min in every macro step;
#  - mr1d-blast100-3b (a strong blast wave in three blocks) with nonlocal-cfl = redo and force-redo: 1 and 3
#    processes identical, with at least one attempt that violates condition (4) (the redo across the processes);
#  - mr1d-blast-wide-3b (the blast tube with far coarse cells, in three blocks) with max-level and nonlocal-cfl left out
#    (10 levels, redo): states that are not physical inside attempts, which are rejected and redone; 1x1, 3x1 and 3x2
#    identical, with at least one attempt rejected for a state that is not physical (the decision to redo is collective,
#    and the process that found the state sends the others the states it holds from then on);
#  - mr1d-blast-hot-3b (the hot gas at 9990 K, the contact in block 2) with nonlocal-cfl = redo on 3 processes: every
#    attempt of the third macro step takes a face that is not physical on process 1, and the macro step falls below 2**-10
#    of its first attempt; every process prints what it holds and waits for the others before the run ends: exit code 1,
#    the lines of the stop and of the face (process 1), the line of the stop of each of the 3 processes, and no line of
#    the launcher about a process ended by a signal;
#  - mr1d-blast-wide-3b with the keys above and max-attempts = 3 (its first macro step has 4 attempts) on 1x1, 3x1 and
#    3x2: the run stops at the third rejection, with exit code 1, the same lines of the limit, of the causes (1 / 0 / 2)
#    and of the last rejection, one line of the stop per process, and no line of the launcher about a process ended by
#    a signal;
#  - mr1d-sod-pneg-3b (a cell at a negative pressure, in block 2) on 3 processes: the run stops at the start of the first
#    macro step with exit code 1, the line of the local steps once, the line of the cell (block 2, cell 5) from process 1,
#    one line of the stop per process, and no line of a process ended by a signal. The local step of that cell comes from
#    an invalid operation (the speed of sound of a negative temperature), which a build that traps floating-point
#    exceptions (DEBUG) stops there with SIGFPE, before the check: this part holds in the builds that do not trap them
#    (TESTING, RELEASE);
#  - PROBE (Linux only, where the tests build the probe library; elsewhere this part is left out): the multirate run of
#    mr2d-plate-4b on 2 x 2 again with the probe library of MPIThreadFunnel preloaded (MOSE_MPI_THREAD_PROBE): no MPI
#    call from a thread other than the main one, the exchanges seen on every rank, and
#    the largest tag of the point-to-point calls at most 32767 on every rank, the least MPI_TAG_UB the MPI standard
#    allows (MPI 4.1, Sec. 3.2.3);
#  - MAP: mr2d-plate-4b with the block-to-rank assignment imposed by the map that decomposition-map names,
#    INPUT/decomposition.map (written by
#    ../../mpi/write_map.py in the format of the ATLAS decomposer): the diagonal blocks together on 2 ranks (0 1 1 0;
#    MOSE alone gives 0 1 0 1), 1 and 2 threads, and the reverse order on 4 ranks (3 2 1 0); the log must name the map.
#  - TIMING: the timers of the multirate, on mr2d-plate-4b with blocks 1-3 on rank 0 and block 4 on rank 1 (map 0 0 0 1),
#    1 and 2 threads, so that rank 1 waits for rank 0 in the exchanges of the stages; every macro step prints its lines
#    (shell-diter 1), and the certificate of nonlocal-cfl stays on (the key is left out), so that the stage wall time per
#    level of process 0 comes in its end-of-run line. verify_timing.py on each: the run ends; the exchange wait of the
#    Timing block is above zero and at least the exchange time per level, mean over the processes and summed over the
#    levels, over Solver; that mean is above zero at some level; the exchange start of the Timing block is at least zero,
#    and with the exchange wait and the collective wait at most 100 %; the six end-of-run lines of the times per level over the processes have lm + 1 values each,
#    minimum <= mean <= maximum, and the values of process 0 lie between the minimum and the maximum. On one process (the
#    1x1 multirate run of mr2d-plate-4b): none of these lines, and no exchange, exchange start or collective line in the
#    Timing block.
#  Every run writes its residual at every step (res-diter = 1) and must run on the threads asked; the 1x1 multirate runs of
#  the comparisons with P x T have cells on more than one level; the runs of TIMING must name the map in their log.
#  The launcher is the one CMake found (MOSE_MPIEXEC, MOSE_MPIEXEC_NUMPROC_FLAG; mpirun -n by default).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
MPIEXEC=${MOSE_MPIEXEC:-mpirun}
NPFLAG=${MOSE_MPIEXEC_NUMPROC_FLAG:--n}
PROBE=${MOSE_MPI_THREAD_PROBE:-}
command -v "$MPIEXEC" > /dev/null || { echo "mr-hybrid: FAIL - MPI launcher $MPIEXEC not found"; exit 1; }
[[ -z "$PROBE" || -f "$PROBE" ]] || { echo "mr-hybrid: FAIL - probe library not found (MOSE_MPI_THREAD_PROBE='$PROBE')"; exit 1; }
mr_cases cases mr2d-plate-4b mr2d-plate-sa-4b mr2d-plate-sst-4b mr3d-hot-8b mr2d-box-4b mr1d-acu-3b mr2d-box-4bx \
  mr1d-blast100-3b mr2d-plate-2br \
  mr2d-plate-2b mr1d-deep12-3b mr1d-blast-wide-3b mr1d-blast-hot-3b mr1d-sod-pneg-3b \
  || { echo "mr-hybrid: FAIL - cannot build the cases"; exit 1; }
status=0

# run_pt RUNDIR CASEDIR P T [ 'Section|key|value' ... ] -> runs MOSE on P processes with T threads
run_pt() {
  local run=$1 case=$2 np=$3 nt=$4; shift 4
  mr_prepare "$run" "$case" "$@" || return 1
  ( cd "$run" && OMP_NUM_THREADS=$nt "$MPIEXEC" "$NPFLAG" "$np" "$MR_MASTER" > solve.log 2>&1 )
  local rc=$?
  [[ $rc -eq 0 && -f "$run/OUTPUT/field.vtm" ]] || { echo "  run $run: exit $rc"; tail -3 "$run/solve.log"; return 1; }
  mr_threads_ok "$run" "$nt" || return 1
  return 0
}

# mr_levels RUNDIR -> 0 if the first line of the cells per level of the run has cells on more than one level
mr_levels() {
  awk '/Multirate: cells per level/ {for (i = 5; i <= NF; i++) if ($i > 0) n++; f = 1; exit} END {exit !(f && n > 1)}' "$1/solve.log" \
    || { echo "  run $1: fewer than two levels with cells"; return 1; }
}

MR=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0"
    "MOSE-Parameters|iter-threshold|100" "MOSE-IO|res-diter|1")
SR=("MOSE-Parameters|iter-threshold|300" "MOSE-IO|res-diter|1")
for spec in "mr2d-plate-4b 2x1 4x1 1x2 1x4 2x2 4x2" "mr2d-plate-sa-4b 2x1 4x1 2x2" "mr2d-plate-sst-4b 2x1 4x1 2x2" \
            "mr3d-hot-8b 2x1 4x1 2x2" "mr2d-box-4b 2x1 4x1" \
            "mr1d-acu-3b 2x1 3x1" "mr2d-box-4bx 2x1 4x1" "mr2d-plate-2br 2x1 2x2"; do
  set -- $spec; c=$1; shift; confs=("$@")
  for kind in sr mr; do
    if [[ $kind == mr ]]; then keys=("${MR[@]}"); else keys=("${SR[@]}"); fi
    run_pt "runs/$c-$kind-1x1" "cases/$c" 1 1 "${keys[@]}" || { status=1; continue; }
    [[ $kind == sr ]] || mr_levels "runs/$c-$kind-1x1" || status=1
    for pt in "${confs[@]}"; do
      np=${pt%x*}; nt=${pt#*x}
      run_pt "runs/$c-$kind-$pt" "cases/$c" "$np" "$nt" "${keys[@]}" || { status=1; continue; }
      if mr_same "runs/$c-$kind-1x1" "runs/$c-$kind-$pt"; then
        echo "  $c $kind: 1x1 and $pt identical"
      else
        echo "  $c $kind: 1x1 and $pt DIFFERENT"; status=1
      fi
    done
  done
done
# the turned connection: the 1x1 runs of mr2d-plate-2br against those of mr2d-plate-2b, the same plate without the turn
for kind in sr mr; do
  if [[ $kind == mr ]]; then keys=("${MR[@]}"); else keys=("${SR[@]}"); fi
  run_pt "runs/mr2d-plate-2b-$kind-1x1" cases/mr2d-plate-2b 1 1 "${keys[@]}" || { status=1; continue; }
  echo "  mr2d-plate-2br $kind against mr2d-plate-2b, block 2 turned:"
  python3 -B verify_turned.py "runs/mr2d-plate-2b-$kind-1x1" "runs/mr2d-plate-2br-$kind-1x1" 2 || status=1
done
# SA-RC with max-level 0 on several processes: the stage calls the block routine of the single-rate stage, whose rotation
# and curvature terms take the velocity gradients of the neighbouring processes (exchange_velocity_gradient)
RC=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0"
    "MOSE-Parameters|iter-threshold|100" "MOSE-IO|res-diter|1" "MOSE-Physics|turbulence|SA-RC")
if run_pt runs/mr2d-plate-sa-4b-rc-1x1 cases/mr2d-plate-sa-4b 1 1 "${RC[@]}" && \
   run_pt runs/mr2d-plate-sa-4b-rc-2x1 cases/mr2d-plate-sa-4b 2 1 "${RC[@]}"; then
  if mr_same runs/mr2d-plate-sa-4b-rc-1x1 runs/mr2d-plate-sa-4b-rc-2x1; then
    echo "  mr2d-plate-sa-4b SA-RC, max-level 0: 1x1 and 2x1 identical"
  else
    echo "  mr2d-plate-sa-4b SA-RC, max-level 0: 1x1 and 2x1 DIFFERENT"; status=1
  fi
else
  status=1
fi
echo "  mr2d-box-4bx, conservation on 4 processes:"
python3 -B ../mr-common/verify_conservation.py cases/mr2d-box-4bx runs/mr2d-box-4bx-mr-4x1 --quantities mass,energy \
  --max 1e-12 || status=1

# the depth of the map below max-level on several processes: mr1d-acu-3b with max-level 10, whose local steps span 3.97
# octaves (depth 3), with cells of the slowest octave in every block; every macro step printed (shell-diter 1)
MR10=("${MR[@]}" "MOSE-Multirate|max-level|10" "MOSE-IO|shell-diter|1")
if run_pt runs/mr1d-acu-3b-cap-1x1 cases/mr1d-acu-3b 1 1 "${MR10[@]}" && \
   run_pt runs/mr1d-acu-3b-cap-3x1 cases/mr1d-acu-3b 3 1 "${MR10[@]}"; then
  nl=$(cat runs/mr1d-acu-3b-cap-1x1/solve.log runs/mr1d-acu-3b-cap-3x1/solve.log | command grep -a -c 'Multirate: H =')
  bad=$(cat runs/mr1d-acu-3b-cap-1x1/solve.log runs/mr1d-acu-3b-cap-3x1/solve.log | command grep -a 'Multirate: H =' \
        | command grep -a -v -E 'dt_min, L = +3,' | wc -l)
  if [[ $nl -eq 200 && $bad -eq 0 ]] && mr_same runs/mr1d-acu-3b-cap-1x1 runs/mr1d-acu-3b-cap-3x1; then
    echo "  mr1d-acu-3b max-level 10: 1x1 and 3x1 identical, depth 3 in the $nl macro steps"
  else
    echo "  mr1d-acu-3b max-level 10: 1x1 and 3x1 DIFFERENT, or $bad of $nl macro steps (200 required) without L = 3"
    status=1
  fi
else
  status=1
fi
# the largest local step is that of all the processes: mr1d-deep12-3b with max-level left out (10 by default), whose local
# steps span 12.85 octaves (depth 10); its middle block, alone on a process with 3 processes, spans 8.56 octaves (the bands
# of the levels 4 to 12 of the mesh), so a process that took the largest local step of its own blocks would build its map
# for a depth of 8 and a macro step 4 times shorter. Every macro step printed (shell-diter 1)
DEF=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|safety|2.0" "MOSE-Parameters|iter-threshold|20" "MOSE-IO|res-diter|1"
     "MOSE-IO|shell-diter|1")
if run_pt runs/mr1d-deep12-3b-def-1x1 cases/mr1d-deep12-3b 1 1 "${DEF[@]}" && \
   run_pt runs/mr1d-deep12-3b-def-3x1 cases/mr1d-deep12-3b 3 1 "${DEF[@]}"; then
  nl=$(cat runs/mr1d-deep12-3b-def-1x1/solve.log runs/mr1d-deep12-3b-def-3x1/solve.log | command grep -a -c 'Multirate: H =')
  bad=$(cat runs/mr1d-deep12-3b-def-1x1/solve.log runs/mr1d-deep12-3b-def-3x1/solve.log | command grep -a 'Multirate: H =' \
        | command grep -a -v -E '= +5\.120E\+02 dt_min, L = +10,' | wc -l)
  if [[ $nl -eq 40 && $bad -eq 0 ]] && mr_same runs/mr1d-deep12-3b-def-1x1 runs/mr1d-deep12-3b-def-3x1; then
    echo "  mr1d-deep12-3b max-level left out: 1x1 and 3x1 identical, L = 10 and H = 512 dt_min in the $nl macro steps"
  else
    echo "  mr1d-deep12-3b max-level left out: 1x1 and 3x1 DIFFERENT, or $bad of $nl macro steps (40 required) without" \
         "L = 10 and H = 512 dt_min: $(command grep -a 'Multirate: H =' runs/mr1d-deep12-3b-def-3x1/solve.log | head -1 \
         | sed 's/^ *//')"
    status=1
  fi
else
  status=1
fi

# condition (4) imposed (nonlocal-cfl = redo) on several processes: mr1d-blast100-3b, whose front crosses a connection,
# with force-redo; the decision to redo is collective and the state goes back to t_n on every process
RD=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-IO|res-diter|1"
    "MOSE-Parameters|time-threshold|4.0e-4" "MOSE-Parameters|iter-threshold|1000000"
    "MOSE-Multirate|nonlocal-cfl|redo" "MOSE-Multirate|force-redo|.true.")
run_pt runs/mr1d-blast100-3b-redo-1x1 cases/mr1d-blast100-3b 1 1 "${RD[@]}" || status=1
run_pt runs/mr1d-blast100-3b-redo-3x1 cases/mr1d-blast100-3b 3 1 "${RD[@]}" || status=1
nv=$(command grep -a 'Multirate run: certificate (redo)' runs/mr1d-blast100-3b-redo-3x1/solve.log \
     | sed -n 's/.*attempts with a violation \([0-9]*\).*/\1/p' | head -1)
if [[ -n $nv && $nv -ge 1 ]] && mr_same runs/mr1d-blast100-3b-redo-1x1 runs/mr1d-blast100-3b-redo-3x1; then
  echo "  mr1d-blast100-3b redo: 1x1 and 3x1 identical ($nv attempts with a violation)"
else
  echo "  mr1d-blast100-3b redo: 1x1 and 3x1 DIFFERENT, or no attempt with a violation (${nv:-?})"; status=1
fi

# a state that is not physical inside an attempt on several processes: mr1d-blast-wide-3b with max-level and nonlocal-cfl
# left out (10 levels, redo); the front starts in block 1
WD=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|recompute"
    "MOSE-Multirate|recompute-every|1" "MOSE-IO|res-diter|1" "MOSE-Parameters|time-threshold|6.0e-4"
    "MOSE-Parameters|iter-threshold|1000000")
if run_pt runs/mr1d-blast-wide-3b-1x1 cases/mr1d-blast-wide-3b 1 1 "${WD[@]}" && \
   run_pt runs/mr1d-blast-wide-3b-3x1 cases/mr1d-blast-wide-3b 3 1 "${WD[@]}" && \
   run_pt runs/mr1d-blast-wide-3b-3x2 cases/mr1d-blast-wide-3b 3 2 "${WD[@]}"; then
  ns=$(command grep -a 'Multirate run: redo, attempts rejected for a state that is not physical' \
       runs/mr1d-blast-wide-3b-3x1/solve.log | sed -n 's/.*not physical \([0-9]*\).*/\1/p' | head -1)
  if [[ -n $ns && $ns -ge 1 ]] && mr_same runs/mr1d-blast-wide-3b-1x1 runs/mr1d-blast-wide-3b-3x1 && \
     mr_same runs/mr1d-blast-wide-3b-1x1 runs/mr1d-blast-wide-3b-3x2; then
    echo "  mr1d-blast-wide-3b: 1x1, 3x1 and 3x2 identical ($ns attempts rejected for a state that is not physical)"
  else
    echo "  mr1d-blast-wide-3b: 1x1, 3x1 and 3x2 DIFFERENT, or no attempt rejected for a state that is not physical" \
         "(${ns:-?})"
    status=1
  fi
else
  status=1
fi

# the stop below 2**-10 of the first attempt on several processes: mr1d-blast-hot-3b with the keys of the part blast-hot of
# mr-cfl4 on 3 processes; the first face that is not physical of the last attempt is on process 1
HB=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|recompute"
    "MOSE-Multirate|recompute-every|1" "MOSE-Multirate|nonlocal-cfl|redo" "MOSE-IO|res-diter|1"
    "MOSE-Parameters|time-threshold|4.0e-4" "MOSE-Parameters|iter-threshold|1000000")
run=runs/mr1d-blast-hot-3b-3x1
if mr_prepare "$run" cases/mr1d-blast-hot-3b "${HB[@]}"; then
  ( cd "$run" && OMP_NUM_THREADS=1 "$MPIEXEC" "$NPFLAG" 3 "$MR_MASTER" > solve.log 2>&1 )
  rc=$?
  HL=$run/solve.log
  nr=$(command grep -a -c -E '\[RANK +[0-9]+\] STOP: MOSE-Multirate: macro step below 2\*\*-10' "$HL")
  nk=$(command grep -a -c -E 'BAD TERMINATION|KILLED BY SIGNAL' "$HL")
  if [[ $rc -eq 1 && $nr -eq 3 && $nk -eq 0 ]] \
     && command grep -a -q -E 'fell below 2\*\*-10 of its first attempt .* after 11 attempts' "$HL" \
     && command grep -a -q -E '\[ERROR\]   stage [0-9]+ of the attempt: .*\(process 1\)' "$HL" \
     && command grep -a -q 'face between cell' "$HL"; then
    echo "  mr1d-blast-hot-3b on 3 processes: the stop below 2**-10 with the face of process 1, every process stopped" \
         "with the others (exit 1)"
  else
    echo "  mr1d-blast-hot-3b on 3 processes: exit $rc (1 expected), lines of the stop of the processes $nr (3 expected)," \
         "lines of a process ended by a signal $nk (0 expected), or the lines of the stop or of the face missing"
    status=1
  fi
else
  status=1
fi

# max-attempts on several processes: mr1d-blast-wide-3b with the keys of the part of the states that are not physical and
# max-attempts = 3; its first macro step has 4 attempts, so the run stops at the third rejection on every process together
CAP=("${WD[@]}" "MOSE-Multirate|max-attempts|3")
first=""
for pt in 1x1 3x1 3x2; do
  np=${pt%x*}; nt=${pt#*x}; run=runs/mr1d-blast-wide-3b-cap-$pt
  mr_prepare "$run" cases/mr1d-blast-wide-3b "${CAP[@]}" || { status=1; continue; }
  ( cd "$run" && OMP_NUM_THREADS=$nt "$MPIEXEC" "$NPFLAG" "$np" "$MR_MASTER" > solve.log 2>&1 )
  rc=$?
  L=$run/solve.log
  nr=$(command grep -a -c -E '\[RANK +[0-9]+\] STOP: MOSE-Multirate: macro step rejected in all its attempts, max-attempts = 3 ' "$L")
  nk=$(command grep -a -c -E 'BAD TERMINATION|KILLED BY SIGNAL' "$L")
  lines=$(command grep -a -E '\[ERROR\] \[MOSE-Multirate\] nonlocal-cfl: the macro step was rejected|\[ERROR\]   attempts of the macro step rejected|\[ERROR\]   the last attempt' "$L")
  [[ -z $first ]] && first=$lines
  if [[ $rc -eq 1 && $nr -eq $np && $nk -eq 0 && -n $lines && $lines == "$first" ]] \
     && command grep -a -q -F '(max-attempts = 3;' "$L" \
     && command grep -a -q -F '1 for a violation of condition (4), 0 for a local step that is not a positive number, 2 for a state that is not physical' "$L"; then
    echo "  mr1d-blast-wide-3b max-attempts 3, $pt: stopped at the third rejection (causes 1 / 0 / 2), every process with the others"
  else
    echo "  mr1d-blast-wide-3b max-attempts 3, $pt: exit $rc (1 expected), lines of the stop of the processes $nr ($np expected), lines" \
         "of a process ended by a signal $nk (0 expected), or the lines of the limit missing or not those of 1x1"
    status=1
  fi
done

# a local step that is not a positive finite number on several processes: mr1d-sod-pneg-3b, the cell in block 2 (process 1)
run=runs/mr1d-sod-pneg-3b-3x1
if mr_prepare "$run" cases/mr1d-sod-pneg-3b "${MR[@]}"; then
  ( cd "$run" && OMP_NUM_THREADS=1 "$MPIEXEC" "$NPFLAG" 3 "$MR_MASTER" > solve.log 2>&1 )
  rc=$?
  L=$run/solve.log
  nr=$(command grep -a -c -E '\[RANK +[0-9]+\] STOP: MOSE-Multirate: a local step that is not a positive finite number' "$L")
  nk=$(command grep -a -c -E 'BAD TERMINATION|KILLED BY SIGNAL' "$L")
  nl=$(command grep -a -c -E '\[ERROR\] \[MOSE-Multirate\] at the start of the macro step \(time +0\.0+E\+00 s\) the local step is not a positive finite number in 1 cell;' "$L")
  if [[ $rc -eq 1 && $nr -eq 3 && $nk -eq 0 && $nl -eq 1 ]] \
     && command grep -a -q -E '\[ERROR\]   cell b, i, j, k: 2 5 1 1: local step +NaN s \(process 1\)$' "$L"; then
    echo "  mr1d-sod-pneg-3b on 3 processes: stopped at the start of the first macro step with the cell of process 1, every" \
         "process with the others (exit 1)"
  else
    echo "  mr1d-sod-pneg-3b on 3 processes: exit $rc (1 expected), lines of the stop of the processes $nr (3 expected), lines" \
         "of a process ended by a signal $nk (0 expected), lines of the local steps $nl (1 expected), or the line of the cell missing"
    status=1
  fi
else
  status=1
fi

# PROBE: MPI calls of the multirate only from the main thread (MPI_THREAD_FUNNELED), with the probe library of
# MPIThreadFunnel (MOSE_MPI_THREAD_PROBE), which the tests build on Linux only
if [[ -z $PROBE ]]; then
  echo "  PROBE: left out, no probe library (MOSE_MPI_THREAD_PROBE not set; the tests build it on Linux only)"
else
  run=runs/mr2d-plate-4b-mr-2x2-probe
  mr_prepare "$run" cases/mr2d-plate-4b "${MR[@]}" || status=1
  if ( cd "$run" && OMP_NUM_THREADS=2 "$MPIEXEC" "$NPFLAG" 2 env LD_PRELOAD="$PROBE" "$MR_MASTER" > solve.log 2> probe.log ); then
    reports=$(command grep -a -c '^MPI-THREAD-PROBE ' "$run/probe.log")
    seen=$(command grep -a '^MPI-THREAD-PROBE ' "$run/probe.log" | command grep -a -c ' mpi_isend_:[1-9]')
    bad=$(command grep -a '^MPI-THREAD-PROBE ' "$run/probe.log" | command grep -a -vc ' other=0 ')
    maxtag=$(command grep -a '^MPI-THREAD-PROBE ' "$run/probe.log" | command grep -a -o ' maxtag=-\?[0-9]*' | cut -d= -f2 \
             | sort -n | tail -1)
    if [[ $reports -eq 2 && $seen -eq 2 && $bad -eq 0 && -n $maxtag ]] && [[ $maxtag -le 32767 ]] && \
       command grep -a -qE 'Number of threads --> *2$' "$run/solve.log"; then
      echo "  PROBE: 2x2, all MPI calls from the main thread on every rank (exchanges seen on both), largest MPI tag $maxtag"
    else
      echo "  PROBE: 2x2, probe reports $reports, exchanges seen on $seen ranks, ranks with other-thread calls $bad," \
           "largest MPI tag ${maxtag:-not reported} (at most 32767): FAIL"
      command grep -a '^MPI-THREAD-PROBE ' "$run/probe.log" | sed 's/^/    /'
      status=1
    fi
  else
    echo "  PROBE: the 2x2 run with the probe failed"; tail -3 "$run/solve.log"; status=1
  fi
fi
# MAP: block-to-rank assignment imposed by a decomposition map, single rate and multirate
for spec in "list:0,1,1,0 2 1" "list:0,1,1,0 2 2" "reverse:4 4 1"; do
  set -- $spec; own=$1; np=$2; nt=$3
  for kind in sr mr; do
    if [[ $kind == mr ]]; then keys=("${MR[@]}"); else keys=("${SR[@]}"); fi
    run=runs/mr2d-plate-4b-$kind-map${np}x$nt
    mr_prepare "$run" cases/mr2d-plate-4b "${keys[@]}" "MOSE-IO|decomposition-map|decomposition.map" \
      || { status=1; continue; }
    python3 -B ../../mpi/write_map.py --layout 2x2 --cells 32 32 --owners "$own" --out "$run/INPUT/decomposition.map" \
      > /dev/null || { echo "  MAP: cannot write the map"; status=1; continue; }
    ( cd "$run" && OMP_NUM_THREADS=$nt "$MPIEXEC" "$NPFLAG" "$np" "$MR_MASTER" > solve.log 2>&1 )
    rc=$?
    if [[ $rc -ne 0 || ! -f $run/OUTPUT/field.vtm ]]; then echo "  run $run: exit $rc"; tail -3 "$run/solve.log"; status=1; continue; fi
    mr_threads_ok "$run" "$nt" || { status=1; continue; }
    command grep -a -q "MPI partition read from INPUT/decomposition.map" "$run/solve.log" \
      || { echo "  MAP: $run did not read the map"; status=1; continue; }
    if mr_same "runs/mr2d-plate-4b-$kind-1x1" "$run"; then
      echo "  MAP mr2d-plate-4b $kind: 1x1 and ${np}x$nt with the map $own identical"
    else
      echo "  MAP mr2d-plate-4b $kind: 1x1 and ${np}x$nt with the map $own DIFFERENT"; status=1
    fi
  done
done
# TIMING: the exchange and collective waits and the times per level of the multirate, on 2 ranks with unequal work
for nt in 1 2; do
  run=runs/mr2d-plate-4b-mr-timing2x$nt
  if mr_prepare "$run" cases/mr2d-plate-4b "${MR[@]}" "MOSE-IO|shell-diter|1" "MOSE-IO|decomposition-map|decomposition.map" && \
     python3 -B ../../mpi/write_map.py --layout 2x2 --cells 32 32 --owners list:0,0,0,1 --out "$run/INPUT/decomposition.map" \
       > /dev/null; then
    ( cd "$run" && OMP_NUM_THREADS=$nt "$MPIEXEC" "$NPFLAG" 2 "$MR_MASTER" > solve.log 2>&1 )
    rc=$?
    if [[ $rc -ne 0 || ! -f $run/OUTPUT/field.vtm ]]; then
      echo "  TIMING: run $run: exit $rc"; tail -3 "$run/solve.log"; status=1
    else
      mr_threads_ok "$run" "$nt" || status=1
      command grep -a -q "MPI partition read from INPUT/decomposition.map" "$run/solve.log" \
        || { echo "  TIMING: $run did not read the map"; status=1; }
      python3 -B verify_timing.py "$run/solve.log" 2 || status=1
    fi
  else
    echo "  TIMING: cannot prepare $run"; status=1
  fi
done
python3 -B verify_timing.py runs/mr2d-plate-4b-mr-1x1/solve.log 1 || status=1
if [[ $status -eq 0 ]]; then echo "mr-hybrid: ok"; else echo "mr-hybrid: FAIL"; fi
exit $status
