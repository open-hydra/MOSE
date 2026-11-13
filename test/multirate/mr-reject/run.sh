#!/usr/bin/env bash
#===============================================================================
#  Multirate set-up rejections: with
#  [MOSE-Multirate] enabled = .true., every configuration outside the scope of the multirate integrator stops at set-up
#  or at the first macro step with a non-zero exit code and a message "[MOSE-Multirate] <item>". The control run (the
#  same case with the multirate on and the accepted keys) must complete. Several blocks are supported, and IRS is
#  checked on the keys right after input.ini (the set-up of MOSE also refuses
#  IRS in a time-accurate run on its own, before the first macro step). The laminar
#  Navier-Stokes equations are supported; item 28: a boundary type that BC_Flux_Record does not know
#  stops the run (only 0 has no flux, as in BC_Fluxes; the records of type 409 and 501 stop every run earlier, at
#  the reading of bc.txt, so they have no item here). The Spalart-Allmaras and SST
#  models are supported (control: the wall channel with Spalart-Allmaras runs; SA-R, k-energy-coupling and
#  sst-production = compressible run in the test mr-conservation), and the other models have their own item (29:
#  turbulence = SAcomp, whose compressibility correction the stages do not evaluate (the single-rate integrator runs it);
#  Wilcox2006, whose sources the stages do not evaluate, on the wall channel; SA-RC and SST-RC with
#  max-level 3, whose rotation and curvature terms the cell-wise sources of a level do not compute, while with
#  max-level 0 they give the single-rate step bit for bit (test mr-onelevel)). Item 13: an SRM grain in mode 1 (a
#  wall that blows the products in; in mode 2, an inflow, the grain is the inlet 404, test mr-onelevel).
#  Items 31 and 32: a restart (newrun = false) with levels = fixed (the level map is not written with the
#  solution) or with fixed-macro-step (a multiple of dt_min of the first macro step of the run) stops at set-up.
#  Item 33: cfl > 1 with max-level > 0 and nonlocal-cfl = redo, levels = recompute and no fixed-macro-step (within a
#  macro step the certificate of condition (4) of Gnedin et al. accepts a sub-step up to the step that cfl admits, and
#  Heun is linearly stable up to cfl = 1), with van Leer (33), first-order (33b) and minmod (33c). Item 36: vnn > 0.5
#  there in a Navier-Stokes run, the same for the diffusion; vnn 0.5 there, vnn 0.6 with max-level 0 and in an Euler
#  run, cfl 1.2 and vnn 0.6 with nonlocal-cfl = off and count, and cfl 1.2 with redo and levels = fixed or a
#  fixed-macro-step (the certificate only counts there) run: the set-up checks neither bound there (controls). A value of
#  nonlocal-cfl outside off, count, redo is refused by the key registry itself (allowed values).
#  Runs 34 and 34b: levels = fixed with fixed-macro-step above the step of the map (2**L/safety dt_min, L the depth of the
#  map), with max-level 3 and with max-level 10, where the depth is 3 because the local steps of the case span fewer than
#  4 octaves: the runs go on, with the warning of the cells beyond their limit once. Item 35: max-attempts = 0, below
#  the rule >= 1 of the key registry, refused by the registry itself when input.ini is read ("[ERROR] max-attempts must be
#  >= 1"). Item 6 has three forms: the shock detector implied by the Riemann solver (HLLC+Tramel) and set by the key
#  shock-detector (Tramel, Chen), stopped at the first macro step.
#  input.ini read again during the run (runs 50-53, the laminar plate, read again at iteration 30 of 40): the keys of
#  [MOSE-Multirate] are fixed for the run (the class of the key registry), so a reading that changes one is an error of
#  the reading. With ini-error left out (warn) a reading that writes max-level = 2 and cfl = 0.4 is not taken: the run ends at
#  iteration 40 with the warning that names the key, no [ERROR], and the residual history of the run without the change,
#  byte for byte; with ini-error = stop the change stops the run at iteration 30 with the message; enabled = .true. written
#  during a single-rate run is not taken either (residual history of the single-rate run without the change).
#  Runs 54-59: the keys of other sections that the multirate integrator checks at set-up, read again: a reading that
#  takes one of them out of the scope is an error of the reading too: with nonlocal-cfl = redo cfl above 1 (54; with
#  first-order, 54b; with ini-error = stop, 59) and vnn above 0.5 in the laminar plate (54c), taken with nonlocal-cfl = off
#  (54d, 54e), and
#  cfl-rise-threshold > 0 (56), whose checks the multirate integrator makes again at every reading, time-accurate = .false. (55) and, in a RANS run (the wall channel with SST), k-energy-coupling = .true.
#  (57), keys fixed for every run; a cfl within the bound is taken, and the interface margin
#  kappa = max(cfl (1 + theta), 1.1) follows it: the run with that reading gives the outputs of the run made of 30 macro
#  steps and a restart with the new cfl, byte for byte, above the floor of kappa (58) and at the floor (58b).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs                     # outputs of a previous run of this test
mr_cases cases mr1d-acu mr3d-sector mr1d-acu-sa mr2d-chan-sa mr2d-chan-sst \
  || { echo "mr-reject: FAIL - cannot build the cases"; exit 1; }
ON=("MOSE-Multirate|enabled|.true." "MOSE-Parameters|iter-threshold|2")
# viscous set-up: a laminar Schmidt number of 0.7 (with one species the mass diffusion flux is zero whatever Sc)
VISC=("MOSE-Physics|equations|navier-stokes" "MOSE-Physics|schmidt|0.7")
status=0

# expect_reject LABEL KEYWORD CASE [--bin BIN] [--bc FACE TYPE DATA...] [-- settings...]: a non-zero exit code and the
# line "[MOSE-Multirate] KEYWORD" in the log; a KEYWORD that starts with "[ERROR]" is a message of the set-up of MOSE,
# which stops the run before the checks of the multirate
expect_reject() {
  local lab=$1 key=$2 case=$3; shift 3
  local bin=$MR_MASTER run=runs/$lab text="[MOSE-Multirate] $key"
  local bc=()
  [[ $key == "[ERROR]"* ]] && text=$key
  if [[ ${1:-} == --bin ]]; then bin=$2; shift 2; fi
  if [[ ${1:-} == --bc ]]; then shift; while [[ $# -gt 0 && $1 != -- ]]; do bc+=("$1"); shift; done; fi
  [[ ${1:-} == -- ]] && shift
  mr_prepare "$run" "cases/$case" "${ON[@]}" "$@" || { echo "  $lab: cannot prepare"; status=1; return; }
  if [[ ${#bc[@]} -gt 0 ]]; then python3 -B ../mr-common/set_bc.py "$run" "${bc[@]}" || { status=1; return; }; fi
  ( cd "$run" && "$bin" > solve.log 2>&1 ); local rc=$?
  if [[ $rc -ne 0 ]] && command grep -a -q -F "$text" "$run/solve.log"; then
    echo "  $lab: rejected (exit $rc): $(command grep -a -F "$text" "$run/solve.log" | head -1 | cut -c1-110)"
  else
    echo "  $lab: NOT rejected as expected (exit $rc)"; tail -3 "$run/solve.log"; status=1
  fi
}

expect_reject 01-rk3        "time-scheme = RK3"            mr1d-acu -- "MOSE-Numerics|time-scheme|RK3"
expect_reject 02-euler      "time-scheme = euler"          mr1d-acu -- "MOSE-Numerics|time-scheme|euler"
expect_reject 03-prim       "integration-variables = prim" mr1d-acu -- "MOSE-Numerics|integration-variables|prim"
expect_reject 04-prec       "integration-variables = prec" mr1d-acu -- "MOSE-Numerics|integration-variables|prec"
expect_reject 05-irs        "IRS"                          mr1d-acu -- "MOSE-Numerics|irs-beta|0.5"
expect_reject 06-shock      "shock detector"               mr1d-acu -- "MOSE-Numerics|riemann-solver|HLLC+Tramel"
expect_reject 06-shock-tr   "shock detector"               mr1d-acu -- "MOSE-Numerics|shock-detector|Tramel"
expect_reject 06-shock-ch   "shock detector"               mr1d-acu -- "MOSE-Numerics|shock-detector|Chen"
expect_reject 07-chemistry  "chemistry"                    mr1d-acu -- "MOSE-Physics|chemistry|finite-rate"
expect_reject 08-soot       "soot"                         mr1d-acu -- "MOSE-Physics|soot-generation|LL91"
expect_reject 09-rotating   "rotating frame"               mr1d-acu -- "MOSE-Physics|rotational-frame|rigid-body"
expect_reject 10-chimera    "chimera"      mr1d-acu --bc 3 102 "1 1" "1 1 1 1 1.0" "1 2 1 1 1.0"
expect_reject 11-coupled    "coupled wall" mr1d-acu --bc 3 103 "1 1 1 1 4 1 0 0 1"
# property line of the grain: mdot_ref, n, p_ref, T_af, mode, mass fractions
expect_reject 13-srm        "SRM grain"    mr1d-acu --bc 3 502 "1.0 0.4 1.0e6 3000.0 1.0 1.0"
expect_reject 14-gsi503     "GSI wall"     mr1d-acu --bc 3 503 "1000.0 300.0 300.0 1.0e5 0.0 0.0 1.0" -- "${VISC[@]}"
expect_reject 15-gsi504     "GSI wall"     mr1d-acu --bc 3 504 "1 0.0 0.0 1.0" -- "${VISC[@]}"
# 505 on the air of the case stops at set-up, before the checks of the multirate: its surface reactions need their eight
# species (Setup_BC); the multirate rejection of the gas-surface interaction walls is the branch of 503 and 504 (14, 15)
expect_reject 16-gsi505     "[ERROR] bc.txt: the surface reactions of the walls 505 and 506 (Bradley) need" mr1d-acu --bc 3 505 "1 0.0" -- "${VISC[@]}"
# 506 (pyro_id, reac_id, phi, h_v0, q_rad, eps_wall) is read now and stops there too, on the same air; the multirate branch
# of the gas-surface interaction walls, 503-506, is that of items 14 and 15
expect_reject 17-gsi506     "[ERROR] bc.txt: the surface reactions of the walls 505 and 506 (Bradley) need" mr1d-acu --bc 3 506 "0 1 0.0 0.0 0.0 0.0" -- "${VISC[@]}"
expect_reject 18-rot3d      "3D rotational periodicity"    mr3d-sector
expect_reject 19-steady     "time-accurate"                mr1d-acu -- "MOSE-Numerics|time-accurate|.false."
expect_reject 20-cflramp    "cfl-rise-threshold"           mr1d-acu -- "MOSE-Numerics|cfl-rise-threshold|10"
expect_reject 21-multigrid  "multigrid"                    mr1d-acu -- "MOSE-Multigrid|levels|2"
expect_reject 22-extfun     "external function"            mr1d-acu --bin "$MR_MMS"
# SAcomp: the single-rate integrator runs it; the stages of the multirate do not evaluate its correction
expect_reject 29-sacomp    "turbulence model: turbulence = SAcomp is not supported" mr1d-acu-sa -- "${VISC[@]}" \
  "MOSE-Physics|turbulence|SAcomp"
# Wilcox2006, whose sources the stages do not evaluate, and the rotation and curvature correction with levels, whose
# terms are a field of the whole domain that the cell-wise sources of a level do not compute: stopped at the first
# macro step
expect_reject 29-wilcox    "turbulence model: turbulence = Wilcox2006 is not supported" mr2d-chan-sst -- \
  "MOSE-Physics|turbulence|Wilcox2006"
expect_reject 29-sarc      "turbulence model: turbulence = SA-RC is not supported with max-level > 0" mr2d-chan-sa -- \
  "MOSE-Physics|turbulence|SA-RC" "MOSE-Multirate|max-level|3"
expect_reject 29-sstrc     "turbulence model: turbulence = SST-RC is not supported with max-level > 0" mr2d-chan-sst -- \
  "MOSE-Physics|turbulence|SST-RC" "MOSE-Multirate|max-level|3"
# restart (newrun = false) with levels = fixed or with fixed-macro-step: rejected on the keys, before
# the restart solution is read
expect_reject 31-restart-fixed "levels fixed"             mr1d-acu -- "MOSE-Multirate|levels|fixed" "MOSE-Parameters|newrun|false"
expect_reject 32-restart-fstep "fixed macro step"         mr1d-acu -- "MOSE-Multirate|fixed-macro-step|2.0" "MOSE-Parameters|newrun|false"
# cfl above 1 with max-level > 0 and nonlocal-cfl = redo (the default): within a macro step the certificate accepts a
# sub-step up to the step that cfl admits, and Heun is linearly stable up to cfl = 1 (van Leer here, first order and
# minmod below: the bound is the same for every limiter)
expect_reject 33-cfl      "cfl"                          mr1d-acu -- "MOSE-Multirate|max-level|3" "MOSE-Numerics|cfl|1.2"
expect_reject 33b-cfl-first "cfl: cfl must be <= 1 when max-level > 0" mr1d-acu -- "MOSE-Multirate|max-level|3" \
  "MOSE-Numerics|space-reconstruction|first-order" "MOSE-Numerics|cfl|1.5"
expect_reject 33c-cfl-minmod "cfl: cfl must be <= 1 when max-level > 0" mr1d-acu -- "MOSE-Multirate|max-level|3" \
  "MOSE-Numerics|flux-limiter|minmod" "MOSE-Numerics|cfl|1.2"
# vnn above 0.5 with max-level > 0 and nonlocal-cfl = redo in a Navier-Stokes run: the same for the diffusion, whose step
# vnn sets (the real stability interval of Heun is [-2, 0]); controls below run vnn 0.5 there, vnn 0.6 in an Euler run,
# where vnn does not enter the step, and both bounds exceeded with nonlocal-cfl = off and count
expect_reject 36-vnn      "vnn: vnn must be <= 0.5 when max-level > 0" mr1d-acu -- "MOSE-Multirate|max-level|3" \
  "${VISC[@]}" "MOSE-Numerics|vnn|0.6"
# levels = fixed with a fixed macro step above the step of the map (2**L dt_min/safety = 4 dt_min, L = 3 the depth of the
# map; with max-level 10 the same, the local steps of mr1d-acu spanning 3.97 octaves): the finest cells run beyond their
# limit in every macro step and the map never rises; the run goes on with the warning of the largest h/(dtcap/sigma),
# once, and its run line at the end, as with levels = recompute where max-level levels are not enough
for spec in "34-fixed-step 3" "34b-fixed-step-depth 10"; do
  read -r lab ml <<< "$spec"
  if mr_run runs/$lab cases/mr1d-acu "${ON[@]}" "MOSE-Multirate|max-level|$ml" "MOSE-Multirate|levels|fixed" \
       "MOSE-Multirate|fixed-macro-step|8.0"; then
    nw=$(command grep -a -c -F '[WARNING] [MOSE-Multirate] the macro step exceeds the limit of some cell' runs/$lab/solve.log)
    nl=$(command grep -a -c -F 'Multirate run: max h/(dtcap/sigma)' runs/$lab/solve.log)
    if [[ $nw -eq 1 && $nl -ge 1 ]]; then
      echo "  $lab: runs (exit 0), the warning of the cells beyond their limit once and its run line"
    else
      echo "  $lab: FAIL - warnings $nw (required 1), run lines $nl (required at least 1)"; status=1
    fi
  else
    echo "  $lab: FAIL - the run with a fixed macro step above the step of the map does not complete"; status=1
  fi
done
# max-attempts below 1: the rule >= 1 of the key registry stops the run when input.ini is read
expect_reject 35-max-attempts "[ERROR] max-attempts must be >= 1" mr1d-acu -- "MOSE-Multirate|max-attempts|0"
# a boundary type without a branch in BC_Flux_Record (MOSE reads any type: an unknown one gets no flux, silently)
expect_reject 28-bctype     "boundary type"                mr1d-acu --bc 3 600

# control: the accepted configuration runs
mr_run runs/00-control cases/mr1d-acu "${ON[@]}" || { echo "  00-control: the accepted configuration does not run"; status=1; }
[[ -f runs/00-control/OUTPUT/field.vtm ]] && echo "  00-control: runs (exit 0)"
# control of the viscous scope: the laminar Navier-Stokes equations run
mr_run runs/00-control-visc cases/mr1d-acu "${ON[@]}" "${VISC[@]}" || { echo "  00-control-visc: the laminar run does not run"; status=1; }
[[ -f runs/00-control-visc/OUTPUT/field.vtm ]] && echo "  00-control-visc: runs (exit 0)"
# control of the RANS scope: the wall channel with Spalart-Allmaras runs
mr_run runs/00-control-rans cases/mr2d-chan-sa "${ON[@]}" || { echo "  00-control-rans: the RANS run does not run"; status=1; }
[[ -f runs/00-control-rans/OUTPUT/field.vtm ]] && echo "  00-control-rans: runs (exit 0)"
# controls of the bound of vnn, which both run: vnn 0.5, the bound itself, in a Navier-Stokes run, and vnn 0.6 with
# max-level 0 (one level, the single-rate step, whose limits the set-up does not check)
mr_run runs/00-control-vnn-edge cases/mr1d-acu "${ON[@]}" "MOSE-Multirate|max-level|3" "${VISC[@]}" "MOSE-Numerics|vnn|0.5" \
  || { echo "  00-control-vnn-edge: the laminar run with vnn 0.5 does not run"; status=1; }
[[ -f runs/00-control-vnn-edge/OUTPUT/field.vtm ]] && echo "  00-control-vnn-edge: runs (exit 0)"
mr_run runs/00-control-vnn-l0 cases/mr1d-acu "${ON[@]}" "MOSE-Multirate|max-level|0" "MOSE-Multirate|safety|1.0" \
  "${VISC[@]}" "MOSE-Numerics|vnn|0.6" \
  || { echo "  00-control-vnn-l0: the laminar run with vnn 0.6 and max-level 0 does not run"; status=1; }
[[ -f runs/00-control-vnn-l0/OUTPUT/field.vtm ]] && echo "  00-control-vnn-l0: runs (exit 0)"
# controls of the mode: vnn 0.6 in an Euler run (vnn does not enter the step), and cfl 1.2 and, in a Navier-Stokes run,
# vnn 0.6 with nonlocal-cfl = off and count, where nothing bounds the sub-steps after the start of the macro step and the
# set-up checks neither bound: these runs complete
for spec in "00-control-vnn-euler redo -- MOSE-Numerics|vnn|0.6" "00-control-cfl-off off -- MOSE-Numerics|cfl|1.2" \
            "00-control-cfl-count count -- MOSE-Numerics|cfl|1.2" "00-control-vnn-off off visc MOSE-Numerics|vnn|0.6" \
            "00-control-vnn-count count visc MOSE-Numerics|vnn|0.6"; do
  read -r lab mode eq kv <<< "$spec"
  keys=("${ON[@]}" "MOSE-Multirate|max-level|3" "MOSE-Multirate|nonlocal-cfl|$mode" "$kv")
  [[ $eq == visc ]] && keys+=("${VISC[@]}")
  if mr_run runs/$lab cases/mr1d-acu "${keys[@]}" \
     && ! command grep -a -q -F '[MOSE-Multirate]' <(command grep -a -F '[ERROR]' runs/$lab/solve.log); then
    echo "  $lab: runs (exit 0)"
  else
    echo "  $lab: FAIL - the run with nonlocal-cfl = $mode and $kv does not complete"; status=1
  fi
done
# controls of the redo with levels = fixed or a fixed macro step, where the certificate only counts and does not redo
# the sub-steps: the set-up checks neither bound there, and cfl 1.2 runs
for spec in "00-control-cfl-redo-fixed MOSE-Multirate|levels|fixed" \
            "00-control-cfl-redo-fstep MOSE-Multirate|fixed-macro-step|2.0"; do
  read -r lab kv <<< "$spec"
  if mr_run runs/$lab cases/mr1d-acu "${ON[@]}" "MOSE-Multirate|max-level|3" "MOSE-Multirate|nonlocal-cfl|redo" \
       "MOSE-Numerics|cfl|1.2" "$kv" \
     && ! command grep -a -q -F '[MOSE-Multirate]' <(command grep -a -F '[ERROR]' runs/$lab/solve.log); then
    echo "  $lab: runs (exit 0)"
  else
    echo "  $lab: FAIL - the run with nonlocal-cfl = redo, $kv and cfl 1.2 does not complete"; status=1
  fi
done

# input.ini read again during the run: the keys of [MOSE-Multirate] are fixed for the run, so a reading that changes one
# is an error of the reading, which ini-error decides. The laminar plate with the multirate
# (max-level 3), 40 macro steps, input.ini read again at iteration 30, numbered solutions every 5; the change is written
# whole (a new file moved in place) once the solution of iteration 5 is written, and must come before the run reads
# input.ini at iteration 30 (else the part says so and fails: a host too slow to poll)
mr_cases cases mr2d-plate || { echo "mr-reject: FAIL - cannot build the case of the reading"; exit 1; }
RD=("MOSE-Parameters|iter-threshold|40" "MOSE-IO|sol-diter|5" "MOSE-IO|sol-overwrite|false" "MOSE-IO|ini-diter|30"
    "MOSE-IO|res-diter|1")
MR3=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0")
# reading RUNDIR 'Section|key|value'... -- CHANGE...: the run with the keys, the keys CHANGE written into its input.ini once
# its first numbered solution exists; sets first, late, rc and last
reading() {
  local run=$1 n pid keys=() chg=()
  shift
  while [[ $# -gt 0 && $1 != -- ]]; do keys+=("$1"); shift; done
  shift
  chg=("$@")
  mr_prepare "$run" "${RCASE:-cases/mr2d-plate}" "${RD[@]}" "${keys[@]}" || { first=no; late=no; rc=1; last=; return; }
  ( cd "$run" && exec "$MR_MASTER" > solve.log 2>&1 ) &
  pid=$!
  for ((n = 0; n < 6000; n++)); do
    [[ -f $run/OUTPUT/field1.vtm ]] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.05
  done
  first=no; late=no
  [[ -f $run/OUTPUT/field1.vtm ]] && first=yes
  cp "$run/input.ini" "$run/input.ini.new"
  python3 -B ../mr-common/set_ini.py "$run/input.ini.new" "${chg[@]}"
  mv "$run/input.ini.new" "$run/input.ini"
  # too late if the run had already written its row of iteration RIT, after which it reads input.ini
  [[ $(awk 'NF {n = $1} END {print n + 0}' "$run/OUTPUT/residual-history.dat" 2>/dev/null) -ge ${RIT:-30} ]] && late=yes
  wait "$pid"
  rc=$?
  last=$(awk 'NF {n = $1} END {print n}' "$run/OUTPUT/residual-history.dat" 2>/dev/null)
}
# not_taken RUNDIR LABEL REFDIR ERROR: last iteration RLAST (40 by default), exit code 0, the warning of the reading of
# iteration RIT (30 by default) once with the error, no [ERROR] in the log, and the residual history of REFDIR, byte for byte
not_taken() {
  local run=$1 lab=$2 ref=$3 err=$4 nw ne same=no
  nw=$(command grep -a -c -F "[WARNING] input.ini read again at iteration ${RIT:-30} is not taken: $err" "$run/solve.log")
  ne=$(command grep -a -c -F '[ERROR]' "$run/solve.log")
  cmp -s "$ref/OUTPUT/residual-history.dat" "$run/OUTPUT/residual-history.dat" && same=yes
  if [[ $first != yes || $late == yes ]]; then
    echo "  $lab: FAIL - the change did not come between its first numbered solution and iteration ${RIT:-30} (host too" \
         "slow to poll)"; status=1
  elif [[ $rc -eq 0 && $last == "${RLAST:-40}" && $nw -eq 1 && $ne -eq 0 && $same == yes ]]; then
    echo "  $lab: the reading of iteration ${RIT:-30} is not taken, the run goes on with the values of its set-up"
  else
    echo "  $lab: FAIL - exit code $rc, last iteration ${last:-none} (required ${RLAST:-40}), warnings $nw (required 1), [ERROR] lines" \
         "$ne, residual history of $ref: $same"; status=1
  fi
}
mr_run runs/50-read-ref cases/mr2d-plate "${RD[@]}" "${MR3[@]}" \
  || { echo "  50-read-ref: the run without a change does not run"; status=1; }
mr_run runs/50-read-ref-sr cases/mr2d-plate "${RD[@]}" \
  || { echo "  50-read-ref-sr: the single-rate run without a change does not run"; status=1; }
# warn (ini-error left out): max-level = 2 and cfl = 0.4 in the same reading; cfl = 0.4 taken would change the history
reading runs/51-read-warn "${MR3[@]}" -- "MOSE-Multirate|max-level|2" "MOSE-Numerics|cfl|0.4"
not_taken runs/51-read-warn "51-read-warn" runs/50-read-ref \
  "input.ini, [MOSE-Multirate] max-level = 2: the run has max-level = 3 from its set-up, and max-level cannot change during the run"
# stop (ini-error = stop from the start): the run stops at the reading of iteration 30
reading runs/52-read-stop "${MR3[@]}" "MOSE-IO|ini-error|stop" -- "MOSE-Multirate|max-level|2"
err="[ERROR] input.ini, [MOSE-Multirate] max-level = 2: the run has max-level = 3 from its set-up, and max-level cannot"
err="$err change during the run"
if [[ $first != yes || $late == yes ]]; then
  echo "  52-read-stop: FAIL - the change did not come between iterations 5 and 30 (host too slow to poll)"; status=1
elif [[ $rc -ne 0 && $last == 30 ]] && command grep -a -q -F "$err" runs/52-read-stop/solve.log \
     && command grep -a -q -E '(ABORT|STOP): ' runs/52-read-stop/solve.log; then
  echo "  52-read-stop: the run stopped at the reading of iteration 30 (exit $rc)"
else
  echo "  52-read-stop: FAIL - exit code $rc, last iteration ${last:-none} (required 30), message present:" \
       "$(command grep -a -c -F "$err" runs/52-read-stop/solve.log)"; status=1
fi
# a single-rate run (the multirate disabled at the set-up) with enabled = .true. written: the reading is not taken
reading runs/53-read-on -- "MOSE-Multirate|enabled|.true."
m53="input.ini, [MOSE-Multirate] enabled = .true.: the run started without enabled in input.ini (its default: .false.),"
m53="$m53 and enabled cannot change during the run"
not_taken runs/53-read-on "53-read-on" runs/50-read-ref-sr "$m53"

# the keys of other sections that the multirate integrator checks at set-up: taken out of the scope, not taken
reading runs/54-read-cfl "${MR3[@]}" -- "MOSE-Numerics|cfl|1.2"
not_taken runs/54-read-cfl "54-read-cfl" runs/50-read-ref "input.ini, [MOSE-Numerics] cfl = 1.2:"
# cfl 1.2 in a run with first-order (theta = 0), where the total-variation bound 1/(1 + theta) alone would admit cfl up
# to safety = 2 and only the linear bound 1 refuses it, and vnn 0.6 in the
# laminar plate
FO=("MOSE-Numerics|space-reconstruction|first-order")
mr_run runs/50-read-ref-fo cases/mr2d-plate "${RD[@]}" "${MR3[@]}" "${FO[@]}" \
  || { echo "  50-read-ref-fo: the first-order run without a change does not run"; status=1; }
reading runs/54b-read-cfl-fo "${MR3[@]}" "${FO[@]}" -- "MOSE-Numerics|cfl|1.2"
not_taken runs/54b-read-cfl-fo "54b-read-cfl-fo" runs/50-read-ref-fo "input.ini, [MOSE-Numerics] cfl = 1.2:"
reading runs/54c-read-vnn "${MR3[@]}" -- "MOSE-Numerics|vnn|0.6"
not_taken runs/54c-read-vnn "54c-read-vnn" runs/50-read-ref "input.ini, [MOSE-Numerics] vnn = 0.6:"
# with nonlocal-cfl = off the readings of cfl 1.2 and vnn 0.6 are taken: the run ends at iteration 40 without the
# warning of a reading that is not taken
for spec in "54d-read-cfl-off MOSE-Numerics|cfl|1.2" "54e-read-vnn-off MOSE-Numerics|vnn|0.6"; do
  read -r lab kv <<< "$spec"
  reading runs/$lab "${MR3[@]}" "MOSE-Multirate|nonlocal-cfl|off" -- "$kv"
  nw=$(command grep -a -c -F 'is not taken' runs/$lab/solve.log)
  if [[ $first != yes || $late == yes ]]; then
    echo "  $lab: FAIL - the change did not come between its first numbered solution and iteration 30 (host too" \
         "slow to poll)"; status=1
  elif [[ $rc -eq 0 && $last == 40 && $nw -eq 0 ]]; then
    echo "  $lab: with nonlocal-cfl = off the reading of iteration 30 is taken"
  else
    echo "  $lab: FAIL - exit code $rc, last iteration ${last:-none} (required 40), warnings $nw (required 0)"; status=1
  fi
done
reading runs/55-read-ta "${MR3[@]}" -- "MOSE-Numerics|time-accurate|.false."
not_taken runs/55-read-ta "55-read-ta" runs/50-read-ref "input.ini, [MOSE-Numerics] time-accurate = .false.:"
reading runs/56-read-ramp "${MR3[@]}" -- "MOSE-Numerics|cfl-rise-threshold|1000"
not_taken runs/56-read-ramp "56-read-ramp" runs/50-read-ref "input.ini, [MOSE-Numerics] cfl-rise-threshold = 1000:"
# the wall channel with SST runs 40 macro steps in a tenth of a second: 1000 macro steps, solutions every 100, the reading
# at iteration 800, so that the change comes in time
RSST=("MOSE-Parameters|iter-threshold|1000" "MOSE-IO|sol-diter|100" "MOSE-IO|ini-diter|800")
mr_run runs/50-read-ref-sst cases/mr2d-chan-sst "${RD[@]}" "${MR3[@]}" "${RSST[@]}" \
  || { echo "  50-read-ref-sst: the SST run without a change does not run"; status=1; }
RCASE=cases/mr2d-chan-sst RIT=800 RLAST=1000
reading runs/57-read-kec "${MR3[@]}" "${RSST[@]}" -- "MOSE-Turbulence|k-energy-coupling|.true."
not_taken runs/57-read-kec "57-read-kec" runs/50-read-ref-sst "input.ini, [MOSE-Turbulence] k-energy-coupling = .true.:"
RCASE=cases/mr2d-plate RIT=30 RLAST=40
# a cfl within the bound, 1, is taken: the run equals 30 macro steps + a restart with that cfl;
# 1.0 is above the floor of kappa and 0.4 at the floor (cfl (1 + theta) = 0.8 < 1.1)
for spec in "58 58-read-cfl-ok 1.0" "58b 58b-read-cfl-low 0.4"; do
  read -r nn lab cv <<< "$spec"
  reading runs/$lab "${MR3[@]}" -- "MOSE-Numerics|cfl|$cv"
  nw=$(command grep -a -c -F 'is not taken' runs/$lab/solve.log)
  split_ok=no
  if mr_prepare runs/$nn-split cases/mr2d-plate "${RD[@]}" "${MR3[@]}" "MOSE-Parameters|iter-threshold|30" \
     && ( cd runs/$nn-split && "$MR_MASTER" > solve.log 2>&1 ) \
     && python3 -B ../mr-common/set_ini.py runs/$nn-split/input.ini "MOSE-Parameters|newrun|false" \
          "MOSE-Parameters|iter-threshold|10" "MOSE-Numerics|cfl|$cv" \
     && ( cd runs/$nn-split && "$MR_MASTER" >> solve.log 2>&1 ); then split_ok=yes; fi
  if [[ $first != yes || $late == yes ]]; then
    echo "  $lab: FAIL - the change did not come between iterations 5 and 30 (host too slow to poll)"; status=1
  elif [[ $rc -eq 0 && $last == 40 && $nw -eq 0 && $split_ok == yes ]] \
       && python3 -B ../../checks/restart/check_restart_vtk.py runs/$lab/OUTPUT runs/$nn-split/OUTPUT > runs/$nn-compare.log 2>&1; then
    echo "  $lab: the reading of cfl = $cv is taken, and the run equals 30 macro steps + a restart with cfl = $cv"
  else
    echo "  $lab: FAIL - exit code $rc, last iteration ${last:-none} (required 40), warnings $nw (required 0)," \
         "split run $split_ok, outputs against 30 + restart with cfl $cv: $(tail -1 runs/$nn-compare.log 2>/dev/null)"; status=1
  fi
done
# ini-error = stop: a cfl above the bound stops the run at the reading
reading runs/59-read-stop-cfl "${MR3[@]}" "MOSE-IO|ini-error|stop" -- "MOSE-Numerics|cfl|1.2"
err="[ERROR] input.ini, [MOSE-Numerics] cfl = 1.2:"
if [[ $first != yes || $late == yes ]]; then
  echo "  59-read-stop-cfl: FAIL - the change did not come between iterations 5 and 30 (host too slow to poll)"; status=1
elif [[ $rc -ne 0 && $last == 30 ]] && command grep -a -q -F "$err" runs/59-read-stop-cfl/solve.log \
     && command grep -a -q -E '(ABORT|STOP): ' runs/59-read-stop-cfl/solve.log; then
  echo "  59-read-stop-cfl: the run stopped at the reading of iteration 30 (exit $rc)"
else
  echo "  59-read-stop-cfl: FAIL - exit code $rc, last iteration ${last:-none} (required 30), message present:" \
       "$(command grep -a -c -F "$err" runs/59-read-stop-cfl/solve.log)"; status=1
fi

if [[ $status -eq 0 ]]; then echo "mr-reject: ok"; else echo "mr-reject: FAIL"; fi
exit $status
