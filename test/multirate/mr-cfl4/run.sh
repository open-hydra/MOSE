#!/usr/bin/env bash
#===============================================================================
#  Multirate, condition (4) of Gnedin, Semenov and Kravtsov (2018) imposed by the certificate and the redo
#  (nonlocal-cfl, interface-margin, force-redo, max-attempts; docs/theory/time-integration.md, Non-local CFL condition). Levels
#  recomputed at every macro step, max-level 3, safety 2, residual at every macro step:
#   acu3b     mr1d-acu-3b, 100 macro steps. The cell that fixes H runs exactly at its limit, so without the interface
#             margin its coarser neighbour has no margin and the certificate fails (count, interface-margin off: at least
#             one attempt with a violation; positive control); with the margin no violation (count, interface-margin
#             on); with redo no macro step redone, and the outputs identical byte by byte to those of count with the
#             margin (the redo path does not change a run that it does not redo); with force-redo (the first attempt of
#             every macro step rejected and redone with the same H and map) every macro step redone once and the outputs
#             identical byte by byte again (the state goes back to t_n to the bit).
#   uniform   mr1d-uni-3b (the mesh of mr1d-acu-3b, a uniform flow), cfl 0.5, 100 macro steps, count and redo: the local steps
#             change only by the rounding, so the cell that fixes H starts with hhat/dcap = 1 within some units in the last
#             place: no violation (the limit up to 1e-12 of it), no macro step redone.
#   floor     mr1d-acu-3b, 100 macro steps, redo, where kappa = max(cfl (1 + theta), 1.1) is at the floor: van Leer with cfl 0.5
#             (max-level 10 and 3) and 0.25, minmod with cfl 0.6, first order with cfl 0.8; and above it, van Leer with cfl
#             0.6. The acoustic wave makes the local steps fall by about 4e-6 within a macro step, far less than 1 - 1/1.1:
#             no macro step redone (docs/theory/time-integration.md, Appendix A, Proposition), the kappa printed at the
#             first macro step is the expected one, and at least one cell raised by the margin. Without the floor these
#             runs have all their 100 macro steps redone, with H halved at max-level 3.
#   blast100  mr1d-blast100 up to 0.4 ms: the front of a strong blast wave (shock Mach 5.9) crosses about 2.5 cells per
#             macro step, beyond the face neighbours covered by the margin. count with the margin: at least one attempt
#             with a violation (positive control); redo: at least one macro step redone, the attempts with a violation
#             are exactly the rejected ones (attempts - macro steps: every accepted attempt satisfies the certificate),
#             mass and energy drift within 1e-12, 1 and 3 OpenMP threads identical byte by byte; force-redo: outputs
#             identical byte by byte to those of redo, with violations + forced rejections = attempts - macro steps.
#             With max-level 10 (the depth of the map from the octaves of the local steps, 3 at the start): at least one
#             macro step whose map a redo made deeper than that of its first attempt (the run line of the redo), the same
#             balance of the attempts, mass and energy drift within 1e-12.
#             With solutions every sol-dtime = 4e-6 s (below the natural macro step, so nearly every macro step is cut to
#             land on an output instant) and redo: at least one landed macro step whose H a redo halved (it no
#             longer lands, a later one reaches the instant); the 100 numbered solutions at k*sol-dtime and the final one
#             at time-threshold bit for bit; one residual row per macro step. With time-threshold = 10 sol-dtime the last
#             macro step lands on an output instant and on the end of the run: the solution due there is written at the
#             last iteration with its number, so the 10 numbered solutions at k*sol-dtime and the final one at
#             time-threshold bit for bit (10 is the first k for which (k-1) sol-dtime + sol-dtime rounds above
#             k sol-dtime in double precision, with sol-dtime = 4e-6 s: the time of the last output plus sol-dtime is
#             then above the time of the run, and a comparison of times would miss that solution).
#   sod       mr1d-sod up to 1e-3 s, redo: the same balance of the attempts; at most 2 attempts per macro step.
#   blast-wide mr1d-blast-wide up to 3 ms with max-level and nonlocal-cfl left out (10 levels, redo): a macro step of the
#             front takes a state that is not physical inside an attempt (a cell or a face), so the attempt is rejected
#             and redone with H halved on the same level map. At least one attempt rejected for it (the run line of the
#             redo), the balance of the attempts with these rejections added, mass and energy drift within 1e-12, 1 and
#             3 OpenMP threads identical byte by byte; and in every macro step whose rejected attempts were all rejected
#             for a state that is not physical, the halvings of H of the macro step equal those rejections (one halving
#             each: the line of the macro step, shell-diter 1).
#   sod-hot   mr1d-sod-hot (a cell outside the tables of the gas at the start), redo: the faces of that cell are not
#             physical at the first stage of the first macro step, which does not depend on H, so the run stops there as
#             the single-rate step does: exit code 1, the lines of the face that is not physical, no macro step done and
#             no line of a rejected attempt.
#   blast-hot mr1d-blast-hot (the hot gas at 9990 K, within the tables up to 10000 K) up to 0.4 ms, redo: at the contact the
#             right state of the face of the last hot cell, reconstructed at the second stage of a sub-step, leaves the
#             tables, so the attempts are rejected with a face that is not physical and redone with H halved; the first two
#             macro steps are accepted that way, and in the third every attempt takes it, so after 11 attempts the macro step
#             falls below 2**-10 of its first attempt and the run stops: exit code not 0, at least one line of a macro step
#             with an attempt rejected with a face, the line of the stop with 11 attempts, the causes 0 / 0 / 11, and the
#             first state that is not physical of the last attempt, a face, with its lines.
#   max-attempts  the attempts of a macro step rejected for a cause, at most. mr1d-blast-wide as in blast-wide, whose
#             first macro step has 4 attempts (the first three rejected), with max-attempts 1, 2 and 3: each run stops in
#             that macro step with exit code 1 and the line of the limit, the causes of the line that follows add up to the
#             limit and do not decrease from one limit to the next (the same attempts), the line after them is that of the
#             cause of the last rejection, and with 3 the causes are 1 / 0 / 2 (two rejections for a state that is not
#             physical, which halve H, and one for a violation of condition (4)); no line of the floor and no macro step
#             done. With max-attempts 4 the fourth attempt is accepted: the outputs and the run line of the redo are those
#             of blast-wide byte by byte. mr1d-blast-hot as in blast-hot with max-attempts 11: the limit is checked after
#             the floor, so the macro step accepted at its 11th attempt and the stop below 2**-10 after 11 attempts give
#             the lines of blast-hot, character by character; with max-attempts 10 the second macro step stops at its 10th
#             rejection with the line of the limit, the causes 0 / 0 / 10 and the state that is not physical of the last
#             attempt, after one macro step done. mr1d-acu-3b as acu3b with count and max-attempts 3: the warning that the key
#             acts only with redo, once, and the outputs of the run without the key byte by byte; the runs without the key
#             print no such warning. mr1d-acu-3b as acu3b with force-redo and max-attempts 1: the attempt rejected by
#             force-redo is not counted, so every macro step is accepted at its second attempt, and the outputs and the run
#             line of force-redo are those of acu3b with force-redo byte by byte.
#   local steps  mr1d-sod-pneg (a cell at a negative pressure: its local step is not a number) and mr1d-sod-p0 (a cell at
#             zero pressure, at rest: its local step is infinite), redo; mr1d-sod-nneg (a cell at a negative density and
#             pressure, about 300 K: with the Navier-Stokes equations its diffusivity mu/rho is negative and so is its local
#             step, finite, which only the comparison with 0 finds) and mr1d-sod-vac (cell 100, in the finest band, at a
#             density of 1e-307 kg/m3, about 300 K: with the Navier-Stokes equations its diffusion rate overflows and its
#             local step is zero), redo: the run stops at the start of the first macro step, before its map and its first
#             stage: exit code 1, the line of the local steps with one cell, the line of the cell (block 1, cell 10 or 100)
#             with its local step and that of its primitive variables, the line of the stop, no macro step done, no line
#             of a rejected attempt and no line of a face. mr1d-sod-pneg, mr1d-sod-p0 and mr1d-sod-vac reach their local
#             steps through an invalid operation (the speed of sound of a negative temperature), a division by zero and an
#             overflow: a build that traps these floating-point exceptions (DEBUG) stops there with SIGFPE, before the
#             check, so these three parts hold in the builds that do not trap them (TESTING, RELEASE).
#   theta     mr1d-acu-3b, one macro step with each limiter (minmod, vanalbada, vanleer, mc, superbee) and with first-order:
#             the theta of the limiter in the interface margin kappa = max(cfl (1 + theta), 1.1), reported at the first macro step,
#             is that of the limiter of input.ini (1/2, sqrt(2)/2, 1, 1, 1, 0), the value that Initialise takes at the first
#             macro step (Theta_Limiter).
#   In every run with the certificate, its kernel (Dtcap_Cell, a copy of compute of Mod_dt.f90) gives dtlocal back at t_n
#   (no mismatch above 1e-12).
#   The run lines of the certificate (attempts with a violation, mismatches), of the redo (macro steps redone, macro steps
#   cut to land whose H a redo halved, attempts rejected by force-redo) come only when their count is not zero: a missing
#   line counts zero, and a run whose count is zero must print none (count with the margin and redo on acu3b); the lines
#   of the certificate in the macro steps come only with a violation or a macro step redone (none in redo on acu3b).
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-acu-3b mr1d-uni-3b mr1d-blast100 mr1d-sod mr1d-blast-wide mr1d-sod-hot mr1d-blast-hot mr1d-sod-pneg mr1d-sod-p0 \
  mr1d-sod-nneg mr1d-sod-vac || { echo "mr-cfl4: FAIL - cannot build the cases"; exit 1; }
status=0
BASE=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" "MOSE-Multirate|safety|2.0"
      "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" "MOSE-IO|res-diter|1")
FORCE=("MOSE-Multirate|force-redo|.true.")

# field RUN TEXT -> the integer after TEXT in the run lines of the certificate or of the redo
field() { command grep -a 'Multirate run: \(certificate\|redo\)' "$1/solve.log" | sed -n "s/.*$2 \([0-9][0-9]*\).*/\1/p" | head -1; }
# forced RUN -> the attempts rejected by force-redo
forced() { command grep -a 'Multirate run: redo, attempts rejected by force-redo' "$1/solve.log" | sed -n 's/.*(test device) \([0-9]*\).*/\1/p' | head -1; }

# balance RUN [FORCED] -> attempts with a violation (+ forced rejections + rejections for a state that is not physical) =
# attempts - macro steps
balance() {
  local line nm na nr ns nf=${2:-0}
  line=$(command grep -a 'Multirate run: redo, macro steps' "$1/solve.log")
  nr=$(field "$1" 'attempts with a violation')
  ns=$(field "$1" 'attempts rejected for a state that is not physical')
  if [[ -z $line ]]; then   # no macro step redone: no rejected attempt, so no attempt with a violation and no forced one
    echo "    no macro step redone (no run line of the redo)"
    [[ -z $nr && -z $ns && $nf -eq 0 ]] \
      || { echo "    attempts with a violation ${nr:-0} + forced $nf + not physical ${ns:-0}, but no macro step redone"; return 1; }
    return 0
  fi
  nm=$(sed -n 's/.*redo, macro steps \([0-9]*\),.*/\1/p' <<< "$line"); na=$(sed -n 's/.*, attempts \([0-9]*\),.*/\1/p' <<< "$line")
  echo "    ${line#*Multirate run: }"
  echo "    $(command grep -a 'Multirate run: redo, macro steps with' "$1/solve.log" | sed 's/.*Multirate run: //')"
  [[ -n $nm && -n $na && $((${nr:-0} + nf + ${ns:-0})) -eq $((na - nm)) ]] \
    || { echo "    attempts with a violation ${nr:-0} + forced $nf + not physical ${ns:-0}, attempts - macro steps $((na - nm)):" \
              "not equal"; return 1; }
}

echo "  acu3b:"
N=("MOSE-Parameters|iter-threshold|100")
mr_run runs/acu3b-count cases/mr1d-acu-3b "${BASE[@]}" "${N[@]}" "MOSE-Multirate|nonlocal-cfl|count" || status=1
mr_run runs/acu3b-count-K cases/mr1d-acu-3b "${BASE[@]}" "${N[@]}" "MOSE-Multirate|nonlocal-cfl|count" \
  "MOSE-Multirate|interface-margin|.true." || status=1
mr_run runs/acu3b-redo cases/mr1d-acu-3b "${BASE[@]}" "${N[@]}" "MOSE-Multirate|nonlocal-cfl|redo" || status=1
mr_run runs/acu3b-force cases/mr1d-acu-3b "${BASE[@]}" "${N[@]}" "MOSE-Multirate|nonlocal-cfl|redo" "${FORCE[@]}" || status=1
v0=$(field runs/acu3b-count 'attempts with a violation'); v1=$(field runs/acu3b-count-K 'attempts with a violation')
rd=$(command grep -a 'Multirate run: redo, macro steps' runs/acu3b-redo/solve.log | sed -n 's/.*macro steps redone \([0-9]*\),.*/\1/p' | head -1)
nf=$(forced runs/acu3b-force)
echo "    count without the margin: attempts with a violation ${v0:-0 (no line)}; with the margin: ${v1:-0 (no line)}; redo: macro steps redone ${rd:-0 (no line)}; force-redo: forced ${nf:-0 (no line)}"
[[ -n $v0 && $v0 -ge 1 ]] || { echo "    no violation without the margin (positive control failed)"; status=1; }
[[ -z $v1 ]] || { echo "    a run line of the certificate with the margin ($v1 attempts with a violation): none expected"; status=1; }
[[ -z $rd ]] || { echo "    a run line of the redo ($rd macro steps redone): none expected"; status=1; }
# the lines of the certificate in the macro steps come only with a violation or a macro step redone: none here
nc=$(command grep -a -c -E 'Multirate: certificate, violations per level|Multirate: attempts of the macro step' runs/acu3b-redo/solve.log)
[[ $nc -eq 0 ]] || { echo "    redo: $nc lines of the certificate in the macro steps, none expected"; status=1; }
[[ -n $nf && $nf -eq 100 ]] || { echo "    force-redo: forced rejections ${nf:-?}, 100 expected"; status=1; }
if mr_same runs/acu3b-count-K runs/acu3b-redo; then echo "    redo and count with the margin: identical"; else echo "    redo and count with the margin: DIFFERENT"; status=1; fi
if mr_same runs/acu3b-redo runs/acu3b-force; then echo "    force-redo and redo: identical"; else echo "    force-redo and redo: DIFFERENT"; status=1; fi

echo "  uniform flow, a tie with the limit at the rounding:"
U=("${BASE[@]}" "${N[@]}" "MOSE-Numerics|cfl|0.5")
mr_run runs/uni-count cases/mr1d-uni-3b "${U[@]}" "MOSE-Multirate|nonlocal-cfl|count" || status=1
mr_run runs/uni-redo cases/mr1d-uni-3b "${U[@]}" "MOSE-Multirate|nonlocal-cfl|redo" || status=1
v0=$(field runs/uni-count 'attempts with a violation'); v1=$(field runs/uni-redo 'attempts with a violation')
echo "    attempts with a violation, count: ${v0:-0 (no line)}; redo: ${v1:-0 (no line)}"
[[ -z $v0 && -z $v1 ]] || { echo "    a violation on a uniform flow: none expected"; status=1; }

echo "  the floor of the interface margin:"
for spec in "vanleer 0.5 10 1.10000" "vanleer 0.5 3 1.10000" "vanleer 0.25 3 1.10000" "minmod 0.6 3 1.10000" \
            "first-order 0.8 3 1.10000" "vanleer 0.6 3 1.20000"; do
  read -r lim cf ml ka <<< "$spec"
  if [[ $lim == first-order ]]; then K="MOSE-Numerics|space-reconstruction|first-order"; else K="MOSE-Numerics|flux-limiter|$lim"; fi
  r=runs/floor-$lim-$cf-$ml
  mr_run "$r" cases/mr1d-acu-3b "${BASE[@]}" "${N[@]}" "MOSE-Multirate|nonlocal-cfl|redo" "MOSE-Numerics|cfl|$cf" \
    "MOSE-Multirate|max-level|$ml" "$K" || { status=1; continue; }
  got=$(command grep -a 'Multirate: theta of the limiter' "$r/solve.log" | sed -n 's/.*kappa = max(cfl (1 + theta), 1.1) = *\([0-9.]*\),.*/\1/p' | head -1)
  rd=$(command grep -a 'Multirate run: redo, macro steps' "$r/solve.log" | sed -n 's/.*macro steps redone \([0-9]*\),.*/\1/p' | head -1)
  ring=$(command grep -a 'Multirate run: cells raised by kappa' "$r/solve.log" | sed 's/.*mean over the maps//' | awk '{s = 0; for (i = 1; i <= NF; i++) s += $i; print (s > 0) ? "yes" : "no"}')
  echo "    $lim, cfl $cf, max-level $ml: kappa ${got:-?}, macro steps redone ${rd:-0 (no line)}, cells raised by the margin ${ring:-no line}"
  [[ $got == "$ka" ]] || { echo "    kappa $ka expected"; status=1; }
  [[ -z $rd ]] || { echo "    macro steps redone: none expected"; status=1; }
  [[ $ring == yes ]] || { echo "    no cell raised by the margin: the margin does not act"; status=1; }
done

echo "  blast100:"
T=("MOSE-Parameters|time-threshold|4.0e-4" "MOSE-Parameters|iter-threshold|1000000")
mr_run runs/blast100-count-K cases/mr1d-blast100 "${BASE[@]}" "${T[@]}" "MOSE-Multirate|nonlocal-cfl|count" \
  "MOSE-Multirate|interface-margin|.true." || status=1
vk=$(field runs/blast100-count-K 'attempts with a violation')
echo "    count with the margin: attempts with a violation ${vk:-?}"
[[ -n $vk && $vk -ge 1 ]] || { echo "    no violation with the margin alone (positive control failed)"; status=1; }
R=("${BASE[@]}" "${T[@]}" "MOSE-Multirate|nonlocal-cfl|redo")
mr_run runs/blast100-redo cases/mr1d-blast100 "${R[@]}" || status=1
balance runs/blast100-redo || status=1
rd=$(command grep -a 'Multirate run: redo, macro steps' runs/blast100-redo/solve.log | sed -n 's/.*macro steps redone \([0-9]*\),.*/\1/p' | head -1)
[[ -n $rd && $rd -ge 1 ]] || { echo "    no macro step redone: the redo path is not exercised"; status=1; }
python3 -B ../mr-common/verify_conservation.py cases/mr1d-blast100 runs/blast100-redo --quantities mass,energy --max 1e-12 || status=1
OMP_NUM_THREADS=3 mr_run runs/blast100-redo-t3 cases/mr1d-blast100 "${R[@]}" || status=1
if mr_same runs/blast100-redo runs/blast100-redo-t3; then echo "    1 and 3 threads: identical"; else echo "    1 and 3 threads: DIFFERENT"; status=1; fi
mr_run runs/blast100-force cases/mr1d-blast100 "${R[@]}" "${FORCE[@]}" || status=1
nf=$(forced runs/blast100-force)
[[ -n $nf && $nf -ge 1 ]] || { echo "    force-redo: no forced rejection"; status=1; }
balance runs/blast100-force "${nf:-0}" || status=1
if mr_same runs/blast100-redo runs/blast100-force; then echo "    force-redo ($nf forced) and redo: identical"; else echo "    force-redo and redo: DIFFERENT"; status=1; fi

# the redo below the depth of the map: max-level 10, where the depth of the map built at each macro step follows the
# octaves of the local steps (3 at the start); a rejected attempt whose cells need a finer level is redone with a deeper
# map (up to max-level) instead of a halved H
echo "  blast100 with max-level 10, the redo below the depth of the map:"
if mr_run runs/blast100-deep cases/mr1d-blast100 "${R[@]}" "MOSE-Multirate|max-level|10"; then
  balance runs/blast100-deep || status=1
  nd=$(command grep -a 'Multirate run: redo, macro steps [0-9]' runs/blast100-deep/solve.log \
       | sed -n 's/.*maps made deeper by a redo \([0-9]*\).*/\1/p' | head -1)
  echo "    maps made deeper by a redo: ${nd:-none (no field)}"
  [[ -n $nd && $nd -ge 1 ]] || { echo "    no map made deeper by a redo: the redo below the depth of the map is not exercised"; status=1; }
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-blast100 runs/blast100-deep --quantities mass,energy --max 1e-12 \
    || status=1
else
  status=1
fi

# landing on the output instants together with the redo: solutions every DTO, below the natural macro step, so that
# nearly every macro step is cut to land; a macro step whose H a redo halves does not land, and a later one reaches the instant
echo "  blast100, output instants with redo:"
DTO=4.0e-6
TT=$(python3 -c "print(repr(100 * $DTO + $DTO / 3))")
mr_run runs/blast100-land cases/mr1d-blast100 "${BASE[@]}" "MOSE-Parameters|time-threshold|$TT" \
  "MOSE-Parameters|iter-threshold|1000000" "MOSE-Multirate|nonlocal-cfl|redo" "MOSE-IO|sol-dtime|$DTO" \
  "MOSE-IO|sol-overwrite|false" || status=1
nu=$(command grep -a 'Multirate run: redo, macro steps cut to land' runs/blast100-land/solve.log | sed -n 's/.* halved \([0-9]*\).*/\1/p' | head -1)
echo "    macro steps cut to land whose H a redo halved: ${nu:-?}"
[[ -n $nu && $nu -ge 1 ]] || { echo "    no landed macro step halved: the halving of a landed macro step by a redo is not exercised"; status=1; }
python3 -B - runs/blast100-land/OUTPUT "$DTO" "$TT" <<'EOF' || status=1
import sys
from pathlib import Path
sys.path.insert(0, '../../fast/restart/restart-counters')
from vtk_time import time_of
out, dto, tt = Path(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
num = sorted(out.glob('field[0-9]*.vtm'), key=lambda p: int(p.stem[5:]))
bad = [p.name for p in num if time_of(p) != float(int(p.stem[5:])) * dto]
tf = time_of(out / 'field.vtm')
rows = [l.split() for l in (out / 'residual-history.dat').read_text().splitlines() if l.strip()]
ok = len(num) == 100 and not bad and tf == tt and len(rows) == int(rows[-1][0])
print('    %d numbered solutions, TIME = k*sol-dtime bit for bit: %s; final TIME = time-threshold bit for bit: %s; residual '
      'rows %d, iteration of the last row %s' % (len(num), 'yes' if not bad else 'NO ' + ' '.join(bad[:5]),
      'yes' if tf == tt else 'NO (%.17g)' % tf, len(rows), rows[-1][0]))
sys.exit(0 if ok else 1)
EOF

# an output instant at the end of the run: the solution due at the last iteration, written with its number
echo "  blast100, an output instant at the end of the run:"
TE=$(python3 -c "print(repr(10 * $DTO))")
mr_run runs/blast100-end cases/mr1d-blast100 "${BASE[@]}" "MOSE-Parameters|time-threshold|$TE" \
  "MOSE-Parameters|iter-threshold|1000000" "MOSE-Multirate|nonlocal-cfl|redo" "MOSE-IO|sol-dtime|$DTO" \
  "MOSE-IO|sol-overwrite|false" || status=1
python3 -B - runs/blast100-end/OUTPUT "$DTO" "$TE" <<'EOF' || status=1
import sys
from pathlib import Path
sys.path.insert(0, '../../fast/restart/restart-counters')
from vtk_time import time_of
out, dto, te = Path(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
num = sorted(out.glob('field[0-9]*.vtm'), key=lambda p: int(p.stem[5:]))
ks = [int(p.stem[5:]) for p in num]
bad = [p.name for p in num if time_of(p) != float(int(p.stem[5:])) * dto]
tf = time_of(out / 'field.vtm') if (out / 'field.vtm').exists() else None
ok = ks == list(range(1, 11)) and not bad and tf == te
print('    numbered solutions %s; TIME = k*sol-dtime bit for bit: %s; final TIME = time-threshold bit for bit: %s'
      % (' '.join(map(str, ks)) or 'none', 'yes' if not bad else 'NO ' + ' '.join(bad[:5]),
         'yes' if tf == te else 'NO (%r)' % tf))
sys.exit(0 if ok else 1)
EOF

echo "  sod:"
mr_run runs/sod-redo cases/mr1d-sod "${BASE[@]}" "MOSE-Parameters|time-threshold|1.0e-3" "MOSE-Parameters|iter-threshold|1000000" \
  "MOSE-Multirate|nonlocal-cfl|redo" || status=1
balance runs/sod-redo || status=1
amax=$(command grep -a 'Multirate run: redo, macro steps ' runs/sod-redo/solve.log | sed -n 's/.*, largest \([0-9]*\),.*/\1/p' | head -1)
amax=${amax:-1}   # no run line of the redo: no macro step redone, one attempt each
[[ $amax -le 2 ]] || { echo "    largest number of attempts $amax > 2"; status=1; }

# a state that is not physical inside an attempt: mr1d-blast-wide with max-level and nonlocal-cfl left out (10 levels, redo)
echo "  blast-wide, a state that is not physical inside an attempt:"
W=("MOSE-Multirate|enabled|.true." "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|recompute"
   "MOSE-Multirate|recompute-every|1" "MOSE-IO|res-diter|1" "MOSE-IO|shell-diter|1"
   "MOSE-Parameters|time-threshold|3.0e-3" "MOSE-Parameters|iter-threshold|1000000")
if mr_run runs/blast-wide cases/mr1d-blast-wide "${W[@]}"; then
  balance runs/blast-wide || status=1
  ns=$(field runs/blast-wide 'attempts rejected for a state that is not physical')
  echo "    $(command grep -a 'Multirate run: redo, attempts rejected for a state' runs/blast-wide/solve.log | sed 's/.*Multirate run: //')"
  [[ -n $ns && $ns -ge 1 ]] || { echo "    no attempt rejected for a state that is not physical: the path is not exercised"; status=1; }
  # every macro step whose rejected attempts were all rejected for a state that is not physical: one halving of H for each
  python3 -B - runs/blast-wide/solve.log <<'EOF' || status=1
import re
import sys
att, n_all, bad = None, 0, []
for ln in open(sys.argv[1], errors='replace'):
    m = re.search(r'Multirate: attempts of the macro step (\d+),', ln)
    if m:
        att = int(m.group(1))
    m = re.search(r'Multirate: attempts of the macro step rejected for a state that is not physical (\d+) .*halvings of H in '
                  r'the macro step (\d+)', ln)
    if m:
        ns, nh = int(m.group(1)), int(m.group(2))
        if att is not None and att - 1 == ns:
            n_all += 1
            if nh != ns:
                bad.append((att, ns, nh))
        att = None
print('    macro steps whose rejected attempts were all for a state that is not physical: %d; with halvings of H other than '
      'those rejections: %d %s' % (n_all, len(bad), bad[:3] if bad else ''))
sys.exit(0 if n_all >= 1 and not bad else 1)
EOF
  python3 -B ../mr-common/verify_conservation.py cases/mr1d-blast-wide runs/blast-wide --quantities mass,energy --max 1e-12 \
    || status=1
  OMP_NUM_THREADS=3 mr_run runs/blast-wide-t3 cases/mr1d-blast-wide "${W[@]}" || status=1
  if mr_same runs/blast-wide runs/blast-wide-t3; then echo "    1 and 3 threads: identical"; else echo "    1 and 3 threads: DIFFERENT"; status=1; fi
else
  status=1
fi

# a face that is not physical at the first stage of the macro step: it does not depend on H, and the run stops there as the
# single-rate step does
echo "  sod-hot, a face that is not physical at the start of the first macro step:"
mr_prepare runs/sod-hot cases/mr1d-sod-hot "${BASE[@]}" "MOSE-Parameters|iter-threshold|10" "MOSE-Multirate|nonlocal-cfl|redo" \
  || status=1
( cd runs/sod-hot && "$MR_MASTER" > solve.log 2>&1 )
rc=$?
SL=runs/sod-hot/solve.log
nh=$(command grep -a -c 'Multirate: H =' "$SL")
nb=$(command grep -a -c -E 'fell below 2\*\*-10|rejected for a state that is not physical' "$SL")
echo "    exit $rc; $(command grep -a 'face between cell' "$SL" | head -1 | sed 's/^ *//'); $(command grep -a 'error flags L, R' "$SL" \
     | head -1 | sed 's/^ *//')"
if [[ $rc -eq 1 && $nh -eq 0 && $nb -eq 0 ]] \
   && command grep -a -q -F '[ERROR] non-physical state detected at the interface.' "$SL" \
   && command grep -a -q 'Out-of-bounds temperature' "$SL" && command grep -a -q 'face between cell' "$SL" \
   && command grep -a -q -E 'error flags L, R: (2 [0-9]|[0-9] 2)$' "$SL"; then
  echo "    stopped at the first stage with the lines of the face"
else
  echo "    exit $rc (1 expected), macro steps done $nh (0 expected), lines of a rejected attempt or of the stop below 2**-10 $nb" \
       "(0 expected), or the lines of the face missing"
  status=1
fi

# a face that is not physical after the first stage: the attempts are rejected and redone; in the third macro step every
# attempt takes one, and the run stops below 2**-10 of its first attempt with the causes and the face
echo "  blast-hot, faces that are not physical after the first stage:"
mr_prepare runs/blast-hot cases/mr1d-blast-hot "${BASE[@]}" "MOSE-Parameters|time-threshold|4.0e-4" \
  "MOSE-Parameters|iter-threshold|1000000" "MOSE-Multirate|nonlocal-cfl|redo" "MOSE-IO|shell-diter|1" || status=1
( cd runs/blast-hot && "$MR_MASTER" > solve.log 2>&1 )
rc=$?
HL=runs/blast-hot/solve.log
nf=$(command grep -a 'Multirate: attempts of the macro step rejected for a state' "$HL" | sed -n 's/.*(with a face \([0-9]*\)).*/\1/p' \
     | awk '$1 >= 1' | wc -l)
command grep -a -E 'Multirate: attempts of the macro step rejected|\[ERROR\]' "$HL" | sed 's/^ */    /' | head -12
if [[ $rc -ne 0 && $nf -ge 1 ]] && command grep -a -q -E 'fell below 2\*\*-10 of its first attempt \(H = .* s after 11 attempts\)' "$HL" \
   && command grep -a -q -F '0 for a violation of condition (4), 0 for a local step that is not a positive number, 11 for a state that is not physical' "$HL" \
   && command grep -a -q 'the last attempt took a state that is not physical' "$HL" \
   && command grep -a -q -E '\[ERROR\]   stage [0-9]+ of the attempt' "$HL" \
   && command grep -a -q -F '[ERROR] non-physical state detected at the interface.' "$HL" && command grep -a -q 'face between cell' "$HL"; then
  echo "    $nf macro steps with an attempt rejected with a face; stopped below 2**-10 after 11 attempts with the face"
else
  echo "    exit $rc (not 0 expected), macro steps with an attempt rejected with a face $nf (at least 1 expected), or the lines" \
       "of the stop missing"
  status=1
fi

# max-attempts: the attempts of a macro step rejected for a cause, at most; the run stops when they are all rejected
echo "  max-attempts:"
rc_hot=$rc   # exit code of blast-hot, above
# causes RUN -> "violations local-steps states" of the line of the causes of the stop
causes() {
  command grep -a '\[ERROR\]   attempts of the macro step rejected:' "$1/solve.log" | head -1 \
    | sed -n 's/.*rejected: \([0-9]*\) for a violation of condition (4), \([0-9]*\) for a local step .*, \([0-9]*\) for a state.*/\1 \2 \3/p'
}
# the first macro step of blast-wide has 4 attempts: with max-attempts 1, 2, 3 the run stops in it, at the k-th rejection
pv=0; pn=0; pe=0
for k in 1 2 3; do
  run=runs/blast-wide-cap$k
  mr_prepare "$run" cases/mr1d-blast-wide "${W[@]}" "MOSE-Multirate|max-attempts|$k" || { status=1; continue; }
  ( cd "$run" && "$MR_MASTER" > solve.log 2>&1 )
  rc=$?
  L=$run/solve.log
  read -r cv cn ce <<< "$(causes "$run")"
  nh=$(command grep -a -c 'Multirate: H =' "$L")
  nf=$(command grep -a -c 'fell below 2\*\*-10' "$L")
  ns=$(command grep -a -c -F "STOP: MOSE-Multirate: macro step rejected in all its attempts, max-attempts = $k (nonlocal-cfl = redo)" "$L")
  # the line after the causes: that of the cause of the k-th rejection, the one that grew from k-1 to k
  last=bad
  if [[ -n ${ce:-} ]]; then
    if [[ $ce -gt $pe ]] && command grep -a -q 'the last attempt took a state that is not physical' "$L"; then last=ok; fi
    if [[ $cv -gt $pv ]] && command grep -a -q -E 'the last attempt violated condition \(4\): violations [1-9][0-9]*, cells with a violation [1-9]' "$L"; then last=ok; fi
    if [[ $cn -gt $pn ]] && command grep -a -q 'the last attempt gave a local step that is not a positive number' "$L"; then last=ok; fi
  fi
  echo "    max-attempts $k: exit $rc, causes ${cv:-?} / ${cn:-?} / ${ce:-?}, line of the last rejection $last, macro steps done $nh"
  if [[ $rc -eq 1 && -n ${ce:-} && $((cv + cn + ce)) -eq $k && $cv -ge $pv && $cn -ge $pn && $ce -ge $pe && $last == ok && $nh -eq 0 \
        && $nf -eq 0 && $ns -eq 1 ]] \
     && command grep -a -q -F "the macro step was rejected in all its attempts (max-attempts = $k;" "$L"; then
    pv=$cv; pn=$cn; pe=$ce
  else
    echo "    max-attempts $k: exit 1, the line of the limit, causes adding up to $k and not below those of $((k - 1)), the line of" \
         "the last cause, no macro step done, no line of the floor and one line of the stop expected"
    status=1
  fi
done
[[ $pv -eq 1 && $pn -eq 0 && $pe -eq 2 ]] || { echo "    max-attempts 3: causes 1 / 0 / 2 expected"; status=1; }
# the fourth attempt of that macro step is accepted: max-attempts 4 changes nothing
if mr_run runs/blast-wide-cap4 cases/mr1d-blast-wide "${W[@]}" "MOSE-Multirate|max-attempts|4"; then
  a=$(command grep -a 'Multirate run: redo, macro steps [0-9]' runs/blast-wide/solve.log)
  b=$(command grep -a 'Multirate run: redo, macro steps [0-9]' runs/blast-wide-cap4/solve.log)
  if [[ -n $a && $a == "$b" ]] && mr_same runs/blast-wide runs/blast-wide-cap4; then
    echo "    max-attempts 4: the outputs and the run line of the redo of blast-wide"
  else
    echo "    max-attempts 4: the outputs or the run line of the redo differ from blast-wide"; status=1
  fi
else
  status=1
fi
# the limit after the floor: blast-hot with max-attempts 11 gives the lines of blast-hot; with 10 the second macro step,
# accepted at its 11th attempt without the limit, stops at its 10th rejection
HK=("${BASE[@]}" "MOSE-Parameters|time-threshold|4.0e-4" "MOSE-Parameters|iter-threshold|1000000" "MOSE-Multirate|nonlocal-cfl|redo"
    "MOSE-IO|shell-diter|1")
sel='^ \[ERROR\]|Multirate: attempts of the macro step|STOP: '
mr_prepare runs/blast-hot-cap11 cases/mr1d-blast-hot "${HK[@]}" "MOSE-Multirate|max-attempts|11" || status=1
( cd runs/blast-hot-cap11 && "$MR_MASTER" > solve.log 2>&1 )
rc11=$?
a=$(command grep -a -E "$sel" runs/blast-hot/solve.log); b=$(command grep -a -E "$sel" runs/blast-hot-cap11/solve.log)
if [[ $rc11 -eq $rc_hot && -n $a && $a == "$b" ]]; then
  echo "    blast-hot with max-attempts 11: exit $rc11 and the $(wc -l <<< "$a") lines of the attempts and of the stop of blast-hot"
else
  echo "    blast-hot with max-attempts 11: exit $rc11 ($rc_hot expected) or lines of the attempts and of the stop not those of blast-hot"
  status=1
fi
mr_prepare runs/blast-hot-cap10 cases/mr1d-blast-hot "${HK[@]}" "MOSE-Multirate|max-attempts|10" || status=1
( cd runs/blast-hot-cap10 && "$MR_MASTER" > solve.log 2>&1 )
rc10=$?
L=runs/blast-hot-cap10/solve.log
nh=$(command grep -a -c 'Multirate: H =' "$L")
if [[ $rc10 -eq 1 && $nh -eq 1 && $(causes runs/blast-hot-cap10) == "0 0 10" ]] \
   && command grep -a -q -F 'the macro step was rejected in all its attempts (max-attempts = 10;' "$L" \
   && command grep -a -q 'the last attempt took a state that is not physical' "$L" && ! command grep -a -q 'fell below 2\*\*-10' "$L"; then
  echo "    blast-hot with max-attempts 10: the second macro step stopped at its 10th rejection (causes 0 / 0 / 10)"
else
  echo "    blast-hot with max-attempts 10: exit $rc10 (1 expected), macro steps done $nh (1 expected), causes $(causes runs/blast-hot-cap10)" \
       "(0 0 10 expected), or the lines of the limit missing"
  status=1
fi
# the key with the redo off: a warning, and nothing else changes
WK='[WARNING] [MOSE-Multirate] max-attempts acts only with nonlocal-cfl = redo'
mr_run runs/acu3b-count-cap cases/mr1d-acu-3b "${BASE[@]}" "${N[@]}" "MOSE-Multirate|nonlocal-cfl|count" \
  "MOSE-Multirate|max-attempts|3" || status=1
nw=$(command grep -a -c -F "$WK" runs/acu3b-count-cap/solve.log)
n0=$(cat runs/acu3b-count/solve.log runs/acu3b-redo/solve.log | command grep -a -c -F "$WK")
if [[ $nw -eq 1 && $n0 -eq 0 ]] && mr_same runs/acu3b-count runs/acu3b-count-cap; then
  echo "    acu3b count with max-attempts: the warning once, the outputs of the run without the key"
else
  echo "    acu3b count with max-attempts: warnings $nw (1 expected), warnings in the runs without the key $n0 (0 expected), or" \
       "outputs not those of the run without the key"
  status=1
fi
# the attempt rejected by force-redo is not counted: with max-attempts 1 every macro step is accepted at its second attempt
FL='Multirate run: redo, attempts rejected by force-redo'
if mr_run runs/acu3b-force-cap1 cases/mr1d-acu-3b "${BASE[@]}" "${N[@]}" "MOSE-Multirate|nonlocal-cfl|redo" "${FORCE[@]}" \
     "MOSE-Multirate|max-attempts|1"; then
  a=$(command grep -a "$FL" runs/acu3b-force/solve.log); b=$(command grep -a "$FL" runs/acu3b-force-cap1/solve.log)
  nf=$(forced runs/acu3b-force-cap1)
  if [[ -n $a && $a == "$b" && ${nf:-0} -eq 100 ]] && mr_same runs/acu3b-force runs/acu3b-force-cap1; then
    echo "    acu3b force-redo with max-attempts 1: forced $nf, the outputs and the run line of force-redo of the run without the key"
  else
    echo "    acu3b force-redo with max-attempts 1: forced ${nf:-?} (100 expected), or the outputs or the run line of force-redo" \
         "not those of the run without the key"
    status=1
  fi
else
  status=1
fi

# a local step that is not a positive finite number at the start of the macro step stops the run there
echo "  local steps that are not a positive finite number:"
NS=("MOSE-Physics|equations|navier-stokes" "MOSE-Numerics|vnn|0.5")
# label, case, cell, equations, the local step printed (an extended regular expression)
for spec in "sod-pneg mr1d-sod-pneg 10 euler NaN" "sod-p0 mr1d-sod-p0 10 euler Infinity" \
            "sod-nneg mr1d-sod-nneg 10 ns -[0-9]\.[0-9]{7}E[-+][0-9]{3}" "sod-vac mr1d-sod-vac 100 ns 0\.0000000E\+000"; do
  read -r lab c cell eq val <<< "$spec"
  E=(); [[ $eq == ns ]] && E=("${NS[@]}")
  mr_prepare runs/$lab cases/$c "${BASE[@]}" "MOSE-Parameters|iter-threshold|10" "MOSE-Multirate|nonlocal-cfl|redo" "${E[@]}" \
    || { status=1; continue; }
  ( cd runs/$lab && "$MR_MASTER" > solve.log 2>&1 )
  rc=$?
  L=runs/$lab/solve.log
  nh=$(command grep -a -c 'Multirate: H =' "$L")
  nb=$(command grep -a -c -E 'fell below 2\*\*-10|rejected in all its attempts|rejected for a state that is not physical|face between cell' "$L")
  nl=$(command grep -a -c -E '\[ERROR\] \[MOSE-Multirate\] at the start of the macro step \(time +0\.0+E\+00 s\) the local step is not a positive finite number in 1 cell;' "$L")
  ns=$(command grep -a -c 'STOP: MOSE-Multirate: a local step that is not a positive finite number at the start of the macro step' "$L")
  if [[ $rc -eq 1 && $nh -eq 0 && $nb -eq 0 && $nl -eq 1 && $ns -eq 1 ]] \
     && command grep -a -q -E "\[ERROR\]   cell b, i, j, k: 1 $cell 1 1: local step +\+?$val s\$" "$L" \
     && command grep -a -q -F '[ERROR]   primitive variables:' "$L"; then
    echo "    $lab: stopped at the start of the first macro step, cell $cell with the local step" \
         "$(command grep -a 'cell b, i, j, k:' "$L" | sed 's/.*local step *//')"
  else
    echo "    $lab: exit $rc (1 expected), macro steps done $nh (0 expected), lines of a rejected attempt or of a face $nb (0)," \
         "lines of the local steps $nl (1), lines of the stop $ns (1), or the line of cell $cell with the local step $val missing"
    status=1
  fi
done

# theta of the limiter in the interface margin at the first macro step, after the set-up, against the value that the set-up
# check of cfl takes from the keys of input.ini before it
echo "  theta of the limiter at the first macro step:"
for lt in minmod:0.50000 vanalbada:0.70711 vanleer:1.00000 mc:1.00000 superbee:1.00000 first-order:0.00000; do
  lim=${lt%%:*}; th=${lt#*:}
  if [[ $lim == first-order ]]; then K="MOSE-Numerics|space-reconstruction|first-order"; else K="MOSE-Numerics|flux-limiter|$lim"; fi
  mr_run runs/theta-$lim cases/mr1d-acu-3b "${BASE[@]}" "MOSE-Parameters|iter-threshold|1" "MOSE-Multirate|nonlocal-cfl|count" \
    "$K" || status=1
  got=$(command grep -a 'Multirate: theta of the limiter' runs/theta-$lim/solve.log | sed -n 's/.*theta of the limiter *\([0-9.]*\),.*/\1/p' | head -1)
  echo "    $lim: theta ${got:-?} (expected $th)"
  [[ $got == "$th" ]] || status=1
done

# the kernel of the certificate (Dtcap_Cell, a copy of compute of Mod_dt.f90) gives dtlocal back at t_n in every run
echo "  certificate against dtlocal at t_n:"
for r in uni-count uni-redo floor-vanleer-0.5-10 floor-vanleer-0.25-3 acu3b-count acu3b-count-K acu3b-redo acu3b-force blast100-count-K blast100-redo blast100-force blast100-deep \
         blast100-land blast100-end sod-redo blast-wide blast-wide-t3 blast-wide-cap4 acu3b-count-cap acu3b-force-cap1; do
  mm=$(command grep -a 'Multirate run: certificate (' runs/$r/solve.log | sed -n 's/.*mismatches with dtlocal at t_n: \([0-9]*\).*/\1/p' | head -1)
  [[ ${mm:-0} -eq 0 ]] || { echo "    $r: mismatches $mm"; status=1; }   # no run line of the certificate: no mismatch
done

if [[ $status -eq 0 ]]; then echo "mr-cfl4: ok"; else echo "mr-cfl4: FAIL"; fi
exit $status
