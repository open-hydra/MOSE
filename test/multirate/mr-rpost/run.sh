#!/usr/bin/env bash
#===============================================================================
#  Multirate, step taken against the step admitted on the new state: r_post =
#  h_c(n)/dtcap_c(t_n+1), and r_nb = h_c(n)/dtcap_v(t_n+1) for the face neighbours v of c in its block (condition (4)
#  of Gnedin, Semenov and Kravtsov 2018, section 2; docs/theory/time-integration.md, Non-local CFL condition), with
#  the level map and the macro step of the last macro step and the local step that Compute_dt gives at the start of
#  the next one. Sod tube with the levels recomputed at every macro step, up to t >= 1e-3 s:
#   max-level 0, safety 1  the single-rate step on the multirate path, H = the smallest local step at the start of the
#             macro step: when the wave speeds grow within a step the cell that sets it exceeds the step admitted on
#             the state it reached, r_post > 1 (+ 1e-12) in at least one macro step (positive control). With one level
#             every cell has the step H and every cell of the block is a face neighbour of another, so the largest r_nb
#             is the largest r_post (to the rounding of the division: relative 1e-14 here) and the macro steps with a
#             cell above 1 are the same for r_post, r_nb and both together. With max-level 3 and safety 1 the Sod
#             tube breaks down under nonlocal-cfl = off (r_post up to 2.16) and completes under redo: this test keeps
#             safety 1 for the one-level run only;
#   max-level 3, safety 8  the cells run at one eighth of their limit: r_post and r_nb <= 1 in every macro step
#             (negative control).
#  Both runs end with exit 0. The run lines "Multirate run: max r_post" and "Multirate run: max r_nb" (15 digits) come only
#  when a macro step had a cell above 1: max-level 0 prints them, safety 8 prints none (the line of the largest r_post and
#  r_nb of the macro steps, every shell-diter = 100 iterations, shows that it evaluated both; the line of their counts comes
#  only when a cell is above 1 or above 1/cfl, so safety 8 prints none). The run line "Multirate run: max h/(dtcap/sigma)" comes only when a macro
#  step started with a cell above its limit: neither run prints it; fixed2, the acoustic wave with one level and a fixed
#  macro step of twice the step of the first macro step (max-level 0, safety 1), prints it with every macro step counted
#  (positive control).
#  Positive control with several levels: mr1d-blast (hot gas at 64 times the temperature of the cold gas, 5 times its
#  pressure), safety 2, max-level 3, up to 6e-4 s. The blast front enters cells of coarser levels, whose admissible step
#  drops by up to 2.4 behind the shock, while on two cells one level apart (3) holds at t_n with a margin between 1 and 2
#  only: r_nb > 1 in at least one macro step, and in at least one macro step a cell with r_nb > 1 >= r_post (more cells
#  with r_post or r_nb above 1 than with r_post above 1, from the lines of the counts of the macro steps). Both are partial measures:
#  the state at the end of the macro step only, and for r_nb the face neighbours in the block only.
#  Every run sets nonlocal-cfl = off: the test checks the measures without the certificate, and with the default (redo) the
#  interface margin raises the coarser cells ahead of the front, which removes the r_nb > 1 that the positive control needs.
#===============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$HERE"
source ../mr-common/lib.sh
rm -rf cases runs
mr_cases cases mr1d-sod mr1d-blast mr1d-acu || { echo "mr-rpost: FAIL - cannot build the cases"; exit 1; }
status=0
# runs/sigma<s>: the number of macro steps with a cell above 1 (and above 1/cfl, informative) from the run line
for s in 1 8; do
  L=3; [[ $s == 1 ]] && L=0
  mr_run "runs/sigma$s" cases/mr1d-sod "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|$L" \
    "MOSE-Multirate|safety|$s.0" "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" \
    "MOSE-Parameters|time-threshold|1.0e-3" "MOSE-Parameters|iter-threshold|1000000" "MOSE-IO|res-diter|1" \
    "MOSE-IO|shell-diter|100" "MOSE-Multirate|nonlocal-cfl|off" || { status=1; continue; }
  line=$(command grep -a 'Multirate run: max r_post' "runs/sigma$s/solve.log")
  lnb=$(command grep -a 'Multirate run: max r_nb' "runs/sigma$s/solve.log")
  if command grep -a -q 'Multirate run: max h/(dtcap/sigma)' "runs/sigma$s/solve.log"; then
    echo "  safety $s: a run line of h/(dtcap/sigma), but no macro step starts with a cell above its limit"; status=1
  fi
  if [[ $s == 8 ]]; then
    ne=$(command grep -a -c 'Multirate: r_post = ' "runs/sigma$s/solve.log")
    nc=$(command grep -a -c -E 'Multirate: .*(cells with r_post > 1|cells > 1: )' "runs/sigma$s/solve.log")
    if [[ -z $line && -z $lnb && $ne -ge 1 && $nc -eq 0 ]]; then
      echo "  safety 8: no run line of r_post or r_nb (no macro step above 1); $ne lines of the macro steps with the" \
           "largest of both, none with their counts"
    else
      echo "  safety 8: run lines of r_post or r_nb, no line of the macro steps with the largest of both, or $nc lines" \
           "with their counts (negative control failed)"
      [[ -n $line ]] && echo "    ${line#*Multirate run: }"
      [[ -n $lnb ]] && echo "    ${lnb#*Multirate run: }"
      status=1
    fi
    continue
  fi
  [[ -n $line && -n $lnb ]] || { echo "  safety $s: no run line of r_post or r_nb (positive control failed)"; status=1; continue; }
  n1=$(sed -E 's/.*a cell above 1: ([0-9]+);.*/\1/' <<< "$line")
  rp=$(sed -E 's/.*h\/dtcap\(t_n\+1\) +([-+.0-9E]+) \(.*/\1/' <<< "$line")
  m1=$(sed -E 's/.*a cell above 1: ([0-9]+);.*/\1/' <<< "$lnb")
  ma=$(sed -E 's/.*with r_post or r_nb above 1: ([0-9]+);.*/\1/' <<< "$lnb")
  rn=$(sed -E 's/.*face neighbour +([-+.0-9E]+) \(.*/\1/' <<< "$lnb")
  echo "  safety $s: ${line#*Multirate run: }"
  echo "  safety $s: ${lnb#*Multirate run: }"
  if [[ $s == 1 ]]; then
    [[ $n1 -ge 1 ]] || { echo "  max-level 0, safety 1: r_post never above 1 (positive control failed)"; status=1; }
    [[ $m1 -eq $n1 && $ma -eq $n1 ]] \
      || { echo "  max-level 0, safety 1: macro steps above 1: r_post $n1, r_nb $m1, either $ma (should be equal)"; status=1; }
    awk -v a="$rn" -v b="$rp" 'BEGIN { d = a - b; if (d < 0) d = -d; exit !(b > 0 && d <= 1e-14 * b) }' \
      || { echo "  max-level 0, safety 1: largest r_nb $rn differs from largest r_post $rp"; status=1; }
  fi
done
# runs/fixed2: positive control of the run line of h/(dtcap/sigma): one level and a fixed macro step of twice the step of
# the first macro step, so that every macro step starts with the smallest cells at about twice their limit (3 macro steps;
# the warning that the fixed macro step exceeds the limit of some cell comes too)
if mr_run runs/fixed2 cases/mr1d-acu "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|0" \
    "MOSE-Multirate|safety|1.0" "MOSE-Multirate|fixed-macro-step|2.0" "MOSE-Multirate|nonlocal-cfl|off" \
    "MOSE-Parameters|iter-threshold|3" "MOSE-IO|res-diter|1"; then
  lr=$(command grep -a 'Multirate run: max h/(dtcap/sigma)' runs/fixed2/solve.log)
  nr=$(sed -n -E 's/.*macro steps above 1 \+ 1e-12: ([0-9]+) of ([0-9]+)\).*/\1 \2/p' <<< "$lr")
  if [[ $nr == "3 3" ]] && command grep -a -q 'exceeds the limit of some cell' runs/fixed2/solve.log; then
    echo "  fixed2: ${lr#*Multirate run: }"
  else
    echo "  fixed2: no run line of h/(dtcap/sigma) with the 3 macro steps above the limit (positive control failed): ${lr:-none}"
    status=1
  fi
else
  status=1
fi
# runs/blast: positive control with several levels (mr1d-blast, safety 2, max-level 3, up to 6e-4 s, lines of every macro
# step): r_nb above 1 in at least one macro step, and in at least one macro step more cells with r_post or r_nb above 1 than
# with r_post above 1, i.e. a cell with r_nb > 1 >= r_post (the coarser cell ahead of the front, whose own step is admitted)
if mr_run runs/blast cases/mr1d-blast "MOSE-Multirate|enabled|.true." "MOSE-Multirate|max-level|3" \
    "MOSE-Multirate|safety|2.0" "MOSE-Multirate|levels|recompute" "MOSE-Multirate|recompute-every|1" \
    "MOSE-Parameters|time-threshold|6.0e-4" "MOSE-Parameters|iter-threshold|1000000" "MOSE-IO|res-diter|1" \
    "MOSE-IO|shell-diter|1" "MOSE-Multirate|nonlocal-cfl|off"; then
  lnb=$(command grep -a 'Multirate run: max r_nb' runs/blast/solve.log)
  if [[ -z $lnb ]]; then
    echo "  blast: no run line of r_nb"; status=1
  else
    lp=$(command grep -a 'Multirate run: max r_post' runs/blast/solve.log)
    echo "  blast: ${lp:+${lp#*Multirate run: }}${lp:-no run line of r_post (no macro step with r_post above 1)}"
    echo "  blast: ${lnb#*Multirate run: }"
    m1=$(sed -E 's/.*a cell above 1: ([0-9]+);.*/\1/' <<< "$lnb")
    # the counts of a macro step come only when a cell is above 1 or above 1/cfl (a missing line counts zero)
    nx=$(awk '/Multirate: cells with r_post > 1: / { split($0, a, "cells with r_post > 1: "); split(a[2], b, ","); np = b[1] + 0
                                                    split($0, c, "with r_post or r_nb > 1: "); if (c[2] + 0 > np) n++ }
              END { print n + 0 }' runs/blast/solve.log)
    echo "  blast: macro steps with a cell where r_nb > 1 >= r_post: $nx"
    [[ $m1 -ge 1 ]] || { echo "  blast: r_nb never above 1 (positive control failed)"; status=1; }
    [[ $nx -ge 1 ]] || { echo "  blast: no cell with r_nb > 1 >= r_post (positive control failed)"; status=1; }
  fi
else
  status=1
fi
if [[ $status -eq 0 ]]; then echo "mr-rpost: ok"; else echo "mr-rpost: FAIL"; fi
exit $status
