#!/usr/bin/env python3
"""The timers of the multirate in the log of a run (solve.log), on NPROC processes.

NPROC > 1 (a run whose every macro step prints its lines, shell-diter 1, with the certificate of nonlocal-cfl on, as by
default, so that the stage wall time per level of process 0 comes in its end-of-run line):
  - the Timing block reports the exchange wait, the collective wait and the compute spread over ranks;
  - the exchange wait is above zero and at least 100 times the exchange time per level, mean over the processes and
    summed over the levels, over Solver, less 0.1 for the print: the exchange wait of every process contains the time
    of the exchanges of its stages, and the sum of the iteration times over the processes is at most NPROC times Solver,
    their largest;
  - the exchange time per level, mean over the processes, is above zero at some level;
  - the Timing block reports the exchange start, the start of the exchanges of the multirate, a part of the compute: at
    least zero, and with the exchange wait and the collective wait at most 100 %, 0.15 for the print of the three (on
    every process the three times are disjoint parts of the time of its iterations). These checks pass also if one of the
    five exchanges of the multirate were left out of the exchange start: that each is timed, the code shows;
  - the six lines of the times per level over the processes, the stage wall time without the exchanges and the exchange
    time, each its minimum, mean and maximum, have lm + 1 values, those of the end-of-run line of the stage wall time
    per level; at every level minimum <= mean <= maximum; and the values of process 0, its exchange time per level from
    its last line of the macro steps and its stage wall time per level from the end-of-run line less that exchange time,
    lie within [minimum, maximum] of their level, within the precision of the print (4 significant digits of the lines
    of the macro steps, 6 of the end-of-run lines).
NPROC = 1: none of these six lines, and no exchange, exchange start or collective line in the Timing block.
Usage: verify_timing.py SOLVE_LOG NPROC   (exit 0 when every check passes)"""
import re
import sys

NUM = r'[-+]?\d+\.\d*(?:[eE][-+]?\d+)?'
STAGE = 'Multirate run: stage wall time per level without the exchanges [s], '
XCH = 'Multirate run: exchange time per level [s], '
KINDS = ('min over the processes', 'mean over the processes', 'max over the processes')


def values(line, prefix):
    return [float(x) for x in re.findall(NUM, line.split(prefix, 1)[1])]


def last(lines, prefix):
    found = [ln for ln in lines if prefix in ln]
    return values(found[-1], prefix) if found else None


def block_value(lines, name):
    for ln in lines:
        m = re.match(r'^\s+' + re.escape(name) + r'\s+(' + NUM + ')', ln)
        if m:
            return float(m.group(1))
    return None


def main():
    log, nproc = sys.argv[1], int(sys.argv[2])
    lines = open(log, errors='replace').read().splitlines()
    fails = []
    six = {(q, k): last(lines, q + k) for q in (STAGE, XCH) for k in KINDS}
    exch = block_value(lines, 'Exchange wait')
    start = block_value(lines, 'Exchange start')
    coll = block_value(lines, 'Collective wait')
    spread = block_value(lines, 'Compute spread over ranks')
    solver = block_value(lines, 'Solver')
    if solver is None:
        fails.append('no Solver line in the Timing block')
    if nproc == 1:
        for key, v in six.items():
            if v is not None:
                fails.append('one process, but the line "%s%s"' % key)
        if exch is not None or start is not None or coll is not None:
            fails.append('one process, but an exchange, exchange start or collective line in the Timing block')
        print('  timing, 1 process: %s' % ('ok' if not fails else 'FAIL: ' + '; '.join(fails)))
        return 0 if not fails else 1
    stage0 = last(lines, 'Multirate run: stage wall time per level [s], all attempts')
    xch0 = last(lines, 'Multirate: exchange time per level [s]')
    if None in (exch, start, coll, spread, solver):
        fails.append('Timing block without Exchange wait, Exchange start, Collective wait, Compute spread or Solver')
    if stage0 is None or xch0 is None:
        fails.append('no stage wall time per level of process 0 (end-of-run line) or no exchange time per level line')
    for key, v in six.items():
        if v is None:
            fails.append('no line "%s%s"' % key)
    if not fails:
        n = len(stage0)
        for key, v in six.items():
            if len(v) != n:
                fails.append('"%s%s": %d values, %d expected (lm + 1)' % (key[0], key[1], len(v), n))
        if len(xch0) != n:
            fails.append('exchange time per level of process 0: %d values, %d expected' % (len(xch0), n))
    if not fails:
        for q in (STAGE, XCH):
            lo, me, hi = (six[(q, k)] for k in KINDS)
            for lev in range(n):
                tol = 1e-5 * max(abs(hi[lev]), 1e-300)
                if not (lo[lev] <= me[lev] + tol and me[lev] <= hi[lev] + tol):
                    fails.append('%slevel %d: min %g, mean %g, max %g' % (q, lev, lo[lev], me[lev], hi[lev]))
        for lev in range(n):
            x_0, stage_x_0 = xch0[lev], stage0[lev] - xch0[lev]
            lo, hi = six[(XCH, KINDS[0])][lev], six[(XCH, KINDS[2])][lev]
            tol = 5e-4 * abs(x_0) + 5e-6 * abs(hi) + 1e-12
            if not (lo - tol <= x_0 <= hi + tol):
                fails.append('exchange time of process 0 at level %d, %g, outside [%g, %g]' % (lev, x_0, lo, hi))
            lo, hi = six[(STAGE, KINDS[0])][lev], six[(STAGE, KINDS[2])][lev]
            tol = 5e-6 * abs(stage0[lev]) + 5e-4 * abs(x_0) + 5e-6 * abs(hi) + 1e-12
            if not (lo - tol <= stage_x_0 <= hi + tol):
                fails.append('stage time without the exchanges of process 0 at level %d, %g, outside [%g, %g]'
                             % (lev, stage_x_0, lo, hi))
        bound = 100.0 * sum(six[(XCH, KINDS[1])]) / solver
        if not exch >= bound - 0.1:
            fails.append('exchange wait %.1f %% below 100 sum_l mean exchange time / Solver = %.3f %%' % (exch, bound))
        if not exch > 0.0:
            fails.append('exchange wait %.1f %%: the exchanges of the stages were not counted' % exch)
        if not any(v > 0.0 for v in six[(XCH, KINDS[1])]):
            fails.append('exchange time per level, mean over the processes: zero at every level')
        if not start >= 0.0:
            fails.append('exchange start %.1f %% below zero' % start)
        if not start + exch + coll <= 100.0 + 0.15:
            fails.append('exchange start %.1f %%, exchange wait %.1f %% and collective wait %.1f %%: above 100 %% together'
                         % (start, exch, coll))
    if fails:
        print('  timing, %d processes: FAIL: %s' % (nproc, '; '.join(fails)))
        return 1
    print('  timing, %d processes: exchange wait %.1f %% (at least %.3f %%), exchange start %.1f %%, collective wait %.1f %%, '
          'compute spread %.1f %%, %d levels, values of process 0 within the range of the processes'
          % (nproc, exch, bound, start, coll, spread, n))
    return 0


if __name__ == '__main__':
    sys.exit(main())
