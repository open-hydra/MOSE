"""Checks of the TimeControl fast test (exit code 0 = pass, 1 = fail).

  python3 -B check.py order MIN SCHEME [SCHEME ...]
      OUTPUT/<scheme>_<k>.tec are runs with dt0/2^k over the same time, and
      OUTPUT/<scheme>_ref.tec a reference with a much smaller step. The order
      observed between the last two runs, for p and u in the max norm, must
      be >= MIN.
  python3 -B check.py same TAG_A TAG_B TOL
      OUTPUT/<TAG_A>.tec and OUTPUT/<TAG_B>.tec agree within TOL (relative,
      max norm, rho, u, p).
  python3 -B check.py reached TAG TIME
      the solution time of OUTPUT/<TAG>.tec is >= TIME.
  python3 -B check.py threshold EPS SPLIT
      OUTPUT/residual-history-{bar,tiny,split-restart}.dat are the steady runs
      of part 4 with res-threshold = EPS: at 1 bar, at 1 bar x 2^-24, and at
      1 bar in two parts with the restart after SPLIT iterations. The run at
      low pressure ends after its first iteration with its last density
      residual at most EPS times the first; the runs at the two pressures
      end at the same iteration, after SPLIT; the run in two parts ends at
      that iteration too, with the same ratio of the last density residual
      to the first.
"""
import math
import sys

from make_case import read_field


def max_diff(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def order(min_order, schemes):
    ok_all = True
    for scheme in schemes:
        ref = read_field(f"OUTPUT/{scheme}_ref.tec")
        errors, k = [], 0
        while True:
            try:
                q = read_field(f"OUTPUT/{scheme}_{k}.tec")
            except FileNotFoundError:
                break
            errors.append((max_diff(q["p"], ref["p"]), max_diff(q["u"], ref["u"])))
            k += 1
        print(f"  {scheme}: max-norm distance from the reference, steps dt0/2^k")
        for k, (ep, eu) in enumerate(errors):
            line = f"    k = {k}:  p {ep:.3e} Pa   u {eu:.3e} m/s"
            if k > 0 and ep > 0 and eu > 0:
                line += (f"   order p {math.log2(errors[k - 1][0] / ep):.2f}"
                         f"  u {math.log2(errors[k - 1][1] / eu):.2f}")
            print(line)
        if len(errors) < 2 or min(errors[-1]) <= 0:
            print(f"  [FAIL] {scheme}: runs missing")
            ok_all = False
            continue
        op = math.log2(errors[-2][0] / errors[-1][0])
        ou = math.log2(errors[-2][1] / errors[-1][1])
        ok = op >= min_order and ou >= min_order
        print(f"  [{'ok  ' if ok else 'FAIL'}] {scheme}: order of the last pair p {op:.2f}, "
              f"u {ou:.2f} (required >= {min_order})")
        ok_all = ok_all and ok
    return ok_all


def same(tag_a, tag_b, tol):
    try:
        a, b = read_field(f"OUTPUT/{tag_a}.tec"), read_field(f"OUTPUT/{tag_b}.tec")
    except FileNotFoundError as err:
        print(f"  [FAIL] {tag_b} vs {tag_a}: missing field ({err.filename})")
        return False
    worst = max(max_diff(a[v], b[v]) / max(abs(x) for x in a[v] + [1e-300])
                for v in ("rho(1)", "u", "p"))
    ok = worst <= tol
    print(f"  [{'ok  ' if ok else 'FAIL'}] {tag_b} vs {tag_a}: max relative difference "
          f"{worst:.2e} (required <= {tol:.0e})")
    return ok


def reached(tag, time):
    try:
        t = read_field(f"OUTPUT/{tag}.tec")["time"]
    except FileNotFoundError:
        print(f"  [FAIL] {tag}: no field")
        return False
    ok = t >= time
    print(f"  [{'ok  ' if ok else 'FAIL'}] {tag}: final time {t:.6e} s (required >= {time:.6e} s)")
    return ok


def history(path):
    """Iterations and density residuals of a residual-history.dat (columns: iteration, time, density, ...)."""
    rows = []
    with open(path) as f:
        for line in f:
            v = line.split()
            if len(v) >= 3:
                rows.append((int(v[0]), float(v[2].replace("D", "E"))))
    return rows


def threshold(eps, n_split):
    runs = {}
    for tag in ("bar", "tiny", "split-restart"):
        try:
            runs[tag] = history(f"OUTPUT/residual-history-{tag}.dat")
        except FileNotFoundError:
            runs[tag] = []
        if not runs[tag]:
            print(f"  [FAIL] {tag}: no residual history")
            return False
        h = runs[tag]
        print(f"    {tag}: last iteration {h[-1][0]}, first density residual {h[0][1]:.6e}, "
              f"last {h[-1][1]:.6e}, last/first {h[-1][1] / h[0][1]:.3e}")
    hb, ht, hs = runs["bar"], runs["tiny"], runs["split-restart"]
    checks = [
        (ht[-1][0] > 1 and ht[-1][1] <= eps * ht[0][1],
         f"the run at 1 bar x 2^-24 ends once its residual has fallen by {eps:.0e} from the first iteration"),
        (hb[-1][0] == ht[-1][0], "the runs at 1 bar and at 1 bar x 2^-24 end at the same iteration"),
        (hb[-1][0] > n_split, f"the run at 1 bar ends after iteration {n_split}, where the run in two parts restarts"),
        (hs[-1][0] == hb[-1][0] and hs[-1][1] / hs[0][1] == hb[-1][1] / hb[0][1],
         "the run in two parts ends where the run in one go ends, with the same last/first residual"),
    ]
    for ok, text in checks:
        print(f"  [{'ok  ' if ok else 'FAIL'}] {text}")
    return all(ok for ok, _ in checks)


def main():
    what = sys.argv[1]
    if what == "order":
        ok = order(float(sys.argv[2]), sys.argv[3:])
    elif what == "same":
        ok = same(sys.argv[2], sys.argv[3], float(sys.argv[4]))
    elif what == "reached":
        ok = reached(sys.argv[2], float(sys.argv[3]))
    elif what == "threshold":
        ok = threshold(float(sys.argv[2]), int(sys.argv[3]))
    else:
        sys.exit(f"unknown check {what}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
