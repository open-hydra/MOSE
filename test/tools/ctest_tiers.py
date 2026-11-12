#!/usr/bin/env python3
"""Keep the tiers of test/tiers.txt honest.

A test's tier follows from its time on the CI runner (one thread):
    light   under 5 s      medium   5 s to 100 s      heavy   100 s and more

    ctest ... --output-junit junit.xml
    ctest_tiers.py report junit.xml [--markdown]     tests whose time puts them in another tier, and failures
    ctest_tiers.py update junit.xml [--scale S]      rewrite tiers.txt from the times (times * S, S = 1 on the CI runner)

`report` only warns (exit code 0) unless --strict: a time varies by tens of percent from run to run, so a test is
reported only when it is outside its tier by more than a factor of TOL (1.5); `update` moves a test only for the
same reason, and fills the tests that have no time yet strictly.  Standard library only.
"""
import argparse
import os
import sys
import xml.etree.ElementTree as ET

LIGHT_MAX = 5.0
MEDIUM_MAX = 100.0
TOL = 1.5
HERE = os.path.dirname(os.path.abspath(__file__))
TIERS = os.path.join(HERE, "..", "tiers.txt")
ORDER = ("light", "medium", "heavy")


def tier_of(seconds):
    return "light" if seconds < LIGHT_MAX else ("medium" if seconds < MEDIUM_MAX else "heavy")


def inside(tier, seconds):
    """True if `seconds` is within the tier, with the tolerance band at both ends."""
    lo = {"light": 0.0, "medium": LIGHT_MAX / TOL, "heavy": MEDIUM_MAX / TOL}[tier]
    hi = {"light": LIGHT_MAX * TOL, "medium": MEDIUM_MAX * TOL, "heavy": float("inf")}[tier]
    return lo <= seconds <= hi


def read_tiers(path=TIERS):
    """(header lines, {name: (tier, seconds or None)})"""
    header, rows = [], {}
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            if not rows:
                header.append(line.rstrip("\n"))
            continue
        name, tier, sec = s.split()
        rows[name] = (tier, None if sec == "-" else float(sec))
    return header, rows


def write_tiers(header, rows, path=TIERS):
    out = list(header) + [""]
    for tier in ORDER:
        for name in sorted(n for n, (t, _) in rows.items() if t == tier):
            sec = rows[name][1]
            out.append("%-28s %-7s %s" % (name, tier, "-" if sec is None else ("%.1f" % sec if sec < 10 else "%.0f" % sec)))
    open(path, "w").write("\n".join(out) + "\n")


def read_junit(path):
    """{name: (seconds, passed)} from the JUnit file of ctest (--output-junit)."""
    res = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        failed = any(c.tag in ("failure", "error") for c in case)
        skipped = case.get("status") == "notrun" or any(c.tag == "skipped" for c in case)
        if not skipped:
            res[case.get("name")] = (float(case.get("time", "0")), not failed)
    return res


def cmd_report(a):
    _, rows = read_tiers()
    times = read_junit(a.junit)
    bad = [(n, rows[n][0], t) for n, (t, ok) in sorted(times.items()) if n in rows and not inside(rows[n][0], t * a.scale)]
    failed = sorted(n for n, (t, ok) in times.items() if not ok)
    unknown = sorted(n for n in times if n not in rows)
    total = {k: sum(t * a.scale for n, (t, _) in times.items() if n in rows and rows[n][0] == k) for k in ORDER}
    if a.markdown:
        print("### Tests by tier\n")
        print("| tier | tests | time |\n|---|--:|--:|")
        for k in ORDER:
            n = sum(1 for x in times if x in rows and rows[x][0] == k)
            print("| %s | %d | %.0f s |" % (k, n, total[k]))
        if failed:
            print("\n**Failed:** " + ", ".join("`%s`" % n for n in failed))
        print("\n### Slowest tests\n\n| test | tier | s |\n|---|---|--:|")
        for n, (t, _) in sorted(times.items(), key=lambda x: -x[1][0])[:15]:
            print("| `%s` | %s | %.1f |" % (n, rows.get(n, ("?",))[0], t * a.scale))
        if bad:
            print("\n### Outside their tier (by more than a factor of %.1f)\n\n| test | tier | s |\n|---|---|--:|" % TOL)
            for n, tier, t in bad:
                print("| `%s` | %s | %.1f |" % (n, tier, t * a.scale))
    else:
        for n, tier, t in bad:
            print("%-28s listed %-6s, takes %7.1f s -> %s" % (n, tier, t * a.scale, tier_of(t * a.scale)))
        for n in unknown:
            print("%-28s has no line in tiers.txt" % n)
        print("%d tests, %d outside their tier, %d failed" % (len(times), len(bad), len(failed)))
    return 1 if a.strict and (bad or unknown or failed) else 0


def cmd_update(a):
    header, rows = read_tiers()
    times = read_junit(a.junit)
    moved = []
    for n, (t, ok) in sorted(times.items()):
        if n not in rows or not ok:
            continue
        sec = t * a.scale
        old_tier, old_sec = rows[n]
        new_tier = old_tier if (old_sec is not None and inside(old_tier, sec)) else tier_of(sec)
        if new_tier != old_tier:
            moved.append((n, old_tier, new_tier, sec))
        rows[n] = (new_tier, sec)
    write_tiers(header, rows)
    for n, o, k, s in moved:
        print("%-28s %-6s -> %-6s (%.1f s)" % (n, o, k, s))
    print("tiers.txt: %d times written, %d tests moved" % (sum(1 for n in times if n in rows), len(moved)))
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd")
    sub.required = True
    for name, fn in (("report", cmd_report), ("update", cmd_update)):
        s = sub.add_parser(name)
        s.add_argument("junit")
        s.add_argument("--scale", type=float, default=1.0, help="factor from the measured time to the time on the CI runner")
        if name == "report":
            s.add_argument("--markdown", action="store_true")
            s.add_argument("--strict", action="store_true")
        s.set_defaults(fn=fn)
    a = p.parse_args()
    sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
