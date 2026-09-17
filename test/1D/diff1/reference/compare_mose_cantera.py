"""Validate MOSE's diffusion against an independent Cantera unity-Lewis solution.

Run MOSE with  Sc = 1.0,  Prl = 1.0  (unity Lewis: D = mu/rho, Pr = 1) so that MOSE
and the Cantera reference use the *same* transport closure -- any difference is then
purely numerical (discretization, grid, the compressible convection MOSE carries).

    python compare_mose_cantera.py            # uses OUTPUT/field.tec
    python compare_mose_cantera.py path.tec   # explicit snapshot

Produces mose_vs_cantera.png and prints L2 / Linf errors for CH4 and T.
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import reference_unity_lewis as ref

I_CH4, I_T, N_SPECIES = 0, 8, 4


def _setup_orion_import():
    root = Path(__file__).resolve().parent
    for parent in [root, *root.parents]:
        cand = parent / "lib" / "ORION" / "src" / "python"
        if cand.exists():
            sys.path.insert(0, str(cand)); return


_setup_orion_import()
from ORION import read_TEC  # noqa: E402


def solution_time(path):
    with open(path) as fh:
        for line in fh:
            m = re.search(r"SOLUTIONTIME\s*=\s*([0-9.+\-eEdD]+)", line)
            if m:
                return float(m.group(1).replace("D", "E").replace("d", "e"))
            if line.lstrip().startswith(("0", "-", "+")):
                break
    return None


def load_mose(path):
    xa, _y, _z, var, _names = read_TEC(str(path))
    xn = np.asarray(xa[0])[:, 0, 0]
    xc = 0.5 * (xn[:-1] + xn[1:])
    rho = [np.asarray(var[0][k]).reshape(-1) for k in range(N_SPECIES)]
    Y = rho[I_CH4] / np.sum(rho, axis=0)
    T = np.asarray(var[0][I_T]).reshape(-1)
    n = min(xc.size, T.size, Y.size)
    return xc[:n], T[:n], Y[:n]


def norms(a, b):
    """Relative L2 and Linf of (a-b) normalised by the range of b."""
    scale = np.ptp(b) or 1.0
    e = a - b
    return np.sqrt(np.mean(e**2)) / scale, np.max(np.abs(e)) / scale


def main():
    here = Path(__file__).resolve().parent
    field = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "OUTPUT" / "field.tec"
    if not field.exists():
        sys.exit(f"MOSE field not found: {field}\n"
                 f"Run MOSE with Sc=1, Prl=1 first (OUTPUT/field.tec).")

    t = solution_time(field) or 50e-3
    xm, Tm, Ym = load_mose(field)
    xr, Tr, Yr = ref.solve(t, N=len(xm))          # same grid as MOSE
    Yr_i = np.interp(xm, xr, Yr)                   # guard against grid mismatch
    Tr_i = np.interp(xm, xr, Tr)

    l2_Y, li_Y = norms(Ym, Yr_i)
    l2_T, li_T = norms(Tm, Tr_i)

    print(f"MOSE vs Cantera unity-Lewis @ {t*1e3:.2f} ms")
    print(f"  CH4 center: MOSE {Ym[len(Ym)//2]:.4f}  Cantera {Yr_i[len(Ym)//2]:.4f}")
    print(f"  T   center: MOSE {Tm[len(Tm)//2]:.1f} K  Cantera {Tr_i[len(Tm)//2]:.1f} K")
    print(f"  CH4  error: L2 {l2_Y*100:.2f}%  Linf {li_Y*100:.2f}%  (of profile range)")
    print(f"  T    error: L2 {l2_T*100:.2f}%  Linf {li_T*100:.2f}%  (of profile range)")

    fig, axY = plt.subplots(figsize=(9, 5.5))
    axT = axY.twinx()
    xmm = xm * 1e3
    axY.plot(xmm, Ym,   '-',  color='C0', lw=2.4, label=r"MOSE  Y$_{CH_4}$")
    axY.plot(xmm, Yr_i, '--', color='C0', lw=1.6, label=r"Cantera  Y$_{CH_4}$")
    axT.plot(xmm, Tm,   '-',  color='C3', lw=2.4, label="MOSE  T")
    axT.plot(xmm, Tr_i, '--', color='C3', lw=1.6, label="Cantera  T")
    axY.set_xlabel("x [mm]"); axY.set_ylabel(r"CH$_4$ mass fraction")
    axT.set_ylabel("Temperature [K]"); axY.grid(True, alpha=0.3)
    lines = axY.get_lines() + axT.get_lines()
    axY.legend(lines, [l.get_label() for l in lines], frameon=False, loc="center", ncol=2)
    fig.suptitle(f"MOSE (Sc=1, Prl=1) vs Cantera unity-Lewis @ {t*1e3:.0f} ms")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = here / "mose_vs_cantera.png"
    fig.savefig(out, dpi=150)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
