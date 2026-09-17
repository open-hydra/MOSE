#!/usr/bin/env python3
"""Verify the ablating-wall case against its analytic solution.

The case (see build_case.py) is a gas column blown from below by an ablating
wall and vented at the top, initialised isothermal at the wall temperature.
That state is an exact steady solution of the governing equations:

  * the column is isothermal, so the conductive flux at the wall vanishes and
    the surface energy balance reduces to  q_pyro(Tw) = qrad  -- one scalar
    equation whose root is computed here, independently of the solver, by
    bisection (melting instead prescribes Tw, and mdot follows in closed form);
  * with p and T uniform the density is uniform, so continuity makes the
    blowing velocity uniform at  v = mdot / rho  all the way up the column.

The solver is therefore checked against a reference it never sees, on both the
wall quantities it reports and the field it produces.

Exit status: 0 = pass, 1 = fail.
"""

import re
import sys

import numpy as np

from gsi_models import (MELT, PYRO_ID, mdot_melting, mdot_pyrolysis,
                        q_pyrolysis, wall_temperature_radiative_equilibrium)

BC_MELTING, BC_PYROLYSIS, BC_SRM = 503, 504, 502

failures = []
checks = 0


def check(what, ok, detail=''):
    global checks
    checks += 1
    status = 'ok  ' if ok else 'FAIL'
    print(f'   [{status}] {what}{("  " + detail) if detail else ""}')
    if not ok:
        failures.append(what)


def close(what, got, want, rtol, unit=''):
    denom = max(abs(got), abs(want))
    rel = abs(got - want) / denom if denom > 0 else 0.0
    check(what, rel <= rtol,
          f'got {got:.8g}{unit}  want {want:.8g}{unit}  rel {rel:.2e} (tol {rtol:g})')


def read_tec(path):
    """Minimal reader for the Tecplot BLOCK files MOSE writes."""
    with open(path) as fh:
        lines = fh.read().split('\n')
    names = re.findall(r'"([^"]+)"', lines[0])
    zone = lines[1]
    I, J, K = (int(re.search(rf'{c}=\s*(\d+)', zone).group(1)) for c in 'IJK')
    vals = [float(v) for v in ' '.join(lines[2:]).split()]
    n_node = I * J * K
    n_cell = max(I - 1, 1) * max(J - 1, 1) * max(K - 1, 1)
    out, pos = {}, 0
    for idx, name in enumerate(names):
        n = n_node if idx < 3 else n_cell
        out[name] = np.array(vals[pos:pos + n])
        pos += n
    return out, (I, J, K)


def read_case():
    cfg = {}
    for line in open('INPUT/case.txt'):
        key, val = line.split(None, 1)
        cfg[key] = val.strip()
    return cfg


def verify_inert(cfg, wall, field, dims, wall_face, T_init):
    """An inert surface: the flow cannot sustain the ablation, so the boundary
    condition must fall back to an impermeable wall and inject nothing."""
    I, J, K = dims
    Tw_ref = float(cfg['Tw_prescribed'])

    check('no mass is injected (mdot is exactly zero)',
          np.all(wall['mdot'] == 0.0), f'max |mdot| = {np.max(np.abs(wall["mdot"])):.3e}')
    close('wall held at the prescribed temperature',
          float(np.mean(wall['Tw'])), Tw_ref, 1e-12, ' K')
    check('surface is losing heat to the gas, which is why it cannot melt',
          np.all(wall['qw'] < 0.0), f'qw = {float(np.mean(wall["qw"])):.4g} W/m2')

    shape = (max(K - 1, 1), max(J - 1, 1), max(I - 1, 1))
    T = field['T'].reshape(shape)
    v = field['v'].reshape(shape)
    rho = field['rho(1)'].reshape(shape)

    near = 0 if wall_face == 3 else -1
    check('the hot wall conducts heat into the gas',
          float(np.mean(T[:, near, :])) > T_init,
          f'T(wall cell) = {float(np.mean(T[:, near, :])):.2f} K  vs  T_init = {T_init:.2f} K')
    check('gas never exceeds the wall temperature',
          float(np.max(T)) <= Tw_ref * (1.0 + 1e-6), f'max T = {float(np.max(T)):.2f} K')

    # Impermeable wall: whatever motion is left is thermal expansion venting, an
    # order of magnitude below the mass flux an ablating surface of this size blows.
    mass_flux = float(np.max(np.abs(rho * v)))
    check('no blowing, only weak expansion-driven flow',
          mass_flux < 5e-2, f'max |rho*v| = {mass_flux:.3e} kg/m2/s')

    print(f'   -> {checks - len(failures)}/{checks} checks passed')
    return 1 if failures else 0


def verify_srm(cfg, field, dims, wall_face):
    """Burning propellant grain (BC 502).

    The rate is set by the pressure through Saint-Robert's law rather than by an
    energy balance, and the wall is held at the flame temperature, but the blown
    column is the same exact solution as for the ablative walls.

    BC 502 is not a viscous wall, so it contributes nothing to OUTPUT/wall.tec and
    everything here is read from the field.  That is the stronger check anyway: the
    column can only stay at Taf if the injected products carry their enthalpy, so
    this case fails loudly if the surface injects cold gas.
    """
    I, J, K = dims
    Taf = float(cfg['srm_Taf'])
    a = float(cfg['srm_a'])
    n = float(cfg['srm_n'])
    p_ref = float(cfg['p'])
    R_gas = float(cfg['R_gas'])

    shape = (max(K - 1, 1), max(J - 1, 1), max(I - 1, 1))
    T = field['T'].reshape(shape)
    v = field['v'].reshape(shape)
    rho = field['rho(1)'].reshape(shape)
    p = field['p'].reshape(shape)

    # Saint-Robert's law evaluated at the pressure the solver actually reached.
    mdot_ref = a * (float(np.mean(p)) / p_ref) ** n

    dT = np.max(np.abs(T - Taf)) / Taf
    check('column stays at the flame temperature',
          dT < 5e-3, f'max |T - Taf| / Taf = {dT:.2e}  (Taf = {Taf:g} K)')

    dp = np.ptp(p) / p_ref
    check('pressure stays uniform', dp < 5e-3, f'spread {dp:.2e}')

    flux = np.abs(rho * v)
    err = np.max(np.abs(flux - mdot_ref)) / mdot_ref
    check('|rho*v| = a (p/pRef)^n through the column',
          err < 5e-3, f'max rel deviation {err:.2e}  (mdot = {mdot_ref:.6g} kg/m2/s)')

    expected_sign = 1.0 if wall_face == 3 else -1.0
    check('blowing is directed away from the grain',
          np.all(np.sign(v) == expected_sign),
          f'mean v = {float(np.mean(v)):+.4g} m/s')

    rho_ref = p_ref / (R_gas * Taf)
    close('blowing speed = mdot / rho', float(np.mean(np.abs(v))),
          mdot_ref / rho_ref, 5e-3, ' m/s')

    print(f'   -> {checks - len(failures)}/{checks} checks passed')
    return 1 if failures else 0


def main():
    cfg = read_case()
    name = cfg['config']
    bc = int(cfg['bc'])
    qrad = float(cfg['qrad'])
    gsi_id = int(cfg['gsi_id'])
    p_ref = float(cfg['p'])
    R_gas = float(cfg['R_gas'])
    wall_face = int(cfg['wall_face'])
    T_init = float(cfg['T_init'])
    inert = bool(int(cfg['inert']))

    print(f'\n== ablating wall: {name} (BC {bc}, face {wall_face}, '
          f'qrad = {qrad:g} W/m2) ==')

    field, (I, J, K) = read_tec('OUTPUT/field.tec')

    # A burning grain is not a viscous wall, so MOSE writes no wall file for it.
    if bc == BC_SRM:
        check('solution is finite everywhere',
              all(np.all(np.isfinite(v)) for v in field.values()))
        return verify_srm(cfg, field, (I, J, K), wall_face)

    wall, _ = read_tec('OUTPUT/wall.tec')

    # The BC reports mdot and qw signed by the face orientation: the factor is
    # +1 on a low face and -1 on a high one. Undo it, so every identity below is
    # checked on the physical (unsigned) quantities.
    orient = 1.0 if wall_face == 3 else -1.0
    Tw_s = wall['Tw']
    qw_s = orient * wall['qw']
    md_s = orient * wall['mdot']

    # ---- 0. nothing went non-finite -------------------------------------
    finite = all(np.all(np.isfinite(v)) for v in field.values()) and \
             all(np.all(np.isfinite(v)) for v in wall.values())
    check('solution is finite everywhere', finite)

    if inert:
        return verify_inert(cfg, wall, field, (I, J, K), wall_face, T_init)

    # ---- 1. the wall state is one-dimensional ---------------------------
    for label, arr in (('Tw', Tw_s), ('mdot', md_s), ('qw', qw_s)):
        spread = np.ptp(arr) / max(abs(np.mean(arr)), 1e-30)
        check(f'{label} uniform along the wall',
              spread < 1e-8, f'spread {spread:.2e}')

    Tw, qw, mdot = float(np.mean(Tw_s)), float(np.mean(qw_s)), float(np.mean(md_s))
    check('blown mass flux is positive (mass leaves the solid)', mdot > 0.0,
          f'mdot = {mdot:+.6g} kg/m2/s')

    # ---- 2. analytic reference, derived here and not by the generator ---
    if bc == BC_MELTING:
        Tw_ref = float(cfg['Tw_prescribed'])
        mdot_ref = mdot_melting(qrad, Tw=Tw_ref)
        close('wall temperature = prescribed value', Tw, Tw_ref, 1e-12, ' K')
    else:
        model = PYRO_ID[gsi_id]
        Tw_ref = wall_temperature_radiative_equilibrium(model, qrad)
        mdot_ref = mdot_pyrolysis(model, Tw_ref)
        close('wall temperature = radiative-equilibrium root',
              Tw, Tw_ref, 1e-3, ' K')

    close('blown mass flux = analytic value', mdot, mdot_ref, 1e-3, ' kg/m2/s')

    # ---- 3. exact algebraic relations of the wall model ------------------
    # These hold for the state the solver actually reached, whatever it is, so
    # they test the boundary condition rather than the convergence.
    if bc == BC_MELTING:
        absorbed = MELT['dh'] + MELT['cp'] * (Tw - MELT['Ti'])
        close('mdot * heat_absorbed = heat into the wall',
              mdot * absorbed, qw, 1e-9, ' W/m2')
    else:
        model = PYRO_ID[gsi_id]
        close('mdot = Arrhenius law at the reported Tw',
              mdot, mdot_pyrolysis(model, Tw), 1e-9, ' kg/m2/s')
        # Surface energy balance: with the column isothermal the conductive
        # flux the BC reports is what has to balance the radiative input.
        close('q_pyro(Tw) = qrad + q_conv  (surface energy balance)',
              q_pyrolysis(model, Tw), qrad + qw, 2e-3, ' W/m2')

    # ---- 4. the isothermal column is preserved ---------------------------
    shape = (max(K - 1, 1), max(J - 1, 1), max(I - 1, 1))
    T = field['T'].reshape(shape)
    v = field['v'].reshape(shape)
    rho = field['rho(1)'].reshape(shape)
    p = field['p'].reshape(shape)

    dT = np.max(np.abs(T - Tw)) / Tw
    check('column stays isothermal at the wall temperature',
          dT < 5e-3, f'max |T - Tw| / Tw = {dT:.2e}')

    dp = np.ptp(p) / p_ref
    check('pressure stays uniform', dp < 5e-3, f'spread {dp:.2e}')

    # ---- 5. continuity: the injected mass is what flows along the column ---
    # Every cell counts, the wall-adjacent one included: the boundary condition
    # carries the wall pressure itself, so the blown mass flux is delivered right
    # from the first cell.
    interior = slice(None)
    flux = np.abs((rho * v)[:, interior, :])
    err = np.max(np.abs(flux - mdot)) / mdot
    check('|rho*v| = mdot through the column',
          err < 5e-3, f'max rel deviation {err:.2e}')

    # Mass must leave the wall, i.e. travel away from it: +y off face 3, -y off face 4.
    v_int = v[:, interior, :]
    expected_sign = 1.0 if wall_face == 3 else -1.0
    check('blowing is directed away from the wall',
          np.all(np.sign(v_int) == expected_sign),
          f'mean v = {float(np.mean(v_int)):+.4g} m/s')

    rho_ref = p_ref / (R_gas * Tw)
    close('blowing speed = mdot / rho', float(np.mean(np.abs(v_int))),
          mdot / rho_ref, 5e-3, ' m/s')

    print(f'   -> {checks - len(failures)}/{checks} checks passed')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
