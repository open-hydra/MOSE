#!/usr/bin/env python3
"""Write gsi_ablator_table.dat, the pyrolysis-gas table of the unit check GSIAblator (gsi_ablator_check.f90): one table,
gas-only, of CO, H2, H2O and CO2 on 2 pressures x 3 temperatures, every row a state of the shift CO + H2O = CO2 + H2 of the
mixture 1 CO + 1 H2O + 0.5 H2 (moles) at the extent xi = 0.05 + 0.10 (iT - 1) + 0.05 (ip - 1): the elements C, H, O of
every row are those of that mixture, 1, 3 and 2 moles, i.e. the makeup C 12.011/M, H 3 1.008/M, O 2 15.999/M,
M = 47.033 g (the element weights of Cantera, those of phase.txt of test/common/CHO-N2-viscous).
Usage: python3 gsi_ablator_table.py > gsi_ablator_table.dat"""
W = {'CO': 28.010, 'H2': 2.016, 'H2O': 18.015, 'CO2': 44.009}
WC, WH, WO = 12.011, 1.008, 15.999
M = WC + 3 * WH + 2 * WO
e = lambda x: '%24.16E' % x
print('# pyrolysis-gas table of the unit check GSIAblator (test/unit/gsi_ablator_table.py)')
print('format 1')
print('tables 1')
print('table shift')
print('treatment gas-only')
print(f'elements 3 C {WC / M:.16E} H {3 * WH / M:.16E} O {2 * WO / M:.16E}')
print('condensed 0')
print('species 4')
for s in W:
    print(s)
print('pressures 2')
print(' '.join(e(p) for p in (5.0e4, 2.0e5)))
print('temperatures 3')
print(' '.join(e(T) for T in (1000.0, 2000.0, 3000.0)))
print('data')
for ip in range(2):
    for iT in range(3):
        xi = 0.05 + 0.10 * iT + 0.05 * ip
        n = {'CO': 1 - xi, 'H2': 0.5 + xi, 'H2O': 1 - xi, 'CO2': xi}
        m = {s: n[s] * W[s] for s in W}
        tot = sum(m.values())
        print(' '.join(e(m[s] / tot) for s in W))
print('end table')
