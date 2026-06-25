"""Independent low-Mach unity-Lewis reference for the 1D diffusion test, with Cantera.

Apples-to-apples ground truth for MOSE run at Sc=1, Prl=1 (unity Lewis, Le=1):
D_k = D = mu/rho for every species and thermal diffusivity alpha = mu/rho (Pr=1).
Under Le=1 both the species mass fractions and the specific enthalpy obey the *same*
scalar diffusion equation with coefficient mu:

    d(rho Y_k)/dt + d(rho u Y_k)/dx = d/dx ( mu  dY_k/dx )
    d(rho h)/dt   + d(rho u h)/dx   = d/dx ( mu  dh/dx   )

(standard Le=1 collapse: k.grad T - sum h_k j_k = rho D grad h = mu grad h).

Crucially this version is *not* zero-velocity.  The T contrast makes the gas
expand/contract as heat diffuses, inducing a dilatational velocity that advects the
scalars -- exactly what MOSE's compressible solver carries.  At constant
thermodynamic pressure p0, continuity + the ideal-gas EOS fix the velocity
divergence (low-Mach constraint, derived from D ln(rho)/Dt = -div u and
rho = p0 W / (Ru T)):

    div u = D ln T/Dt - D ln W/Dt
          = [ Lh - sum_k h_k LY_k ] / (rho cp T)  +  W sum_k (LY_k / W_k) / rho

with LY_k = d/dx(mu dY_k/dx), Lh = d/dx(mu dh/dx) the diffusion operators and
h_k the species specific enthalpies.  In 1D periodic this is integrated directly
for the face velocity (the mean of div u is removed, i.e. p0 is held fixed).
Because rho = sum_k (rho Y_k), continuity is satisfied automatically by
conservative advection of the partial densities.

Solved with Cantera (GRI30) thermo/transport on the same grid as MOSE
(200 cells, L=0.05 m, periodic), constant p=101325 Pa.

Importable: solve(t_end, N, advection=True) -> (xc, T, Y_CH4).  advection=False
recovers the old zero-velocity diffusion limit.  Run as a script to also dump CSVs.
"""
import numpy as np
import cantera as ct
from pathlib import Path

# ---- setup (matches build_ic.py / input.ini) ------------------------------
L      = 0.05
p0     = 101325.0
x0     = 25e-3
d      = 2.5e-3
SP     = ['CH4', 'H2O', 'N2', 'O2']
Y_ox   = {'CH4': 0.0,   'H2O': 0.1, 'N2': 0.758, 'O2': 0.142}   # f=0 (center)
Y_fu   = {'CH4': 0.214, 'H2O': 0.0, 'N2': 0.591, 'O2': 0.195}   # f=1 (edges)
T_ox, T_fu = 1350.0, 320.0


def solve(t_end, N=200, cfl=0.2, advection=True):
    """Integrate the low-Mach unity-Lewis problem to t_end.

    Returns xc[m], T[K], Y_CH4 (N,).  advection=True adds the dilatational velocity
    induced by thermal expansion (apples-to-apples with MOSE); advection=False is
    the pure zero-velocity diffusion limit.
    """
    gas = ct.Solution('gri30.yaml')
    ksp = [gas.species_index(s) for s in SP]
    Wk  = gas.molecular_weights[ksp]            # kg/kmol, per tracked species
    nsp = len(SP)
    dx  = L / N
    xc  = (np.arange(N) + 0.5) * dx
    f   = 1.0 - np.exp(-((xc - x0) ** 2) / d ** 2)

    Y = np.zeros((N, nsp))
    for j, s in enumerate(SP):
        Y[:, j] = Y_ox[s] + (Y_fu[s] - Y_ox[s]) * f
    T = T_ox + (T_fu - T_ox) * f

    full = np.zeros(gas.n_species)

    def lap_flux(phi, mu):
        """d/dx( mu dphi/dx ) on the periodic grid (per unit volume)."""
        muf = 0.5 * (mu + np.roll(mu, -1))          # face i+1/2
        Fr  = muf * (np.roll(phi, -1) - phi) / dx   # flux at i+1/2
        Fl  = np.roll(Fr, 1)                         # flux at i-1/2
        return (Fr - Fl) / dx

    def state_from_TY(T, Y):
        """Properties at (T, p0, Y): rho, mu, cp, W, h, h_k(N,nsp)."""
        rho = np.empty(N); mu = np.empty(N); cp = np.empty(N)
        W = np.empty(N); h = np.empty(N); hk = np.empty((N, nsp))
        for i in range(N):
            full[:] = 0.0; full[ksp] = Y[i]
            gas.TPY = T[i], p0, full
            rho[i] = gas.density; mu[i] = gas.viscosity; cp[i] = gas.cp_mass
            W[i] = gas.mean_molecular_weight; h[i] = gas.enthalpy_mass
            hk[i] = gas.partial_molar_enthalpies[ksp] / Wk
        return rho, mu, cp, W, h, hk

    def upwind_div(q, uf):
        """d/dx( u q ) with first-order upwind face values (q is a conserved density)."""
        flux = np.where(uf >= 0.0, uf * q, uf * np.roll(q, -1))  # face i+1/2
        return (flux - np.roll(flux, 1)) / dx                    # cell i

    # conservative state: partial densities U_k = rho Y_k, total enthalpy E = rho h
    rho, mu, cp, W, h, hk = state_from_TY(T, Y)
    U = rho[:, None] * Y
    E = rho * h

    t = 0.0
    while t < t_end:
        rho = U.sum(axis=1)
        Y   = U / rho[:, None]
        h   = E / rho
        # recover T from (h, Y) at p0, and read transport/thermo at that state
        mu = np.empty(N); cp = np.empty(N); Wm = np.empty(N); hk = np.empty((N, nsp))
        for i in range(N):
            full[:] = 0.0; full[ksp] = Y[i]
            gas.HPY = h[i], p0, full
            T[i] = gas.T; mu[i] = gas.viscosity; cp[i] = gas.cp_mass
            Wm[i] = gas.mean_molecular_weight
            hk[i] = gas.partial_molar_enthalpies[ksp] / Wk

        dt = min(cfl * dx * dx / np.max(mu / rho), t_end - t)

        # diffusion operators
        LY = np.column_stack([lap_flux(Y[:, j], mu) for j in range(nsp)])  # (N,nsp)
        Lh = lap_flux(h, mu)

        # advance diffusion (conservative; sum_k LY_k == 0 so continuity is kept)
        U += dt * LY
        E += dt * Lh

        if advection:
            # low-Mach velocity divergence from thermal expansion
            S = (Lh - (hk * LY).sum(axis=1)) / (rho * cp * T) \
                + Wm * (LY / Wk).sum(axis=1) / rho
            S -= S.mean()                                # enforce periodicity (p0 fixed)
            uf = np.cumsum(S) * dx                        # u at face i+1/2, u_{-1/2}=0
            uf -= uf.mean()                               # remove arbitrary mean velocity
            for j in range(nsp):
                U[:, j] -= dt * upwind_div(U[:, j], uf)
            E -= dt * upwind_div(E, uf)

        t += dt

    rho = U.sum(axis=1)
    Y   = U / rho[:, None]
    h   = E / rho
    for i in range(N):
        full[:] = 0.0; full[ksp] = Y[i]
        gas.HPY = h[i], p0, full
        T[i] = gas.T
    return xc, T, Y[:, 0]


if __name__ == '__main__':
    import sys
    t_end = float(sys.argv[1]) if len(sys.argv) > 1 else 50e-3
    xc, T, Y_CH4 = solve(t_end)
    ref = Path(__file__).resolve().parent / 'reference'
    ref.mkdir(exist_ok=True)
    lab = f"{t_end*1e3:g}ms"
    np.savetxt(ref / f'CH4_{lab}.csv', np.column_stack([xc, Y_CH4]), delimiter=',')
    np.savetxt(ref / f'T_{lab}.csv',   np.column_stack([xc, T]),     delimiter=',')
    N = len(xc)
    print(f"Cantera low-Mach unity-Lewis reference @ {t_end*1e3:.1f} ms")
    print(f"  Y_CH4 center={Y_CH4[N//2]:.4f}  edge={Y_CH4[0]:.4f}")
    print(f"  T     center={T[N//2]:.1f} K  edge={T[0]:.1f} K")
    print(f"  wrote {ref}/CH4_{lab}.csv, T_{lab}.csv")
