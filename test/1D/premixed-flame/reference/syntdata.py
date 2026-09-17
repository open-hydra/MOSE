import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

summary_rows = []

# --------------------------------------------------
# Configuration
# --------------------------------------------------

mechanism = "ZK.yaml"

phis = [1.0]

Tin = 300.0
Pin = ct.one_atm

width = 0.05

species_to_save = [
    "CH4",
    "O2",
    "CO",
    "CO2",
    "H2O",
    "OH"
]

flame_speeds = []

# --------------------------------------------------
# Loop over equivalence ratios
# --------------------------------------------------

for phi in phis:

    print(f"\nRunning phi = {phi}")

    gas = ct.Solution(mechanism)

    gas.set_equivalence_ratio(
        phi,
        fuel="CH4",
        oxidizer="O2:1.0,N2:3.76"
    )

    gas.TP = Tin, Pin

    yCH4 = gas["CH4"].Y[0]
    yO2  = gas["O2"].Y[0]
    yN2  = gas["N2"].Y[0]

    rho = gas.density

    flame = ct.FreeFlame(gas, width=width)
    flame.transport_model = "unity-Lewis-number"
    flame.set_refine_criteria(ratio=3,slope=0.05,curve=0.1)
    flame.solve(loglevel=1, auto=True)

    Su = flame.velocity[0]

    flame_speeds.append(Su)

    print(f"Flame speed = {Su:.4f} m/s")

    # ------------------------------------------
    # Derive additional quantities
    # ------------------------------------------

    hrr = np.zeros(flame.flame.n_points)
    mach = np.zeros(flame.flame.n_points)

    for i in range(flame.flame.n_points):

        gas.TPY = (flame.T[i],flame.P,flame.Y[:, i])

        hrr[i] = -np.dot(gas.net_production_rates,gas.partial_molar_enthalpies)
        mach[i] = flame.velocity[i] / gas.sound_speed

    # ------------------------------------------
    # Export profile
    # ------------------------------------------

    data = {
        "x": flame.grid,
        "T": flame.T,
        "u": flame.velocity,
        "rho": flame.density,
        "Mach": mach,
        "HRR": hrr
    }

    for sp in species_to_save:
        data[f"Y_{sp}"] = flame.Y[
            gas.species_index(sp)
        ]

    df = pd.DataFrame(data)

    df.to_csv(
        f"flame_phi_{phi:.1f}.csv",
        index=False
    )

    mdot = rho * Su

    summary_rows.append({
    "phi": phi,
    "Pin": gas.P,
    "Tin": gas.T,
    "rho": rho,
    "Su": Su,
    "mdot": mdot,
    "Y_CH4": yCH4,
    "Y_O2": yO2,
    "Y_N2": yN2
    })

# --------------------------------------------------
# Flame speed summary
# --------------------------------------------------

summary = pd.DataFrame({
    "phi": phis,
    "Su_mps": flame_speeds
})

summary.to_csv(
    "flame_speeds.csv",
    index=False
)

plt.figure(figsize=(8,5))

for phi in phis:

    df = pd.read_csv(
        f"flame_phi_{phi:.1f}.csv"
    )

    plt.plot(
        df["x"]*1000,
        df["T"],
        label=f"phi={phi}"
    )

plt.xlabel("x [mm]")
plt.ylabel("Temperature [K]")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    "temperature_profiles.png",
    dpi=300
)

plt.show()

plt.figure(figsize=(8,5))

for phi in phis:

    df = pd.read_csv(
        f"flame_phi_{phi:.1f}.csv"
    )

    plt.plot(
        df["x"]*1000,
        df["Y_OH"],
        label=f"phi={phi}"
    )

plt.xlabel("x [mm]")
plt.ylabel("OH mass fraction")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    "OH_profiles.png",
    dpi=300
)

plt.show()

plt.figure(figsize=(8,5))

for phi in phis:

    df = pd.read_csv(
        f"flame_phi_{phi:.1f}.csv"
    )

    plt.plot(
        df["x"]*1000,
        df["u"],
        label=f"phi={phi}"
    )

plt.xlabel("x [mm]")
plt.ylabel("u [m/s]")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    "u_profiles.png",
    dpi=300
)

plt.show()

plt.figure(figsize=(8,5))

for phi in phis:

    df = pd.read_csv(
        f"flame_phi_{phi:.1f}.csv"
    )

    plt.plot(
        df["x"]*1000,
        df["Mach"],
        label=f"phi={phi}"
    )

plt.xlabel("x [mm]")
plt.ylabel("Mach number")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    "Mach_profiles.png",
    dpi=300
)

plt.show()

pd.DataFrame(summary_rows).to_csv(
    "verification_cases.csv",
    index=False
)