import numpy as np
import math
import sys

def isentropicp2p1(M1,M2):
    delta = 0.5*(gamma-1)
    return ((1+delta*M1*M1)/(1+delta*M2*M2))**(0.5*gamma/delta)

def isentropicT2T1(M1,M2):
    delta = 0.5*(gamma-1)
    return ((1+delta*M1*M1)/(1+delta*M2*M2))

def nufunc(M, gamma):

    return math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(
        math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - math.atan(math.sqrt(M ** 2 - 1))

def mach_from_nu(nu, gamma, a=1.01, b=30.0, tol=1e-6, max_iter=1000000):
    """
    Calculate the Mach number from the Prandtl-Meyer angle using the bisection method.

    Parameters:
    - nu: Prandtl-Meyer angle in degrees
    - gamma: Specific heat ratio (ratio of specific heats, typically for air ~ 1.4)
    - a: Lower bound for the Mach number
    - b: Upper bound for the Mach number
    - tol: Tolerance for convergence
    - max_iter: Maximum number of iterations

    Returns:
    - M: Mach number
    """
    nu_rad = nu

    for _ in range(max_iter):
        M_mid = (a + b) / 2
        nu_err = nufunc(M_mid, gamma) - nu_rad

        if abs(nu_err) < tol:
            return M_mid

        if np.sign(nu_err)<0:
            a = M_mid
        else:
            b = M_mid

    raise RuntimeError("Bisection method did not converge within the specified number of iterations.")

def prandtl_meyer_expansion(theta, M):
    # Calculate the Mach angle (ν) using the Prandtl-Meyer function
    nu1 = nufunc(M,gamma)

    if (theta<0):
        nu2 = nu1-theta*np.pi/180
    else:
        nu2 = nu1+theta*np.pi/180

    # Calculate the Mach number (M) using the inverse of the Prandtl-Meyer function
    M2 = mach_from_nu(nu2, gamma)

    return M2

# Given parameters
gamma = 1.4  # Specific heat ratio for air

# Input values
M1 = 2.0          #float(input("Enter the initial Mach number (M1): "))
theta = 34.971    #float(input("Enter the expansion angle in degrees: "))

# Calling the function
M2 = prandtl_meyer_expansion(theta, M1)

# Displaying the results
# print("\nResults:")
# print(f"Initial Mach number (M1): {M1}")
# print(f"Expansion angle (theta): {theta} degrees")
# print(f"Final Mach number (M2): {M2}")
# print(f"Temperature ratio (T2/T1): {isentropicT2T1(M1,M2)}")
# print(f"Pressure ratio (p2/p1): {isentropicp2p1(M1,M2)}")
# print(f"Density ratio (rho2/rho1): {isentropicp2p1(M1,M2)/isentropicT2T1(M1,M2)}")

with open("OUTPUT/exit.txt", 'r') as file:
        lines = file.readlines()  # Read all lines into a list
        last_line = lines[-1]  # Get the last line
        columns = last_line.strip().split()  # Split the last line into columns
        CFDval = columns[1]  # Return the second column

# print(f"CFD Final Mach number: {CFDval}")

perr = abs(float(CFDval)-M2)/M2*100
if (perr > 0.5):
    print('Mach number mismatch (CFD vs analytical)')
    print(CFDval, M2)
    print('Error = ',perr, '%, threshold = ',0.5,'%')
    sys.exit(1)
else:
    print('Correct!')


