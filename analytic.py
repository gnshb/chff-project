import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import argparse

# A1 = 100.0  # m^2 (cross-sectional area at point 1)
# P1 = 1e17  # Pa (100 kPa)
# N = 20  # meters

# plot = True
# SR = False

parser = argparse.ArgumentParser(description="Relativistic pipe solver")

parser.add_argument("-A", type=float, default=100.0,
                    help="Cross-sectional area at point 1 (default: 100.0 m^2)")
parser.add_argument("-P", type=float, default=1e17,
                    help="Pressure at point 1 (default: 1e17 Pa)")
parser.add_argument("-N", type=int, default=20,
                    help="Length of the pipe (default: 20 meters)")

parser.add_argument("--plot", dest="plot", action="store_true", help="Enable plotting (default: enabled)")
parser.add_argument("--no-plot", dest="plot", action="store_false", help="Disable plotting")
parser.set_defaults(plot=False)

parser.add_argument("--SR", dest="SR", action="store_true", help="Use Special Relativity (default: disabled)")
parser.set_defaults(SR=False)

args = parser.parse_args()

# Display the parsed values
A1=args.A
P1=args.P
N = args.N - 1
plot = args.plot
SR = args.SR


# Given values
A2 = 1  # m^2 (cross-sectional area at point 2)
P2 = 0   # Pa (50 kPa)
rho = 1.0    # kg/m^3 (density)
c = 3e8      # m/s (speed of light)


# Function for the relativistic Bernoulli equation
def bernoulli_normal(v1):
    # Lorentz factor for point 1

    # Velocity at point 2 from continuity equation
    v2 = (A1 / A2) * v1

    # Lorentz factor for point 2
    # Relativistic Bernoulli equation
    left = P1 + 0.5 * rho * v1**2
    right = P2 +  0.5 *rho * (v2**2)

    return left - right


def bernoulli_rel(v1):
    # Lorentz factor for point 1
    gamma1 = 1 / np.sqrt(1 - v1**2 / c**2)

    # Velocity at point 2 from continuity equation
    temp = (A1 / A2) * v1 * gamma1
    v2 = temp/np.sqrt(1+temp**2 / c**2)

    # Lorentz factor for point 2
    gamma2 = 1 / np.sqrt(1 - v2**2 / c**2)

    # Relativistic Bernoulli equation
    left = (4 * P1 / rho + 1) * gamma1 + 0.5 * v1**2 * gamma1**2
    right = (4 * P2 / rho + 1) * gamma2 + 0.5 * (v2**2) * gamma2**2

    return left - right


# Function to calculate velocity profile along the pipe
def velocity_profile_rel(L, A1, A2, v1_solution):
    # Generate positions along the pipe
    x = np.linspace(0, L, L)  # 100 points along the pipe

    # Interpolate cross-sectional area linearly along the pipe
    A_x = A1 + (A2 - A1) * x / L

    # Calculate the velocity at each point using the continuity equation
    gamma1 = 1 / np.sqrt(1 - v1_solution[0]**2 / c**2)
    temp = (A1 / A_x) * v1_solution[0] * gamma1
    v_x = temp/np.sqrt(1+temp**2 / c**2)

    return x, v_x

# Function to calculate velocity profile along the pipe
def velocity_profile_normal(L, A1, A2, v1_solution):
    # Generate positions along the pipe
    x = np.linspace(0, L, L)  # 100 points along the pipe

    # Interpolate cross-sectional area linearly along the pipe
    A_x = A1 + (A2 - A1) * x / L

    # Calculate the velocity at each point using the continuity equation
    v_x = (A1 / A_x) * v1_solution[0]

    return x, v_x

# Function to calculate the pressure profile along the pipe
def pressure_profile_rel(x, v_x, A1, A2, P1, P2, rho, c, v1_solution):
    # Pressure array
    p_x = np.zeros_like(x)

    # Use Bernoulli equation to calculate the pressure at each position
    for i, v in enumerate(v_x):
        # Calculate the Lorentz factor at the current velocity
        gamma = 1 / np.sqrt(1 - v**2 / c**2)
        gamma1 = 1 / np.sqrt(1 - v1_solution[0]**2 / c**2)
        # Calculate the corresponding pressure at point x[i]
        # From Bernoulli's equation:
        k1 = v1_solution[0]**2 * gamma1**2 /2
        k2 = v**2 * gamma**2
        s1 = (1+4*P1/rho)*gamma1

        p_x[i] = (rho/4)* ((k1 - k2 + s1)/gamma - 1)

    return p_x


def pressure_profile_normal(x, v_x, A1, A2, P1, P2, rho, c, v1_solution):
    # Pressure array
    p_x = np.zeros_like(x)

    # Use Bernoulli equation to calculate the pressure at each position
    for i, v in enumerate(v_x):
        # Calculate the Lorentz factor at the current velocity

        # Calculate the corresponding pressure at point x[i]
        # From Bernoulli's equation:
        k1 = rho * v1_solution[0]**2  / 2
        k2 = rho * v**2 /2
        s1 = P1

        p_x[i] =  ((k1 - k2 + s1))

    return p_x


if SR:
    bernoulli = bernoulli_rel
    velocity_profile = velocity_profile_rel
    pressure_profile = pressure_profile_rel
    text = "Relativistic"
else:
    bernoulli = bernoulli_normal
    velocity_profile = velocity_profile_normal
    pressure_profile = pressure_profile_normal
    text = "Non Reltivistic"

# Initial guess for velocity at point 1
v1_guess = 1e5  # m/s

# Solve for v1 using fsolve (find the root of the equation)
v1_solution = fsolve(bernoulli, v1_guess)

# Calculate the velocity profile
x, v_x = velocity_profile(N, A1, A2, v1_solution)

if plot:
    # Plot the velocity profile
    plt.figure(figsize=(10, 6))
    plt.plot(x, v_x, label="Velocity Profile", color='b')
    plt.xlabel('Position along pipe (m)')
    plt.ylabel('Velocity (m/s)')
    plt.title(f'{text} Velocity Profile in a Pipe with Varying Cross Section')
    plt.grid(True)
    plt.legend()
    plt.show()

p_x = pressure_profile(x, v_x, A1, A2, P1, P2, rho, c, v1_solution)

if plot:
    # Plot the pressure profile
    plt.figure(figsize=(10, 6))
    plt.plot(x, p_x, label="Pressure Profile", color='r')
    plt.xlabel('Position along pipe (m)')
    plt.ylabel('Pressure (kPa)')
    plt.title(f'{text} Pressure Profile in a Pipe with Varying Cross Section')
    plt.grid(True)
    plt.legend()
    plt.show()


if not plot:
    print("u:")
    print(*v_x)

    print("p:")
    print(*p_x, 0)
