# December 2020 Past Paper
#------------------------------------------------------------------------------
#                                 Q1 
# optimize for 2 values given initial conditions and one equation
print('\n------ Q1. a)------\n')

import sympy as sym

# Define the variables
x, y, w, TA, y0 = sym.symbols('x y w TA y0')

# Define the equation C
C = (TA / w) * sym.cosh((w / TA) * x) + y0 - (TA / w)

# Differentiate equation_C with respect to x
dydx = sym.diff(C, x)

# Calculate the second derivative
d2yd2x = sym.diff(dydx, x)

# Express cosh in terms of sinh
d2yd2x_sub = d2yd2x.subs(sym.cosh(w / TA * x), sym.sqrt(1 + sym.sinh(w / TA * x)**2))

# Define the RHS of Equation B
B = (w / TA) * sym.sqrt(1 + (dydx**2))

# Display the results
print("Equation C:")
print(C)
print("\n1st Derivative of Equation C with respect to x:")
print(dydx)
print("\n2nd Derivative of Equation C with respect to x:")
print(d2yd2x_sub)
print("\nRight-hand side of Equation B:")
print(B)

# Simplify the expressions for comparison
d2yd2x_simp = sym.simplify(d2yd2x_sub)
B_simp = sym.simplify(B)

# Check if the second derrivative of C is equal to the RHS of B
is_equal = sym.simplify(d2yd2x_simp - B_simp) == 0

print("\nIs Equation C the general solution to Equation B?")
print(is_equal)

print('\n------ Q1. b)------\n')
#use a root finding techique to calc tension with initial conditions
import math

def bisection(f, a, b, N):

    if f(a) * f(b) >= 0:
        print("Initial values 'a' and 'b' must have opposite signs.")
        return None

    a_n = a
    b_n = b

    for n in range(1, N + 1):
        m_n = (a_n + b_n) / 2
        f_m_n = f(m_n)

        if f(a_n) * f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif f(b_n) * f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            print("Found exact solution.")
            return m_n
        else:
            print("Bisection method fails.")
            return None

    return (a_n + b_n) / 2

# constants provided in the question
y = 15
x = 50
y0 = 5
w = 10

#function set to = 0 so we can find the root
f = lambda TA: (TA / w) * math.cosh((w / TA) * x) + y0 - (TA / w) - y

# Attempt to find the root using the bisection method
# guesses were found using desmos
TA = bisection(f, 1250, 1300, 100)
print("Approximate Tension:", TA)

#------------------------------------------------------------------------------
#                                 Q2 
# solve the non-linear systems and linear systems
# report the sol for theta and w at t: 10, 20 ,30 for both cases
# inital cond: t = 0, w0 = 0 and theta0 = pi/4
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

# Constants
g = 9.81  # acceleration due to gravity (m/s^2)
l = 9.81  # length of the pendulum (m)
h = 0.001  # time step for Euler method (s)
t_max = 30  # maximum time (s)

# Function to compute derivatives for the non-linear system
def derivatives_nonlinear(theta, omega):
    dtheta_dt = omega
    domega_dt = -g / l * np.sin(theta)
    return dtheta_dt, domega_dt

# Function to compute derivatives for the linearized system
def derivatives_linear(theta, omega):
    dtheta_dt = omega
    domega_dt = -g / l * theta
    return dtheta_dt, domega_dt

# Function to solve the system using the Euler method
def solve_system(derivatives, initial_conditions):
    t_values = np.arange(0, t_max + h, h)
    theta_values = np.zeros_like(t_values)
    omega_values = np.zeros_like(t_values)

    # Set initial conditions
    theta_values[0] = initial_conditions['theta']
    omega_values[0] = initial_conditions['omega']

    # Euler method
    for i in range(1, len(t_values)):
        dtheta_dt, domega_dt = derivatives(theta_values[i - 1], omega_values[i - 1])
        theta_values[i] = theta_values[i - 1] + h * dtheta_dt
        omega_values[i] = omega_values[i - 1] + h * domega_dt

    return t_values, theta_values, omega_values

# Function to calculate the period of oscillation
def calculate_period(t_values, theta_values):
    # Find the indices where the system crosses theta = 0 from positive to negative
    zero_crossings = np.where(np.diff(np.sign(theta_values)) < 0)[0]

    # Calculate the period based on the time values at the zero crossings
    periods = np.diff(t_values[zero_crossings])
    average_period = np.mean(periods)

    return average_period

# Function to print the values of theta and omega at specified times
def print_values_at_times(t_values, theta_values, omega_values, times):
    header = ["Time (s)", "Theta", "Omega"]

    rows = []
    for time in times:
        index = int(time / h)
        rows.append([time, theta_values[index], omega_values[index]])

    print(tabulate(rows, headers=header, tablefmt="grid"))

# Function to find 𝜃0 for a given difference range
def find_theta0_for_difference_range(difference_range):
    theta0_range = np.arange(np.pi / 4, np.pi / 2, 0.01)  # Adjust the step size as needed

    for theta0 in theta0_range:
        initial_conditions = {'theta': theta0, 'omega': 0}

        t_values_nonlinear, _, _ = solve_system(
            derivatives_nonlinear, initial_conditions
        )

        t_values_linear, _, _ = solve_system(
            derivatives_linear, initial_conditions
        )

        period_nonlinear = calculate_period(t_values_nonlinear, _)
        period_linear = calculate_period(t_values_linear, _)

        difference_percentage = abs((period_linear - period_nonlinear) / period_nonlinear) * 100

        if difference_range[0] <= difference_percentage <= difference_range[1]:
            return theta0

    return None  # Return None if no suitable 𝜃0 is found within the range

# Evaluate periods and plot oscillations for different initial conditions
initial_conditions_list = [{'theta': np.pi / 4, 'omega': 0},
                            {'theta': np.pi / 2, 'omega': 0},
                            {'theta': np.pi / 8, 'omega': 0}]

for initial_conditions in initial_conditions_list:
    t_values_nonlinear, theta_values_nonlinear, omega_values_nonlinear = solve_system(
        derivatives_nonlinear, initial_conditions
    )

    t_values_linear, theta_values_linear, omega_values_linear = solve_system(
        derivatives_linear, initial_conditions
    )

    period_nonlinear = calculate_period(t_values_nonlinear, theta_values_nonlinear)
    period_linear = calculate_period(t_values_linear, theta_values_linear)

    print(f"Initial conditions: Theta0 = {initial_conditions['theta']:.4f}")
    print(f"Non-linear system period: {period_nonlinear:.4f} seconds")
    print(f"Linearized system period: {period_linear:.4f} seconds")

    # Print values of theta and omega at specified times
    print("\nValues at specified times:")
    print("Non-linear system:")
    print_values_at_times(t_values_nonlinear, theta_values_nonlinear, omega_values_nonlinear, [10, 20, 30])
    print("\nLinearized system:")
    print_values_at_times(t_values_linear, theta_values_linear, omega_values_linear, [10, 20, 30])
    print("\n")

    # Plot oscillations for each initial condition
    plt.figure(figsize=(10, 5))
    plt.plot(t_values_nonlinear, theta_values_nonlinear, label="Non-linearized")
    plt.plot(t_values_linear, theta_values_linear, label="Linearized")
    plt.xlabel("Time (seconds)")
    plt.ylabel("Theta (radians)")
    plt.title(f"Theta Motion: Theta0={initial_conditions['theta']:.4f}")
    plt.legend()
    plt.show()

# Find 𝜃0 for a specified difference range
difference_range = [0.95, 1.05]  # Adjust the range as needed
found_theta0 = find_theta0_for_difference_range(difference_range)

if found_theta0 is not None:
    print(f"The value of 𝜃0 for a difference between linear and non-linear periods within {difference_range[0]:.2f}% and {difference_range[1]:.2f}% is approximately: {found_theta0:.4f}")
else:
    print("No suitable 𝜃0 found within the specified difference range.")


#------------------------------------------------------------------------------
#                                 Q3 
print('\n------ Q3.------\n')
# calculate the inerest rate using a given equation

from scipy.optimize import newton

# Given values 
present_worth = 115000
annual_payments = 25500
number_of_years = 6

# Define the equation to find the interest rate, make it = 0
equation = lambda i: present_worth * (i * (1 + i)**number_of_years) / ((1 + i)**number_of_years - 1) - annual_payments

# Use Newton-Raphson method to find the root (interest rate in % )
interest_rate = newton(equation, x0=0.05) * 100  # Provide an initial guess (e.g., 0.05)

# Display the calculated interest rate
print(f"The interest rate is approximately: {interest_rate:.6f} %")

#------------------------------------------------------------------------------
#                                 Q4 
#find the displacement of the rocket after 30s using a given eq
print('\n------ Q4 ------\n')
import numpy as np

#using simpsons rule from lecure
def simps(f, a, b, N):

    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    
    dx = (b - a) / N
    x = np.linspace(a, b, N + 1)
    y = f(x)
    S = dx / 3 * np.sum(y[0:-1:2] + 4 * y[1::2] + y[2::2])
    return S

# given constants
u = 1.8 * (10**3)
m0 = 160 * (10**3)
q = 2.5 * (10**3)

#define the equation tht needs to be integrated (sub for x)
f = lambda x: u * np.log( m0 / (m0 - q * x))

#integrate f between 0 and 3 for N = 100
numerical_integration_result = simps(f, 0, 30, 100)
print(f"The distance traveled by the rocket is: {numerical_integration_result:.2f} m")

#------------------------------------------------------------------------------
#                                 Q5 
#find the slope of the tank s that will double the time it takes to drain the water
print('\n------ Q5 ------\n')    
                                                                                            
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

# Constants
A_orifice = 0.01  # Area of the orifice in m^2
h_initial = 1  # Initial height of the water in the tank
b_base = 1  # Base radius of the cylindrical tank
g = 9.81  # Acceleration due to gravity in m/s^2

# Calculate the flow rate using Torricelli's law
Q = A_orifice * np.sqrt(2 * g * h_initial)
print('The flow rate is', Q)

# Function for the integrand of the cylindrical tank
def integrand_cylinder(x):
    return (np.pi * b_base**2) / (A_orifice * np.sqrt(2 * g * x))

# Perform the integration for the cylindrical tank
integral_result_cylinder, error_cylinder = quad(integrand_cylinder, 0, h_initial)
print('Time to drain the cylindrical tank:', integral_result_cylinder)


# Function for the integrand of the conical tank
def integrand_conical(x, s):
    return (np.pi * (b_base + s * x)**2) / (A_orifice * np.sqrt(2 * g * x))

# Function to calculate the time to drain the conical tank args it to let the code know that s is only a constant and not to be integrated w.r.t
def time_to_drain_conical(h, s):
    return quad(integrand_conical, 0, h, args=(s))[0]

# Function to find the slope s that doubles the time to drain the tank
def find_slope(s):
    return time_to_drain_conical(h_initial, s) - 2 * integral_result_cylinder

# Initial guess for the slope s
initial_guess = 0.1
# Solve for s
s_solution = fsolve(find_slope, initial_guess)

# Output the solution for s
print ('value of s: ', s_solution)
