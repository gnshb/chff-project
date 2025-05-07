import numpy as np
import matplotlib.pyplot as plt
def parse_cpp_output(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    data = {}
    u_count = 0
    p_count = 0
    current = None

    for line in lines:
        if line == 'u:':
            u_count += 1
            current = f'u{u_count}'
            data[current] = []
        elif line == 'p:':
            p_count += 1
            current = f'p{p_count}'
            data[current] = []
        elif current:
            numbers = list(map(float, line.replace('\t', ' ').split()))
            data[current].extend(numbers)

    for key in data:
        data[key] = np.array(data[key])

    return data

# Usage
data = parse_cpp_output('ex.txt')

# Example access
u_sr_simple, p_sr_simple = data['u1'], data['p1'][:-1]
u_no_sr_simple, p_no_sr_simple = data['u2'], data['p2'][:-1]
u_sr_analytic, p_sr_analytic = data['u3'], data['p3'][:-1]
u_no_sr_analytic, p_no_sr_analytic = data['u4'], data['p4'][:-1]

x = np.linspace(0, 2, 19)

# Velocity comparison
plt.plot(x, u_sr_analytic, label="Relativistic Exact", color='blue')
plt.plot(x, u_no_sr_analytic, label="Non-Relativistic Exact", color='green')
plt.scatter(x, u_sr_simple, label="Relativistic Simulation", color='purple')
plt.scatter(x, u_no_sr_simple, label="Non-Relativistic Simulation", color='orange')
plt.title("Velocity Comparison")
plt.xlabel("x axis")
plt.ylabel("Velocity (m/s)")
plt.grid(True)
plt.legend()
plt.savefig("velocity_comparison_plot.png", dpi=300)
plt.close()

# Pressure comparison
plt.plot(x, p_sr_analytic, label="Relativistic Exact", color='blue')
plt.plot(x, p_no_sr_analytic, label="Non-Relativistic Exact", color='green')
plt.scatter(x, p_sr_simple, label="Relativistic Simulation", color='purple')
plt.scatter(x, p_no_sr_simple, label="Non-Relativistic Simulation", color='orange')
plt.title("Pressure Comparison")
plt.xlabel("x axis")
plt.ylabel("Pressure (Pa)")
plt.grid(True)
plt.legend()
plt.savefig("pressure_comparison_plot.png", dpi=300)
