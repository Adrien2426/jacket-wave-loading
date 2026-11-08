import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rho = 1025.0
CD = 1.2
CM = 0.8
D = 4.0           # Column diameter (m)
a = 6.0           # Wave amplitude (m)
d = 40.0          # Water depth (m)
k = 0.043         # Wave number (1/m)
omega = 0.628     # Angular frequency (rad/s)
x_column = 0.0    # Column position (m)

# FIXED INTEGRATION LIMITS: seabed to crest
z = np.linspace(-40.0, 6.0, 20001)

# Trapezoidal integration weights
dz = z[1] - z[0]
weights = np.full(z.size, dz)
weights[0] *= 0.5
weights[-1] *= 0.5

# HORIZONTAL VELOCITY AMPLITUDE
# u(z,t) = U(z) * sin(theta)
U = a * omega * np.cosh(k * (z + d)) / np.sinh(k * d)

# INTEGRATE THE DEPTH-DEPENDENT TERMS
drag_amplitude = (
    0.5 * rho * CD * D * np.sum(weights * U**2)
)

inertia_amplitude = (
    rho * CM * (np.pi * D**2 / 4)
    * omega * np.sum(weights * U)
)

# FORCES OVER ONE WAVE PERIOD
T = 2 * np.pi / omega
time = np.linspace(0, T, 20001)
theta = omega * time - k * x_column

# Signed drag: u * abs(u)
Fd = drag_amplitude * np.sin(theta) * np.abs(np.sin(theta))

# Signed inertia
Fi = inertia_amplitude * np.cos(theta)

# Add forces AT THE SAME TIME
Ftotal = Fd + Fi

# FIND MAXIMUM POSITIVE TOTAL FORCE
idx = np.argmax(Ftotal)
t_max = time[idx]
phase_max = np.degrees(theta[idx]) % 360

# RESULTS
print("\nVERTICAL COLUMN — HORIZONTAL FORCES")
print("Fixed integration limits: -40 m to +6 m")
print("All forces below are in kN.\n")

print(f"Maximum drag alone:     {np.max(Fd) / 1000:10.2f}")
print(f"Maximum inertia alone:  {np.max(Fi) / 1000:10.2f}")

print("\nAT THE MAXIMUM TOTAL FORCE:")
print(f"Drag contribution:      {Fd[idx] / 1000:10.2f}")
print(f"Inertia contribution:   {Fi[idx] / 1000:10.2f}")
print(f"TOTAL FORCE:            {Ftotal[idx] / 1000:10.2f}")
print(f"\nTime:                   {t_max:.3f} s")
print(f"Local phase:            {phase_max:.2f} degrees")

# PLOT
plt.figure(figsize=(10, 5))

plt.plot(time, Fd / 1000, label="Horizontal drag")
plt.plot(time, Fi / 1000, label="Horizontal inertia")
plt.plot(time, Ftotal / 1000, label="Horizontal total", linewidth=2)

plt.scatter(
    t_max, Ftotal[idx] / 1000,
    color="red", zorder=5,
    label=f"Maximum total = {Ftotal[idx] / 1000:.2f} kN"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Time (s)")
plt.ylabel("Horizontal force (kN)")
plt.title("Vertical column: drag + inertia")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()