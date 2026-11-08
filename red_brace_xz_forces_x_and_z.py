import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho = 1025.0
Cd = 1.2
Cm = 0.8
D = 1.0

a = 6.0
k = 0.043
omega = 0.628
d = 40.0

z_min = -35.0
z_max = 6.0              # Fixed-crest approximation

rise = 42.0
span = 30.0

alpha = np.arctan2(rise, span)  
s = np.sin(alpha)
c = np.cos(alpha)
m = span / rise

T = 2*np.pi / omega

# Brace rises in +x.
# x = 0 at the SWL crossing, so x = m*z.


# ---------- DEPTH AND TIME GRIDS ----------
z = np.linspace(z_min, z_max, 4001)
t_list = np.linspace(0, T, 8000, endpoint=False)

dz = z[1] - z[0]

# Trapezoidal integration weights:
# integral(f dz) ≈ dz * (f0/2 + f1 + ... + fn/2)
weights = np.full(z.size, dz)
weights[0] *= 0.5
weights[-1] *= 0.5

B = a*omega / np.sinh(k*d)

cosh_R = np.cosh(k*(z + d))
sinh_R = np.sinh(k*(z + d))

Fd_x = np.zeros(t_list.size)
Fi_x = np.zeros(t_list.size)
Fd_z = np.zeros(t_list.size)
Fi_z = np.zeros(t_list.size)


# ---------- CALCULATE FORCES ----------
for i, t in enumerate(t_list):

    theta = omega*t - k*m*z

    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    # Horizontal and vertical velocities
    u = B*cosh_R*sin_theta
    w = B*sinh_R*cos_theta

    # Horizontal and vertical accelerations
    ax = B*omega*cosh_R*cos_theta
    az = -B*omega*sinh_R*sin_theta

    # Normal components for a brace rising in +x
    un = u*s - w*c
    an = ax*s - az*c

    # NumPy integration with respect to z
    I_drag = np.sum(weights * un * np.abs(un))
    I_inertia = np.sum(weights * an)

    # Normal forces: dl = dz/sin(alpha)
    Fd_normal = 0.5*rho*Cd*D * I_drag / s
    Fi_normal = rho*Cm*np.pi*D**2/4 * I_inertia / s

    # Horizontal components: positive in +x
    Fd_x[i] = Fd_normal*s
    Fi_x[i] = Fi_normal*s

    # Vertical components: positive upwards
    Fd_z[i] = -Fd_normal*c
    Fi_z[i] = -Fi_normal*c


# Add drag and inertia AT THE SAME TIME
Ft_x = Fd_x + Fi_x
Ft_z = Fd_z + Fi_z


# ---------- CLEAR RESULTS ----------
def print_results(title, force_arrays, positive=True):

    print("\n" + "="*72)
    print(title)
    print("="*72)

    print(
        f"{'Force':<15}"
        f"{'Magnitude (kN)':>17}"
        f"{'Time (s)':>16}"
        f"{'Phase (deg)':>18}"
    )
    print("-"*72)

    for name, force in zip(
        ["Drag", "Inertia", "TOTAL"], force_arrays
    ):
        if positive:
            idx = np.argmax(force)
        else:
            idx = np.argmin(force)

        magnitude = abs(force[idx])/1000
        time = t_list[idx]
        phase = np.degrees(omega*time)

        print(
            f"{name:<15}"
            f"{magnitude:>17.2f}"
            f"{time:>16.3f}"
            f"{phase:>18.2f}"
        )


print("\nDIAGONAL BRACE — WAVE FORCE RESULTS")
print(f"Brace angle: {np.degrees(alpha):.2f} degrees")
print(f"Water depth: {d:.1f} m")
print(f"Fixed integration limits: {z_min:.1f} to +{z_max:.1f} m")
print("All table entries are positive force MAGNITUDES.")

print_results(
    "HORIZONTAL — PEAK FORCE TOWARDS +x",
    [Fd_x, Fi_x, Ft_x],
    positive=True
)

print_results(
    "HORIZONTAL — PEAK FORCE TOWARDS -x",
    [Fd_x, Fi_x, Ft_x],
    positive=False
)

print_results(
    "VERTICAL — PEAK UPWARD FORCE",
    [Fd_z, Fi_z, Ft_z],
    positive=True
)

print_results(
    "VERTICAL — PEAK DOWNWARD FORCE",
    [Fd_z, Fi_z, Ft_z],
    positive=False
)

print("\nTOTAL is maximised after adding the forces at each time.")
print("It is NOT the sum of the separate drag and inertia maxima.")
print("Phase means omega*t at the reference position x = 0.")


# ---------- PLOTS ----------
fig, axes = plt.subplots(
    2, 1, figsize=(11, 8), sharex=True
)

# Include t=T to close the periodic curves
plot_times = np.append(t_list, T)

plot_settings = [
    (axes[0], [Fd_x, Fi_x, Ft_x], "Horizontal forces (+x)"),
    (axes[1], [Fd_z, Fi_z, Ft_z], "Vertical forces (+z upwards)")
]

for ax, force_arrays, title in plot_settings:

    for force, name, color in zip(
        force_arrays,
        ["Drag", "Inertia", "Total"],
        ["tab:blue", "tab:purple", "black"]
    ):
        plot_force = np.append(force, force[0])/1000

        ax.plot(
            plot_times,
            plot_force,
            label=name,
            color=color,
            linewidth=2.2 if name == "Total" else 1.4
        )

    # Mark maximum and minimum total force
    total = force_arrays[2]

    for idx in [np.argmax(total), np.argmin(total)]:
        ax.scatter(
            t_list[idx],
            total[idx]/1000,
            color="red",
            zorder=5
        )

    ax.axhline(0, color="grey", linewidth=0.7)
    ax.set_title(title)
    ax.set_ylabel("Force (kN)")
    ax.grid(alpha=0.3)
    ax.legend()

axes[1].set_xlabel("Time (s)")

fig.suptitle(
    "Diagonal brace — fixed-crest integration: -35 to +6 m"
)

plt.tight_layout()
plt.show()
