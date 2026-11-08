import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628

alpha = np.arctan2(42.0, 30.0)  # Angle from horizontal
dl_dz = 1.23
x_brace = 0.0

# FIXED-CREST DEPTH INTEGRATION
z = np.linspace(-35.0, 6.0, 20001)
weights = np.full(z.size, z[1] - z[0])
weights[[0, -1]] *= 0.5

# Vertical velocity amplitude and its normal projection
W = a * omega * np.sinh(k * (z + d)) / np.sinh(k * d)
Wn = W * np.cos(alpha)

# Normal drag and inertia amplitudes
AD = 0.5 * rho * CD * D * dl_dz * np.sum(weights * Wn * np.abs(Wn))
AI = rho * CM * np.pi * D**2 / 4 * dl_dz * omega * np.sum(weights * Wn)


def forces(x, t):
    """Horizontal y forces caused by w: drag, inertia, total, in N."""
    theta = omega * np.asarray(t) - k * np.asarray(x)

    Fd_n = AD * np.cos(theta) * np.abs(np.cos(theta))
    Fi_n = -AI * np.sin(theta)

    # Front brace falls in +y: positive normal points towards +y and +z.
    Fd_y = Fd_n * np.sin(alpha)
    Fi_y = Fi_n * np.sin(alpha)

    return Fd_y, Fi_y, Fd_y + Fi_y


def main():
    T = 2 * np.pi / omega
    time = np.linspace(0, T, 100001, endpoint=False)
    Fd, Fi, Ft = forces(x_brace, time)
    phase = np.mod(np.degrees(omega * time - k * x_brace), 360)

    print("\nYELLOW BRACE: HORIZONTAL y FORCE CAUSED BY w")
    print("Positive direction: +y. Fixed limits: z = -35 to +6 m.")
    print("Phase = omega*t - k*x, measured locally at the brace.")

    print("\nMAXIMUM POSITIVE FORCES — each row has its own peak")
    print(f'{"Force":<12}{"Maximum (kN)":>16}{"Time (s)":>12}{"Phase (deg)":>16}')

    for name, values in [("Drag", Fd), ("Inertia", Fi), ("Total", Ft)]:
        j = np.argmax(values)
        print(f"{name:<12}{values[j]/1000:>16.2f}"
              f"{time[j]:>12.3f}{phase[j]:>16.2f}")

    j = np.argmax(Ft)
    print("\nAT MAXIMUM TOTAL — values at the SAME instant")
    print(f"Drag    = {Fd[j]/1000:+.2f} kN")
    print(f"Inertia = {Fi[j]/1000:+.2f} kN")
    print(f"Total   = {Ft[j]/1000:+.2f} kN")

    plt.figure(figsize=(10, 5))
    for name, values in [("Drag from w", Fd),
                         ("Inertia from w", Fi),
                         ("Total from w", Ft)]:
        plt.plot(time, values/1000, label=name)

    plt.scatter(time[j], Ft[j]/1000, color="red", zorder=5)
    plt.axhline(0, color="black", linewidth=0.7)
    plt.xlabel("Time (s)")
    plt.ylabel("Horizontal force Fy (kN)")
    plt.title("Yellow brace: horizontal y force caused by w")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()