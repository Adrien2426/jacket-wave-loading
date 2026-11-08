"""Beam 2 along x: vertical wave forces using the lecture's w|w| formula.
Source: 6_Morisons (1).pdf, printed page 83 / PDF page 24.
Uniform depth z=-35 m; x=0 to 30 m; waves along +x; +z is upward.
"""
import numpy as np
import matplotlib.pyplot as plt

rho, CD, CM = 1025.0, 1.2, 0.8
D, L = 2.0, 30.0
z_beam, d = -35.0, 40.0
a, k, omega = 6.0, 0.043, 0.628

x = np.linspace(0.0, L, 4001)
weights = np.full(x.size, x[1] - x[0])
weights[[0, -1]] *= 0.5
T = 2*np.pi / omega
time = np.linspace(0.0, T, 8001)
phase = omega * time
W = a*omega*np.sinh(k*(z_beam+d)) / np.sinh(k*d)

Fd_z = np.zeros_like(time)
Fi_z = np.zeros_like(time)
for j, phi in enumerate(phase):
    theta = phi - k*x
    w = W*np.cos(theta)
    az = -omega*W*np.sin(theta)
    Fd_z[j] = 0.5*rho*CD*D*np.sum(weights*w*np.abs(w))
    Fi_z[j] = rho*CM*np.pi*D**2/4*np.sum(weights*az)

Ftotal_z = Fd_z + Fi_z
print('\nBEAM 2 ALONG x: VERTICAL FORCES')
print('Integration: x=0 to 30 m, at z=-35 m')
print('Lecture formulation: signed drag w*abs(w)')
print('Positive force is upward; phase is omega*t at x=0.')
for heading, upward in [('MAXIMUM UPWARD', True), ('MAXIMUM DOWNWARD', False)]:
    print(f'\n{heading}')
    print(f'{"Force":<12}{"Force (kN)":>14}{"Phase (deg)":>16}{"Time (s)":>12}')
    for name, force in [('Drag', Fd_z), ('Inertia', Fi_z), ('Total', Ftotal_z)]:
        j = np.argmax(force) if upward else np.argmin(force)
        print(f'{name:<12}{force[j]/1000:>14.2f}'
              f'{np.degrees(phase[j]):>16.2f}{time[j]:>12.3f}')

idx = np.argmax(Ftotal_z)
print('\nAT MAXIMUM UPWARD TOTAL:')
print(f'Drag: {Fd_z[idx]/1000:+.2f} kN')
print(f'Inertia: {Fi_z[idx]/1000:+.2f} kN')
print(f'Total: {Ftotal_z[idx]/1000:+.2f} kN')

plt.figure(figsize=(10, 5))
for name, force in [('Drag', Fd_z), ('Inertia', Fi_z), ('Total', Ftotal_z)]:
    plt.plot(np.degrees(phase), force/1000, label=name)
plt.scatter(np.degrees(phase[idx]), Ftotal_z[idx]/1000, color='red', zorder=5)
plt.axhline(0, color='black', linewidth=0.7)
plt.xlabel('Reference phase at x=0 (degrees)')
plt.ylabel('Vertical force (kN), positive upward')
plt.title('Beam 2 along x: lecture-note formulation')
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
