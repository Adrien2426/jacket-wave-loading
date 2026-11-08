"""NumPy implementation of IMG_4392's simplified handwritten expression.

Assumptions: common wave phase along the member, horizontal velocity only,
fixed z=-35 to +6 m integration. This is NOT the earlier full x-z brace
model with spatially varying phase and vertical wave velocity.
The requested rounded length conversion dl/dz = 1.23 is used.
"""
import numpy as np
import matplotlib.pyplot as plt

rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628
rise, span = 42.0, 30.0
angle = np.arctan2(rise, span)
length_factor = 1.23

z = np.linspace(-35.0, 6.0, 20001)
weights = np.full(z.size, z[1] - z[0])
weights[[0, -1]] *= 0.5
R = k * (z + d)
I1 = np.sum(weights * np.cosh(R)**2)
I2 = np.sum(weights * np.cosh(R))

B = a * omega / np.sinh(k * d)
C1 = 0.5 * rho * CD * D * B**2
C2 = rho * CM * np.pi * D**2 / 4 * omega * B
AD = length_factor * C1 * I1
AI = length_factor * C2 * I2

phase = np.linspace(0, 2*np.pi, 20001)
time = phase / omega
# Signed extension of the sheet's sin^2 expression to the whole cycle.
Fd = AD * np.sin(phase) * np.abs(np.sin(phase))
Fi = AI * np.cos(phase)
Ft = Fd + Fi
idx = np.argmax(Ft)

print('\nHANDWRITTEN MODEL: COMMON PHASE, HORIZONTAL VELOCITY ONLY')
print('Fixed integration limits: z = -35 to +6 m')
print(f'Angle = {np.degrees(angle):.4f} degrees')
print(f'dl/dz used = {length_factor:.2f}; exact = {1/np.sin(angle):.6f}')
print(f'Depth integrals: I1 = {I1:.4f} m; I2 = {I2:.4f} m')
print('\nForce          Maximum positive (kN)   Reference phase (deg)')
for name, values in [('Drag', Fd), ('Inertia', Fi), ('Total', Ft)]:
    j = np.argmax(values)
    print(f'{name:<15}{values[j]/1000:>19.2f}{np.degrees(phase[j]):>24.2f}')
print(f'\nAt total maximum: drag = {Fd[idx]/1000:.2f} kN, '
      f'inertia = {Fi[idx]/1000:.2f} kN, time = {time[idx]:.3f} s')

plt.figure(figsize=(10, 5))
for values, name in [(Fd, 'Drag'), (Fi, 'Inertia'), (Ft, 'Total')]:
    plt.plot(np.degrees(phase), values/1000, label=name)
plt.scatter(np.degrees(phase[idx]), Ft[idx]/1000, color='red', zorder=5)
plt.axhline(0, color='black', linewidth=0.7)
plt.xlabel('Common reference phase (degrees)')
plt.ylabel('Horizontal force (kN)')
plt.title('Handwritten approximation: fixed crest, horizontal velocity only')
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
