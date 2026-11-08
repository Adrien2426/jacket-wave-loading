"""Yellow y-z brace: forces caused by vertical wave velocity w.

Lecture scalar Morison law applied separately to w_n = w*cos(alpha),
then geometric resolution into y and z. This is not the coupled-speed
drag model. Total means drag + inertia caused by w only.
Fixed crest approximation: integrate z=-35 to +6 m at every phase.
dl/dz=1.23 as requested. Brace rises in +y at constant x=30 m, matching
the viewer. Positive normal direction is (0, -sin(alpha), cos(alpha)).
"""
import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628
alpha = np.arctan2(42.0, 30.0)
s, c = np.sin(alpha), np.cos(alpha)
length_factor = 1.23
x_brace = 30.0  # Constant x in the structure viewer

# DEPTH INTEGRATION: fixed crest limits for both force terms
z = np.linspace(-35.0, 6.0, 20001)
weights = np.full(z.size, z[1] - z[0])
weights[[0, -1]] *= 0.5

# w(z,t) = W(z)*cos(theta); theta is constant along this y-z brace.
W = a*omega*np.sinh(k*(z+d)) / np.sinh(k*d)
Wn = W*c

# Integrate the depth-dependent amplitudes once.
AD = 0.5*rho*CD*D*length_factor*np.sum(weights*Wn*np.abs(Wn))
AI = rho*CM*(np.pi*D**2/4)*length_factor*omega*np.sum(weights*Wn)

# LOCAL PHASE AT THE BRACE, covering a full cycle
theta = np.linspace(0.0, 2*np.pi, 100001)
local_deg = np.degrees(theta)
reference_phase = np.mod(theta + k*x_brace, 2*np.pi)
reference_deg = np.degrees(reference_phase)
time = reference_phase / omega

# SIGNED NORMAL FORCES CAUSED BY w
Fd_n = AD*np.cos(theta)*np.abs(np.cos(theta))
Fi_n = -AI*np.sin(theta)
Ft_n = Fd_n + Fi_n

# RESOLVE EACH FORCE INTO THE USER'S y AND z AXES
Fd_y, Fi_y, Ft_y = -s*Fd_n, -s*Fi_n, -s*Ft_n
Fd_z, Fi_z, Ft_z = c*Fd_n, c*Fi_n, c*Ft_n

def print_peaks(title, forces):
    print('\n' + title)
    print(f'{"Force":<12}{"Max positive (kN)":>20}'
          f'{"Local phase (deg)":>21}{"omega*t (deg)":>18}')
    for name, force in zip(['Drag', 'Inertia', 'Total'], forces):
        j = np.argmax(force)
        print(f'{name:<12}{force[j]/1000:>20.2f}'
              f'{local_deg[j]:>21.2f}{reference_deg[j]:>18.2f}')

print('\nYELLOW BRACE — FORCES CAUSED BY VERTICAL VELOCITY w')
print('Total = drag + inertia from w at the same instant.')
print('Fixed depth limits: -35 to +6 m; dl/dz = 1.23')
print(f'Brace angle = {np.degrees(alpha):.4f} degrees')
print(f'Constant brace x = {x_brace:.1f} m')
print('Local phase theta = omega*t - k*x_brace.')
print('Positive normal direction: towards -y and +z.')

print_peaks('NORMAL FORCE FROM w', [Fd_n, Fi_n, Ft_n])
print_peaks('FORCE IN +y', [Fd_y, Fi_y, Ft_y])
print_peaks('UPWARD FORCE (+z)', [Fd_z, Fi_z, Ft_z])
print('\nNegative extrema have the same magnitudes, half a cycle apart.')

j = np.argmax(Ft_n)
print('\nAT MAXIMUM POSITIVE NORMAL TOTAL:')
print(f'Normal drag:    {Fd_n[j]/1000:+.2f} kN')
print(f'Normal inertia: {Fi_n[j]/1000:+.2f} kN')
print(f'Normal total:   {Ft_n[j]/1000:+.2f} kN')
print(f'y total:        {Ft_y[j]/1000:+.2f} kN')
print(f'z total:        {Ft_z[j]/1000:+.2f} kN')
print(f'Time in the reference cycle: {time[j]:.3f} s')

# PLOTS: LOCAL PHASE, NOT TIME, ON THE HORIZONTAL AXIS
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
groups = [
    ('Normal force from w (+normal towards -y and +z)', [Fd_n, Fi_n, Ft_n]),
    ('y component (+y positive)', [Fd_y, Fi_y, Ft_y]),
    ('z component (upward positive)', [Fd_z, Fi_z, Ft_z])
]
for ax, (title, forces) in zip(axes, groups):
    for name, force in zip(['Drag', 'Inertia', 'Total'], forces):
        ax.plot(local_deg, force/1000, label=name)
    jmax = np.argmax(forces[2])
    ax.scatter(local_deg[jmax], forces[2][jmax]/1000, color='red', zorder=5)
    ax.axhline(0, color='black', linewidth=0.7)
    ax.set_title(title)
    ax.set_ylabel('Force (kN)')
    ax.grid(alpha=0.3)
    ax.legend()
axes[-1].set_xlabel('Local wave phase at the yellow brace (degrees)')
plt.tight_layout()
plt.show()
