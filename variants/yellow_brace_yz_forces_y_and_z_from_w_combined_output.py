"""Yellow y-z brace: normal, y and z forces from vertical velocity w.
Lecture separate-component Morison law; fixed crest limits -35 to +6 m.
x is constant along each brace. Total = drag + inertia from w only.
Front brace falls in +y: (0, 0, +7) to (0, 30, -35).
alpha is measured from horizontal; positive normal n=(0, sin(alpha), cos(alpha)).
Resolve each already-calculated normal force: Fy=Fn*sin(alpha), Fz=Fn*cos(alpha).
This keeps the established w-only scalar Morison model and rounded dl/dz=1.23.
"""
import numpy as np
import matplotlib.pyplot as plt

rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628
alpha = np.arctan2(42.0, 30.0)
s, c = np.sin(alpha), np.cos(alpha)
dl_dz = 1.23

z = np.linspace(-35.0, 6.0, 20001)
weights = np.full(z.size, z[1] - z[0])
weights[[0, -1]] *= 0.5
Wn = c*a*omega*np.sinh(k*(z+d)) / np.sinh(k*d)
AD = 0.5*rho*CD*D*dl_dz*np.sum(weights*Wn*np.abs(Wn))
AI = rho*CM*np.pi*D**2/4*dl_dz*omega*np.sum(weights*Wn)

def forces(x, t):
    """Return signed forces in N; x in m, t in s (scalars or arrays).
    The x coordinate must stay constant along the y-z brace.
    """
    theta = omega*np.asarray(t) - k*np.asarray(x)
    drag = AD*np.cos(theta)*np.abs(np.cos(theta))
    inertia = -AI*np.sin(theta)
    total = drag + inertia
    return {
        'Dn': drag, 'In': inertia, 'Tn': total,
        'Dy': s*drag, 'Iy': s*inertia, 'Ty': s*total,
        'Dz': c*drag, 'Iz': c*inertia, 'Tz': c*total,
    }

def main():
    x_brace = 0.0  # Front y-z brace, matching the corrected viewer
    T = 2*np.pi/omega
    time = np.linspace(0, T, 100001, endpoint=False)
    F = forces(x_brace, time)
    local_phase = np.mod(np.degrees(omega*time-k*x_brace), 360)
    reference_phase = np.mod(np.degrees(omega*time), 360)
    rows = [
        ('Drag', 'Dn', 'Dy', 'Dz'),
        ('Inertia', 'In', 'Iy', 'Iz'),
        ('Total', 'Tn', 'Ty', 'Tz'),
    ]
    groups = [
        ('HORIZONTAL FORCE IN y CAUSED BY w', ['Dy', 'Iy', 'Ty']),
        ('VERTICAL FORCE IN z CAUSED BY w', ['Dz', 'Iz', 'Tz']),
    ]
    print('\nYELLOW y-z BRACE: FORCES CAUSED BY VERTICAL VELOCITY w')
    print(f'x = {x_brace:.1f} m; fixed limits z = -35 to +6 m')
    print(f'Angle from horizontal = {np.degrees(alpha):.4f} degrees')
    print('Front brace falls in +y; positive normal points towards +y and +z.')
    print('Horizontal here means y. The w-only normal force has no x component.')
    print('Total = drag + inertia from w at the same instant.')
    print('\nMAXIMUM POSITIVE FORCES')
    print('Each row has its own peak. Its y and z components peak together.')
    print(f'{"Force":<10}{"y (kN)":>12}{"z (kN)":>12}{"Time (s)":>12}'
          f'{"theta (deg)":>16}{"omega*t (deg)":>18}')
    for name, normal, horizontal, vertical in rows:
        j = np.argmax(F[normal])
        print(f'{name:<10}{F[horizontal][j]/1000:>12.2f}'
              f'{F[vertical][j]/1000:>12.2f}{time[j]:>12.3f}'
              f'{local_phase[j]:>16.2f}{reference_phase[j]:>18.2f}')
    print('theta = omega*t - k*x: local wave phase at this brace.')
    print('Negative peaks have equal magnitudes, half a cycle later.')
    j = np.argmax(F['Tn'])
    print('\nAT THE MAXIMUM POSITIVE TOTAL (same instant for both directions)')
    print(f'{"Direction":<15}{"Drag (kN)":>14}{"Inertia (kN)":>16}{"Total (kN)":>14}')
    for label, keys in [('Horizontal y', ['Dy', 'Iy', 'Ty']),
                        ('Vertical z', ['Dz', 'Iz', 'Tz'])]:
        print(f'{label:<15}' + ''.join(f'{F[key][j]/1000:>{width}.2f}'
              for key, width in zip(keys, [14, 16, 14])))
    print('\nNormal force is the parent vector: do not add it to its y/z components.')
    print('The positive-normal peak magnitudes, for reference:')
    for name, normal, _, _ in rows:
        print(f'{name}: {np.max(F[normal])/1000:.2f} kN')

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for ax, (heading, keys) in zip(axes, groups):
        for name, key in zip(['Drag from w', 'Inertia from w', 'Total from w'], keys):
            ax.plot(time, F[key]/1000, label=name)
        j = np.argmax(F[keys[2]])
        ax.scatter(time[j], F[keys[2]][j]/1000, color='red', zorder=5)
        ax.axhline(0, color='black', linewidth=0.7)
        ax.set_title(heading)
        ax.set_ylabel('Force (kN)')
        ax.grid(alpha=0.3)
        ax.legend()
    axes[-1].set_xlabel('Time t (s)')
    fig.suptitle('Yellow brace: horizontal y and vertical z components from w')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
