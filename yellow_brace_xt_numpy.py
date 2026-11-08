"""Yellow y-z brace: normal, y and z forces from vertical velocity w.
Lecture separate-component Morison law; fixed crest limits -35 to +6 m.
x is constant along each brace. Total = drag + inertia from w only.
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
        'Dy': -s*drag, 'Iy': -s*inertia, 'Ty': -s*total,
        'Dz': c*drag, 'Iz': c*inertia, 'Tz': c*total,
    }

def main():
    x_brace = 0.0  # Opposite y-z brace, selected in the viewer
    T = 2*np.pi/omega
    time = np.linspace(0, T, 100001)
    F = forces(x_brace, time)
    local_phase = np.mod(np.degrees(omega*time-k*x_brace), 360)
    reference_phase = np.mod(np.degrees(omega*time), 360)
    groups = [
        ('NORMAL FROM w: positive towards -y and +z', ['Dn', 'In', 'Tn']),
        ('y COMPONENT: positive towards +y', ['Dy', 'Iy', 'Ty']),
        ('z COMPONENT: positive upward', ['Dz', 'Iz', 'Tz']),
    ]
    print('\nYELLOW BRACE: FORCES AS FUNCTIONS OF x AND t')
    print(f'x = {x_brace:.1f} m; fixed limits z = -35 to +6 m')
    print('Total = drag + inertia caused by w only.')
    print('theta = omega*t - k*x; normal force is resolved into y and z.')
    for heading, keys in groups:
        print('\n' + heading)
        print(f'{"Force":<10}{"Max + (kN)":>14}{"Time (s)":>12}'
              f'{"theta (deg)":>16}{"omega*t (deg)":>18}')
        for name, key in zip(['Drag', 'Inertia', 'Total'], keys):
            j = np.argmax(F[key])
            print(f'{name:<10}{F[key][j]/1000:>14.2f}{time[j]:>12.3f}'
                  f'{local_phase[j]:>16.2f}{reference_phase[j]:>18.2f}')
    j = np.argmax(F['Tn'])
    print('\nSIMULTANEOUS FORCES AT MAXIMUM POSITIVE NORMAL TOTAL:')
    for key in ['Dn', 'In', 'Tn', 'Ty', 'Tz']:
        print(f'{key}: {F[key][j]/1000:+.2f} kN')

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    for ax, (heading, keys) in zip(axes, groups):
        for name, key in zip(['Drag', 'Inertia', 'Total'], keys):
            ax.plot(time, F[key]/1000, label=name)
        j = np.argmax(F[keys[2]])
        ax.scatter(time[j], F[keys[2]][j]/1000, color='red', zorder=5)
        ax.axhline(0, color='black', linewidth=0.7)
        ax.set_title(heading)
        ax.set_ylabel('Force (kN)')
        ax.grid(alpha=0.3)
        ax.legend()
    axes[-1].set_xlabel('Time t (s)')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
