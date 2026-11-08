"""Run directly in Spyder.

Wave-only Morison forces on a brace rising in +x, with a 42 m rise
and 30 m horizontal span. Origin: x=0 at the brace's SWL crossing.
Fixed integration limits z=-35 to +6 m for drag AND inertia throughout
one cycle. This is the requested fixed-crest approximation.
Positive x is the wave direction; positive z is upwards.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt

# ---------- Parameters ----------
rho = 1025.0
Cd = 1.2
Cm = 0.8
D = 1.0
a = 6.0
k = 0.043
omega = 0.628
d = 40.0
z_min = -35.0
z_max = 6.0                 # Fixed crest elevation for BOTH force terms
rise = 42.0
span = 30.0

alpha = np.arctan2(rise, span)
s = np.sin(alpha)
c = np.cos(alpha)
m = span / rise
T = 2 * np.pi / omega
B = a * omega / np.sinh(k * d)


def normal_kinematics(z, t):
    theta = omega*t - k*m*z
    R = k*(z + d)
    u = B*np.cosh(R)*np.sin(theta)
    w = B*np.sinh(R)*np.cos(theta)
    ax = B*omega*np.cosh(R)*np.cos(theta)
    az = -B*omega*np.sinh(R)*np.sin(theta)
    # Unit normal n = (sin(alpha), 0, -cos(alpha)).
    return u*s - w*c, ax*s - az*c


def force_components(t):
    def drag_integrand(z):
        un, _ = normal_kinematics(z, t)
        return un*abs(un)

    I_drag = quad(drag_integrand, z_min, z_max,
                  epsabs=1e-7, epsrel=1e-8, limit=150)[0]
    I_inertia = quad(lambda z: normal_kinematics(z, t)[1],
                     z_min, z_max, epsabs=1e-7, epsrel=1e-8)[0]

    # dl = dz/sin(alpha): signed normal forces on the entire brace.
    Fd_n = 0.5*rho*Cd*D*I_drag/s
    Fi_n = rho*Cm*np.pi*D**2/4*I_inertia/s

    # Project using n: horizontal positive right, vertical positive up.
    Fd_x = Fd_n*s
    Fi_x = Fi_n*s
    Fd_z = -Fd_n*c
    Fi_z = -Fi_n*c

    # Add forces AT THE SAME TIME, then find the total's maximum.
    Ft_x = Fd_x + Fi_x
    Ft_z = Fd_z + Fi_z

    # This order is used in the results and plotting below.
    return np.array([Fd_x, Fi_x, Ft_x, Fd_z, Fi_z, Ft_z])


def find_extreme(t_list, values, component, maximum=True):
    sign = 1 if maximum else -1
    idx = np.argmax(sign*values[:, component])
    guess = t_list[idx]
    dt = t_list[1] - t_list[0]
    result = minimize_scalar(
        lambda t: -sign*force_components(t)[component],
        bounds=(guess-dt, guess+dt), method='bounded',
        options={'xatol': 1e-10})
    time = result.x % T
    force = force_components(time)[component]
    phase = np.degrees(omega*time)
    return force, time, phase


def main():
    t_list = np.linspace(0, T, 2000, endpoint=False)
    values = np.array([force_components(t) for t in t_list])
    names = ['Drag x', 'Inertia x', 'Total x',
             'Drag z', 'Inertia z', 'Total z']
    print(f'Brace angle = {np.degrees(alpha):.3f} degrees')
    print(f'Wave period = {T:.5f} s')
    print(f'Fixed integration limits: {z_min} to {z_max} m')
    print('Positive x: wave direction. Positive z: upwards.')
    print('Total = drag + inertia at the same instant.\n')
    print(f'{"Force":<12}{"Extreme":<10}{"kN":>12}{"Time (s)":>13}{"omega*t (deg)":>17}')
    extremes = {}
    for j, name in enumerate(names):
        for maximum, label in [(True, 'Maximum'), (False, 'Minimum')]:
            f, time, phase = find_extreme(t_list, values, j, maximum)
            extremes[j, maximum] = (f, time, phase)
            print(f'{name:<12}{label:<10}{f/1000:>12.3f}{time:>13.5f}{phase:>17.3f}')

    # Close the periodic curves at t=T.
    ts = np.r_[t_list, T]
    ys = np.vstack([values, values[0]])/1000
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for ax, columns, direction in [
            (axes[0], [0, 1, 2], 'Horizontal (+x)'),
            (axes[1], [3, 4, 5], 'Vertical (+z upwards)')]:
        for j, label, color in zip(columns, ['Drag', 'Inertia', 'Total'],
                                    ['tab:blue', 'tab:purple', 'black']):
            ax.plot(ts, ys[:, j], label=label, color=color,
                    linewidth=2.2 if label == 'Total' else 1.4)
        total_column = columns[-1]
        for maximum in [True, False]:
            f, time, _ = extremes[total_column, maximum]
            ax.scatter(time, f/1000, color='tab:red', zorder=5)
        ax.axhline(0, color='grey', linewidth=0.7)
        ax.set_ylabel('Force (kN)')
        ax.set_title(direction)
        ax.grid(alpha=0.3)
        ax.legend()
    axes[1].set_xlabel('Time (s)')
    fig.suptitle('Diagonal brace: fixed-crest integration, z = -35 to +6 m')
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
