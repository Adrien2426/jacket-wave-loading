"""Wave-only horizontal drag and inertia on a straight x-z brace.

Run: python3 variants/red_brace_xz_force_x_surface_options.py
Optional: python3 variants/red_brace_xz_force_x_surface_options.py --surface instantaneous

Default: fixed-crest integration, z=-35 to +6 m, as requested.
This fixed wetted extent is an approximation throughout the cycle.
Use --surface swl for the lecture baseline, or --surface instantaneous
for the moving water-surface intersection.
Geometry: 42 m rise, 30 m horizontal span; z=-35 to +7 m; d=40 m.
The default brace rises in +x, matching the posted equations. For the
falling-in-+x brace shown on the front/rear face of the coursework figure,
use --slope falling. x=0 is the brace's SWL crossing in either case.

Above-SWL kinematics in instantaneous mode use direct Airy extrapolation.
This is an approximation discussed in the lecture, not nonlinear theory.
"""
import argparse
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

RHO = 1025.0
CD = 1.2
CM = 0.8  # Total inertia coefficient specified in the coursework.
D = 1.0
A = 6.0
DEPTH = 40.0  # Water depth; not the brace's lower endpoint depth.
K = 0.043
OMEGA = 0.628  # Rounded values retained from the posted code.
RISE = 42.0
SPAN = 30.0
Z_BOTTOM = -35.0
Z_TOP = 7.0
ALPHA = np.arctan2(RISE, SPAN)
S = np.sin(ALPHA)
C = np.cos(ALPHA)
M = SPAN / RISE  # dx/dz for the rising brace.
PERIOD = 2 * np.pi / OMEGA
B = A * OMEGA / np.sinh(K * DEPTH)


def kinematics(z, phase, slope=1):
    """phase = omega*t at x=0; u rightwards, w upwards."""
    x = slope * M * z
    theta = phase - K * x
    r = K * (z + DEPTH)
    u = B * np.cosh(r) * np.sin(theta)
    w = B * np.sinh(r) * np.cos(theta)
    ax = OMEGA * B * np.cosh(r) * np.cos(theta)
    az = -OMEGA * B * np.sinh(r) * np.sin(theta)
    # Unit normal n=(sin(alpha), 0, -slope*cos(alpha)).
    un = S * u - slope * C * w
    an = S * ax - slope * C * az
    return un, an


def upper_limit(phase, surface='swl', slope=1):
    if surface == 'swl':
        return 0.0
    if surface == 'crest':
        return min(A, Z_TOP)
    # Intersection of z=eta(x,t) with x=slope*M*z.
    zs = brentq(lambda z: z - A*np.sin(phase-K*slope*M*z), -A, A)
    return min(zs, Z_TOP)


def forces(phase, surface='swl', slope=1):
    """Return signed horizontal drag and inertia in N.

    ds=dz/sin(alpha); normal-to-horizontal projection multiplies by
    sin(alpha). These cancel, leaving the z-integrals below. This does
    NOT remove the sin/cos factors already present in un and an.
    """
    top = upper_limit(phase, surface, slope)
    if top <= Z_BOTTOM:
        return np.array([0.0, 0.0])
    def drag_integrand(z):
        un, _ = kinematics(z, phase, slope)
        return un * abs(un)  # Retain the sign; do not replace globally by un**2.
    fd = 0.5*RHO*CD*D*quad(drag_integrand, Z_BOTTOM, top,
                          epsabs=1e-7, epsrel=1e-8, limit=150)[0]
    fi = RHO*CM*np.pi*D**2/4*quad(
        lambda z: kinematics(z, phase, slope)[1], Z_BOTTOM, top,
        epsabs=1e-7, epsrel=1e-8)[0]
    return np.array([fd, fi])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--surface', choices=['crest', 'swl', 'instantaneous'], default='crest')
    parser.add_argument('--slope', choices=['rising', 'falling'], default='rising')
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args()
    slope = 1 if args.slope == 'rising' else -1
    phases = np.linspace(0, 2*np.pi, 2000, endpoint=False)
    values = np.array([forces(p, args.surface, slope) for p in phases])
    step = phases[1] - phases[0]
    print(f'Brace angle = {np.degrees(ALPHA):.5f} deg; T = {PERIOD:.5f} s')
    print(f'Model: {args.surface}; brace {args.slope} in +x; horizontal forces')
    for j, name in enumerate(['Drag', 'Inertia']):
        for sign, label in [(1, 'maximum'), (-1, 'minimum')]:
            guess = phases[np.argmax(sign*values[:, j])]
            result = minimize_scalar(
                lambda p: -sign*forces(p, args.surface, slope)[j],
                bounds=(guess-step, guess+step), method='bounded',
                options={'xatol': 1e-10})
            phase = result.x % (2*np.pi)
            f = forces(phase, args.surface, slope)[j]
            top = upper_limit(phase, args.surface, slope)
            print(f'{name:7s} {label}: {f/1000:10.4f} kN; '
                  f't={phase/OMEGA:.5f} s; omega*t={np.degrees(phase):.5f} deg; '
                  f'z upper={top:.5f} m')
    if not args.no_plot:
        import matplotlib.pyplot as plt
        ts = np.r_[phases, 2*np.pi] / OMEGA
        ys = np.vstack([values, values[0]]) / 1000
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(ts, ys[:, 0], label='Horizontal drag')
        ax.plot(ts, ys[:, 1], label='Horizontal inertia')
        ax.axhline(0, color='grey', linewidth=.7)
        ax.set(xlabel='Time (s)', ylabel='Horizontal force (kN)',
               title=f'Diagonal brace: {args.surface} integration')
        ax.grid(alpha=.25)
        ax.legend()
        fig.tight_layout()
        plt.show()


if __name__ == '__main__':
    main()
