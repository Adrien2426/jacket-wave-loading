# Jacket wave loading

Python calculations and an interactive jacket structure viewer developed for our fluid-loading coursework.

**[Open the interactive structure viewer](https://jacket-structure-viewer.y655375436.chatgpt.site)**

## What is included

- Seven `.py` files: saved scripts from the calculation development.
- `index.html`: the standalone interactive viewer. Download and open it in a browser; it also works offline.
- `requirements.txt`: Python libraries needed to run the calculations.

These scripts record several models and stages of the work. They are not a single consolidated final solution. Read the assumptions and coordinate notes below before comparing results.

## Run a calculation

Install the required libraries once:

```sh
python -m pip install -r requirements.txt
```

Open a script in Spyder and press Run, or run it from this folder:

```sh
python yellow_brace_xt_numpy.py
```

Results appear in the console and plots open in a window. Forces returned by functions are in N; result tables and plots generally use kN. Total force means drag plus inertia at the same instant, followed by a search for the maximum.

## Script guide

| Script | Purpose and assumptions |
| --- | --- |
| `beam2_vertical_numpy.py` | Vertical force on the horizontal beam along x, at z = -35 m; signed scalar drag from w. |
| `brace_forces_spyder.py` | Drag, inertia and total x/z forces on a brace rising in +x; normal velocity includes u and w; fixed crest limits. Its x origin is the brace's SWL crossing. |
| `brace_forces_checked.py` | x-z brace calculation with crest, SWL and instantaneous-surface options; rising/falling options. Its x origin is also the brace's SWL crossing. The filename does not indicate final coursework approval. |
| `diagonal_handwritten_numpy.py` | Common-phase, u-only horizontal-force calculation for a y-z brace, using the requested rounded length conversion 1.23. |
| `yellow_brace_xt_numpy.py` | Normal/y/z forces caused by w, as functions of constant brace position x and time t; default x = 0. Uses the earlier rising-y sign convention noted below. |
| `yellow_brace_vertical_velocity_numpy.py` | Earlier w-only normal/y/z calculation at x = 30; local and reference phases reported separately. Uses the earlier rising-y sign convention. |
| `yellow_brace_numpy.py` | Earlier alternative using combined normal speed from u and w in the drag law. This is not the separate-component lecture-style model used in the later calculations. |

## Geometry, coordinates and model status

The viewer follows a 30 m by 30 m plan, 40 m water depth, brace elevations -35 m to +7 m, and a 42 m rise over a 30 m span. The origin is on the front-left column at SWL. +x points towards the back-left column and is the wave direction; +y points towards the front-right column; +z points upward.

The front yellow brace is at x = 0 and falls in +y, from `(0, 0, 7)` to `(0, 30, -35)`. For this geometry, the normal contribution from upward w points towards +y and +z.

**Coordinate reconciliation still needed:** the two saved w-only yellow-brace scripts retain the earlier rising-y convention, with the positive normal towards -y and +z. Their y-component signs do not match the corrected viewer. The files are preserved as developed; reverse the y-component signs and update their normal-direction labels when adopting the viewer's falling-y geometry. The normal-force and z-force magnitudes are unchanged by that orientation reversal. The x-z scripts use a different horizontal origin from the viewer, so reference phases also need a consistent origin before being compared.

The brace scripts generally use the requested fixed integration limits z = -35 to +6 m for both drag and inertia throughout the cycle. This is a fixed-crest approximation, not an instantaneous wetted-length model. Above-SWL evaluations extend the linear-wave kinematics. The checked script exposes alternative surface choices; select them explicitly rather than comparing different choices as if they were equivalent.

Most scripts use rho = 1025 kg/m³, Cd = 1.2, Cm = 0.8, a = 6 m, k = 0.043 m⁻¹ and omega = 0.628 rad/s. The saved scripts themselves are the source of their parameter values.

The viewer demonstrates geometry and illustrative velocity/force components. Its sliders do not calculate wave-load maxima.

## Working together

Create a branch for a change, explain which assumptions or parameters changed, and compare results using the same coordinates, surface treatment and force component before merging.
