# Jacket wave loading

Python calculations and an interactive jacket structure viewer developed for our fluid-loading coursework.

**[Open the interactive structure viewer](https://jacket-structure-viewer.y655375436.chatgpt.site)**

## Naming convention

- `along_x` or `along_y` identifies the beam's orientation.
- `xz` or `yz` identifies the brace's plane.
- `force_x`, `force_y` or `force_z` identifies the force component.
- `from_w` means the contribution caused by vertical wave velocity; `from_u` means horizontal wave velocity.
- `from_u_and_w` means the saved calculation uses both velocities, including a combined normal speed in its drag term.

Every main script includes drag, inertia and their simultaneous total. For the red brace, `forces_x_and_z` means both coordinate directions.

## The eight reviewed scripts

These are the files supplied from the coursework folder, with clearer names. Their contents and calculations are unchanged by this organisation update.

| Script | Member | Force output | Velocity dependence |
| --- | --- | --- | --- |
| [beam_along_y_force_x_from_u_and_w.py](beam_along_y_force_x_from_u_and_w.py) | Lower beam along y | Horizontal x | Drag uses u and w; inertia uses horizontal acceleration. |
| [beam_along_x_force_z_from_w.py](beam_along_x_force_z_from_w.py) | Lower beam along x | Vertical z | w and its acceleration; phase varies along x. |
| [beam_along_y_force_z_from_w.py](beam_along_y_force_z_from_w.py) | Lower beam along y | Vertical z | w and its acceleration; constant x. |
| [red_brace_xz_forces_x_and_z.py](red_brace_xz_forces_x_and_z.py) | Brace rising in +x, in the x-z plane | Horizontal x and vertical z | Normal velocity and acceleration projected from u and w. |
| [yellow_brace_yz_force_y_from_w.py](yellow_brace_yz_force_y_from_w.py) | Front yellow brace in the y-z plane | Horizontal y | w-only normal force, resolved into y. |
| [yellow_brace_yz_force_z_from_w.py](yellow_brace_yz_force_z_from_w.py) | Front yellow brace in the y-z plane | Vertical z | w-only normal force, resolved into z. |
| [yellow_brace_yz_force_x_from_u_and_w.py](yellow_brace_yz_force_x_from_u_and_w.py) | Yellow brace in the y-z plane | Horizontal x | Drag uses u and w; inertia uses horizontal acceleration. |
| [vertical_column_force_x_from_u.py](vertical_column_force_x_from_u.py) | Vertical column along z | Horizontal x | u and its acceleration. |

The [variants folder](variants/README.md) contains earlier implementations, surface options and the combined-output yellow-brace script. These scripts record different models and stages of the work; they are not a single consolidated final solution.

## Run a calculation

Install the required libraries once:

```sh
python -m pip install -r requirements.txt
```

Open a script in Spyder and press Run, or run it from this folder:

```sh
python yellow_brace_yz_force_y_from_w.py
python yellow_brace_yz_force_z_from_w.py
```

Results appear in the console and plots open in a window. Functions return N; tables and plots generally use kN. Total force means drag plus inertia at the same instant, followed by a search for the maximum.

`index.html` is the standalone viewer: download it and open it in a browser. Its sliders illustrate velocity and force projections; they do not calculate wave-load maxima.

## Geometry, coordinates and model status

The viewer follows a 30 m by 30 m plan, 40 m water depth, brace elevations -35 m to +7 m, and a 42 m rise over a 30 m span. The origin is on the front-left column at SWL. +x points towards the back-left column and is the wave direction; +y points towards the front-right column; +z points upward.

The front yellow brace is at x = 0 and falls in +y, from `(0, 0, 7)` to `(0, 30, -35)`. For this geometry, the normal contribution from upward w points towards +y and +z.

### Yellow brace: separate the force caused by w into y and z

`yellow_brace_yz_force_y_from_w.py` reports the horizontal y component; `yellow_brace_yz_force_z_from_w.py` reports the vertical z component. Both drag and inertia are first calculated in the normal direction using the existing w-only scalar model, then resolved geometrically. This does not introduce a new hydrodynamic model.

With alpha = atan(42/30) = 54.4623 degrees, the positive unit normal for the front brace is `(0, sin(alpha), cos(alpha))`. For each of drag, inertia and their simultaneous total:

- `Fy_from_w = Fn_from_w * sin(alpha)` (horizontal, along y).
- `Fz_from_w = Fn_from_w * cos(alpha)` (vertical, along z).
- `Fx_from_w = 0` for this normal-force contribution.

The positive y and z peaks of each force term occur together because they are positive constant multiples of the same signed normal force. Drag, inertia and total have different peak phases. Using the saved fixed-crest parameters at x = 0:

| Force caused by w | Maximum +y (kN) | Maximum +z (kN) | Local phase (degrees) |
| --- | ---: | ---: | ---: |
| Drag | 52.04 | 37.17 | 0.00 |
| Inertia | 20.29 | 14.49 | 270.00 |
| Total at the same instant | 54.02 | 38.58 | 348.76 |

The separate scripts plot y and z individually and print the simultaneous drag/inertia values at maximum total. The y signs now match the viewer; normal and z magnitudes and integration limits are unchanged.

**Earlier versions:** the rising-y legacy script is retained in [variants](variants/README.md), with its positive normal towards -y and +z. The two root w-only scripts match the corrected falling-y front brace. The red x-z script still sets x=0 at the brace SWL crossing, so its reference phase differs from the viewer's origin.

**Model differences retained:** the horizontal beam script uses drag proportional to `u*sqrt(u**2+w**2)`. The yellow horizontal-x script uses `u*sqrt(u**2+(w*cos(alpha))**2)`. Their names therefore say `from_u_and_w`; they are not the separate u-only `u*abs(u)` model. The two yellow w-only scripts retain the separate-component calculation discussed in the coursework work. This naming update does not change or reconcile these hydrodynamic models.

The brace scripts generally use the requested fixed integration limits z = -35 to +6 m for both drag and inertia throughout the cycle. This is a fixed-crest approximation, not an instantaneous wetted-length model. Above-SWL evaluations extend the linear-wave kinematics. The surface-options script in `variants/` exposes alternative surface choices; select them explicitly rather than comparing different choices as if they were equivalent.

Most scripts use rho = 1025 kg/m³, Cd = 1.2, Cm = 0.8, a = 6 m, k = 0.043 m⁻¹ and omega = 0.628 rad/s. The saved scripts themselves are the source of their parameter values.

The viewer demonstrates geometry and illustrative velocity/force components. Its sliders do not calculate wave-load maxima.

## Working together

Create a branch for a change, explain which assumptions or parameters changed, and compare results using the same coordinates, surface treatment and force component before merging.
